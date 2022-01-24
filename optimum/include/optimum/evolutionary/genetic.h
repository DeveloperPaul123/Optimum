#pragma once

#include <algorithm>
#include <concepts>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <vector>

#include "optimum/evolutionary/selection.h"

namespace dp {
    namespace detail {
        template <class Fn, class T>
        concept mutation_operator =
            std::invocable<Fn, T> && std::is_same_v<std::invoke_result_t<Fn, T>, T>;

        template <class Fn, class T, class Result = std::invoke_result_t<Fn, T>>
        concept fitness_operator = std::invocable<Fn, T> &&
            (std::integral<Result> || std::floating_point<Result>);

        template <class Fn, class T, class Result = std::invoke_result_t<Fn, T, T>>
        concept crossover_operator = std::invocable<Fn, T, T> && std::is_same_v<Result, T>;

        template <class Fn, class T, class Numeric,
                  class Result = std::invoke_result_t<Fn, T, Numeric>>
        concept termination_operator = std::invocable<Fn, T, Numeric> &&
            std::is_same_v<Result, bool> &&
            ((std::integral<Numeric> || std::floating_point<Numeric>));
    }  // namespace detail

    template <typename T>
    struct default_mutator {
        T operator()(T&& t) { return t; }
    };

    template <typename T>
    struct default_crossover {
        T operator()(T&& first, T&& second) { return first + second; }
    };

    template <typename T, typename U = double>
    requires std::floating_point<U> || std::integral<U>
    struct default_evaluation {
        constexpr U operator()(T&& t) { return std::accumulate(t.begin(), t.end(), {}); }
    };

    namespace ga {
        /// @brief Settings type for probabilities
        struct algorithm_settings {
            double elitism = 0.0;
            double mutation_rate = 0.5;
            double crossover_rate = 0.2;
        };
    }  // namespace ga

    template <class ChromosomeType>
    requires std::default_initializable<ChromosomeType>
    class genetic_algorithm {
      public:
        using value_type = ChromosomeType;
        using list_type = std::vector<value_type>;

        /// @brief Result type that holds the best result
        struct results {
            value_type best;
            double fitness{};
        };

        struct iteration_statistics {
            results current_best;
            std::size_t current_generation_count{};
            std::size_t population_size{};
        };

        template <class FitnessOperator, class MutationOperator, class CrossoverOperator,
                  class TerminationOperator>
        requires detail::mutation_operator<MutationOperator, ChromosomeType> &&
            detail::fitness_operator<FitnessOperator, ChromosomeType> &&
            detail::crossover_operator<CrossoverOperator, ChromosomeType> &&
            detail::termination_operator<TerminationOperator, ChromosomeType,
                                         std::invoke_result_t<FitnessOperator, ChromosomeType>>
            genetic_algorithm(ga::algorithm_settings settings, list_type initial_population,
                              MutationOperator&& mutator, CrossoverOperator&& crossover_operator,
                              FitnessOperator&& fitness_operator,
                              TerminationOperator&& termination_operator)
            : mutator_(std::forward<MutationOperator>(mutator)),
              crossover_(std::forward<CrossoverOperator>(crossover_operator)),
              fitness_(std::forward<FitnessOperator>(fitness_operator)),
              termination_(std::forward<TerminationOperator>(termination_operator)),
              settings_(settings) {
            // initialize our population
            std::ranges::transform(initial_population, std::back_inserter(population_),
                                   [this](ChromosomeType value) {
                                       return chromosome_metadata{value, fitness_(value)};
                                   });
            // sort by fitness
            std::ranges::sort(population_, fitness_sort_op{});
        }

        template <typename IterationCallback = std::function<void(const iteration_statistics&)>>
        requires std::invocable<IterationCallback, const iteration_statistics&>
        [[nodiscard]] results solve(const IterationCallback& callback =
                                        [](const iteration_statistics&) {}) {
            auto best_element = std::ranges::max_element(
                population_, [](chromosome_metadata first, chromosome_metadata second) {
                    return first.fitness < second.fitness;
                });

            // TODO: Make selector configurable
            // default selector
            dp::ga::selection::rank_selection selector{};

            iteration_statistics stats{};
            stats.current_best.best = best_element->value;

            while (!termination_(best_element->value, best_element->fitness)) {
                // perform elitism
                auto number_elitism = static_cast<std::size_t>(
                    std::round(static_cast<double>(population_.size()) * settings_.elitism));
                // ensure we do elitism selection if it is enabled
                if (number_elitism == 0 && settings_.elitism > 0.0) number_elitism = 2;

                auto elite_population = elitism(population_, number_elitism);

                // cross over
                auto crossover_number = static_cast<std::size_t>(
                    std::round(static_cast<double>(population_.size()) * settings_.crossover_rate));
                if (crossover_number <= 1) crossover_number = 4;

                // create a new "generation", we will also insert the elite population into this
                // one.
                population crossover_population;
                crossover_population.reserve(crossover_number * 2 + elite_population.size());

                for (std::size_t i = 0; i < crossover_number; i++) {
                    // randomly select 2 parents
                    const auto& [parent1, parent2] = selector(
                        population_,
                        [this](const chromosome_metadata& value) { return fitness_(value.value); });

                    // generate two children from each parent sets
                    auto child1 = crossover_(parent1.value, parent2.value);
                    auto child2 = crossover_(parent2.value, parent1.value);

                    child1 = mutator_(child1);
                    child2 = mutator_(child2);

                    // insert them into the new population
                    // do not evaluate fitness yet
                    crossover_population.push_back({child1, fitness_(child1)});
                    crossover_population.push_back({child2, fitness_(child2)});
                }

                // add elite population
                crossover_population.insert(crossover_population.end(), elite_population.begin(),
                                            elite_population.end());

                // sort crossover population by fitness (lowest first)
                std::ranges::sort(crossover_population, fitness_sort_op{});

                // reset the current population
                population_.clear();
                // assign/move the new generation
                population_ = std::move(crossover_population);

                // update the best element
                best_element = std::ranges::max_element(population_, fitness_sort_op{});

                // send callback stats for each generation
                stats.current_best.best = best_element->value;
                stats.current_best.fitness = best_element->fitness;
                stats.population_size = population_.size();
                ++stats.current_generation_count;
                callback(std::add_const_t<iteration_statistics>(stats));
            }

            return {best_element->value, best_element->fitness};
        }

      private:
        using chromosome = ChromosomeType;
        using mutator_operator = std::function<ChromosomeType(ChromosomeType)>;
        using crossover_operator = std::function<ChromosomeType(ChromosomeType, ChromosomeType)>;
        using fitness_evaluator = std::function<double(ChromosomeType)>;
        using termination_criterion = std::function<bool(ChromosomeType, double)>;

        struct chromosome_metadata {
            chromosome value{};
            double fitness{};
            auto operator<=>(const chromosome_metadata&) const = default;
        };

        template <typename Comparator = std::less<>>
        struct fitness_sort_op {
            bool operator()(const chromosome_metadata& first, const chromosome_metadata& second) {
                Comparator cmp;
                return cmp(first.fitness, second.fitness);
            }
        };

        using population = std::vector<chromosome_metadata>;
        population population_;
        mutator_operator mutator_;
        crossover_operator crossover_;
        fitness_evaluator fitness_;
        termination_criterion termination_;

        ga::algorithm_settings settings_{};

        static auto elitism(population& current_population, std::size_t number_elitism) {
            // perform elitism selection

            // sort so that largest fitness item is at front
            std::ranges::partial_sort(current_population,
                                      current_population.begin() +
                                          std::min(current_population.size(), number_elitism + 1),
                                      fitness_sort_op<std::greater<>>{});

            // select the first n in the current population
            return std::ranges::views::take(current_population, number_elitism);
        }
    };
}  // namespace dp
