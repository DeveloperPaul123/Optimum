#pragma once

#include <algorithm>
#include <concepts>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <ranges>
#include <vector>

namespace dp {
    namespace detail {
        template <class Fn, class T>
        concept mutation_operator =
            std::invocable<Fn, T> && std::is_same_v<std::invoke_result_t<Fn, T>, T>;

        template <class Fn, class T, class Result = std::invoke_result_t<Fn, T>>
        concept fitness_operator = std::invocable<Fn, T> && std::is_convertible_v<Result, double>;

        template <class Fn, class T, class Result = std::invoke_result_t<Fn, T, T>>
        concept crossover_operator = std::invocable<Fn, T, T> && std::is_same_v<Result, T>;
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
    class genetic_algorithm {
      public:
        using value_type = ChromosomeType;
        using list_type = std::vector<value_type>;

        /// @brief Result type that holds the best result
        struct results {
            genetic_algorithm::value_type best;
            double fitness{};
        };

        template <class FitnessOperator, class MutationOperator, class CrossoverOperator,
                  class TerminationOperator>
        requires detail::mutation_operator<MutationOperator, ChromosomeType> &&
            detail::fitness_operator<FitnessOperator, ChromosomeType> &&
            detail::crossover_operator<CrossoverOperator, ChromosomeType>
            genetic_algorithm(list_type initial_population, MutationOperator&& mutator,
                              CrossoverOperator&& crossover_operator,
                              FitnessOperator&& fitness_operator,
                              TerminationOperator&& termination_operator,
                              ga::algorithm_settings settings = ga::algorithm_settings{})
            : mutator_(std::forward<MutationOperator>(mutator)),
              crossover_(std::forward<CrossoverOperator>(crossover_operator)),
              fitness_(std::forward<FitnessOperator>(fitness_operator)),
              termination_(std::forward<TerminationOperator>(termination_operator)),
              generator_(device_()), settings_(settings) {
            // initialize our population
            std::ranges::transform(initial_population, std::back_inserter(population_),
                                   [this](ChromosomeType value) {
                                       return chromosome_metadata{value, fitness_(value)};
                                   });
        }

        [[nodiscard]] results solve() {
            auto best_element = std::ranges::max_element(
                population_, [](chromosome_metadata first, chromosome_metadata second) {
                    return first.fitness < second.fitness;
                });

            while (!termination_(best_element->value, best_element->fitness)) {
                // perform elitism
                std::size_t number_elitism =
                    (std::size_t)std::round(population_.size() * settings_.elitism);
                // ensure we do elitism selection if it is enabled
                if (number_elitism == 0 && settings_.elitism > 0.0) number_elitism = 1;

                auto elite_population = elitism(population_, number_elitism);

                // cross over
                std::size_t crossover_number =
                    (std::size_t)std::round(population_.size() * settings_.crossover_rate);
                if (crossover_number == 0) crossover_number = 2;

                population crossover_population;
                crossover_population.reserve(crossover_number * 2);

                population parents(2);

                for (std::size_t i = 0; i < crossover_number; i++) {
                    // randomly select 2 parents
                    // TODO: Introduce selection criteria
                    std::ranges::sample(population_, parents.begin(), 2, generator_);
                    // generate two children from each parent sets
                    auto child1 = crossover_(parents[0].value, parents[1].value);
                    auto child2 = crossover_(parents[1].value, parents[0].value);
                    // insert them into the new population
                    // do not evaluate fitness yet
                    crossover_population.push_back({child1, 0.0});
                    crossover_population.push_back({child2, 0.0});
                }

                // mutation in crossover population
                for (auto& [chrm, fitness] : crossover_population) {
                    // mutate based on mutation rate
                    chrm = mutator_(chrm);
                    fitness = fitness_(chrm);
                }

                population_.clear();
                population_.reserve(elite_population.size() + crossover_population.size());
                // update the full population
                population_.insert(population_.end(), elite_population.begin(),
                                   elite_population.end());
                population_.insert(population_.end(), crossover_population.begin(),
                                   crossover_population.end());

                // update the best element
                best_element = std::ranges::max_element(
                    population_, [](chromosome_metadata first, chromosome_metadata second) {
                        return first.fitness < second.fitness;
                    });
            }

            return {best_element->value, best_element->fitness};
        }

      private:
        using chromosome = ChromosomeType;
        using mutator_operator = std::function<ChromosomeType(ChromosomeType)>;
        using crossover_operator = std::function<ChromosomeType(ChromosomeType, ChromosomeType)>;
        using fitness_evaluator = std::function<double(ChromosomeType)>;
        using termination_critereon = std::function<bool(ChromosomeType, double)>;
        struct chromosome_metadata {
            chromosome value;
            double fitness;
        };

        using population = std::vector<chromosome_metadata>;
        population population_;
        mutator_operator mutator_;
        crossover_operator crossover_;
        fitness_evaluator fitness_;
        termination_critereon termination_;

        ga::algorithm_settings settings_{};
        std::random_device device_;
        std::mt19937 generator_;

        population elitism(population& current_population, std::size_t number_elitism) {
            // perform elitism selection

            // sort so that largest fitness item is at front
            std::ranges::sort(current_population,
                              [](chromosome_metadata first, chromosome_metadata second) {
                                  return first.fitness > second.fitness;
                              });

            // select the first n in the current population
            population out;
            out.reserve(number_elitism);
            std::ranges::copy_n(current_population.begin(), number_elitism,
                                std::back_inserter(out));
            return out;
        }
    };
}  // namespace dp
