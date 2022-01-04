#pragma once

#include <algorithm>
#include <iterator>
#include <random>
#include <ranges>
#include <utility>

namespace dp::ga::selection {
    struct roulette_selection {
        roulette_selection() : generator_(device_()), distribution_(0.0, 1.0) {}

        template <typename Container, typename UnaryOperator,
                  typename T = std::ranges::range_value_t<Container>,
                  typename FitnessResult = std::invoke_result_t<UnaryOperator, T>>
        std::pair<T, T> operator()(const Container& population, UnaryOperator fitness_op) {
            // generate sum
            FitnessResult sum =
                std::accumulate(population.begin(), population.end(), FitnessResult{},
                                [&](FitnessResult current_sum, const T& value) {
                                    return current_sum + fitness_op(value);
                                });

            auto first_value = distribution_(generator_);
            auto second_value = distribution_(generator_);

            auto threshold1 = first_value * sum;
            auto threshold2 = second_value * sum;

            std::pair<T, T> return_pair{};

            auto first_found = false;
            auto second_found = false;

            FitnessResult accumulator{};
            for (const auto& value : population) {
                accumulator += fitness_op(value);
                if (accumulator >= threshold1 && !first_found) {
                    return_pair.first = value;
                    first_found = true;
                }
                if (accumulator >= threshold2 && !second_found) {
                    return_pair.second = value;
                    second_found = true;
                }
                if (first_found && second_found) break;
            }

            // pick 2 parents and return them
            return return_pair;
        }

      private:
        std::random_device device_;
        std::mt19937 generator_;
        std::uniform_real_distribution<double> distribution_;

    };

    struct rank_selection {
        template <typename Container, typename UnaryOperator,
                  typename T = std::ranges::range_value_t<Container>,
                  typename FitnessResult = std::invoke_result_t<UnaryOperator, T>>
        std::pair<T, T> operator()(const Container& population, UnaryOperator fitness_op) {
            // assume population is sorted already
            std::ranges::reverse_view reverse_view(population);

            // fitness evaluator/operator
            auto rank_fitness_op = [&](const T& value) -> FitnessResult {
                auto location = std::find(reverse_view.begin(), reverse_view.end(), value);
                return static_cast<FitnessResult>(std::distance(reverse_view.begin(), location) +
                                                  1.0);
            };

            // use roulette selection
            return selection_(reverse_view, rank_fitness_op);
        }

      private:
        roulette_selection selection_{};
    };
}  // namespace dp::ga::selection
