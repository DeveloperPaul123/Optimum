#pragma once

#include <algorithm>
#include <iterator>
#include <random>
#include <ranges>
#include <utility>
namespace dp::ga::selection {
    struct roulette_selection {
        roulette_selection() : generator_(device_()) {}

        template <typename Container, typename UnaryOperator,
                  typename T = typename Container::value_type,
                  typename FitnessResult = std::invoke_result_t<UnaryOperator, T>>
        std::pair<T, T> operator()(const Container& population, UnaryOperator fitness_op) {
            // generate sum
            FitnessResult sum =
                std::accumulate(population.begin(), population.end(), FitnessResult{},
                                [&](FitnessResult current_sum, const T& value) {
                                    return current_sum + fitness_op(value);
                                });
            // pick 2 parents and return them
            return std::make_pair(pick_one(population, fitness_op, sum),
                                  pick_one(population, fitness_op, sum));
        }

      private:
        std::random_device device_;
        std::mt19937 generator_;
        template <typename Container, typename UnaryOperator,
                  typename T = typename Container::value_type,
                  typename FitnessResult = std::invoke_result_t<UnaryOperator, T>>
        T pick_one(const Container& population, UnaryOperator fitness_op, FitnessResult sum) {
            const std::uniform_real_distribution dis(0.0, 1.0);
            auto random_value = dis(generator_);

            auto threshold = random_value * sum;
            FitnessResult accumulator{};
            for (auto& value : population) {
                accumulator += fitness_op(value);
                if (accumulator >= threshold) {
                    return value;
                }
            }
            return population.at(0);
        }
    };

    struct rank_selection {
        template <typename Container, typename UnaryOperator,
                  typename T = typename Container::value_type,
                  typename FitnessResult = std::invoke_result_t<UnaryOperator, T>>
        std::pair<T, T> operator()(const Container& population, UnaryOperator fitness_op) {
            // TODO: Avoid this copy if possible
            Container copy{};
            copy.reserve(population.size());
            std::copy(population.begin(), population.end(), std::back_inserter(copy));
            // sort by fitness, highest to lowest
            std::sort(copy.begin(), copy.end(), [&](const T& first, const T& second) {
                return fitness_op(first) < fitness_op(second);
            });

            // fitness evaluator/operator
            auto rank_fitness_op = [&](const T& value) -> FitnessResult {
                auto location = std::find(copy.begin(), copy.end(), value);
                return static_cast<FitnessResult>(std::distance(copy.begin(), location) + 1.0);
            };

            // use roulette selection
            return selection_(copy, rank_fitness_op);
        }

      private:
        roulette_selection selection_{};
    };
}  // namespace dp::ga::selection
