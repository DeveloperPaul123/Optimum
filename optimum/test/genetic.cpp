#include "optimum/evolutionary/genetic.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <random>
#include <string>

#include "optimum/evolutionary/selection.h"
#include "optimum/evolutionary/termination.h"

struct random_word_generator {
    random_word_generator() : generator_(device_()) {}
    [[nodiscard]] std::string operator()(std::string_view char_set, std::size_t max_length) {
        std::string out(max_length, '\0');
        std::generate_n(out.begin(), max_length, [&]() {
            std::uniform_int_distribution dist{{}, char_set.size()-1};

            return char_set[dist(generator_)];
        });
        return out;
    }

  private:
    std::random_device device_;
    std::mt19937 generator_;
};

std::string get_half_string(std::string str) {
    auto mid_point = str.size() / 2;
    return str.substr(0, mid_point);
}

std::string get_second_half_string(std::string str) {
    auto mid_point = str.size() / 2;
    return str.substr(mid_point + 1);
}

TEST(GeneticAlgorithm, HelloWorld) {
    const std::string solution = "hello world";
    const std::string alphabet = "abcdefghijklmnopqrstuvwxyz";
    const std::string available_chars = alphabet + " ";

    const auto word_length = solution.size();
    random_word_generator word_generator{};
    // generate initial population
    constexpr auto initial_pop_size = 1000;
    std::vector<std::string> initial_population(initial_pop_size);
    std::generate_n(initial_population.begin(), initial_pop_size,
                    [&]() { return word_generator(available_chars, word_length); });

    auto fitness_op = [solution, available_chars](std::string value) -> double {
        double score = 0.0;
        const auto length = std::min(value.size(), solution.size());
        for (auto i = 0; i < length; i++) {
            score += solution[i] == value[i];
        }
        long diff = solution.size() - value.size();
        if (solution.size() != value.size()) score -= std::abs(diff) * 2;
        return score;
    };

    std::random_device device;
    std::mt19937 engine(device());

    auto string_mutator = [available_chars, &engine](std::string value) -> std::string {
        std::uniform_int_distribution<> distribution(0, value.size() - 1);
        auto random_index = distribution(engine);
        auto return_string = value;
        // now get random character
        std::uniform_int_distribution<> char_dist(0, available_chars.size() - 1);
        auto random_char_index = char_dist(engine);

        return_string[random_index] = available_chars[random_char_index];
        return return_string;
    };

    auto cross_over = [&engine](std::string first, std::string second) {
        std::string first_half = get_half_string(first);
        std::string second_half = get_second_half_string(second);
        return first_half + second_half;
    };

    auto termination =
        dp::ga::termination::fitness_termination_criteria<std::string>(fitness_op(solution));
    dp::ga::algorithm_settings settings{0.1, 0.5, 0.25};
    dp::genetic_algorithm<std::string> genetics(initial_population, string_mutator, cross_over,
                                                fitness_op, termination, settings);

    auto results = genetics.solve();
    EXPECT_EQ(results.best, solution);
}

TEST(GeneticSelection, RouletteSelection) {
    const std::string test_value = "test";
    const std::string alphabet = "abcdefghijklmnopqrstuvwxyz";
    const std::string available_chars = alphabet + " ";

    const auto word_length = test_value.size();
    random_word_generator word_generator{};
    // generate initial population
    constexpr auto initial_pop_size = 1000;
    std::vector<std::string> initial_population(initial_pop_size);
    std::generate_n(initial_population.begin(), initial_pop_size,
                    [&]() { return word_generator(available_chars, word_length); });

    auto fitness_op = [test_value, available_chars](std::string value) -> double {
        double score = 0.0;
        for (auto i = 0; i < test_value.size(); i++) {
            score += test_value[i] == value[i];
        }
        return score;
    };

    dp::ga::selection::roulette_selection selection{};
    const auto& [parent1, parent2] = selection(initial_population, fitness_op);
    EXPECT_NE(parent1, parent2);
}