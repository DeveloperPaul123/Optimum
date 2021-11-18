#include "optimum/evolutionary/genetic.h"
#include "optimum/evolutionary/termination.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <random>
#include <string>

struct random_word_generator {
    random_word_generator() : generator_(device_()) {}
    [[nodiscard]] std::string operator()(std::string_view char_set, std::size_t max_length) {
        std::string out;
        std::ranges::sample(char_set, std::back_inserter(out), max_length, generator_);
        return out;
    }

  private:
    std::random_device device_;
    std::mt19937 generator_;
};

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
        for (auto i = 0; i < solution.size(); i++) {
            score += solution[i] == value[i];
        }
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

    auto cross_over = [](std::string first, std::string second) {
        // assume that string sizes are the same.
        assert(first.size() == second.size());
        std::string::size_type str_size = first.size();
        auto is_odd = str_size % 2 == 1;
        std::string::size_type half_size = str_size / 2;

        std::string first_half = first.substr(0, half_size);
        std::string second_half = second.substr(second.size() - half_size - 1);
        return first_half + second_half;
    };

    auto termination =
        dp::ga::termination::fitness_termination_criteria<std::string>(fitness_op(solution));
    dp::genetic_algorithm<std::string>::algorithm_settings settings{0.40, 0.5, 0.2};
    dp::genetic_algorithm<std::string> ga(initial_population, string_mutator, cross_over,
                                          fitness_op, termination, settings);

    auto results = ga.solve();
    EXPECT_EQ(results.best, solution);
}
