#include "optimum/evolutionary/genetic.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <random>
#include <string>
#include <unordered_map>

#include "optimum/evolutionary/selection.h"
#include "optimum/evolutionary/termination.h"

struct random_word_generator {
    random_word_generator() : generator_(device_()) {}
    [[nodiscard]] std::string operator()(std::string_view char_set, std::size_t max_length) {
        std::string out(max_length, '\0');
        std::generate_n(out.begin(), max_length, [&]() {
            std::uniform_int_distribution dist{{}, char_set.size() - 1};

            return char_set[dist(generator_)];
        });
        return out;
    }

  private:
    std::random_device device_;
    std::mt19937 generator_;
};

std::string get_half_string(std::string str) {
    const auto mid_point = str.size() / 2;
    return str.substr(0, mid_point);
}

std::string get_second_half_string(std::string str) {
    const auto mid_point = str.size() / 2;
    return str.substr(mid_point + 1);
}

TEST(GeneticAlgorithm, HelloWorld) {
    const std::string solution = "hello world from a genetic algorithm";
    const std::string alphabet = "abcdefghijklmnopqrstuvwxyz";
    const std::string available_chars = alphabet + " ,'";

    const auto word_length = solution.size();
    random_word_generator word_generator{};
    // generate initial population
    constexpr auto initial_pop_size = 1000;
    std::vector<std::string> initial_population(initial_pop_size);
    std::generate_n(initial_population.begin(), initial_pop_size,
                    [&]() { return word_generator(available_chars, word_length); });

    auto fitness_op = [solution](std::string value) -> double {
        double score = 0.0;
        const auto length = std::min(value.size(), solution.size());
        for (auto i = 0; i < length; i++) {
            score += solution[i] == value[i];
        }
        const long diff = solution.size() - value.size();
        if (solution.size() != value.size()) score -= std::abs(diff) * 2;
        return score;
    };

    std::random_device device;
    std::mt19937 engine(device());

    auto string_mutator = [available_chars, &engine, &word_generator,
                           word_length](const std::string& value) -> std::string {
        if (value.empty()) {
            return word_generator(available_chars, word_length);
        }

        const auto value_max = value.size() > 1 ? value.size() - 1 : 1;
        const std::uniform_int_distribution<> distribution(0, value_max);
        const auto random_index = distribution(engine);
        auto return_string = value;
        // now get random character
        const std::uniform_int_distribution<> char_dist(0, available_chars.size() - 1);
        const auto random_char_index = char_dist(engine);

        return_string[random_index] = available_chars[random_char_index];
        return return_string;
    };

    auto cross_over = [](std::string first, std::string second) {
        const auto first_string_midpoint = first.size() / 2;
        const auto max_size = std::max(first.size(), second.size());
        const std::string first_half = first.substr(0, first_string_midpoint);
        const auto remaining_length = max_size - first_half.size();
        const unsigned long long second_half_start = second.size() - remaining_length;
        const std::string second_half = second.substr(second_half_start);
        return first_half + second_half;
    };

    auto termination =
        dp::ga::termination::fitness_termination_criteria<std::string>(fitness_op(solution));
    dp::ga::algorithm_settings settings{0.1, 0.5, 0.25};
    dp::genetic_algorithm genetics(settings, initial_population, string_mutator,
                                                cross_over, fitness_op, termination);
    auto results = genetics.solve();
    EXPECT_EQ(results.best, solution);
}

template <typename Selector>
std::unordered_map<std::string, int> test_selection(
    Selector& selector, const std::string& test_value,
    const std::vector<std::string>& initial_population) {
    auto fitness_op = [test_value](const std::string& value) -> double {
        double score = 0.0;
        score -=
            std::abs(static_cast<double>(test_value.size()) - static_cast<double>(value.size())) *
            10.0;
        for (auto i = 0; i < test_value.size(); i++) {
            score += test_value[i] == value[i];
        }
        return score;
    };

    std::unordered_map<std::string, int> selection_histogram;

    // do 1000 selections and build up a selection histogram
    for (auto i = 0; i < 1000; i++) {
        const auto& [parent1, parent2] = selector(initial_population, fitness_op);
        if (!selection_histogram.contains(parent1)) {
            selection_histogram[parent1] = 0;
        } else {
            ++selection_histogram[parent1];
        }

        if (!selection_histogram.contains(parent2)) {
            selection_histogram[parent2] = 0;
        } else {
            ++selection_histogram[parent2];
        }
    }

    return selection_histogram;
}

TEST(GeneticSelection, RouletteSelection) {
    const std::string test_value = "test";

    // create a sample population where 1 member has a much larger fitness compared to the others.
    const std::vector<std::string> initial_population{"tesa", "aaaa", "bbbb", "aaa", "bbb"};

	dp::ga::selection::roulette_selection selection{};
    const auto selection_histogram = test_selection(selection, test_value, initial_population);

    for (auto& [string, count] : selection_histogram) {
        std::cout << string << " : " << std::to_string(count) << "\n";
    }

    const auto [string_value, count] = *std::ranges::max_element(
        selection_histogram, [](auto first, auto second) { return first.second < second.second; });
    EXPECT_EQ(string_value, "tesa");
}

TEST(GeneticSelection, RankSelection) {
    dp::ga::selection::rank_selection selector{};
    const std::string test_value = "test";

    // create a sample population where 1 member has a much larger fitness compared to the others.
    const std::vector<std::string> initial_population{"tesa", "aaaa", "bbbb", "aaa", "bbb"};

    const auto selection_histogram = test_selection(selector, test_value, initial_population);

    for (auto& [string, count] : selection_histogram) {
        std::cout << string << " : " << std::to_string(count) << "\n";
    }

    constexpr auto comp_op = [](auto first, auto second) { return first.second < second.second; };
    const auto [string_value, count] = *std::ranges::max_element(selection_histogram, comp_op);
    const auto [min_value, min_count] = *std::ranges::min_element(selection_histogram, comp_op);

    EXPECT_EQ(string_value, "tesa");
    EXPECT_EQ(min_value, "bbb");
}
