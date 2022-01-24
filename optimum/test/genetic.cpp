#include "optimum/evolutionary/genetic.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
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

TEST(GeneticAlgorithm, StringSearch) {
    // solution
    const std::string solution = "hello world from a genetic algorithm";

    // define available characters
    const std::string alphabet = "abcdefghijklmnopqrstuvwxyz";
    const std::string available_chars = alphabet + " ,'";

    const auto word_length = solution.size();

    // generate random words as our initial population
    random_word_generator word_generator{};
    // generate initial population
    constexpr auto initial_pop_size = 1000;
    std::vector<std::string> initial_population(initial_pop_size);
    std::generate_n(initial_population.begin(), initial_pop_size,
                    [&]() { return word_generator(available_chars, word_length); });

    // fitness evaluator
    auto fitness_op = [solution](std::string value) -> double {
        double score = 0.0;
        const auto length = std::min(value.size(), solution.size());
        for (auto i = 0; i < length; i++) {
            score += solution[i] == value[i];
        }
        // penalize score if the size of the value is different from the solution
        const long diff = solution.size() - value.size();
        if (solution.size() != value.size()) score -= std::abs(diff) * 2;
        return score;
    };

    std::random_device device;
    std::mt19937 engine(device());

    // random string mutator
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

    // crossover operator (child generator)
    auto cross_over = [](std::string first, std::string second) {
        const auto first_string_midpoint = first.size() / 2;
        const auto max_size = std::max(first.size(), second.size());
        const std::string first_half = first.substr(0, first_string_midpoint);
        const auto remaining_length = max_size - first_half.size();
        const unsigned long long second_half_start = second.size() - remaining_length;
        const std::string second_half = second.substr(second_half_start);
        return first_half + second_half;
    };

    // termination criteria
    auto termination =
        dp::ga::termination::fitness_termination_criteria<std::string>(fitness_op(solution));

    static_assert(
        dp::detail::termination_operator<
            dp::ga::termination::fitness_termination_criteria<std::string>, std::string, double>);

    // algorithm settings
    dp::ga::algorithm_settings settings{0.1, 0.5, 0.25};

    // set up algorithm
    dp::genetic_algorithm genetics(settings, initial_population, string_mutator, cross_over,
                                   fitness_op, termination);

    auto [best, fitness] = genetics.solve([](auto& stats) {
        std::cout << "best: " << stats.current_best.best
                  << " fitness: " << stats.current_best.fitness << "\n";
    });

    EXPECT_EQ(best, solution);
}

// type declaration for knapsack problem
using knapsack = std::array<int, 5>;

// print helper for knapsack 
inline std::ostream& operator<<(std::ostream& out, const knapsack& ks) {
    out << "[ ";
    std::ranges::copy_n(ks.begin(), ks.size() - 1, std::ostream_iterator<int>(out, ", "));
    // print last element
    out << ks[ks.size() - 1];
    out << " ]\n";

    return out;
}

TEST(GeneticAlgorithm, KnapsackProblem) {

    // knapsack problem as described here: https://en.wikipedia.org/wiki/Knapsack_problem

    struct knapsack_box {
        int value;
        int weight;
        auto operator<=>(const knapsack_box&) const = default;
    };

    // weight capacity of our knapsack
    constexpr auto max_weight = 15;

    // available boxes for the knapsack
    std::vector<knapsack_box> available_items = {{4, 12}, {2, 1}, {10, 4}, {1, 1}, {2, 2}};

    // fitness evaluator
    auto fitness = [&available_items, max_weight](const knapsack ks) -> int {
        auto value_sum = 0;
        auto weight_sum = 0;
        for (const auto& index : ks) {
            if (index >= 0 && index < available_items.size()) {
                value_sum += available_items[index].value;
                weight_sum += available_items[index].weight;
            }
        }
        // if the weight is less than the max, it adds to the value
        if (weight_sum > max_weight) value_sum -= 25 * std::abs(weight_sum - max_weight);
        return value_sum;
    };

    std::random_device device;
    std::mt19937 engine(device());

    // random mutation operator
    auto mutator = [&engine, &available_items](const knapsack& ks) {
        knapsack output = ks;

        const std::uniform_int_distribution<int> distribution(0, ks.size() - 1);
        const auto index = distribution(engine);

        if (std::ranges::count(output, -1) > 0) {
            // the output has some empty spaces, so there is the potential for us to add a new
            // number to the output
            const std::uniform_int_distribution<int> item_dist(0, available_items.size() - 1);
            auto new_value = item_dist(engine);
            // only allow unique values
            while (std::ranges::find(output, new_value) != std::end(output)) {
                new_value = item_dist(engine);
            }

            output[index] = new_value;
        } else {
            // the output already has unique numbers in it, so we'll just shuffle it
            std::ranges::shuffle(output, engine);
        }

        return output;
    };

    // crossover operator (i.e. child generator)
    auto crossover = [](const knapsack& first, const knapsack& second) {
        knapsack child{};
        std::ranges::fill(child, -1);

        auto first_copy_end = first.begin() + 3;
        const auto first_negative = std::ranges::find(first, -1);
        if (first_negative != first.end() && first_negative < first_copy_end) {
            first_copy_end = first_negative;
        }

        // copy first elements over.
        std::ranges::copy(first.begin(), first_copy_end, child.begin());

        const auto negative_number_count =
            std::ranges::count_if(child, [](const auto& value) { return value < 0; });

        // find the first negative value in the child knapsack (i.e. it's empty)
        auto child_first_negative = std::ranges::find(child, -1);

        // we need to copy from child_first_negative the first "negative_number_count" numbers that
        // are not already in the child knapsack
        for (const auto& value : second) {
            if (child_first_negative == child.end()) break;
            if (std::ranges::find(child, value) == child.end()) {
                *child_first_negative = value;
                child_first_negative += 1;
            }
        }
        return child;
    };

    // check that the crossover function works correctly
    const knapsack p1 = {1, -1, -1, -1, -1};
    const knapsack p2 = {0, 2, 3, -1, -1};
    const knapsack p3 = {0, 1, -1, -1, -1};
    const knapsack p4 = {0, 1, 2, 3, -1};
    auto c1 = crossover(p1, p2);
    auto c2 = crossover(p2, p1);
    auto c3 = crossover(p2, p4);
    auto c4 = crossover(p3, p4);

    ASSERT_EQ(c1, knapsack({1, 0, 2, 3, -1}));
    ASSERT_EQ(c2, knapsack({0, 2, 3, 1, -1}));
    ASSERT_EQ(c3, knapsack({0, 2, 3, 1, -1}));
    ASSERT_EQ(c4, knapsack({0, 1, 2, 3, -1}));

    // the solution is all the boxes except for the heaviest one.
    const knapsack solution = {-1, 1, 2, 3, 4};
    const knapsack all_items = {0, 1, 2, 3, 4};

    // assert that the actual solution has the highest fitness value.
    // this needs to be true for the algorithm to converge correctly.
    ASSERT_GT(fitness(solution), fitness(all_items));

    // genetic algorithm settings.
    constexpr dp::ga::algorithm_settings settings{0.1, 0.5, 0.25};

    // generate an initial random population
    constexpr auto population_size = 2;
    std::vector<knapsack> initial_population{};  // TODO: Generate initial population
    initial_population.reserve(population_size);

    // random length uniform distribution
    std::uniform_int_distribution<int> length_dist(1, 4);

    // random value uniform distribution
    std::uniform_int_distribution<int> values_dist(0, available_items.size() - 1);

    // generator lambda
    auto knapsack_generator = [&engine, &length_dist, &values_dist]() {
        knapsack basic;
        std::ranges::fill(basic, -1);

        const auto random_length = length_dist(engine);
        for (auto i = 0; i < random_length; ++i) {
            auto value = values_dist(engine);
            // only allow unique values
            while (std::ranges::find(basic, value) != std::end(basic)) {
                value = values_dist(engine);
            }
            basic[i] = value;
        }
        return basic;
    };

    // generate the population
    std::ranges::generate_n(std::back_inserter(initial_population), population_size,
                            knapsack_generator);

    // define the termination criteria
    auto termination =
        dp::ga::termination::fitness_termination_criteria<knapsack>(fitness(solution));

    static_assert(dp::detail::termination_operator<
                  dp::ga::termination::fitness_termination_criteria<knapsack>, knapsack, int>);

    dp::genetic_algorithm genetics(settings, initial_population, mutator, crossover, fitness,
                                   termination);
    auto [best, _] = genetics.solve([](auto& stats) {
        std::cout << "best: " << stats.current_best.best
                  << " fitness: " << stats.current_best.fitness << "\n";
    });

    std::ranges::sort(best);
    EXPECT_EQ(best, solution);
}

// Helper function for testing selectors
template <typename Selector>
std::unordered_map<std::string, int> test_selection(
    Selector& selector, const std::string& test_value,
    const std::vector<std::string>& initial_population, const unsigned& selection_count = 1000) {
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
    for (auto i = 0; i < selection_count; i++) {
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
