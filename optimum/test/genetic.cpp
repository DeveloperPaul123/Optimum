#include "optimum/evolutionary/genetic.h"
#include <gtest/gtest.h>

#include <random>

struct random_word_generator {
	random_word_generator() : generator_(device_()) {

	}
	[[nodiscard]] std::string operator()(std::string_view char_set, std::size_t max_length) {
		std::string copy = char_set.data();
		std::shuffle(copy.begin(), copy.end(), generator_);
		return copy.substr(0, max_length);
	}
private:
	std::random_device device_;
	std::mt19937 generator_;
};

TEST(GeneticAlgorithm, BasicTest) {

	const std::string solution = "hello world";
	const std::string alphabet = "abcdefghijklmnopqrstuvwxyz";
	const std::string available_chars = alphabet + " ";

	const auto word_length = solution.size();
	random_word_generator word_generator{};
	// generate initial population
	std::vector<std::string> initial_population(100);
	std::generate_n(std::back_inserter(initial_population), 100, [&]() {
		return word_generator(available_chars, word_length);
		});

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

	auto cross_over = dp::default_crossover < std::string>();
	dp::genetic_algorithm<std::string> ga(
		initial_population,
		string_mutator,
		cross_over,
		fitness_op
	);

	auto results = ga.solve();
	auto best = results.best;
	// TODO
}
