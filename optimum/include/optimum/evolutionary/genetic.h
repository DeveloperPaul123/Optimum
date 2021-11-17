#pragma once

#include <concepts>
#include <vector>
#include <numeric>
#include <functional>

namespace dp {
	namespace detail {
		template<class Fn, class T>
		concept mutation_operator = std::invocable<Fn, T> && std::is_same_v<std::invoke_result_t<Fn, T>, T>;

		template<class Fn, class T, class Result = std::invoke_result_t<Fn, T>>
		concept fitness_operator = std::invocable<Fn, T> && std::is_arithmetic_v<Result>;

		template<class Fn, class T, class Result = std::invoke_result_t<Fn, T, T>>
		concept crossover_operator = std::invocable<Fn, T> && std::is_same_v<Result, T>;
	}

	template<typename T>
	struct default_mutator {
		T operator()(T&& t) {
			return t;
		}
	};

	template<typename T>
	struct default_crossover {
		T operator()(T&& first, T&& second) {
			return first + second;
		}
	};

	template<typename T, typename U = double>
		requires std::floating_point<U> || std::integral<U>
	struct default_evaluation {
		constexpr U operator()(T&& t) {
			return std::accumulate(t.begin(), t.end(), {});
		}
	};

	template<class ChromosomeType>
	class genetic_algorithm {
	public:

		using value_type = ChromosomeType;
		using list_type = std::vector<value_type>;

		struct results {
			genetic_algorithm::value_type best;

		};

		template<
			class FitnessOperator,
			class MutationOperator,
			class CrossoverOperator
		>
		requires detail::mutation_operator<MutationOperator, ChromosomeType> &&
			detail::fitness_operator<FitnessOperator, ChromosomeType> &&
			detail::crossover_operator<CrossoverOperator, ChromosomeType>

		genetic_algorithm(list_type initial_population,
				MutationOperator&& mutator,
				CrossoverOperator&& crossover_operator,
				FitnessOperator&& fitness_operator)
			: initial_population_(initial_population),
			mutator_(std::forward<MutationOperator>(mutator)), crossover_(std::forward<CrossoverOperator>(crossover_operator)),
			fitness_(std::forward<FitnessOperator>(fitness_operator)) {
		}

		[[nodiscard]] results solve() {
			// TODO
			return {};
		}

	private:
		using mutator = std::function<ChromosomeType(ChromosomeType)>;
		using crossover = std::function<ChromosomeType(ChromosomeType, ChromosomeType)>;
		using fitness = std::function<double(ChromosomeType)>;

		list_type initial_population_{};
		mutator mutator_;
		crossover crossover_;
		fitness fitness_;
	};
}
