#pragma once

namespace dp::ga::termination {
    template <typename T>
    struct generations_termination_criteria {
        generations_termination_criteria(std::size_t max_generations) : count_(max_generations) {}

        bool operator()(T, double) {
            count_--;
            return count_ == 0;
        }

      private:
        std::size_t count_;
    };

    template <typename T>
    struct fitness_termination_criteria {
        fitness_termination_criteria(double target_fitness) : fitness_(target_fitness) {}
        bool operator()(T, double fitness) { return fitness >= fitness_; }

      private:
        double fitness_{};
    };
}  // namespace dp::ga::termination
