<p align="center">
  
<img src="https://socialify.git.ci/DeveloperPaul123/Optimum/image?description=1&descriptionEditable=A%20small%20optimization%20library%20using%20modern%20C%2B%2B&font=Raleway&forks=1&language=1&pattern=Signal&stargazers=1&theme=Dark" alt="optimum" width="640" height="320" />
  
  <br>
  <br>

 <a href="https://www.apache.org/licenses/LICENSE-2.0.html">
    <img src="https://img.shields.io/badge/license-Apache 2.0-blue" alt="License Apache 2.0">
  </a>
  
  <a href="https://github.com/DeveloperPaul123/Optimum/stargazers">
    <img src="https://img.shields.io/badge/Say%20Thanks-👍-1EAEDB.svg" alt="Say thanks">
  </a>
  
  <a href="https://discord.gg/MyGjmfQTFP">
    <img alt="Discord" src="https://img.shields.io/discord/652515194572111872">
  </a>

</p>

`Optimum` is an accessible, optimization algorithms library written in C++.

:warning: *The api is currently under heavy development and is subject to change.* :warning:

## What's in it?

* 🚧 **Nelder-Mead Optimization Algorithm:** This algorithm is an optimization algorithm that uses a simplex to converge on a 'best' point. It takes n+1 points for an n-dimensional problem and creates a simplex that is then expanded, contracted, reflected and so on to move towards the minima.

* ✅ **Genetic Algorithm**: This algorithm is optimization algorithm that introduces random mutations and cross overs (like nature and genetics) to make it more likely to find the true global minima or maxima of a given cost function. This is good to use when you know very little about the cost function or if there are a lot of local minima or maxima.

* 🚧 **Kabsch Algorithm:** This algorithm provides a way to find the optimal translation and rotation for the mapping of two point clouds to each other. In order to use this you need to know the mapping of your point pairs beforehand. 

* 🚧 **Iterative Closest Point** (:construction: Work in progress :construction:): This algorithm provides a means to perform registration between two clouds of points in either 2D or 3D. This can be used when little is know about the mapping between the two point clouds, or if there is no matching points. 

## Why?

I made this library because these types of algorithms are typically only the topic of advanced research papers. I wanted to make a simple, small library that makes some of these algorithms not only accessible to all, but also easy to use. 

## Usage

All classes are in the `dp` namespace. I make judicious use of namespaces but they may be omitted for brevity in any of the following examples.

### Nelder-Mead Optimizer

:construction:

### Genetic Algorithm

At its core, a genetic algorithm proceeds with the following steps:

* Selection
* Crossover
* Mutation
* Fitness evaluation

These steps repeat until a suitable solution, based on given termination criteria, is found. `dp::genetic_algorithm` is designed to be flexible and support any type of problem. As such, you will need to supply operators that perform the basic steps of the algorithm. For this usage example, we will use the algorithm to find the string `hello world`. We will start with selection criterea. The default selection model is [rank selection](https://en.wikipedia.org/wiki/Selection_(genetic_algorithm)). [Roulette selection](https://en.wikipedia.org/wiki/Fitness_proportionate_selection) is also available.

Next we need a way to evaluate the fitness of each string.

```cpp
auto fitness_op = [solution, available_chars](std::string value) -> double {
    double score = 0.0;
    for (auto i = 0; i < solution.size(); i++) {
        score += solution[i] == value[i];
    }
    return score;
};
```

This fitness function compares each character of the input string to the corresponding character of the solution. The boolean result for each comparison is then summed to give a final fitness score. With this, we see that the max score possible is the length of the solution string.

Next we need to create a crossover operator. This will be the operator that generates new "children".

```cpp
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
```

This operator simply takes have of one string (i.e. parent) and half of the other and concatenates them together to form a new string. Note that we create a new string that is the same size as the two parent strings.

Finally, we need a mutation operator. This operator will randomly mutate a string to introduce genetic variance in our population.

```cpp
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
```

This operator simple replaces a character in the string with a random one from the list of used characters.

With these operators defined, we can now set up the `dp::genetic_algorithm` and solve for the solution.

```cpp
auto termination =
    dp::ga::termination::fitness_termination_criteria<std::string>(fitness_op(solution));
dp::ga::algorithm_settings settings{0.40, 0.5, 0.2};
dp::genetic_algorithm<std::string> genetics(initial_population, string_mutator, cross_over,
                                        fitness_op, termination, settings);
auto results = genetics.solve();
// solution is inside results object
```

Note that we define our termination criterea here so that when the target fitness is reached, then the algorithm stops. This works well here because there is only one max to our fitness evaluation and there is only one true solution.

### Kabsch Algorithm

:construction:

## Contributing

Contributions as welcome! See the [contribution guidelines](CONTRIBUTING.md) for information on how to do so.

## License

```txt
Copyright 2015-2021 Paul T

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```

See [LICENSE](LICENSE) for more details.
