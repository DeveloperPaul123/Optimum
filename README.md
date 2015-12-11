# Optimum
A small, lightweight library with various optimization algorithms including the Nelder-Mead algorithm. 

## What's in it? ##
* Nelder-Mead Optimization Algorithm 

#### Usage ####

````cpp
float precision = 0.0001;
int dimension = 2;
NelderMeadOptimizer o(dimension, precision);

// request a simplex to start with
Vector v(0.5, 0.5);
o.initialGuess(v);

while (!o.done()) {
    v = o.step(v, f(v));
}
````

#### In Progress ####

* Genetic algorithm.
