# Optimum
A small, lightweight library with various optimization algorithms including the Nelder-Mead algorithm. 

## What's in it? ##
* Nelder-Mead Optimization Algorithm 
* Genetic Algorithm (Work in progress)

## Why? ##
I made this library because these types of algorithms (for optimization/minimization) are typically only the topic of advanced research papers. I wanted to make this library as a spin-off of a github repository I found [online](https://github.com/blinry/nelder-mead-optimizer). I wanted to improve it and add to it, yet keep it easy to use and incorporate into any c++ project. Future support may include other languages if this does well. 

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

#### License ####
-------
This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

See LICENSE file for a copy of the GNU General Public License.

Copyright (C) 2015  Paul T <developer.paul.123@gmail.com>
