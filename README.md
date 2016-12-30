This code is designed to numerically solve the
 <a href="https://en.wikipedia.org/wiki/Poisson's_equation">Poisson equation</a> 
using the <a href="https://en.wikipedia.org/wiki/Finite_element_method">
finite Element Method (FEM)</a>.


**Requirements**
The requirements for this software is 
<a href="https://www.dealii.org">deal.ii library</a> version 8.3.0 or higher,
<a href="https://www.cmake.org">CMake</a> version 2.8 or higher.

**Installation**
First obtain and install a copy of the dealii
<a href="https://www.dealii.org">deal.ii library</a> version 8.3.0 or higher. 
See the dealii library for installation instructions and help installing trilinos and p4est.

**Compiling**
To generate a makefile for this code using CMake type into the terminal:

cmake . -DDEAL_II_DIR=/path_to_deal.ii

To can compile the code use:

*make release*

**Running**
To run the executable use:

*./main*