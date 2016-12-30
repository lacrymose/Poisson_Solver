#ifndef H_BOUNDARY_FUNCTIONS__
#define H_BOUNDARY_FUNCTIONS__

#include <deal.II/base/function.h>

using namespace dealii;

// Dirichlet BC Class publicly derived from Function class (inheritence)
template <int dim>
class DirichletBoundaryValues : public Function<dim>
{
public:
    DirichletBoundaryValues() : Function<dim>() // base class initializer
    {}

    // NOTE: Only base class' function needs to be declared virtual
    // gives value at point, component for vector value functions
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;
};


// Neumann BC Class publicly derived from Function class (inheritence)
template <int dim>
class NeumannBoundaryValues : public Function<dim>
{
public:
    NeumannBoundaryValues() : Function<dim>() // base class initializer
    {}

    // gives the vector (1-tensor) of values of the neumann condition
    // at the point p
    virtual Tensor<1,dim> vector_value(const Point<dim> &p,
                                       const unsigned int component = 0) const;
};

#endif
