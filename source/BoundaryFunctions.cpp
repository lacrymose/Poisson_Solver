
#include "../include/BoundaryFunctions.hpp"

// Dirichlet BC Value function
template <int dim>
double DirichletBoundaryValues<dim>::value(const Point<dim> &p,
        const unsigned int /*component*/) const
{
    return p.square();

}

// Neumann BC Value function
template <int dim>
Tensor<1,dim> NeumannBoundaryValues<dim>::vector_value(const Point<dim> &p,
        const unsigned int ) const
{
    Tensor<1,dim> return_value = 2.0 * p;
    return return_value;

}
