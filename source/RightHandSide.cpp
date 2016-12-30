#include "../include/RightHandSide.hpp"

// RHS Value function
template <int dim>
double RightHandSide<dim>::value(const Point<dim> &p,
                                 const unsigned int /*component*/) const
{
    /*
    	double return_value = 2*M_PI*M_PI;
    	for( unsigned int i=0; i<dim; i++)
    		return_value *= sin(M_PI*p[i]);
    */
    double return_value = -4.0;
    return return_value;
}


