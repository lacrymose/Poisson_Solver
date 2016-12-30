#ifndef _H_RIGHT_HAND_SIDE__
#define _H_RIGHT_HAND_SIDE__

#include <deal.II/base/function.h>
#include <math.h>
using namespace dealii;

// RHS class publicly derived from Function class (inheritence)
template <int dim>
class RightHandSide : public Function<dim>
{
public:
    RightHandSide (): Function<dim>() // base class initializer
    {}

    // NOTE: Only base class' function needs to be declared virtual
    // gives value at point, component for vector value functions
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;
};

#endif
