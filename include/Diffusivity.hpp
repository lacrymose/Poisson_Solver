#ifndef _H_DIFFUSIVITY__
#define _H_DIFFUSIVITY__

#include<deal.II/base/function.h>

using namespace dealii;

// Diffusivity function class publically derived from Function
// class, use the defualt constructor!!!

template <int dim>
class DiffusivityFunction : public Function<dim>
{
public:
    DiffusivityFunction() : Function<dim>() // base class initializer
    {}

    // gives value of function at a point
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;

    // returns list of values (values) of the function evaluated at
    // the list of points (points)
    virtual void value_list(const std::vector< Point<dim> > &points,
                            std::vector<double>						  &values,
                            const unsigned int						  component=0) const;
};

#endif
