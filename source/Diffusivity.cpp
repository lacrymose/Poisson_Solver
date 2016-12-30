#include "../include/Diffusivity.hpp"

// Actual diffusivity function
template <int dim>
double DiffusivityFunction<dim>::value(const Point<dim> &p,
                                       const unsigned int /*component*/) const
{
    // f(x) = 20 		if (x^2 + y^2 + z^2 ) < 0.25
    // 			= 1  		otherwise
    if (p.square() < 0.5*0.5)
        return 1.0;
    else
        return 1.0;
}



template<int dim>
void DiffusivityFunction<dim>::value_list(const std::vector<Point<dim>> &points,
        std::vector<double>						&values,
        const unsigned int component) const
{
    // Assert that the two list have same size
    Assert(values.size() == points.size(),
           ExcDimensionMismatch (values.size(), points.size() ) );

    // Assert just dealing with scalar
    Assert(component == 0,
           ExcIndexRange(component, 0, 1) );

    // evaluate function at list of points and store in values
    const unsigned int n_points = points.size();
    for(unsigned int i=0; i<n_points; ++i)
    {
        if(points[i].square() < 0.5*0.5)
            values[i] = 1.0;
        else
            values[i] = 1.0;
    }
}


