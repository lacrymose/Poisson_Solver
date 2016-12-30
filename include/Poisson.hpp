#ifndef H_POISSON__
#define H_POISSON__

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/lac/solver_gmres.h>

#include <iostream>
#include <fstream>


using namespace dealii;
using std::cout;
using std::endl;

namespace PoissonProblem
{
template <int dim>
class Poisson
{
    // Member functions
public:
    Poisson(unsigned int degree,
            unsigned int NumberRefinements);
    ~Poisson();
    void run();

private:
    void make_grid();
    void make_boundaries();
    void setup_system();
    void assemble_system();
    void solve();
    void compute_error() const;
    void output_results() const;
    void output_system() const;

    // Object member data
    Triangulation<dim>			triangulation;
    FE_Q<dim>							 	fe;
    DoFHandler<dim>					dof_handler;

    SparsityPattern					sparsity_pattern;
    SparseMatrix<double>		system_matrix;

    Vector<double>					solution;
    Vector<double> 					system_rhs;

    // For adativity
    unsigned 	int 					refinement_number;
};
}

#endif
