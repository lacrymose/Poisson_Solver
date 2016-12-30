/* ---------------------------------------------------------------------
 *
 * $Id: step-1.cc 31349 2013-10-20 19:07:06Z maier $
 *
 * Copyright (C) 1999 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 */

// @sect3{Include files}
///////////////////////////////////////////////////////////////////////////////
// Solves:
// 		 - \nabla \cdot ( k(x) \nabla u(x) ) = f(x)		\in \Omega
// 										u(x) = u_D(x) 	\on \Gamma_D
// 						   -k(x) \nabla u(x) = g(x)		\on \Gamma_N
//
/////////////////////////////////////////////////////////////////////////////
//
//
#include "../include/Poisson.hpp"
#include "BoundaryFunctions.cpp"
#include "Diffusivity.cpp"
#include "RightHandSide.cpp"

namespace PoissonProblem
{

// Empty Constructor
template <int dim>
Poisson<dim>::Poisson(unsigned int degree,
                      unsigned int NumberRefinements)
    :  // member initializers
       fe(degree),
       dof_handler(triangulation),
       refinement_number(NumberRefinements)
{}

// destructor
template<int dim>
Poisson<dim>::~Poisson<dim>()
{
    // FEValues will be destroyed before dof_handler, but dof_handler
    // references FE values though dofs, so clear dofs before destroying
    // FEValues
    dof_handler.clear();
}


// make the grid
template <int dim>
void Poisson<dim>::make_grid ()
{
    // make the triangulation
    GridGenerator::hyper_cube (triangulation,-1,1);
    triangulation.refine_global (refinement_number);


    // print out the triangulation
    std::ofstream out ("grid-1.eps");
    GridOut grid_out;
    grid_out.write_eps (triangulation, out);
}

// assign which portions of the boundary are Neumann boundaries
template<int dim>
void Poisson<dim>::make_boundaries()
{
    typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();
    for(; cell != endc; cell++)
    {
        for(unsigned int face_no=0;
                face_no < GeometryInfo<dim>::faces_per_cell;
                ++face_no)
        {   // test to see if along the bottom boundary
            if( (std::fabs(cell->face(face_no)->center()(0) - (-1) ) < 1e-12)
                    ||	// test to see if along the left boundary
                    (std::fabs(cell->face(face_no)->center()(1) - (-1) ) < 1e-12) )
            {
                // indicate portion of boundary to be neumann (1)
                // default is Dirichlet (0)
                cell->face(face_no)->set_boundary_id(1);
            }
        } // for face
    }	// for cell
}

// distribute dofs and allocate memory
template <int dim>
void Poisson<dim>::setup_system()
{
    // distribute dofs
    dof_handler.distribute_dofs(fe);
    cout << "# of dofs: " << dof_handler.n_dofs() << endl;

    // make sparsity pattern
    DynamicSparsityPattern dsp(dof_handler.n_dofs() );
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    // allocate memory
    system_matrix.reinit(sparsity_pattern);
    solution.reinit(dof_handler.n_dofs() );
    system_rhs.reinit(dof_handler.n_dofs() );
}


// assmble matrix for Ax = b
template <int dim>
void Poisson<dim>::assemble_system()
{
    // quadrature rules depend on the degree of the approximation polynomial
    QGauss<dim> quadrature_formula(fe.degree+1); // body quadrature rule
    QGauss<dim-1> face_quadrature_formula(fe.degree+1); //face quad rules

    // RHS & Neumann  function instantiation
    const RightHandSide<dim> right_hand_side;
    const NeumannBoundaryValues<dim> neumann_condition;

    // create FEValues function
    FEValues<dim> fe_values(fe, quadrature_formula,
                            update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);

    // create FEFaceValues function
    FEFaceValues<dim>	fe_face_values(fe, face_quadrature_formula,
                                       update_values		 | 	update_quadrature_points |
                                       update_normal_vectors |	update_JxW_values);

    const unsigned int		dofs_per_cell 		= fe.dofs_per_cell;
    const unsigned int		n_q_points				= quadrature_formula.size();
    const unsigned int		n_face_q_points		= face_quadrature_formula.size();

    // local cell matrix and local rhs
    FullMatrix<double>		cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>				cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // instantiate diffusivity coefficient object and list to contain values
    // at quadrature points
    const DiffusivityFunction<dim> Diffusivity;
    std::vector<double> diffusivity_values(n_q_points);


    // create iterator over cells
    typename DoFHandler<dim>::active_cell_iterator
    cell	=	dof_handler.begin_active(),
       endc 	= dof_handler.end();

    // iterate over all the cells and construct local matrix and RHS
    for(; cell!=endc; ++cell)
    {
        // reinitialize fe_values for each cell
        fe_values.reinit(cell);


        // get the values of the diffusivity function at quadrature points
        // in this cell
        // NOTE: Get values of diffusivity ONE TIME PER CELL
        Diffusivity.value_list(fe_values.get_quadrature_points(),
                               diffusivity_values);

        cell_matrix = 0;
        cell_rhs = 0;


        // do body integral over cell
        for(unsigned int q_point=0; q_point<n_q_points; ++q_point)
        {
            // over all test functions i on cell
            for(unsigned int i=0; i<dofs_per_cell; ++i)
            {
                // over all test functions j on cell
                for(unsigned int j=0; j<dofs_per_cell; ++j)
                {
                    // A(i,j) = \int_{cell} \nabla \psi_{i} \cdot \nabla \psi_{j}
                    cell_matrix(i,j) += diffusivity_values[q_point] *
                                        fe_values.shape_grad(i, q_point) *
                                        fe_values.shape_grad(j, q_point) *
                                        fe_values.JxW(q_point);
                } // test function j

                // b(i) = \int_{cell} \psi_{i} f(x)
                cell_rhs(i) += ( fe_values.shape_value(i, q_point) *
                                 right_hand_side.value(fe_values.quadrature_point(q_point) ) *
                                 fe_values.JxW(q_point) );
            } // test function i
        }	// q_points


        // do the face integrals over all the faces of this cell that are on the
        // Neumann boundary
        for(unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
            if(cell->face(face_no)->at_boundary() // see if face on boundary
                    &&		// see if on Neumann part of boundary
                    (cell->face(face_no)->boundary_id() == 1) )
            {
                // update fe_face_values on this cell face
                fe_face_values.reinit(cell, face_no);

                for(unsigned int q_point=0; q_point<n_face_q_points; q_point++)
                {
                    // get the neumann BC value at the quadrature point
                    const double neumann_value =
                        (neumann_condition.vector_value(
                             fe_face_values.quadrature_point(q_point) )
                         * fe_face_values.normal_vector(q_point) );

                    // preform quadrature at point over all basis functions with
                    // support includeing that point
                    for(unsigned int i=0; i<dofs_per_cell; i++)
                    {
                        cell_rhs(i) += (neumann_value *
                                        fe_face_values.shape_value(i,q_point) *
                                        fe_face_values.JxW(q_point) );
                    } //for i (dofs)
                }	// for q_point
            }	// if on neumann boundary
        }	// for faces on cell


        // get global dof indices from local dof indices
        cell->get_dof_indices(local_dof_indices);

        // places local matrix and rhs into global
        for(unsigned int i=0; i<dofs_per_cell; ++i)
        {
            for(unsigned int j=0; j<dofs_per_cell; ++j)
            {
                system_matrix.add(local_dof_indices[i],
                                  local_dof_indices[j],
                                  cell_matrix(i,j) );
            }
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }

    } // over all cells

    // interpolate BC and apply them to linear system
    std::map<types::global_dof_index, double> boundary_values;

    VectorTools::interpolate_boundary_values(dof_handler,
            0,
            DirichletBoundaryValues<dim>(),
            boundary_values);

    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
}

template <int dim>
void Poisson<dim>::solve()
{
    // creat the solver object with set max iters and tol
    SolverControl					solver_control(1000,1e-12);
    SolverCG<>						solver(solver_control);
//	SolverGMRES<>						solver(solver_control);
    // solve the sytem
    solver.solve(system_matrix, solution, system_rhs,
                 PreconditionIdentity() );

    cout << "    " << solver_control.last_step()
         << " CG iters needed to obtain conv."
         << endl;

//	SparseDirectUMFPACK direct_solver;
//	direct_solver.initialize(system_matrix);
//	direct_solver.vmult(solution, system_rhs);

}

template<int dim>
void Poisson<dim>::compute_error() const
{
    DirichletBoundaryValues<dim> exact_solution;
    Vector<double> cellwise_errors(triangulation.n_active_cells() );

    QGauss<dim> quadrature_formula(fe.degree+1); // body quadrature rule
    VectorTools::integrate_difference(dof_handler, solution, exact_solution,
                                      cellwise_errors, quadrature_formula,
                                      VectorTools::L2_norm );

    const double l2_error = cellwise_errors.l2_norm();
    cout << "||e||_L2 = " << l2_error << endl;
}

template<int dim>
void Poisson<dim>::output_results() const
{
    // create DataOut object and in attatch solution to it

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    std::ofstream output (dim == 2 ?
                          "solution-2d.vtk" :
                          "solution-3d.vtk");

    data_out.add_data_vector(solution, "solution");

    data_out.build_patches();

    data_out.write_vtk(output);
    output.close();

}

template <int dim>
void Poisson<dim>::output_system() const
{
    std::ofstream output_system("A.mtx");
    system_matrix.print_formatted(output_system);
    output_system.close();
    output_system.open("b.vec");
    system_rhs.print(output_system);
    output_system.close();
}

template <int dim>
void Poisson<dim>::run()
{
    cout << "Solving problem in " << dim << " space dimensions." << endl;

    make_grid();
    make_boundaries();
    setup_system();
    assemble_system();
    solve();
    compute_error();
//	output_system();
    output_results();
}


} // end namespace
