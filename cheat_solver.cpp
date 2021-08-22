/*
OpenCSLSolver: An open-source, parallel, pressure based CFD solver.
Copyright (C) 2021  Bertan ÖZKAN

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
 any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "cheat_solver.h"
#ifdef COMPILE_WITH_OPENMP
#include <omp.h>
#endif

CHeat_solver::CHeat_solver()
{

}


void CHeat_solver::calc_flux(const CFVM::CFVM_face &FVM_face, const CHeat_face &Heat_face, std::shared_ptr<CMatrix_solver> Matrix_solver)
{
    double D = root.configuration->thermal_diffusivity * FVM_face.diffusion_ksi_A;

    solvers_common::Diffusive_flux(Matrix_solver,FVM_face, D);
}

void CHeat_solver::calc_flux(const CFVM::CFVM_BC_face &FVM_BC_face, CHeat_BC_face &Heat_BC_face, std::shared_ptr<CMatrix_solver> Matrix_solver)
{

    // If boundary is adiabatic
    if(Heat_BC_face.ethermal_BC_type == enums::Ethermal_BC_type::adiabatic){
        Heat_BC_face.Temperature = Heat_BC_face.left_cell_ptr->Temperature;
    }

    // IF boundary is constant_temperature
    if(Heat_BC_face.ethermal_BC_type == enums::Ethermal_BC_type::constant_temperature){

        double D = root.configuration->thermal_diffusivity * FVM_BC_face.diffusion_ksi_A;

        solvers_common::Diffusive_flux_dirichlet_BC( Matrix_solver,FVM_BC_face,D,Heat_BC_face.Temperature);

    }
}

void CHeat_solver::cell_loop(const CFVM::CFVM_cell &FVM_cell,const CHeat_cell &Heat_cell, std::shared_ptr<CMatrix_solver> Matrix_solver)
{

    double dum = (FVM_cell.volume / root.configuration->global_delta_t);
    Matrix_solver->matrix.add_insert_element(FVM_cell.index,FVM_cell.index,dum);
    Matrix_solver->b_vector.add_insert_element(FVM_cell.index,Heat_cell.Temperature*dum);
    Matrix_solver->x_vector.insert_element(FVM_cell.index , Heat_cell.Temperature);

}

void CHeat_solver::initialize(const std::shared_ptr<CFVM> FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver)
{

    cells.resize(FVM->cells.size());
    faces.resize(FVM->faces.size());
    BC_faces.resize(FVM->BC_faces.size());
    size_t nnz = 0;

    // Inserting initial conditions to cells
    for(auto& cell : cells)
    {
        cell.Temperature = root.configuration->inital_condition->Temperature;
    }

    // Inserting FVM face pointers to heat_solver faces
    for(size_t i = 0 ; i < faces.size() ; i++)
    {
        faces[i].left_cell_ptr = &cells[FVM->faces[i].left_cell_index];
        faces[i].right_cell_ptr = &cells[FVM->faces[i].right_cell_index];
    }

    // Inserting FVM BC_face pointers to heat_solver BC_faces
    for(size_t i = 0 ; i < BC_faces.size() ; i++)
    {
        BC_faces[i].left_cell_ptr = &cells[FVM->BC_faces[i].left_cell_index];
    }

    // Inserting boundary conditions to faces
    auto BC_face = BC_faces.begin();
    auto FVM_BC_face = FVM->BC_faces.begin();
    for(; BC_face != BC_faces.end(); BC_face++ ,FVM_BC_face++)
    {
        int i_BC = root.geometry->faces[FVM_BC_face->global_face_index]
                .boundary_marker;
        BC_face->Temperature = root.configuration->boundary_conditions[i_BC]
                .Temperature;
        BC_face->ethermal_BC_type = root.configuration->boundary_conditions[i_BC]
                .ethermal_BC_type;
    }
}



void CHeat_solver::iterate(const std::shared_ptr<CFVM> FVM, std::shared_ptr<CMatrix_solver> Matrix_solver)
{

    Matrix_solver->matrix.clear();
    Matrix_solver->b_vector.clear();
    Matrix_solver->x_vector.clear();

    #pragma omp parallel
    {
        size_t i_face;
        #pragma omp for
        for(i_face = 0; i_face < faces.size() ; i_face++)
            calc_flux(FVM->faces[i_face], faces[i_face], Matrix_solver);
        #pragma omp barrier

        size_t i_BC_face;
        #pragma omp for
        for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
            calc_flux(FVM->BC_faces[i_BC_face], BC_faces[i_BC_face], Matrix_solver);
        #pragma omp barrier

        size_t i_cell;
        #pragma omp for
        for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
            cell_loop(FVM->cells[i_cell],cells[i_cell], Matrix_solver);
    }


    temperature_global_residual = 0.0;

    // Solve matrix
    Matrix_solver->solve();

    #pragma omp parallel for reduction (+:temperature_global_residual)
    for(size_t i_cell = 0 ; i_cell < cells.size() ; i_cell++)
    {
        temperature_global_residual += abs(cells[i_cell].Temperature - Matrix_solver->x_vector.get_element(i_cell));
        cells[i_cell].Temperature = Matrix_solver->x_vector.get_element(i_cell);
    }

}


