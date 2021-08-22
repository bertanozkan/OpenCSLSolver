/*
 * OpenCSLSolver: An open-source, parallel, pressure based CFD solver.
 * Copyright (C) 2021  Bertan ÖZKAN
 *
 * Email: ozkanbertan@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SOLVERS_COMMON_H
#define SOLVERS_COMMON_H

#include<croot.h>
#include "cmatrix_solver.h"

namespace solvers_common {


constexpr double alpha = 0.7;

template<class D>
inline D face_average(const CFVM::CFVM_face &FVM_face, const D left_val, const D right_val)
{
    D val;
    val = FVM_face.geo_aver * left_val;
    val += (1.0 - FVM_face.geo_aver) * right_val;
    return val;
}


inline void Diffusive_flux( std::shared_ptr<CMatrix_solver> &Matrix_solver,const CFVM::CFVM_face &FVM_face,const double D)
{
    // Fill matrix for left cell
    Matrix_solver->matrix.add_insert_element(FVM_face.left_cell_index,FVM_face.left_cell_index,D);
    Matrix_solver->matrix.add_insert_element(FVM_face.left_cell_index,FVM_face.right_cell_index,-D);


    // Fill matrix for right cell
    Matrix_solver->matrix.add_insert_element(FVM_face.right_cell_index,FVM_face.right_cell_index,D);
    Matrix_solver->matrix.add_insert_element(FVM_face.right_cell_index,FVM_face.left_cell_index,-D);

}


inline void Diffusive_flux_dirichlet_BC( std::shared_ptr<CMatrix_solver> &Matrix_solver, const CFVM::CFVM_BC_face &FVM_BC_face,
                           const double D, const double BC_scalar_value)
{

    // Fill matrix for left cell
    Matrix_solver->matrix.add_insert_element(FVM_BC_face.left_cell_index,FVM_BC_face.left_cell_index,D);


    // Calc and fill BC value to appropriate places
    double flux = D * BC_scalar_value;
    Matrix_solver->b_vector.add_insert_element(FVM_BC_face.left_cell_index,flux);
}

template<enums::Econv_flux_type econv_flux_type >
inline double Convective_face_interp(const CFVM::CFVM_face &FVM_face, const double Mdot, const double left_scalar, const double right_scalar)
{

    constexpr double upwind_weigth = 1.0 - alpha;
    constexpr double central_weight = alpha;
    double value = 0.0;
    if(Mdot >= 0.0)
    {
        if constexpr (econv_flux_type == enums::Econv_flux_type::Blended_central)
            value = left_scalar*upwind_weigth + central_weight * face_average(FVM_face, left_scalar, right_scalar);
        else
            value = left_scalar;
    }else
    {
        if constexpr (econv_flux_type == enums::Econv_flux_type::Blended_central)
            value = right_scalar*upwind_weigth + central_weight * face_average(FVM_face, left_scalar, right_scalar);
        else
        value = right_scalar;
    }
    return value;
}


template<enums::Econv_flux_type econv_flux_type >
inline void Convective_flux( std::shared_ptr<CMatrix_solver> &Matrix_solver,
                           const CFVM::CFVM_face &FVM_face, const double Mdot, const double left_val, const double right_val)
{

    double val;
    if(Mdot >= 0.0)
    {
        // Fill matrix for left cell
        Matrix_solver->matrix.add_insert_element(FVM_face.left_cell_index,FVM_face.left_cell_index,Mdot);

        // Fill matrix for right cell
        Matrix_solver->matrix.add_insert_element(FVM_face.right_cell_index,FVM_face.left_cell_index,-Mdot);

        if constexpr (econv_flux_type == enums::Econv_flux_type::Blended_central)
        {
            val = alpha * (face_average(FVM_face, left_val, right_val) - left_val);
            // Fill b vector for left cell
            Matrix_solver->b_vector.add_insert_element(FVM_face.left_cell_index, -Mdot * val);

            // Fill b vector for right cell
            Matrix_solver->b_vector.add_insert_element(FVM_face.right_cell_index,Mdot* val);
        }

    }else if(Mdot < 0.0)
    {
        // Fill matrix for left cell
        Matrix_solver->matrix.add_insert_element(FVM_face.left_cell_index,FVM_face.right_cell_index,Mdot);

        // Fill matrix for right cell
        Matrix_solver->matrix.add_insert_element(FVM_face.right_cell_index,FVM_face.right_cell_index,-Mdot);

        if constexpr (econv_flux_type == enums::Econv_flux_type::Blended_central)
        {
            val = alpha * (face_average(FVM_face, left_val, right_val) - right_val);
            // Fill b vector for left cell
            Matrix_solver->b_vector.add_insert_element(FVM_face.left_cell_index,-Mdot*val);

            // Fill b vector for right cell
            Matrix_solver->b_vector.add_insert_element(FVM_face.right_cell_index, Mdot * val);
        }
    }
}


inline void Convective_flux_neumann_zero_BC( std::shared_ptr<CMatrix_solver> &Matrix_solver,
                           const CFVM::CFVM_BC_face &FVM_BC_face, const double Mdot)
{
    if(Mdot >= 0.0)
    {
        // Fill matrix for left cell
        Matrix_solver->matrix.add_insert_element(FVM_BC_face.left_cell_index,FVM_BC_face.left_cell_index,Mdot);

    }else if(Mdot < 0.0)
    {
        // Fill matrix for left cell
        Matrix_solver->matrix.add_insert_element(FVM_BC_face.left_cell_index,FVM_BC_face.left_cell_index,Mdot);
    }
}



}

#endif // SOLVERS_COMMON_H
