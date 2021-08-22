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

#ifndef CMATRIX_SOLVER_SUPPLY_H
#define CMATRIX_SOLVER_SUPPLY_H

#include "croot.h"

namespace unit_tests {


static inline void test_matrix_solver_cg()
{

    std::cout << "Matrix solver GC test " << std::endl;
    CMatrix_solver matrix_solver(10,10,30);

    matrix_solver.matrix.insert_element(0,0,5.0);
    matrix_solver.matrix.insert_element(0,1,1.0);
    matrix_solver.b_vector.insert_element(0,2.0);
    matrix_solver.x_vector.insert_element(0,1.0);
    for(size_t i = 1 ; i < 9 ; i++)
    {
        matrix_solver.matrix.insert_element(i,i,5.0);
        matrix_solver.matrix.insert_element(i,i-1,1.0);

        matrix_solver.matrix.insert_element(i,i+1,1.0);
        matrix_solver.b_vector.insert_element(i,2.0);
        matrix_solver.x_vector.insert_element(i,1.0);
    }
    matrix_solver.matrix.insert_element(matrix_solver.matrix.n_rows-1,matrix_solver.matrix.n_rows-2,1.0);
    matrix_solver.matrix.insert_element(matrix_solver.matrix.n_rows-1,matrix_solver.matrix.n_rows-1,5.0);
    matrix_solver.b_vector.insert_element(matrix_solver.matrix.n_rows-1,2.0);
    matrix_solver.x_vector.insert_element(matrix_solver.matrix.n_rows-1,1.0);

    /*std::cout << " Matrix row_ptr :" << std::endl;
    for(size_t i = 0 ; i < matrix_solver.matrix.row_ptr.size() ; i++)
        std::cout << " " << matrix_solver.matrix.row_ptr[i] << " ";
    std::cout << std::endl;

    std::cout << " Matrix value array :" << std::endl;
    for(size_t i = 0 ; i < matrix_solver.matrix.value.size() ; i++)
        std::cout << " " << matrix_solver.matrix.value[i] << " ";
    std::cout << std::endl;

    std::cout << " Matrix column index array :" << std::endl;
    for(size_t i = 0 ; i < matrix_solver.matrix.column_index.size() ; i++)
        std::cout << " " << matrix_solver.matrix.column_index[i] << " ";
    std::cout << std::endl;

    std::cout << " b vector :" << std::endl;
    for(size_t i = 0 ; i < matrix_solver.b_vector.n_values ; i++)
        std::cout << " " << matrix_solver.b_vector.get_element(i) << " ";
    std::cout << std::endl;

    std::cout << " x vector :" << std::endl;
    for(size_t i = 0 ; i < matrix_solver.x_vector.n_values ; i++)
        std::cout << " " << matrix_solver.x_vector.get_element(i) << " ";
    std::cout << std::endl;*/

    double res = matrix_solver.CG_solve();


    std::cout << "CG solver residual 0.323 to : " << res << std::endl;
}


static inline void test_matrix_solver_ILU_BICGSTAB()
{

    std::cout << "Matrix solver ILU BICGSTAB test " << std::endl;
    CMatrix_solver matrix_solver(10,10,30);

    matrix_solver.matrix.insert_element(0,0,5.0);
    matrix_solver.matrix.insert_element(0,1,1.0);
    matrix_solver.b_vector.insert_element(0,2.0);
    matrix_solver.x_vector.insert_element(0,1.0);
    for(size_t i = 1 ; i < 9 ; i++)
    {
        matrix_solver.matrix.insert_element(i,i,5.0);
        matrix_solver.matrix.insert_element(i,i-1,1.0);

        matrix_solver.matrix.insert_element(i,i+1,1.0);
        matrix_solver.b_vector.insert_element(i,2.0);
        matrix_solver.x_vector.insert_element(i,1.0);
    }
    matrix_solver.matrix.insert_element(matrix_solver.matrix.n_rows-1,matrix_solver.matrix.n_rows-2,1.0);
    matrix_solver.matrix.insert_element(matrix_solver.matrix.n_rows-1,matrix_solver.matrix.n_rows-1,5.0);
    matrix_solver.b_vector.insert_element(matrix_solver.matrix.n_rows-1,2.0);
    matrix_solver.x_vector.insert_element(matrix_solver.matrix.n_rows-1,1.0);


    double res = matrix_solver.ILU_BICGSTAB_solve();


    std::cout << "ILU_BICGSTAB solver residual 0.323 to : " << res << std::endl;
}

}

#endif // CMATRIX_SOLVER_SUPPLY_H
