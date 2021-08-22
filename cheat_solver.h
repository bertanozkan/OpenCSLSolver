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

#ifndef CHEAT_SOLVER_H
#define CHEAT_SOLVER_H

#include "cfvm.h"
#include "enums.h"
#include "memory"
#include "solvers_common.h"

class CHeat_solver
{
public:
    CHeat_solver();

    class CHeat_cell
    {
    public:
        double Temperature = 0.0;
    };


    class CHeat_face
    {
    public:
        CHeat_cell* left_cell_ptr = nullptr;
        CHeat_cell* right_cell_ptr = nullptr;
    };


    class CHeat_BC_face
    {
    public:
        enums::Ethermal_BC_type ethermal_BC_type;
        double Temperature = 0.0;
        CHeat_cell* left_cell_ptr = nullptr;
    };


    std::vector<CHeat_cell> cells;
    std::vector<CHeat_face> faces;
    std::vector<CHeat_BC_face> BC_faces;

    void cell_loop(const CFVM::CFVM_cell &FVM_cell, const CHeat_cell &Heat_cell, std::shared_ptr<CMatrix_solver> Matrix_solver);
    void calc_flux(const CFVM::CFVM_face &FVM_face, const CHeat_face &Heat_face, std::shared_ptr<CMatrix_solver> Matrix_solver);
    void calc_flux(const CFVM::CFVM_BC_face &FVM_BC_face, CHeat_BC_face &Heat_BC_face, std::shared_ptr<CMatrix_solver> Matrix_solver);

    double temperature_global_residual = 0.0;

    void initialize(const std::shared_ptr<CFVM> FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver);

    void iterate(const std::shared_ptr<CFVM> FVM, std::shared_ptr<CMatrix_solver> Matrix_solver);
};




#endif // CHEAT_SOLVER_H
