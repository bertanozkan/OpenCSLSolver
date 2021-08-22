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

#ifndef CROOT_H
#define CROOT_H

#include "vector"
#include "math.h"
#include "string"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <utility>

class Croot;
class Cconfiguration;
class Cpost_processing;
class Cgeometry;
class CFVM;
class Cfluid_properties;
class CHeat_solver;
class CFlow_solver;
class CMatrix_solver;


/*!
 * \brief Root Class for all functions and classes
 */
class Croot
{
public:
    Croot();

    std::shared_ptr <Cconfiguration> configuration = nullptr;
    std::shared_ptr <Cpost_processing> post_processing = nullptr;
    std::shared_ptr <Cgeometry> geometry = nullptr;
    std::shared_ptr <CFVM> FVM = nullptr;
    std::shared_ptr <CHeat_solver> Heat_solver = nullptr;
    std::shared_ptr <CFlow_solver> Flow_solver = nullptr;
    std::shared_ptr <Cfluid_properties> Fluid_properties = nullptr;
    std::shared_ptr <CMatrix_solver> Matrix_solver = nullptr;

    int i_iteration;
    int i_inner_iter;

    // variables for calculating elapsed time
    unsigned long long int total_elapsed_iteration_time = 0.0;
    unsigned long long int total_elapsed_matrix_solution_time = 0.0;

    void run_config(Cconfiguration inp_config);
    void run_config(std::string filename);

    void create_geometry();

    void initialize_solver();

    void start_solver();


};

extern Croot root;



#include "enums.h"
#include "tensor.h"
#include "cconfiguration.h"
#include "cgeometry.h"
#include "cgeometry_supply.h"
#include "cfvm.h"
#include "cpost_processing.h"
#include "cheat_solver.h"
#include "cflow_solver.h"
#include "cfluid_properties.h"
#include "cfluid_properties_supply.h"
#include "cmatrix_solver.h"
#include "cmatrix_solver_supply.h"

#endif // CROOT_H
