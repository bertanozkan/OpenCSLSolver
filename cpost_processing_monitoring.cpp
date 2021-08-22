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

#include "cpost_processing.h"

void Cpost_processing::monitoring()
{
    std::cout.setf(std::ios::left, std::ios::adjustfield);

    if(root.configuration->eproblem_type == enums::Eproblem_type::Heat_conduction)
    {
        if(root.i_iteration % 10 == 0)
        {
            std::cout.width(22);
            std::cout.fill('-');
            std::cout << '-' << std::endl;
            std::cout.fill(' ');
            std::cout << "  Iter  " << "||" << "    Res   " << "||" << std::endl;
            std::cout.width(22);
            std::cout.fill('-');
            std::cout << '-' << std::endl;
            std::cout.fill(' ');
        }
        std::cout << std::setw(8) << root.i_iteration << "||" << std::setw(10)
                  << root.Heat_solver->temperature_global_residual << "||" << std::endl;
    }

    if(root.configuration->eproblem_type == enums::Eproblem_type::Fluid_flow)
    {
        if(root.configuration->etimestepping_type == enums::Etimestepping_type::Steady )
        {
            if(root.i_iteration % 10 == 0)
            {
                std::cout.width(100);
                std::cout.fill('-');
                std::cout << '-' << std::endl;
                std::cout.fill(' ');
                std::cout << "  Iter   " << "||" << " Res Mom X  " << "||" << " Res Mom Y  " << "||" << " Res Mom Z  " << "||" << " Press Res  "<< "||" << " Enthal Res "<< "||" << " UND = "<<
                             root.Flow_solver->implicit_underrelaxation << std::endl;
                std::cout.width(100);
                std::cout.fill('-');
                std::cout << '-' << std::endl;
                std::cout.fill(' ');
            }
            std::cout << std::setw(8) << root.i_iteration;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Momentum_global_residual.val[0] ;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Momentum_global_residual.val[1] ;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Momentum_global_residual.val[2] ;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Pressure_global_residual;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Enthalpy_global_residual << "||" ;
            std::cout << std::endl;
        }
        else if(root.configuration->etimestepping_type == enums::Etimestepping_type::Transient)
        {

            if(root.i_inner_iter == 0)
            {
                std::cout.width(100);
                std::cout.fill('-');
                std::cout << '-' << std::endl;
                std::cout.fill(' ');
                std::cout << "  Iter   " << "||" << " Res Mom X  " << "||" << " Res Mom Y  " << "||" << " Res Mom Z  " << "||" << " Press Res  "<< "||" << " Enthal Res "<< "||" << std::endl;
                std::cout.width(100);
                std::cout.fill('-');
                std::cout << '-' << std::endl;
                std::cout.fill(' ');
            }

            std::cout << std::setw(9) << root.i_inner_iter;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Momentum_global_residual.val[0] ;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Momentum_global_residual.val[1] ;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Momentum_global_residual.val[2] ;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Pressure_global_residual;
            std::cout << "||" << std::setw(12)<< root.Flow_solver->Enthalpy_global_residual << "||" ;
            std::cout << std::endl;

            if(root.i_inner_iter == root.configuration->max_inner_iteration-1)
            {
                std::cout.width(100);
                std::cout.fill('-');
                std::cout << '-' << std::endl;
                std::cout << "Outer iteration : " << root.i_iteration << " Ended " << std::endl;
                std::cout.width(100);
                std::cout.fill('-');
                std::cout << '-' << std::endl;
                std::cout.fill(' ');
            }
        }

    }

}
