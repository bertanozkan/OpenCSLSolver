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

#ifndef CFVM_H
#define CFVM_H

#include "vector"
#include "tensor.h"
#include "cgeometry.h"
#include "memory"
#include "cmatrix_solver.h"



class CFVM
{
public:
    CFVM();
    CFVM(Cgeometry &geo);

    class CFVM_cell    {
    public:
        CFVM_cell();
        size_t index = 0;
        double volume = 0.0;
    };

    class CFVM_face
    {
    public:
        CFVM_face();

        Svector area_vector;
        Svector cell_center_unit_vector;
        Svector alpha;
        double area = 0.0;
        double cell_center_distance = 0.0;
        double diffusion_ksi_A = 0.0;
        double geo_aver = 0.0;
        CFVM_cell* left_cell_ptr = nullptr;
        CFVM_cell* right_cell_ptr = nullptr;
        int global_face_index = -1;
        int left_cell_index = -1;
        int right_cell_index = -1;
    };

    class CFVM_BC_face
    {
    public:
        CFVM_BC_face();

        Svector area_vector;
        Svector cell_center_unit_vector;
        Svector alpha;
        double area = 0;
        double cell_center_distance = 0;
        double diffusion_ksi_A = 0;
        double perp_dist = 0.0;
        CFVM_cell* left_cell_ptr = nullptr;
        int left_cell_index = -1;
        int global_face_index = -1;
        int BC_index;
    };


    std::vector<CFVM_cell> cells;
    std::vector<CFVM_face> faces;
    std::vector<CFVM_BC_face> BC_faces;


    void define_matrix_structure( std::shared_ptr<CMatrix_solver> &Matrix_solver);

    void create_FVM_from_geometry();

};



#include "croot.h"



#endif // CFVM_H
