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

#ifndef CGEOMETRY_H
#define CGEOMETRY_H

#include "vector"
#include "math.h"
#include "string"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>

#include "tensor.h"
#include "enums.h"



class Cgeometry
{
public:
    Cgeometry();

    class Cgeo_face
    {
    public:
        Cgeo_face(enums::Eface_type eface_type);
        Cgeo_face();

        enums::Eface_type eface_type;
        Svector center;
        Svector normal;
        Svector tan1;
        Svector tan2;
        double area = 0;;
        int left_cell_index = -1;
        int right_cell_index = -1;
        int boundary_marker = -1;
        std::vector<int> i_nodes;

        bool is_boundary();

    };

    class Cgeo_cell
    {
    public:
        Cgeo_cell();
        Cgeo_cell(enums::Ecell_type ecell_type);

        enums::Ecell_type ecell_type;
        Svector center;
        double volume = 0;
        std::vector<int> i_nodes;
        std::vector<int> i_faces;
        std::vector<int> i_neighbours;

    };

    int n_boundary = 0;
    int n_BC_face = 0;
    std::vector<Cgeo_cell> cells;
    std::vector<Cgeo_face> faces;
    std::vector<Svector> nodes;
    std::vector<int> BC_face_global_indexes;
    std::vector<std::string> BC_tags;
    std::vector<std::vector<int>> node_cell_indexes;

    double tet_volume(Svector a , Svector b , Svector c, Svector d);
    Cgeo_face return_face(Cgeo_cell &cell , int i_face);
    void add_BC_face(int i_cell, int i_cell_face, int BC_marker);
    void add_BC_face_tri(int nodes_inp[3] , int BC_marker);
    void add_BC_face_quad(int nodes_inp[4] , int BC_marker);

    void calc_face_area_and_normal(Cgeo_face &geo_face);
    void calc_cell_volume_and_center(Cgeo_cell &geo_cell);

    bool if_equal(Cgeometry::Cgeo_face left_face , Cgeometry::Cgeo_face right_face);
    void check_cells_insert_face(Cgeometry::Cgeo_cell &left_cell ,
                                 Cgeometry::Cgeo_cell &right_cell, int common_i_node
                                 , int i_left_cell, int i_right_cell);
    void calc_cell_face_connectivity();
    void calc_geometric_values();
    void extract_BC_faces();

};

#endif // CGEOMETRY_H
