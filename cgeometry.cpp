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

#include "cgeometry.h"

Cgeometry::Cgeo_cell::Cgeo_cell(){}

Cgeometry::Cgeo_face::Cgeo_face(){}

Cgeometry::Cgeometry(){}

Cgeometry::Cgeo_cell::Cgeo_cell(enums::Ecell_type ecell_type)
{
    this->volume = 0.0;
    this->ecell_type = ecell_type;
    if(ecell_type == enums::Ecell_type::hexahedron)
    {
        this->i_nodes.resize(8,-1);
        this->i_faces.resize(6,-1);
        this->i_neighbours.resize(6,-1);
    }else if (ecell_type == enums::Ecell_type::wedge)
    {
        this->i_nodes.resize(6,-1);
        this->i_faces.resize(5,-1);
        this->i_neighbours.resize(5,-1);
    }else
    {
        this->i_nodes.resize(4,-1);
        this->i_faces.resize(4,-1);
        this->i_neighbours.resize(4,-1);
    }
}

Cgeometry::Cgeo_face::Cgeo_face(enums::Eface_type eface_type)
{
    this->eface_type = eface_type;
    if(eface_type == enums::Eface_type::quad)
        this->i_nodes.resize(4,-1);
    else
        this->i_nodes.resize(3,-1);
    this->area = 0.0;
    this->boundary_marker = -1;
    this->left_cell_index = -1;
    this->right_cell_index = -1;
}


bool Cgeometry::Cgeo_face::is_boundary()
{
    if(this->boundary_marker == -1) return false;
    else return true;
}


double Cgeometry::tet_volume(Svector a , Svector b , Svector c , Svector d)
{
    Svector vector0 , vector1 , vector2;
    vector0 = b - a; vector1 = c - a; vector2 = d - a;
    return (1.0/6.0) * dot(vector0,(vector1*vector2));
}


Cgeometry::Cgeo_face Cgeometry::return_face(Cgeo_cell &cell , int i_face)
{
    Cgeo_face face;
    if(cell.ecell_type == enums::Ecell_type::tetra)
    {
        if(i_face == 0)
        {
            face = Cgeo_face(enums::Eface_type::tri);
            face.i_nodes[0] = cell.i_nodes[1];
            face.i_nodes[1] = cell.i_nodes[0];
            face.i_nodes[2] = cell.i_nodes[2];
        }
        else if(i_face == 1)
        {
            face = Cgeo_face(enums::Eface_type::tri);
            face.i_nodes[0] = cell.i_nodes[0];
            face.i_nodes[1] = cell.i_nodes[1];
            face.i_nodes[2] = cell.i_nodes[3];
        }
        else if(i_face == 2)
        {
            face = Cgeo_face(enums::Eface_type::tri);
            face.i_nodes[0] = cell.i_nodes[1];
            face.i_nodes[1] = cell.i_nodes[2];
            face.i_nodes[2] = cell.i_nodes[3];
        }
        else if(i_face == 3)
        {
            face = Cgeo_face(enums::Eface_type::tri);
            face.i_nodes[0] = cell.i_nodes[2];
            face.i_nodes[1] = cell.i_nodes[0];
            face.i_nodes[2] = cell.i_nodes[3];
        }
        else
        {
            std::cout << " Error! Cgeometry::return_face failed " <<
                         "(face connectivity issue or mesh issue) " << std::endl;
            std::terminate();
        }
    }
    else if(cell.ecell_type == enums::Ecell_type::wedge)
    {

        if(i_face == 0)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[0];
            face.i_nodes[1] = cell.i_nodes[1];
            face.i_nodes[2] = cell.i_nodes[4];
            face.i_nodes[3] = cell.i_nodes[3];
        }
        else if(i_face == 1)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[1];
            face.i_nodes[1] = cell.i_nodes[2];
            face.i_nodes[2] = cell.i_nodes[5];
            face.i_nodes[3] = cell.i_nodes[4];
        }
        else if(i_face == 2)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[2];
            face.i_nodes[1] = cell.i_nodes[0];
            face.i_nodes[2] = cell.i_nodes[3];
            face.i_nodes[3] = cell.i_nodes[5];
        }
        else if(i_face == 3)
        {
            face = Cgeo_face(enums::Eface_type::tri);
            face.i_nodes[0] = cell.i_nodes[0];
            face.i_nodes[1] = cell.i_nodes[2];
            face.i_nodes[2] = cell.i_nodes[1];
        }
        else if(i_face == 4)
        {
            face = Cgeo_face(enums::Eface_type::tri);
            face.i_nodes[0] = cell.i_nodes[3];
            face.i_nodes[1] = cell.i_nodes[4];
            face.i_nodes[2] = cell.i_nodes[5];
        }
        else
        {
            std::cout << " Error! Cgeometry::return_face failed " <<
                         "(face connectivity issue or mesh issue) " << std::endl;
            std::terminate();
        }
    }
    else if(cell.ecell_type == enums::Ecell_type::hexahedron)
    {
        if(i_face == 0)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[0];
            face.i_nodes[1] = cell.i_nodes[1];
            face.i_nodes[2] = cell.i_nodes[5];
            face.i_nodes[3] = cell.i_nodes[4];
        }
        else if(i_face == 1)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[1];
            face.i_nodes[1] = cell.i_nodes[3];
            face.i_nodes[2] = cell.i_nodes[7];
            face.i_nodes[3] = cell.i_nodes[5];
        }
        else if(i_face == 2)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[3];
            face.i_nodes[1] = cell.i_nodes[2];
            face.i_nodes[2] = cell.i_nodes[6];
            face.i_nodes[3] = cell.i_nodes[7];
        }
        else if(i_face == 3)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[2];
            face.i_nodes[1] = cell.i_nodes[0];
            face.i_nodes[2] = cell.i_nodes[4];
            face.i_nodes[3] = cell.i_nodes[6];
        }
        else if(i_face == 4)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[1];
            face.i_nodes[1] = cell.i_nodes[0];
            face.i_nodes[2] = cell.i_nodes[2];
            face.i_nodes[3] = cell.i_nodes[3];
        }
        else if(i_face == 5)
        {
            face = Cgeo_face(enums::Eface_type::quad);
            face.i_nodes[0] = cell.i_nodes[4];
            face.i_nodes[1] = cell.i_nodes[5];
            face.i_nodes[2] = cell.i_nodes[7];
            face.i_nodes[3] = cell.i_nodes[6];
        }
        else
        {
            std::cout << " Error! Cgeometry::return_face failed " <<
                         "(face connectivity issue or mesh issue) " << std::endl;
            std::terminate();
        }
    }

    return face;
}


void Cgeometry::add_BC_face(int i_cell, int i_cell_face, int BC_marker)
{
    int i_face = this->faces.size();
    this->faces.push_back(return_face(this->cells[i_cell],i_cell_face));
    this->faces[i_face].boundary_marker = BC_marker;
    this->faces[i_face].left_cell_index = i_cell;
    this->cells[i_cell].i_faces[i_cell_face] = i_face;
}

void Cgeometry::add_BC_face_tri(int nodes_inp[3] , int BC_marker)
{
    // Finding the cell which this face belongs to
    std::vector<int> cell_cand;
    for (size_t i = 0 ; i < node_cell_indexes[nodes_inp[0]].size() ; i++) {
        cell_cand.push_back(node_cell_indexes[nodes_inp[0]][i]);
    }
    for (size_t i = 1 ; i < 3 ; i++) {
        for (size_t j = 0 ; j < cell_cand.size() ; j++) {
            if(std::find(node_cell_indexes[nodes_inp[i]].begin(), node_cell_indexes[nodes_inp[i]].end(), cell_cand[j]) != node_cell_indexes[nodes_inp[i]].end()) {
            } else {
                cell_cand.erase (cell_cand.begin()+j);
                j--;
            }
        }
    }

    // Check if only 1 cell is remained
    if(cell_cand.size() != 1)
        std::cout << " Error! something wrong with the boundaries of mesh file. " << std::endl;

    int i_cell = cell_cand[0];
    //Finding which face of the cell is the right one
    Cgeo_face dumface;
    int i_cell_face = 0;
    for (size_t i = 0 ; i < cells[i_cell].i_faces.size() ; i++) {
        dumface = return_face(cells[i_cell],i);
        bool found = true;
        for (int j = 0 ; j < 3 ; j++) {
            if(std::find(dumface.i_nodes.begin(), dumface.i_nodes.end(), nodes_inp[j]) != dumface.i_nodes.end()) {
            } else {
                found = false;
                break;
            }
        }
        if(found)
        {
            i_cell_face = i;
            break;
        }
    }

    int i_face = this->faces.size();
    this->faces.push_back(dumface);
    this->faces[i_face].boundary_marker = BC_marker;
    this->faces[i_face].left_cell_index = i_cell;
    this->cells[i_cell].i_faces[i_cell_face] = i_face;
}

void Cgeometry::add_BC_face_quad(int nodes_inp[4] , int BC_marker)
{
    // Finding the cell which this face belongs to
    std::vector<int> cell_cand;
    for (size_t i = 0 ; i < node_cell_indexes[nodes_inp[0]].size() ; i++) {
        cell_cand.push_back(node_cell_indexes[nodes_inp[0]][i]);
    }
    for (size_t i = 1 ; i < 4 ; i++) {
        for (size_t j = 0 ; j < cell_cand.size() ; j++) {
            if(std::find(node_cell_indexes[nodes_inp[i]].begin(), node_cell_indexes[nodes_inp[i]].end(), cell_cand[j]) != node_cell_indexes[nodes_inp[i]].end()) {
            } else {
                cell_cand.erase (cell_cand.begin()+j);
                j--;
            }
        }
    }

    // Check if only 1 cell is remained
    if(cell_cand.size() != 1)
        std::cout << " Error! something wrong with the boundaries of mesh file. " << std::endl;

    int i_cell = cell_cand[0];
    //Finding which face of the cell is the right one
    Cgeo_face dumface;
    int i_cell_face = 0;
    for (size_t i = 0 ; i < cells[i_cell].i_faces.size() ; i++) {
        dumface = return_face(cells[i_cell],i);
        bool found = true;
        for (int j = 0 ; j < 4 ; j++) {
            if(std::find(dumface.i_nodes.begin(), dumface.i_nodes.end(), nodes_inp[j]) != dumface.i_nodes.end()) {
            } else {
                found = false;
                break;
            }
        }
        if(found)
        {
            i_cell_face = i;
            break;
        }
    }

    int i_face = this->faces.size();
    this->faces.push_back(dumface);
    this->faces[i_face].boundary_marker = BC_marker;
    this->faces[i_face].left_cell_index = i_cell;
    this->cells[i_cell].i_faces[i_cell_face] = i_face;
}


void Cgeometry::calc_face_area_and_normal(Cgeo_face &geo_face)
{
    Svector vector0, vector1;
    double length;

    vector0 = this->nodes[geo_face.i_nodes[1]]
                - this->nodes[geo_face.i_nodes[0]];
    vector1 = this->nodes[geo_face.i_nodes[2]]
                - this->nodes[geo_face.i_nodes[0]];

    geo_face.normal = vector0 * vector1;
    length = mag(geo_face.normal);
    geo_face.normal = (1.0/length) * geo_face.normal;
    geo_face.tan1 = (1.0/mag(vector0)) * vector0;
    geo_face.tan2 = geo_face.normal * geo_face.tan1;
    geo_face.tan2 = (1.0/mag(geo_face.tan2))*geo_face.tan2;
    geo_face.area = 0.5 * length;

    geo_face.center = this->nodes[geo_face.i_nodes[0]]
            + this->nodes[geo_face.i_nodes[1]]
            + this->nodes[geo_face.i_nodes[2]];

    geo_face.center = (1.0/3.0) * geo_face.center;

    if(geo_face.eface_type == enums::Eface_type::quad)
    {
        geo_face.center = geo_face.area * geo_face.center;
        vector0 = this->nodes[geo_face.i_nodes[2]]
                    - this->nodes[geo_face.i_nodes[0]];
        vector1 = this->nodes[geo_face.i_nodes[3]]
                    - this->nodes[geo_face.i_nodes[0]];
        length = mag(vector0 * vector1);
        geo_face.area = geo_face.area + (0.5 * length);

        geo_face.center = geo_face.center + ((0.5 * length / 3.0) *
                (this->nodes[geo_face.i_nodes[0]]
                + this->nodes[geo_face.i_nodes[2]]
                + this->nodes[geo_face.i_nodes[3]]));

        geo_face.center = (1.0 / geo_face.area) * geo_face.center;
    }
}



void Cgeometry::calc_cell_volume_and_center(Cgeo_cell &geo_cell)
{
    if(geo_cell.ecell_type == enums::Ecell_type::tetra)
    {
        geo_cell.volume = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[1]]
                , this->nodes[geo_cell.i_nodes[2]]
                , this->nodes[geo_cell.i_nodes[3]]);

        geo_cell.center = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[1]]
                + this->nodes[geo_cell.i_nodes[2]]
                + this->nodes[geo_cell.i_nodes[3]];

        geo_cell.center = (1.0/4.0) * geo_cell.center;
    }

    if(geo_cell.ecell_type == enums::Ecell_type::wedge)
    {
        double dumvol;
        Svector center_dum;

        geo_cell.volume = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[4]]
                , this->nodes[geo_cell.i_nodes[5]]
                , this->nodes[geo_cell.i_nodes[3]]);

        geo_cell.center = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[4]]
                + this->nodes[geo_cell.i_nodes[5]]
                + this->nodes[geo_cell.i_nodes[3]];

        geo_cell.center = ( geo_cell.volume/4.0) * geo_cell.center;


        dumvol = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[4]]
                , this->nodes[geo_cell.i_nodes[1]]
                , this->nodes[geo_cell.i_nodes[5]]);

        geo_cell.volume += dumvol;

        center_dum = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[4]]
                + this->nodes[geo_cell.i_nodes[1]]
                + this->nodes[geo_cell.i_nodes[5]];

        center_dum = ( dumvol/4.0) * center_dum;
        geo_cell.center = geo_cell.center + center_dum;


        dumvol = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[1]]
                , this->nodes[geo_cell.i_nodes[2]]
                , this->nodes[geo_cell.i_nodes[5]]);

        geo_cell.volume += dumvol;

        center_dum = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[1]]
                + this->nodes[geo_cell.i_nodes[2]]
                + this->nodes[geo_cell.i_nodes[5]];

        center_dum = ( dumvol/4.0) * center_dum;
        geo_cell.center = geo_cell.center + center_dum;


        geo_cell.center = (1.0/geo_cell.volume) * geo_cell.center;
    }


    if(geo_cell.ecell_type == enums::Ecell_type::hexahedron)
    {
        double dumvol;
        Svector center_dum;

        geo_cell.volume = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[2]]
                , this->nodes[geo_cell.i_nodes[7]]
                , this->nodes[geo_cell.i_nodes[3]]);

        geo_cell.center = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[2]]
                + this->nodes[geo_cell.i_nodes[7]]
                + this->nodes[geo_cell.i_nodes[3]];

        geo_cell.center = ( geo_cell.volume/4.0) * geo_cell.center;


        dumvol = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[2]]
                , this->nodes[geo_cell.i_nodes[6]]
                , this->nodes[geo_cell.i_nodes[7]]);

        geo_cell.volume += dumvol;

        center_dum = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[2]]
                + this->nodes[geo_cell.i_nodes[6]]
                + this->nodes[geo_cell.i_nodes[7]];

        center_dum = ( dumvol/4.0) * center_dum;
        geo_cell.center = geo_cell.center + center_dum;


        dumvol = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[7]]
                , this->nodes[geo_cell.i_nodes[1]]
                , this->nodes[geo_cell.i_nodes[3]]);

        geo_cell.volume += dumvol;

        center_dum = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[7]]
                + this->nodes[geo_cell.i_nodes[1]]
                + this->nodes[geo_cell.i_nodes[3]];

        center_dum = ( dumvol/4.0) * center_dum;
        geo_cell.center = geo_cell.center + center_dum;


        dumvol = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[7]]
                , this->nodes[geo_cell.i_nodes[5]]
                , this->nodes[geo_cell.i_nodes[1]]);

        geo_cell.volume += dumvol;

        center_dum = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[7]]
                + this->nodes[geo_cell.i_nodes[5]]
                + this->nodes[geo_cell.i_nodes[1]];

        center_dum = ( dumvol/4.0) * center_dum;
        geo_cell.center = geo_cell.center + center_dum;


        dumvol = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[7]]
                , this->nodes[geo_cell.i_nodes[4]]
                , this->nodes[geo_cell.i_nodes[5]]);

        geo_cell.volume += dumvol;

        center_dum = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[7]]
                + this->nodes[geo_cell.i_nodes[4]]
                + this->nodes[geo_cell.i_nodes[5]];

        center_dum = ( dumvol/4.0) * center_dum;
        geo_cell.center = geo_cell.center + center_dum;


        dumvol = tet_volume(this->nodes[geo_cell.i_nodes[0]]
                , this->nodes[geo_cell.i_nodes[7]]
                , this->nodes[geo_cell.i_nodes[6]]
                , this->nodes[geo_cell.i_nodes[4]]);

        geo_cell.volume += dumvol;

        center_dum = this->nodes[geo_cell.i_nodes[0]]
                + this->nodes[geo_cell.i_nodes[7]]
                + this->nodes[geo_cell.i_nodes[6]]
                + this->nodes[geo_cell.i_nodes[4]];

        center_dum = ( dumvol/4.0) * center_dum;
        geo_cell.center = geo_cell.center + center_dum;


        geo_cell.center = (1.0/geo_cell.volume) * geo_cell.center;
    }
}


bool Cgeometry::if_equal
(Cgeometry::Cgeo_face left_face , Cgeometry::Cgeo_face right_face)
{
    if(left_face.eface_type != right_face.eface_type)
        return false;

    for(auto& i : left_face.i_nodes)
    {
        auto result = std::find(std::begin(right_face.i_nodes) ,
                            std::end(right_face.i_nodes) , i);
        if(result == std::end(right_face.i_nodes))
            return false;
    }

    return true;
}


void Cgeometry::check_cells_insert_face
(Cgeometry::Cgeo_cell &left_cell , Cgeometry::Cgeo_cell &right_cell
 ,int common_i_node , int i_left_cell , int i_right_cell)
{
    Cgeo_face dumface_left;
    Cgeo_face dumface_right;
    int n_face_left = left_cell.i_faces.size();
    int n_face_right = right_cell.i_faces.size();
    for(int i_face_left = 0 ; i_face_left < n_face_left ; i_face_left++)
    {
        dumface_left = return_face(left_cell,i_face_left);
        auto result = std::find(std::begin(dumface_left.i_nodes)
                                ,std::end(dumface_left.i_nodes),common_i_node);
        if(result == dumface_left.i_nodes.end())
            continue;
        for(int i_face_right = 0 ; i_face_right < n_face_right ; i_face_right++)
        {
            dumface_right = return_face(right_cell,i_face_right);
            auto result = std::find(std::begin(dumface_right.i_nodes)
                                    ,std::end(dumface_right.i_nodes),common_i_node);
            if(result == dumface_right.i_nodes.end())
                continue;
            // If these faces are the same face
            // Insert left face in geometry class
            if(if_equal(dumface_left,dumface_right))
            {
                // If this face is not allready created
                if(left_cell.i_neighbours[i_face_left] == -1)
                {
                    int i_face = this->faces.size();
                    dumface_left.left_cell_index = i_left_cell;
                    dumface_left.right_cell_index = i_right_cell;
                    this->faces.push_back(dumface_left);

                    left_cell.i_faces[i_face_left] = i_face;
                    right_cell.i_faces[i_face_right] = i_face;
                    left_cell.i_neighbours[i_face_left] = i_right_cell;
                    right_cell.i_neighbours[i_face_right] = i_left_cell;
                }
            }
        }
    }
}


void Cgeometry::calc_cell_face_connectivity()
{

    node_cell_indexes.resize(this->nodes.size());
    for(size_t i_cell = 0 ; i_cell < this->cells.size() ; i_cell++)
    {
        for(size_t i_node = 0 ; i_node < this->cells[i_cell].i_nodes.size() ; i_node++)
        {
            node_cell_indexes[this->cells[i_cell].i_nodes[i_node]].push_back(i_cell);
        }
    }

    this->faces.reserve(this->cells.size());

    // Loop over nodes
    for(size_t i_node = 0 ; i_node < node_cell_indexes.size() ; i_node ++)
    {
        // Loop over cells that has "i_node" lets name this cell left cell
        for(size_t i_cell_left = 0 ; i_cell_left < node_cell_indexes[i_node].size()
            ; i_cell_left ++)
        {
            // Over other cells that has "i_node" lets name them right cell
            for(size_t i_cell_right = 0;i_cell_right < node_cell_indexes[i_node].size()
                ; i_cell_right ++)
            {
                // Checking if right cell and left cell is not the same cell
                if(i_cell_left != i_cell_right)
                {
                    // Comparing two cells and adding new face if approprite
                    check_cells_insert_face(
                            this->cells[node_cell_indexes[i_node][i_cell_left]]
                            , this->cells[node_cell_indexes[i_node][i_cell_right]]
                            , i_node , node_cell_indexes[i_node][i_cell_left]
                            , node_cell_indexes[i_node][i_cell_right]);
                }
            }
        }
    }
}



void Cgeometry::calc_geometric_values()
{
    for(auto& it : this->faces)
        calc_face_area_and_normal(it);

    for(auto& it : this->cells)
        calc_cell_volume_and_center(it);
}


void Cgeometry::extract_BC_faces()
{
    this->BC_face_global_indexes.reserve(this->n_BC_face);
    for(size_t i = 0 ; i < this->faces.size() ; i++)
    {
        if(this->faces[i].boundary_marker != -1)
            this->BC_face_global_indexes.push_back(i);
    }
}
