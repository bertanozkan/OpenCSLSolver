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

#include "cfvm.h"
#ifdef COMPILE_WITH_OPENMP
#include <omp.h>
#endif

CFVM::CFVM_cell::CFVM_cell(){}

CFVM::CFVM_face::CFVM_face(){}

CFVM::CFVM_BC_face::CFVM_BC_face(){}

CFVM::CFVM(){}

/*
void CFVM::define_iterator_ranges()
{
    int n_thread = omp_get_num_threads();
    it_cell_start.resize(n_thread);
    it_cell_end.resize(n_thread);
    it_face_start.resize(n_thread);
    it_face_end.resize(n_thread);
    it_BC_face_start.resize(n_thread);
    it_BC_face_end.resize(n_thread);
    size_t cell_chunk_size = cells.size() / n_thread;
    size_t face_chunk_size = faces.size() / n_thread;
    size_t BC_face_chunk_size = BC_faces.size() / n_thread;
    for (int i_thread = 0 ; i_thread < n_thread ; i_thread++)
    {
        it_cell_start[i_thread] = cells.begin();
        it_face_start[i_thread] = faces.begin();
        it_BC_face_start[i_thread] = BC_faces.begin();
        std::advance(it_cell_start[i_thread], i_thread * cell_chunk_size);
        std::advance(it_face_start[i_thread], i_thread * face_chunk_size);
        std::advance(it_BC_face_start[i_thread], i_thread * BC_face_chunk_size);
        it_cell_end[i_thread] = it_cell_start[i_thread];
        it_face_end[i_thread] = it_face_start[i_thread];
        it_BC_face_end[i_thread] = it_BC_face_start[i_thread];
        if(i_thread == (n_thread - 1))
        {
            it_cell_end[i_thread] = cells.end();
            it_face_end[i_thread] = faces.end();
            it_BC_face_end[i_thread] = BC_faces.end();
        }
        else
        {
            std::advance(it_cell_end[i_thread], cell_chunk_size);
            std::advance(it_face_end[i_thread], face_chunk_size);
            std::advance(it_BC_face_end[i_thread], BC_face_chunk_size);
        }
    }
}*/


void CFVM::define_matrix_structure( std::shared_ptr<CMatrix_solver> &Matrix_solver)
{
    // Calculate total number of nonzeros on matrix
    size_t nnz = 0;
    for(size_t i = 0 ; i < cells.size() ; i++)
    {
        nnz += root.geometry->cells[i].i_neighbours.size();
    }

    // Initialize linear equation system and solver
    Matrix_solver = std::shared_ptr<CMatrix_solver> (new CMatrix_solver(cells.size(),cells.size(),nnz));

    // Create matrix structure
    for(size_t i_cell = 0; i_cell < root.geometry->cells.size() ; i_cell ++)
    {
        Matrix_solver->matrix.insert_element(i_cell , i_cell , 0.0);
        for(size_t i_neig = 0 ; i_neig < root.geometry->cells[i_cell].i_neighbours.size() ; i_neig++)
        {
            int neigh_index = root.geometry->cells[i_cell].i_neighbours[i_neig];
            Matrix_solver->matrix.insert_element(i_cell , neigh_index , 0.0);
        }
    }
    Matrix_solver->M = Matrix_solver->matrix;
    Matrix_solver->extract_diag_indexes();
}

void CFVM::create_FVM_from_geometry()
{
    std::cout << "// - Creating FVM grid from geometry - //" << std::endl;
    cells.resize(root.geometry->cells.size());
    faces.reserve(root.geometry->faces.size());
    BC_faces.reserve(root.geometry->BC_face_global_indexes.size());

    // Creating cells
    for(size_t i_cell = 0 ; i_cell < root.geometry->cells.size() ; i_cell++)
    {
        cells[i_cell].volume = root.geometry->cells[i_cell].volume;
        cells[i_cell].index = i_cell;
    }

    // Creating faces
    for(size_t i_face = 0 ; i_face < root.geometry->faces.size() ; i_face++)
    {
        if(!root.geometry->faces[i_face].is_boundary())
        {

            CFVM_face dumface;
            Cgeometry::Cgeo_face geo_face = root.geometry->faces[i_face];
            dumface.global_face_index = i_face;
            dumface.left_cell_index = geo_face.left_cell_index;
            dumface.right_cell_index = geo_face.right_cell_index;
            dumface.left_cell_ptr = &cells[dumface.left_cell_index];
            dumface.right_cell_ptr = &cells[dumface.right_cell_index];
            dumface.area = geo_face.area;
            dumface.area_vector = geo_face.area * geo_face.normal;
            dumface.cell_center_unit_vector = (root.geometry->cells[geo_face.right_cell_index].center
                    - root.geometry->cells[geo_face.left_cell_index].center);
            dumface.alpha = (1.0/dot(dumface.area_vector,dumface.cell_center_unit_vector))*dumface.area_vector;
            dumface.cell_center_distance = mag(dumface.cell_center_unit_vector);
            dumface.cell_center_unit_vector = (1.0/dumface.cell_center_distance) * dumface.cell_center_unit_vector;
            //dumface.diffusion_ksi_A = dumface.area / dumface.cell_center_distance;
            dumface.diffusion_ksi_A = dot(dumface.alpha,dumface.area_vector);

            dumface.geo_aver = mag(root.geometry->cells[geo_face.right_cell_index].center - geo_face.center) /
                    dumface.cell_center_distance;

            faces.push_back(dumface);
        }
    }

    // Creating BC_faces
    for(size_t i_face = 0 ; i_face < root.geometry->BC_face_global_indexes.size() ; i_face++)
    {
        CFVM_BC_face dumBCface;
        Cgeometry::Cgeo_face geo_face =
                root.geometry->faces[root.geometry->BC_face_global_indexes[i_face]];
        dumBCface.global_face_index = root.geometry->BC_face_global_indexes[i_face];
        dumBCface.left_cell_index = geo_face.left_cell_index;
        dumBCface.left_cell_ptr = &cells[dumBCface.left_cell_index];
        dumBCface.area = geo_face.area;
        dumBCface.area_vector = geo_face.area * geo_face.normal;
        dumBCface.cell_center_unit_vector = geo_face.center
                - root.geometry->cells[geo_face.left_cell_index].center;
        dumBCface.cell_center_distance = mag(geo_face.center
                         - root.geometry->cells[geo_face.left_cell_index].center);
        //dumBCface.diffusion_ksi_A = dumBCface.area / dumBCface.cell_center_distance;

        dumBCface.alpha = (1.0/dot(dumBCface.area_vector,geo_face.center
                            - root.geometry->cells[geo_face.left_cell_index].center))*dumBCface.area_vector;

        dumBCface.cell_center_unit_vector = (1.0/dumBCface.cell_center_distance) * dumBCface.cell_center_unit_vector;

        dumBCface.diffusion_ksi_A = dot(dumBCface.alpha,dumBCface.area_vector);

        dumBCface.perp_dist = dot((dumBCface.cell_center_distance * dumBCface.cell_center_unit_vector) , ((1.0/dumBCface.area) * dumBCface.area_vector));

        BC_faces.push_back(dumBCface);
    }
    std::cout << "// - Creating FVM grid from geometry Ended - //" << std::endl;
}


