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

#include "cpost_processing.h"

Cpost_processing::Cpost_processing(){}


void Cpost_processing::write_VTK_data()
{

    this->fstream_VTK<<"SCALARS Volume[m3] double 1\n"
                    <<"LOOKUP_TABLE default"<<std::endl;

    for(auto& i :root.geometry->cells)
        this->fstream_VTK<< i.volume <<std::endl;

    this->fstream_VTK<<"VECTORS Center double" <<std::endl;

    for(auto& i :root.geometry->cells)
        this->fstream_VTK<< i.center.val[0] << " " << i.center.val[1] << " " << i.center.val[2] <<std::endl;

    if(root.configuration->eproblem_type == enums::Eproblem_type::Heat_conduction)
    {
        this->fstream_VTK<<"SCALARS Temperature[K] double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Heat_solver->cells)
            this->fstream_VTK<< i.Temperature <<std::endl;
    }

    if(root.configuration->eproblem_type == enums::Eproblem_type::Fluid_flow)
    {
        this->fstream_VTK<<"SCALARS Density[Kg/m^3] double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.Density <<std::endl;

        this->fstream_VTK<<"SCALARS Pressure[Pa] double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.Pressure <<std::endl;

        /*this->fstream_VTK<<"SCALARS ap double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.ap <<std::endl;

        this->fstream_VTK<<"SCALARS Pcorr double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.Pressure_corr <<std::endl;*/

        this->fstream_VTK<<"SCALARS Temperature[K] double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.Temperature <<std::endl;

        this->fstream_VTK<<"SCALARS viscosity double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< root.Fluid_properties->get_viscosity(i.Temperature) <<std::endl;

       /* this->fstream_VTK<<"SCALARS p_rho_der double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< root.Fluid_properties->get_Pressure_Rho_derivative(i.Temperature) <<std::endl;

        this->fstream_VTK<<"SCALARS total_mdot double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.total_mdot <<std::endl;*/

        this->fstream_VTK<<"SCALARS MACH double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.Mach <<std::endl;

        this->fstream_VTK<<"VECTORS Velocity[m/s] double"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.Velocity.val[0] << " " << i.Velocity.val[1] << " " << i.Velocity.val[2] <<std::endl;

        this->fstream_VTK<<"VECTORS Pres_grad double"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.Pressure_grad.val[0] << " " << i.Pressure_grad.val[1] << " " << i.Pressure_grad.val[2] <<std::endl;

        /*this->fstream_VTK<<"VECTORS Pres_corr_grad double"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.Pressure_corr_grad.val[0] << " " << i.Pressure_corr_grad.val[1] << " " << i.Pressure_corr_grad.val[2] <<std::endl;

        this->fstream_VTK<<"VECTORS ap_vector double"<<std::endl;

        for(auto& i :root.Flow_solver->cells)
            this->fstream_VTK<< i.ap_vec.val[0] << " " << i.ap_vec.val[1] << " " << i.ap_vec.val[2] <<std::endl;*/
    }


}

void Cpost_processing::write_boundary_VTK_data()
{
     this->fstream_VTK<<"SCALARS Area[m2] double 1\n"
                     <<"LOOKUP_TABLE default"<<std::endl;

     for(auto& i_face : root.geometry->BC_face_global_indexes)
         this->fstream_VTK<< root.geometry->faces[i_face].area <<std::endl;

     this->fstream_VTK<<"VECTORS normal double" <<std::endl;

     for(auto& i_face : root.geometry->BC_face_global_indexes)
         this->fstream_VTK<< root.geometry->faces[i_face].normal.val[0] << " " << root.geometry->faces[i_face].normal.val[1] << " " << root.geometry->faces[i_face].normal.val[2] <<std::endl;

    this->fstream_VTK<<"SCALARS Marker double 1\n"
                    <<"LOOKUP_TABLE default"<<std::endl;

    for(auto& i_face : root.geometry->BC_face_global_indexes)
        this->fstream_VTK<< (root.geometry->faces[i_face].boundary_marker)
                         <<std::endl;

    this->fstream_VTK<<"Vectors area_vector double" <<std::endl;

    for(auto& i_face : root.FVM->BC_faces)
        this->fstream_VTK<< (i_face.area_vector.val[0]) << " " << (i_face.area_vector.val[1]) << " " << (i_face.area_vector.val[2])
                         <<std::endl;

    this->fstream_VTK<<"Vectors center double" <<std::endl;

    for(auto& i_face : root.geometry->BC_face_global_indexes)
        this->fstream_VTK<< root.geometry->faces[i_face].center.val[0] << " " << root.geometry->faces[i_face].center.val[1] << " " << root.geometry->faces[i_face].center.val[2] <<std::endl;


    if(root.configuration->eproblem_type == enums::Eproblem_type::Heat_conduction)
    {
        this->fstream_VTK<<"SCALARS Temperature[K] double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i_face : root.Heat_solver->BC_faces)
            this->fstream_VTK<< i_face.Temperature <<std::endl;
    }

    if(root.configuration->eproblem_type == enums::Eproblem_type::Fluid_flow)
    {
        this->fstream_VTK<<"SCALARS MACH double 1\n"
                        <<"LOOKUP_TABLE default"<<std::endl;

        for(auto& i_face : root.Flow_solver->BC_faces)
            this->fstream_VTK<< i_face.left_cell_ptr->Mach <<std::endl;


        this->fstream_VTK<<"Vectors Stress_vector double" <<std::endl;

        auto BC_face = root.Flow_solver->BC_faces.begin();
        auto FVM_BC_face = root.FVM->BC_faces.begin();
        for(; BC_face != root.Flow_solver->BC_faces.end() ; BC_face++, FVM_BC_face++)
        {
            Svector wall_stress;
            wall_stress = root.Flow_solver->calc_wall_stress(*FVM_BC_face,*BC_face,root.Fluid_properties);
            this->fstream_VTK<< wall_stress.val[0] << " " << wall_stress.val[1] << " " << wall_stress.val[2] <<std::endl;
        }
    }


}


/*void Cpost_processing::write_CSV_file(std::string filename)
{
    std::cout<<"// - Writing Cell CSV File - //"<<std::endl;
    std::fstream csv;
    csv.open(filename, std::ios::out);
    if(!csv.is_open())
    {
        std::cout << " Error! file " << filename
                  << " cannot be opened !!!" << std::endl;
        std::terminate();
    }

    if(root.configuration->eproblem_type == enums::Heat_conduction)
    {
        // Writing header
        csv << "x_coord[m] ; y_coord[m] ; z_coord[m] ; Temperature[k]" << std::endl;
        for(size_t i_cell = 0 ; i_cell < root.geometry->cells.size() ; i_cell++)
        {
            csv << root.geometry->cells[i_cell].center[0] << " ; " <<
                   root.geometry->cells[i_cell].center[1] << " ; " <<
                   root.geometry->cells[i_cell].center[2] << " ; " <<
                   root.Heat_solver->cells[i_cell].Temperature << std::endl ;
        }
    }

    csv.close();
    std::cout<<"// - Writing Cell CSV File Ended - //"<<std::endl;
}


void Cpost_processing::write_boundary_CSV_file(std::string filename)
{
    std::cout<<"// - Writing Boundary CSV File - //"<<std::endl;
    std::fstream csv;
    csv.open(filename, std::ios::out);
    if(!csv.is_open())
    {
        std::cout << " Error! file " << filename
                  << " cannot be opened !!!" << std::endl;
        std::terminate();
    }

    if(root.configuration->eproblem_type == enums::Heat_conduction)
    {
        // Writing header
        csv << "Boundary_marker ; x_coord[m] ; y_coord[m] ; z_coord[m] ; area[m2] ; Temperature[k]" << std::endl;
        for(auto& face : root.Heat_solver->BC_faces)
        {
            int i_face = face.FVM_BC_face->global_face_index;
            csv << root.geometry->faces[i_face].boundary_marker << " ; " <<
                   root.geometry->faces[i_face].center[0] << " ; " <<
                   root.geometry->faces[i_face].center[1] << " ; " <<
                   root.geometry->faces[i_face].center[2] << " ; " <<
                   root.geometry->faces[i_face].area << " ; " <<
                   face.Temperature << std::endl ;
        }
    }

    csv.close();
    std::cout<<"// - Writing Boundary CSV File Ended - //"<<std::endl;

}
*/


void Cpost_processing::write_VTK_file(std::string filename)
{
    std::cout<<"// - Writing Cell Paraview File - //"<<std::endl;
    this->fstream_VTK.open(filename, std::ios::out);
    if(!this->fstream_VTK.is_open())
    {
        std::cout << " Error! file " << filename
                  << " cannot be opened !!!" << std::endl;
        std::terminate();
    }
    this->fstream_VTK<<"# vtk DataFile Version 1.0"<<std::endl;
    this->fstream_VTK<<"Solver Data File\nASCII"<<std::endl;

    this->fstream_VTK<<"DATASET UNSTRUCTURED_GRID\nPOINTS "
                    << root.geometry->nodes.size() << " double" << std::endl;

    // Writing nodes
    for(auto& i : root.geometry->nodes)
    {
         this->fstream_VTK << i[0] << " " << i[1] << " " << i[2] << std::endl;
    }

    // Calc cell list size
    int  vtk_cell_list_size = 0;
    for(auto& i : root.geometry->cells)
    {
        if(i.ecell_type == enums::Ecell_type::tetra) vtk_cell_list_size += 5;
        else if(i.ecell_type == enums::Ecell_type::wedge) vtk_cell_list_size += 7;
        if(i.ecell_type == enums::Ecell_type::hexahedron) vtk_cell_list_size += 9;
    }


    // Writing cells
    this->fstream_VTK<<"CELLS "<<(root.geometry->cells.size())
                    <<" "<<vtk_cell_list_size<< std::endl;
    for(auto& i : root.geometry->cells)
    {
        if(i.ecell_type == enums::Ecell_type::tetra)
        {
            this->fstream_VTK<<"4 "<<i.i_nodes[0] <<" "<<i.i_nodes[1]<<" ";
            this->fstream_VTK<<i.i_nodes[2] <<" "<<i.i_nodes[3]<<std::endl;
        }else if(i.ecell_type == enums::Ecell_type::wedge)
        {
            this->fstream_VTK<<"6 "<<i.i_nodes[0] <<" "<<i.i_nodes[2]<<" ";
            this->fstream_VTK<<i.i_nodes[1] <<" "<<i.i_nodes[3]<<" ";
            this->fstream_VTK<<i.i_nodes[5] <<" "<<i.i_nodes[4]<<std::endl;
        }else if(i.ecell_type == enums::Ecell_type::hexahedron)
        {
            this->fstream_VTK<<"8 "<<i.i_nodes[0] <<" "<<i.i_nodes[1]<<" ";
            this->fstream_VTK<<i.i_nodes[3] <<" "<<i.i_nodes[2]<<" ";
            this->fstream_VTK<<i.i_nodes[4] <<" "<<i.i_nodes[5]<<" ";
            this->fstream_VTK<<i.i_nodes[7] <<" "<<i.i_nodes[6]<<std::endl;
        }
    }


    this->fstream_VTK<<"CELL_TYPES "<<(root.geometry->cells.size())<<std::endl;
    // Cell types
    for(auto& i : root.geometry->cells)
    {
        if(i.ecell_type == enums::Ecell_type::tetra) this->fstream_VTK<<"10"<<std::endl;
        else if(i.ecell_type == enums::Ecell_type::wedge) this->fstream_VTK<<"13"<<std::endl;
        else if(i.ecell_type == enums::Ecell_type::hexahedron) this->fstream_VTK<<"12"<<std::endl;
    }


    this->fstream_VTK<<"CELL_DATA "<<(root.geometry->cells.size())<<std::endl;

    this->write_VTK_data();

    this->fstream_VTK.close();
    std::cout<<"// - Writing Cell Paraview File Ended - //"<<std::endl;
}


void Cpost_processing::write_boundary_VTK_file(std::string filename)
{
    std::cout<<"// - Writing Boundary Paraview File - //"<<std::endl;
    this->fstream_VTK.open(filename, std::ios::out);
    if(!this->fstream_VTK.is_open())
    {
        std::cout << " Error! file " << filename
                  << " cannot be opened !!!" << std::endl;
        std::terminate();
    }
    this->fstream_VTK<<"# vtk DataFile Version 1.0"<<std::endl;
    this->fstream_VTK<<"Solver Data File\nASCII"<<std::endl;

    this->fstream_VTK<<"DATASET UNSTRUCTURED_GRID\nPOINTS "
                    << root.geometry->nodes.size() << " double" << std::endl;

    // Writing nodes
    for(auto& i : root.geometry->nodes)
    {
         this->fstream_VTK << i[0] << " " << i[1] << " " << i[2] << std::endl;
    }

    // Calc cell list size
    int  vtk_cell_list_size = 0;
    for(auto& i_face : root.geometry->BC_face_global_indexes)
    {
        auto& face = root.geometry->faces[i_face];
        if(face.eface_type == enums::Eface_type::tri) vtk_cell_list_size += 4;
        if(face.eface_type == enums::Eface_type::quad) vtk_cell_list_size += 5;
    }


    // Writing BC faces
    this->fstream_VTK<<"CELLS "<<(root.geometry->n_BC_face)
                    <<" "<<vtk_cell_list_size<< std::endl;
    for(auto& i_face : root.geometry->BC_face_global_indexes)
    {
        auto face = root.geometry->faces[i_face];
        if(face.eface_type == enums::Eface_type::tri)
        {
            this->fstream_VTK<<"3 "<<face.i_nodes[0] <<" "<<face.i_nodes[1]<<" ";
            this->fstream_VTK<<face.i_nodes[2] <<std::endl;
        }if(face.eface_type == enums::Eface_type::quad)
        {
            this->fstream_VTK<<"4 "<<face.i_nodes[0] <<" "<<face.i_nodes[1]<<" ";
            this->fstream_VTK<<face.i_nodes[2]<< " " <<face.i_nodes[3] <<std::endl;
        }
    }

    this->fstream_VTK<<"CELL_TYPES "<<(root.geometry->n_BC_face)<<std::endl;

    // Bc face type
    for(auto& i_face : root.geometry->BC_face_global_indexes)
    {
        auto& face = root.geometry->faces[i_face];
        if(face.eface_type == enums::Eface_type::tri) this->fstream_VTK<<"5"<<std::endl;
        else if(face.eface_type == enums::Eface_type::quad) this->fstream_VTK<<"9"<<std::endl;
    }


    this->fstream_VTK<<"CELL_DATA "<<(root.geometry->n_BC_face)<<std::endl;

    this->write_boundary_VTK_data();

    this->fstream_VTK.close();
    std::cout<<"// - Writing Boundary Paraview File Ended - //"<<std::endl;
}
