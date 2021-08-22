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

#ifndef CGEOMETRY_SUPPLY_H
#define CGEOMETRY_SUPPLY_H

#include "croot.h"

static inline Cgeometry generate_1d_mesh(double lx, double ly, double lz, int nx)
{
    std::cout << "// - Generating 1D Mesh - //" << std::endl;
    Cgeometry geo;
    double dx = lx/double(nx);
    geo.n_boundary = 3;
    geo.nodes.reserve((nx+1)*4);
    geo.cells.reserve(nx);
    // Create nodes
    for(int i = 0 ; i < (nx+1) ; i++)
    {
        geo.nodes.push_back(Svector((double(i)*dx),0.0,lz));
        geo.nodes.push_back(Svector((double(i)*dx),0.0,0.0));
        geo.nodes.push_back(Svector((double(i)*dx),ly,lz));
        geo.nodes.push_back(Svector((double(i)*dx),ly,0.0));
    }
    //Create cells
    for(int i = 0 ; i < (nx) ; i++)
    {
        geo.cells.push_back(Cgeometry::Cgeo_cell(enums::Ecell_type::hexahedron));
        for(int ii = 0 ; ii < 8 ; ii++)
            geo.cells[i].i_nodes[ii] = 4*i+ii;
    }
    // Calculate connectivity
    geo.calc_cell_face_connectivity();

    // Create boundary faces
    geo.n_BC_face = 0;
    geo.add_BC_face(0,4,0);
    geo.n_BC_face += 1;
    for(int i = 0 ; i < (nx) ; i++)
    {
        geo.add_BC_face(i,0,1);
        geo.add_BC_face(i,1,1);
        geo.add_BC_face(i,2,1);
        geo.add_BC_face(i,3,1);
        geo.n_BC_face += 4;
    }
    geo.add_BC_face((nx-1),5,2);
    geo.n_BC_face += 1;

    geo.calc_geometric_values();
    geo.extract_BC_faces();

    for(size_t i = 0 ; i < geo.node_cell_indexes.size() ; i++)
        geo.node_cell_indexes[i].clear();

    geo.node_cell_indexes.clear();

    std::cout << "// - Generating 1D Mesh Ended - //" << std::endl;
    return geo;
}

namespace unit_tests {

static inline void test_geo_calculations()
{
    std::cout << "Unit test for basic geometry calculations: ";

    bool print_values = false;
    if(print_values) std::cout << std::endl;
    bool passed = true;
    Cgeometry geo;
    geo.nodes.reserve(8);
    geo.nodes.push_back(Svector(0.0,0.0,0.0));
    geo.nodes.push_back(Svector(1.0,0.0,0.0));
    geo.nodes.push_back(Svector(0.0,1.0,0.0));
    geo.nodes.push_back(Svector(1.0,1.0,0.0));
    geo.nodes.push_back(Svector(0.0,0.0,1.0));
    geo.nodes.push_back(Svector(1.0,0.0,1.0));
    geo.nodes.push_back(Svector(0.0,1.0,1.0));
    geo.nodes.push_back(Svector(1.0,1.0,1.0));

    geo.faces.push_back(Cgeometry::Cgeo_face(enums::Eface_type::tri));
    geo.faces[0].i_nodes[0] = 0;
    geo.faces[0].i_nodes[1] = 1;
    geo.faces[0].i_nodes[2] = 2;

    geo.calc_face_area_and_normal(geo.faces[0]);
    if(abs(geo.faces[0].area - 0.5) >= 1.0e-15) passed = false;
    if(abs(geo.faces[0].center[0] - (1.0/3.0)) >= 1.0e-15) passed = false;
    if(abs(geo.faces[0].center[1] - (1.0/3.0)) >= 1.0e-15) passed = false;
    if(abs(geo.faces[0].center[2] - 0.0) >= 1.0e-15) passed = false;
    if(print_values){
        std::cout << "Face 0 area: " << geo.faces[0].area << std::endl;
        std::cout << "Face 0 centroid: x: " << geo.faces[0].center[0] <<
                                 " y: " << geo.faces[0].center[1] <<
                                 " z: " << geo.faces[0].center[2] << std::endl;
    }

    geo.faces.push_back(Cgeometry::Cgeo_face(enums::Eface_type::quad));
    geo.faces[1].i_nodes[0] = 0;
    geo.faces[1].i_nodes[1] = 1;
    geo.faces[1].i_nodes[2] = 3;
    geo.faces[1].i_nodes[3] = 2;

    geo.calc_face_area_and_normal(geo.faces[1]);
    if(abs(geo.faces[1].area - 1.0) >= 1.0e-15) passed = false;
    if(abs(geo.faces[1].center[0] - 0.5) >= 1.0e-15) passed = false;
    if(abs(geo.faces[1].center[1] - 0.5) >= 1.0e-15) passed = false;
    if(abs(geo.faces[1].center[2] - 0.0) >= 1.0e-15) passed = false;
    if(print_values){
        std::cout << "Face 1 area: " << geo.faces[1].area << std::endl;
        std::cout << "Face 1 centroid: x: " << geo.faces[1].center[0] <<
                                 " y: " << geo.faces[1].center[1] <<
                                 " z: " << geo.faces[1].center[2] << std::endl;
    }

    geo.cells.push_back(Cgeometry::Cgeo_cell(enums::Ecell_type::tetra));
    geo.cells[0].i_nodes[0] = 0;
    geo.cells[0].i_nodes[1] = 1;
    geo.cells[0].i_nodes[2] = 2;
    geo.cells[0].i_nodes[3] = 4;

    geo.calc_cell_volume_and_center(geo.cells[0]);
    if(abs(geo.cells[0].volume - (1.0/6.0)) >= 1.0e-15) passed = false;
    if(abs(geo.cells[0].center[0] - 0.25) >= 1.0e-15) passed = false;
    if(abs(geo.cells[0].center[1] - 0.25) >= 1.0e-15) passed = false;
    if(abs(geo.cells[0].center[2] - 0.25) >= 1.0e-15) passed = false;
    if(print_values){
        std::cout << "Cell 0 Volume: " << geo.cells[0].volume << std::endl;
        std::cout << "Cell 0 centroid: x: " << geo.cells[0].center[0] <<
                                 " y: " << geo.cells[0].center[1] <<
                                 " z: " << geo.cells[0].center[2] << std::endl;
    }

    geo.cells.push_back(Cgeometry::Cgeo_cell(enums::Ecell_type::wedge));
    geo.cells[1].i_nodes[0] = 0;
    geo.cells[1].i_nodes[1] = 1;
    geo.cells[1].i_nodes[2] = 2;
    geo.cells[1].i_nodes[3] = 4;
    geo.cells[1].i_nodes[4] = 5;
    geo.cells[1].i_nodes[5] = 6;

    geo.calc_cell_volume_and_center(geo.cells[1]);
    if(abs(geo.cells[1].volume - (1.0/2.0)) >= 1.0e-15) passed = false;
    if(abs(geo.cells[1].center[0] - (1.0/3.0)) >= 1.0e-15) passed = false;
    if(abs(geo.cells[1].center[1] - (1.0/3.0)) >= 1.0e-15) passed = false;
    if(abs(geo.cells[1].center[2] - (1.0/2.0)) >= 1.0e-15) passed = false;
    if(print_values){
        std::cout << "Cell 1 Volume: " << geo.cells[1].volume << std::endl;
        std::cout << "Cell 1 centroid: x: " << geo.cells[1].center[0] <<
                                 " y: " << geo.cells[1].center[1] <<
                                 " z: " << geo.cells[1].center[2] << std::endl;
    }

    geo.cells.push_back(Cgeometry::Cgeo_cell(enums::Ecell_type::hexahedron));
    geo.cells[2].i_nodes[0] = 0;
    geo.cells[2].i_nodes[1] = 1;
    geo.cells[2].i_nodes[2] = 2;
    geo.cells[2].i_nodes[3] = 3;
    geo.cells[2].i_nodes[4] = 4;
    geo.cells[2].i_nodes[5] = 5;
    geo.cells[2].i_nodes[6] = 6;
    geo.cells[2].i_nodes[7] = 7;

    geo.calc_cell_volume_and_center(geo.cells[2]);
    if(abs(geo.cells[2].volume - 1.0) >= 1.0e-15) passed = false;
    if(abs(geo.cells[2].center[0] - (1.0/2.0)) >= 1.0e-15) passed = false;
    if(abs(geo.cells[2].center[1] - (1.0/2.0)) >= 1.0e-15) passed = false;
    if(abs(geo.cells[2].center[2] - (1.0/2.0)) >= 1.0e-15) passed = false;
    if(print_values){
        std::cout << "Cell 2 Volume: " << geo.cells[2].volume << std::endl;
        std::cout << "Cell 2 centroid: x: " << geo.cells[2].center[0] <<
                                 " y: " << geo.cells[2].center[1] <<
                                 " z: " << geo.cells[2].center[2] << std::endl;
    }

    if(passed) std::cout << "PASSED" << std::endl;
    else std::cout << "FAILED" << std::endl;
}

}


static inline Cgeometry Read_SU2_mesh(std::string filename)
{
    std::ifstream inFile;
    std::cout << "// - Generating Mesh from file " << filename << " - //" << std::endl;
    Cgeometry geo;

    inFile.open(filename);
    if (!inFile) {
        std::cout << " ERROR! Unable to open mesh file named " << filename << std::endl;
        exit(1); // terminate with error
    }

    int n_node;
    int n_cell;
    int bc_n_cells;
    int intdum;
    bool breakflag = false ;
    std::string line;
    std::string result;

    while(getline(inFile,line))
    {
        std::istringstream iss(line);
        if(getline(iss,result,'='))
        {
            if(result == "NELEM")
            {
                std::string token;
                while(getline(iss,token))
                {
                    n_cell = ::atof(token.c_str());
                }
                breakflag = true;
                break;
            }
        }
        if(breakflag) break;
    }

    // creating cell objects
    geo.cells.reserve(n_cell);

    for(int i_cell = 0 ;  i_cell < n_cell ; i_cell++)
    {
        inFile >> intdum;
        //cell is a tetrahedra
        if(intdum == 10)
        {
            geo.cells.push_back(Cgeometry::Cgeo_cell(enums::Ecell_type::tetra));
            inFile >> geo.cells[i_cell].i_nodes[0];
            inFile >> geo.cells[i_cell].i_nodes[1];
            inFile >> geo.cells[i_cell].i_nodes[2];
            inFile >> geo.cells[i_cell].i_nodes[3];
        }
        //cell is a wedge
        else if(intdum == 13)
        {
            geo.cells.push_back(Cgeometry::Cgeo_cell(enums::Ecell_type::wedge));
            inFile >> geo.cells[i_cell].i_nodes[0];
            inFile >> geo.cells[i_cell].i_nodes[1];
            inFile >> geo.cells[i_cell].i_nodes[2];
            inFile >> geo.cells[i_cell].i_nodes[3];
            inFile >> geo.cells[i_cell].i_nodes[4];
            inFile >> geo.cells[i_cell].i_nodes[5];
        }
        //cell is a hexahedron
        else if(intdum == 12)
        {
            geo.cells.push_back(Cgeometry::Cgeo_cell(enums::Ecell_type::hexahedron));
            inFile >> geo.cells[i_cell].i_nodes[0];
            inFile >> geo.cells[i_cell].i_nodes[1];
            inFile >> geo.cells[i_cell].i_nodes[3];
            inFile >> geo.cells[i_cell].i_nodes[2];
            inFile >> geo.cells[i_cell].i_nodes[4];
            inFile >> geo.cells[i_cell].i_nodes[5];
            inFile >> geo.cells[i_cell].i_nodes[7];
            inFile >> geo.cells[i_cell].i_nodes[6];
        }
        inFile >> intdum;
    }


    breakflag = false ;
    while(getline(inFile,line))
    {
        std::istringstream iss(line);
        if(getline(iss,result,'='))
        {
            if(result == "NPOIN")
            {
                std::string token;
                while(getline(iss,token))
                {
                    n_node = ::atof(token.c_str());
                }
                breakflag = true;
                break;
            }
        }
        if(breakflag) break;
    }


    // reading nodes
    geo.nodes.reserve(n_node);
    Svector dum;
    for(int i_node = 0 ; i_node < n_node ; i_node++)
    {
        inFile >> dum.val[0];
        inFile >> dum.val[1];
        inFile >> dum.val[2];
        geo.nodes.push_back(dum);
        inFile >> intdum;
    }

    // Calculate connectivity
    geo.calc_cell_face_connectivity();


    // Read boundaries
    breakflag = false ;
    while(getline(inFile,line))
    {
        std::istringstream iss(line);
        if(getline(iss,result,'='))
        {
            if(result == "NMARK")
            {
                std::string token;
                while(getline(iss,token))
                {
                    geo.n_boundary = ::atof(token.c_str());
                }
                breakflag = true;
                break;
            }
        }
        if(breakflag) break;
    }


    geo.BC_tags.resize(geo.n_boundary);
    geo.n_BC_face = 0;
    for(int i_bc = 0 ; i_bc < geo.n_boundary ; i_bc++)
    {
        breakflag = false ;
        while(getline(inFile,line))
        {
            std::istringstream iss(line);
            if(getline(iss,result,'='))
            {
                if(result == "MARKER_TAG")
                {
                    std::string token;
                    while(getline(iss,token))
                    {
                        size_t f = token.find(" ");
                        if(f != std::string::npos)
                        {
                            token.replace(f,std::string(" ").length(),"");
                        }
                        std::cout << "A boundary named " << token << " is found. " << std::endl;
                        geo.BC_tags[i_bc] = token;
                    }
                    breakflag = true;
                    break;
                }
            }
            if(breakflag) break;
        }

        // Read this boundaries elements
        breakflag = false ;
        while(getline(inFile,line))
        {
            std::istringstream iss(line);
            if(getline(iss,result,'='))
            {
                if(result == "MARKER_ELEMS")
                {
                    std::string token;
                    while(getline(iss,token))
                    {
                        bc_n_cells = ::atof(token.c_str());
                        geo.n_BC_face += bc_n_cells;
                    }
                    breakflag = true;
                    break;
                }
            }
            if(breakflag) break;
        }

        for (int i = 0 ; i < bc_n_cells ; i++) {
            inFile >> intdum;
            // Tri case
            if(intdum == 5)
            {
                int nodes[3];
                inFile >> nodes[0];
                inFile >> nodes[1];
                inFile >> nodes[2];
                geo.add_BC_face_tri(nodes,i_bc);
            }

            // Quad case
            else if(intdum == 9)
            {
                int nodes[4];
                inFile >> nodes[0];
                inFile >> nodes[1];
                inFile >> nodes[3];
                inFile >> nodes[2];
                geo.add_BC_face_quad(nodes,i_bc);
            }
        }

    }

    inFile.close();

    geo.calc_geometric_values();
    geo.extract_BC_faces();

    for(size_t i = 0 ; i < geo.node_cell_indexes.size() ; i++)
        geo.node_cell_indexes[i].clear();

    geo.node_cell_indexes.clear();

    std::cout << "// - Generating Mesh from file " << filename << " ended - //" << std::endl;
    return geo;

}

#endif // CGEOMETRY_SUPPLY_H
