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

#include "cconfiguration.h"

#include <iostream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <regex>

Cconfiguration::Cconfiguration()
{
    inital_condition = std::shared_ptr<Cinitial_condition>(new Cinitial_condition());
}

Cconfiguration::Cinitial_condition::Cinitial_condition()
{

}

Cconfiguration::Cboundary_condition::Cboundary_condition()
{

}


void Cconfiguration::set_BC_number(int n_BC)
{
    this->n_BC = n_BC;
    this->boundary_conditions.resize(n_BC);
}

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

void Cconfiguration::read_config_file(std::string filename)
{
    std::ifstream config_file;
    std::smatch match;

    std::cout << "// - Reading configuration from file " << filename << " - //"  << std::endl;
    config_file.open(filename);
    if (!config_file) {
        std::cout << " ERROR! Unable to open configuration file " << std::endl;
        exit(1);
    }

    std::string line;
    while( std::getline(config_file, line) )
    {
        trim(line);
        if(line[0] != '#')
        {
            if(line.size() > 0)
            {

                std::istringstream iss(line);
                std::string result;
                if(getline(iss,result, ' ' ))
                {

                    //---------------------------------------------MESH----------------------------------------------------------

                    // Name of the mesh file.
                    if(result == "MESH_FILE_NAME")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                           emesh_source = enums::Emesh_source::from_file;
                           mesh_filename = token;
                           std::cout << "// - Mesh from file: " << token << " is chosen - //" <<  std::endl;
                        }
                    }


                    // Command to generate 1D mesh
                    // Enter Lx, Ly, Lz, dx with spaces between.
                    else if(result == "GENERATE_1D_MESH")
                    {
                        std::string token;
                        getline(iss,token,' ');
                        trim(token);
                        emesh_source = enums::Emesh_source::generate_1D;
                        geo_parameters_1D[0] = ::atof(token.c_str());

                        getline(iss,token,' ');
                        trim(token);
                        geo_parameters_1D[1] = ::atof(token.c_str());

                        getline(iss,token, ' ');
                        trim(token);
                        geo_parameters_1D[2] = ::atof(token.c_str());

                        getline(iss,token);
                        trim(token);
                        geo_parameters_1D[3] = ::atof(token.c_str());

                        std::cout << "// - 1D Mesh will be generated with Lx: " << geo_parameters_1D[0] <<
                                     " Ly: " << geo_parameters_1D[1] <<
                                     " Lz " << geo_parameters_1D[2] <<
                                     " dx " << geo_parameters_1D[3] << " - //" << std::endl;
                    }

                    //--------------------------------------------EQUATIONS---------------------------------------------------------

                    // Selection of problem type.
                    // "HEAT_CONDUCTION" for solving heat problems
                    // "FLUID_FLOW" for solving fluid flow problems
                    else if(result == "PROBLEM_TYPE")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                            trim(token);
                        } if(token == "HEAT_CONDUCTION")
                        {
                            eproblem_type = enums::Eproblem_type::Heat_conduction;
                            std::cout<<"// - HEAT_CONDUCTION is chosen - //"<< std::endl;
                        } else if(token == "FLUID_FLOW")
                        {
                            eproblem_type = enums::Eproblem_type::Fluid_flow;
                            std::cout<<"// - FLUID_FLOW is chosen - //"<< std::endl;
                        }
                    }


                    // Selection of flow type.
                    // "VISCOUS" for solving viscous flows
                    // "INVISCID" for solving inviscid flows
                    else if(result == "FLOW_TYPE")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                            trim(token);
                        } if(token == "VISCOUS")
                        {
                            eflow_type_visc = enums::Eflow_type_visc::Viscous;
                            std::cout<<"// - VISCOUS is chosen - //"<< std::endl;
                        } else if(token == "INVISCID")
                        {
                            eflow_type_visc = enums::Eflow_type_visc::Inviscid;
                            std::cout<<"// - INVISCID is chosen - //"<< std::endl;
                        }
                    }

                    // Selection of timestepping type.
                    // "STEADY" for steady solution
                    // "TRANSIENT" for transient solution
                    else if(result == "TIMESTEPPING_TYPE")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                            trim(token);
                        } if(token == "STEADY")
                        {
                            etimestepping_type = enums::Etimestepping_type::Steady;
                            std::cout<<"// - STEADY is chosen - //"<< std::endl;
                        } else if(token == "TRANSIENT")
                        {
                            etimestepping_type = enums::Etimestepping_type::Transient;
                            std::cout<<"// - TRANSIENT is chosen - //"<< std::endl;
                        }
                    }


                    //-------------------------------------------NUMERICS---------------------------------------------------------


                    // Selection of matrix solver type.
                    // "CG" for Conjugate gradient matrix solver
                    // "ILU_BICGSTAB" for BICGSTAB solver with ILU preconditioner
                    else if(result == "MATRIX_SOLVER_TYPE")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                            trim(token);
                        } if(token == "CG")
                        {
                            ematrix_solver_type = enums::Ematrix_solver_type::CG;
                            std::cout<<"// - CG is chosen - //"<< std::endl;
                        } else if(token == "ILU_BICGSTAB")
                        {
                            ematrix_solver_type = enums::Ematrix_solver_type::ILU_BICSTAB;
                            std::cout<<"// - ILU_BICGSTAB is chosen - //"<< std::endl;
                        }
                    }


                    // Selection of convective flux type.
                    // "FIRST_ORDER_UPWIND" for first order upwind convective flux
                    // "BLENDED_CENTRAL" for blended central and upwind convective flux
                    else if(result == "CONV_FLUX_TYPE")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                            trim(token);
                        } if(token == "FIRST_ORDER_UPWIND")
                        {
                            econv_flux_type = enums::Econv_flux_type::First_order_upwind;
                            std::cout<<"// - FIRST_ORDER_UPWIND is chosen - //"<< std::endl;
                        } else if(token == "BLENDED_CENTRAL")
                        {
                            econv_flux_type = enums::Econv_flux_type::Blended_central;
                            std::cout<<"// - BLENDED_CENTRAL is chosen - //"<< std::endl;
                        }
                    }

                    // Timestep of unsteady simulation.
                    else if(result == "TIME_STEP")
                    {
                       std::string token;
                       while(getline(iss,token))
                       {
                           trim(token);
                           global_delta_t = ::atof(token.c_str());
                       }
                       std::cout << "// - Transient solution timestep: " << global_delta_t << " - //" << std::endl;
                    }

                    // Ramp parameters for implicit underrelaxation parameter
                    // Ramp start, ramp end, ramp end iteration number
                    else if(result == "UNDERRELAX_RAMP")
                    {
                        std::string token;
                        getline(iss,token,' ');
                        trim(token);
                        implicit_underrelax = ::atof(token.c_str());

                        getline(iss,token, ' ');
                        trim(token);
                        implicit_underrelax_max = ::atof(token.c_str());

                        getline(iss,token);
                        trim(token);
                        underrelax_ramp_iteration = ::atof(token.c_str());
                        std::cout << "// - Underrelax ramp start: " << implicit_underrelax <<
                                     " underrelax end: " << implicit_underrelax_max <<
                                     " underrelax end iteration: " << underrelax_ramp_iteration<< " - //" << std::endl;
                    }


                    //-------------------------------------------INITIAL CONDITIONS---------------------------------------------------------


                    // Initial condition of velocity vector in units of [m/s].
                    // Enter "x" "y" "z" components with spaces between.
                    else if(result == "VELOCITY_IC")
                    {
                        std::string token;
                        getline(iss,token,' ');
                        trim(token);
                        inital_condition->Velocity.val[0] = ::atof(token.c_str());

                        getline(iss,token, ' ');
                        trim(token);
                        inital_condition->Velocity.val[1] = ::atof(token.c_str());

                        getline(iss,token);
                        trim(token);
                        inital_condition->Velocity.val[2] = ::atof(token.c_str());

                        std::cout << "// - Velocity Initial Conditions: " << inital_condition->Velocity.val[0] <<
                                     " ; " << inital_condition->Velocity.val[1] <<
                                     " ; " << inital_condition->Velocity.val[2] << " - //" << std::endl;
                    }


                    // Initial condition of Pressure in unit of [Pa].
                    else if(result == "PRESSURE_IC")
                    {
                       std::string token;
                       while(getline(iss,token))
                       {
                           trim(token);
                           inital_condition->Pressure = ::atof(token.c_str());
                       }
                       std::cout << "// - Pressure Initial Condition: " << inital_condition->Pressure << " - //" << std::endl;
                    }


                    // Initial condition of Temperature in unit of [K].
                    else if(result == "TEMPERATURE_IC")
                    {
                       std::string token;
                       while(getline(iss,token))
                       {
                           trim(token);
                           inital_condition->Temperature = ::atof(token.c_str());
                       }
                       std::cout << "// - Temperature Initial Condition: " << inital_condition->Temperature << " - //" << std::endl;
                    }



                    // Initial condition of right side in 1D mesh, Pressure in unit of [Pa].
                    // Used in shocktube simulation
                    else if(result == "PRESSURE_IC_RIGHT")
                    {
                       std::string token;
                       while(getline(iss,token))
                       {
                           trim(token);
                           inital_condition->Pressure2 = ::atof(token.c_str());
                       }
                       inital_condition->einitial_condition_type = enums::Einitial_condition_type::Shock_tube;
                       std::cout << "// - Pressure Initial Condition of right side: " << inital_condition->Pressure2 << " - //" << std::endl;
                    }


                    // Initial condition of right side in 1D mesh, Temperature in unit of [K].
                    // Used in shocktube simulation
                    else if(result == "TEMPERATURE_IC_RIGHT")
                    {
                       std::string token;
                       while(getline(iss,token))
                       {
                           trim(token);
                           inital_condition->Temperature2 = ::atof(token.c_str());
                       }
                       inital_condition->einitial_condition_type = enums::Einitial_condition_type::Shock_tube;
                       std::cout << "// - Temperature Initial Condition of right side: " << inital_condition->Temperature2 << " - //" << std::endl;
                    }


                    //-------------------------------------------BOUNDARY CONDITIONS---------------------------------------------------------


                    // Command to set Empty boundary condition
                    // Enter only BC index in mesh.
                    else if(result == "SET_EMPTY_BC")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                            trim(token);
                            size_t dum_index = ::atof(token.c_str());
                            if(boundary_conditions.size() < dum_index+1 )
                            {
                                boundary_conditions.resize(dum_index+1);
                            }
                            boundary_conditions[dum_index].eflow_BC_type = enums::Eflow_BC_type::empty;
                            std::cout << "// - BC index " << dum_index << " is set to be Empty BC - //" << std::endl;
                        }
                    }


                    // Command to set velocity inlet boundary condition
                    // Enter bc index, "x", "y", "z" components and temperature with spaces between.
                    else if(result == "SET_VEL_INLET_BC")
                    {
                        std::string token;
                        getline(iss,token,' ');
                        trim(token);
                        size_t dum_index = ::atof(token.c_str());
                        if(boundary_conditions.size() < dum_index+1 )
                        {
                            boundary_conditions.resize(dum_index+1);
                        }

                        boundary_conditions[dum_index].eflow_BC_type = enums::Eflow_BC_type::velocity_inlet;

                        getline(iss,token,' ');
                        trim(token);
                        boundary_conditions[dum_index].Velocity.val[0] = ::atof(token.c_str());

                        getline(iss,token, ' ');
                        trim(token);
                        boundary_conditions[dum_index].Velocity.val[1] = ::atof(token.c_str());

                        getline(iss,token, ' ');
                        trim(token);
                        boundary_conditions[dum_index].Velocity.val[2] = ::atof(token.c_str());

                        getline(iss,token);
                        trim(token);
                        boundary_conditions[dum_index].Temperature = ::atof(token.c_str());

                        std::cout << "// - BC index " << dum_index << " is set to be Velocity inlet BC with velocity: " <<
                                     boundary_conditions[dum_index].Velocity.val[0] <<
                                     " ; " << boundary_conditions[dum_index].Velocity.val[1] <<
                                     " ; " << boundary_conditions[dum_index].Velocity.val[2] <<
                                     " and temperature: " << boundary_conditions[dum_index].Temperature << " - //" << std::endl;
                    }


                    // Command to set supersonic velocity inlet boundary condition
                    // Enter bc index, "x", "y", "z" components, pressure and temperature with spaces between.
                    else if(result == "SET_SUPERSONIC_VEL_INLET_BC")
                    {
                        std::string token;
                        getline(iss,token,' ');
                        trim(token);
                        size_t dum_index = ::atof(token.c_str());
                        if(boundary_conditions.size() < dum_index+1 )
                        {
                            boundary_conditions.resize(dum_index+1);
                        }

                        boundary_conditions[dum_index].eflow_BC_type = enums::Eflow_BC_type::supersonic_velocity_inlet;

                        getline(iss,token,' ');
                        trim(token);
                        boundary_conditions[dum_index].Velocity.val[0] = ::atof(token.c_str());

                        getline(iss,token, ' ');
                        trim(token);
                        boundary_conditions[dum_index].Velocity.val[1] = ::atof(token.c_str());

                        getline(iss,token, ' ');
                        trim(token);
                        boundary_conditions[dum_index].Velocity.val[2] = ::atof(token.c_str());

                        getline(iss,token, ' ');
                        trim(token);
                        boundary_conditions[dum_index].Pressure = ::atof(token.c_str());

                        getline(iss,token);
                        trim(token);
                        boundary_conditions[dum_index].Temperature = ::atof(token.c_str());

                        std::cout << "// - BC index " << dum_index << " is set to be Supersonic Velocity inlet BC with velocity: " <<
                                     boundary_conditions[dum_index].Velocity.val[0] <<
                                     " ; " << boundary_conditions[dum_index].Velocity.val[1] <<
                                     " ; " << boundary_conditions[dum_index].Velocity.val[2] <<
                                     " pressure: " << boundary_conditions[dum_index].Pressure <<
                                     " and temperature: " << boundary_conditions[dum_index].Temperature << " - //" << std::endl;
                    }


                    // Command to set pressure outlet boundary condition
                    // Enter bc index and pressure with spaces between.
                    else if(result == "SET_PRES_OUTLET_BC")
                    {
                        std::string token;
                        getline(iss,token,' ');
                        trim(token);
                        size_t dum_index = ::atof(token.c_str());
                        if(boundary_conditions.size() < dum_index+1 )
                        {
                            boundary_conditions.resize(dum_index+1);
                        }

                        boundary_conditions[dum_index].eflow_BC_type = enums::Eflow_BC_type::pressure_outlet;

                        getline(iss,token);
                        trim(token);
                        boundary_conditions[dum_index].Pressure = ::atof(token.c_str());

                        std::cout << "// - BC index " << dum_index << " is set to be Pressure outlet BC with pressure: " <<
                                     boundary_conditions[dum_index].Pressure << " - //" << std::endl;
                    }


                    // Command to set slip wall boundary condition
                    // Enter only BC index in mesh.
                    else if(result == "SET_SLIP_WALL_BC")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                            trim(token);
                            size_t dum_index = ::atof(token.c_str());
                            if(boundary_conditions.size() < dum_index+1 )
                            {
                                boundary_conditions.resize(dum_index+1);
                            }
                            boundary_conditions[dum_index].eflow_BC_type = enums::Eflow_BC_type::wall;
                            boundary_conditions[dum_index].ewall_BC_type = enums::Ewall_BC_type::slip;
                            std::cout << "// - BC index " << dum_index << " is set to be slip wall BC - //" << std::endl;
                        }
                    }

                    // Command to set no-slip wall boundary condition
                    // Enter only BC index in mesh.
                    else if(result == "SET_NO_SLIP_WALL_BC")
                    {
                        std::string token;
                        while(getline(iss,token))
                        {
                            trim(token);
                            size_t dum_index = ::atof(token.c_str());
                            if(boundary_conditions.size() < dum_index+1 )
                            {
                                boundary_conditions.resize(dum_index+1);
                            }
                            boundary_conditions[dum_index].eflow_BC_type = enums::Eflow_BC_type::wall;
                            boundary_conditions[dum_index].ewall_BC_type = enums::Ewall_BC_type::no_slip;
                            std::cout << "// - BC index " << dum_index << " is set to be no-slip wall BC - //" << std::endl;
                        }
                    }


                    //-------------------------------------------SIMULATION CONTROLS---------------------------------------------------------

                    // Number of  parallel threads to compute on.
                    else if(result == "N_P_THREADS")
                    {
                       std::string token;
                       while(getline(iss,token))
                       {
                          trim(token);
                          openmp_number_of_threads = ::atof(token.c_str());
                       }
                       std::cout << "// - Number of parallel threads: " << openmp_number_of_threads << " - //" << std::endl;
                    }

                    // Maximum iterations to stop the simulation.
                    else if(result == "MAX_ITERATIONS")
                    {
                       std::string token;
                       while(getline(iss,token))
                       {
                           trim(token);
                           max_iteration = ::atof(token.c_str());
                       }
                        std::cout << "// - Maximum number of iterations: " << max_iteration << " - //" << std::endl;
                    }

                    // Maximum number of inner iterations.
                    else if(result == "MAX_INNER_ITER")
                    {
                       std::string token;
                       while(getline(iss,token));
                       {
                          trim(token);
                          max_inner_iteration = ::atof(token.c_str());
                          std::cout<<"// - Maximum number of inner iterations: "<<max_inner_iteration<<" - //"<< std::endl;
                       }
                    }

                    // Frequency of writing Paraview files.
                    else if(result == "VTK_WRITING_FREQUENCY")
                    {
                       std::string token;
                       while(getline(iss,token));
                       {
                           trim(token);
                           VTK_file_frequency = ::atof(token.c_str());
                       }
                       std::cout<<"// - VTK file writing frequency: "<<VTK_file_frequency<<" - //"<< std::endl;
                    }

                }
            }
        }
    }
    config_file.close();
}
