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

#include "croot.h"

#include <chrono>
#ifdef COMPILE_WITH_OPENMP
    #include <omp.h>
#endif

Croot root;

Croot::Croot()
{

}

void Croot::run_config(Cconfiguration inp_config)
{
    configuration = std::shared_ptr<Cconfiguration> (new Cconfiguration(inp_config));

    #ifdef COMPILE_WITH_OPENMP
    omp_set_dynamic( 0 );
    omp_set_num_threads(this->configuration->openmp_number_of_threads);
    #endif

    create_geometry();

    initialize_solver();

    start_solver();
}

void Croot::run_config(std::string filename)
{
    configuration = std::shared_ptr<Cconfiguration> (new Cconfiguration());
    configuration->read_config_file(filename);

    #ifdef COMPILE_WITH_OPENMP
    omp_set_dynamic( 0 );
    omp_set_num_threads(this->configuration->openmp_number_of_threads);
    #endif

    create_geometry();

    initialize_solver();

    start_solver();
}


void Croot::create_geometry()
{
    // If mesh will created from mesh file
    if(configuration->emesh_source == enums::Emesh_source::from_file)
    {
        root.geometry = std::shared_ptr<Cgeometry>
                (new Cgeometry(Read_SU2_mesh(configuration->mesh_filename)));
    }
    else if (configuration->emesh_source == enums::Emesh_source::generate_1D)
    {
        root.geometry = std::shared_ptr<Cgeometry>
                (new Cgeometry(generate_1d_mesh(
                configuration->geo_parameters_1D[0],
                configuration->geo_parameters_1D[1],
                configuration->geo_parameters_1D[2],
                configuration->geo_parameters_1D[3])));
    }
}


void Croot::initialize_solver()
{
    // Config sanity check
    if(configuration->etimestepping_type == enums::Etimestepping_type::Steady || configuration->eproblem_type == enums::Eproblem_type::Heat_conduction)
        configuration->max_inner_iteration = 1;

    //Initialize FVM
    FVM = std::shared_ptr<CFVM> (new CFVM());
    FVM->create_FVM_from_geometry();
    FVM->define_matrix_structure(Matrix_solver);

    // If problem is heat conduction problem Heat solver is initialised
    if(configuration->eproblem_type == enums::Eproblem_type::Heat_conduction)
    {
        Heat_solver = std::shared_ptr<CHeat_solver> (new CHeat_solver());
        Heat_solver->initialize(FVM,Matrix_solver);
    }

    // If problem is fluid flow problem Flow solver is initialised
    else if(configuration->eproblem_type == enums::Eproblem_type::Fluid_flow)
    {
        Flow_solver = std::shared_ptr<CFlow_solver> (new CFlow_solver());
        Fluid_properties = std::shared_ptr<Cfluid_properties> (new Cfluid_properties(get_AIR_properties()));
        Flow_solver->initialize(FVM,Matrix_solver,Fluid_properties);
    }



}

void Croot::start_solver()
{
    post_processing = std::shared_ptr<Cpost_processing> (new Cpost_processing());

    std::chrono::steady_clock::time_point begin,end;

    for(i_iteration = 0 ; i_iteration < configuration->max_iteration ; i_iteration++)
    {
        for(i_inner_iter = 0 ; i_inner_iter < configuration->max_inner_iteration ; i_inner_iter++)
        {
            begin = std::chrono::steady_clock::now();

            // If problem is heat conduction
            if(configuration->eproblem_type == enums::Eproblem_type::Heat_conduction)
            {
                Heat_solver->iterate(FVM,Matrix_solver);
            }

            // If problem is fluid flow
            else if(configuration->eproblem_type == enums::Eproblem_type::Fluid_flow)
            {
                if(configuration->etimestepping_type == enums::Etimestepping_type::Transient)
                {
                    if(configuration->eflow_type_visc == enums::Eflow_type_visc::Inviscid)
                    {
                        if(configuration->econv_flux_type == enums::Econv_flux_type::Blended_central)
                        {
                            Flow_solver->iterate<enums::Etimestepping_type::Transient,
                                    enums::Eflow_type_visc::Inviscid ,
                                    enums::Econv_flux_type::Blended_central >
                                (FVM,Matrix_solver,Fluid_properties);
                        }
                        else
                        {
                            Flow_solver->iterate<enums::Etimestepping_type::Transient,
                                    enums::Eflow_type_visc::Inviscid ,
                                    enums::Econv_flux_type::First_order_upwind>
                                (FVM,Matrix_solver,Fluid_properties);
                        }
                    }
                    else
                    {
                        if(configuration->econv_flux_type == enums::Econv_flux_type::Blended_central)
                        {
                            Flow_solver->iterate<enums::Etimestepping_type::Transient,
                                    enums::Eflow_type_visc::Viscous ,
                                    enums::Econv_flux_type::Blended_central >
                                (FVM,Matrix_solver,Fluid_properties);
                        }
                        else
                        {
                            Flow_solver->iterate<enums::Etimestepping_type::Transient,
                                    enums::Eflow_type_visc::Viscous ,
                                    enums::Econv_flux_type::First_order_upwind>
                                (FVM,Matrix_solver,Fluid_properties);
                        }
                    }
                }
                else
                {
                    if(configuration->eflow_type_visc == enums::Eflow_type_visc::Inviscid)
                    {
                        if(configuration->econv_flux_type == enums::Econv_flux_type::Blended_central)
                        {
                            Flow_solver->iterate<enums::Etimestepping_type::Steady,
                                    enums::Eflow_type_visc::Inviscid ,
                                    enums::Econv_flux_type::Blended_central >
                                (FVM,Matrix_solver,Fluid_properties);
                        }
                        else
                        {
                            Flow_solver->iterate<enums::Etimestepping_type::Steady,
                                    enums::Eflow_type_visc::Inviscid ,
                                    enums::Econv_flux_type::First_order_upwind>
                                (FVM,Matrix_solver,Fluid_properties);
                        }
                    }
                    else
                    {
                        if(configuration->econv_flux_type == enums::Econv_flux_type::Blended_central)
                        {
                            Flow_solver->iterate<enums::Etimestepping_type::Steady,
                                    enums::Eflow_type_visc::Viscous ,
                                    enums::Econv_flux_type::Blended_central >
                                (FVM,Matrix_solver,Fluid_properties);
                        }
                        else
                        {
                            Flow_solver->iterate<enums::Etimestepping_type::Steady,
                                    enums::Eflow_type_visc::Viscous ,
                                    enums::Econv_flux_type::First_order_upwind>
                                (FVM,Matrix_solver,Fluid_properties);
                        }
                    }
                }
            }
            root.post_processing->monitoring();

            end = std::chrono::steady_clock::now();
            this->total_elapsed_iteration_time += std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
        }


        // If problem is fluid flow
        if(configuration->eproblem_type == enums::Eproblem_type::Fluid_flow)
        {
            Flow_solver->update_timestep();
        }

        if(configuration->VTK_file_frequency != 0 && i_iteration != 0 && i_iteration % configuration->VTK_file_frequency == 0)
        {
            std::ostringstream filenmstream;
            filenmstream << "solution_" << i_iteration << ".vtk";
            post_processing->write_VTK_file(filenmstream.str());
            std::ostringstream bcfilenmstream;
            bcfilenmstream << "solution_BC_" << i_iteration << ".vtk";
            post_processing->write_boundary_VTK_file(bcfilenmstream.str());
        }
    }

    post_processing->write_VTK_file("solution.vtk");
    post_processing->write_boundary_VTK_file("solution_BC.vtk");

    long int av_elapsed_iteration_time = this->total_elapsed_iteration_time
            / (configuration->max_iteration * configuration->max_inner_iteration);

    long int av_elapsed_matrix_solution_time = this->total_elapsed_matrix_solution_time
            / (configuration->max_iteration * configuration->max_inner_iteration);


    std::cout<<"// - Average elapsed iteration time : " << av_elapsed_iteration_time <<" microseconds - //"<<std::endl;
    std::cout<<"// - Average elapsed matrix solution time : " << av_elapsed_matrix_solution_time <<" microseconds - //"<<std::endl;

}
