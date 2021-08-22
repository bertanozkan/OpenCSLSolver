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

#ifndef CCONFIGURATION_H
#define CCONFIGURATION_H

#include "enums.h"
#include "croot.h"

class Cconfiguration
{
public:
    Cconfiguration();

    int n_BC;

    class Cinitial_condition;
    std::shared_ptr<Cinitial_condition> inital_condition;

    class Cboundary_condition;
    std::vector<Cboundary_condition> boundary_conditions;

    // Geometry variables
    double geo_parameters_1D[4];
    enums::Emesh_source emesh_source;
    std::string mesh_filename;

    // Common FVM variables
    enums::Eproblem_type eproblem_type;
    enums::Econv_flux_type econv_flux_type = enums::Econv_flux_type::First_order_upwind;
    double global_delta_t;

    // Flow variables
    enums::Eflow_type_visc eflow_type_visc = enums::Eflow_type_visc::Viscous;

    // Heat conduction solver variables
    double thermal_diffusivity;

    // Simulation controls
    enums::Etimestepping_type etimestepping_type;
    int openmp_number_of_threads = 1;
    int max_iteration;
    int max_inner_iteration;
    double implicit_underrelax = 0.4;
    double implicit_underrelax_max = 0.8;
    int underrelax_ramp_iteration = 100;
    int VTK_file_frequency = 0;

    // Matrix solver options
    enums::Ematrix_solver_type ematrix_solver_type;

    void set_BC_number(int n_BC);

    void read_config_file(std::string filename);

};


class Cconfiguration::Cinitial_condition
{
public:
    Cinitial_condition();

    double Temperature;
    double Pressure;
    Svector Velocity;
    enums::Einitial_condition_type einitial_condition_type = enums::Einitial_condition_type::Uniform;
    double Pressure2;
    double Temperature2;



};


class Cconfiguration::Cboundary_condition
{
public:
    Cboundary_condition();

    std::string BC_tag;
    double Temperature;
    double Pressure;
    Svector Velocity;
    enums::Ethermal_BC_type ethermal_BC_type = enums::Ethermal_BC_type::adiabatic;
    enums::Eflow_BC_type eflow_BC_type;
    enums::Ewall_BC_type ewall_BC_type = enums::Ewall_BC_type::no_slip;

};


#endif // CCONFIGURATION_H
