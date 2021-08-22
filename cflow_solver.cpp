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

#include "cflow_solver.h"
#include "cfluid_properties.h"
#ifdef COMPILE_WITH_OPENMP
#include <omp.h>
#endif

CFlow_solver::CFlow_solver()
{

}



void CFlow_solver::initialize(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{

    cells.resize(FVM->cells.size());
    faces.resize(FVM->faces.size());
    BC_faces.resize(FVM->BC_faces.size());
    size_t nnz = 0;

    // Inserting initial conditions to cells
    for(auto& cell : cells)
    {
        cell.Temperature = root.configuration->inital_condition->Temperature;
        cell.Pressure = root.configuration->inital_condition->Pressure;
        cell.Velocity.val[0]= root.configuration->inital_condition->Velocity.val[0];
        cell.Velocity.val[1]= root.configuration->inital_condition->Velocity.val[1];
        cell.Velocity.val[2]= root.configuration->inital_condition->Velocity.val[2];
        cell.Density = Fluid_properties->get_Density(cell.Pressure,cell.Temperature);
        cell.Enthalpy = Fluid_properties->get_Enthalpy(cell.Temperature,cell.Velocity);

        cell.Density_last_t = cell.Density;
        cell.Pressure_last_t = cell.Pressure;
        cell.Velocity_last_t = cell.Velocity;
        cell.Enthalpy_last_t = cell.Enthalpy;
    }
    if(root.configuration->inital_condition->einitial_condition_type == enums::Einitial_condition_type::Shock_tube)
    {
        for(size_t i_cell = round(cells.size()/2)-1 ; i_cell < cells.size() ; i_cell++)
        {
            cells[i_cell].Pressure = root.configuration->inital_condition->Pressure2;
            cells[i_cell].Temperature = root.configuration->inital_condition->Temperature2;
            cells[i_cell].Density = Fluid_properties->get_Density(cells[i_cell].Pressure,cells[i_cell].Temperature);
            cells[i_cell].Enthalpy = Fluid_properties->get_Enthalpy(cells[i_cell].Temperature,cells[i_cell].Velocity);

            cells[i_cell].Density_last_t = cells[i_cell].Density;
            cells[i_cell].Pressure_last_t = cells[i_cell].Pressure;
            cells[i_cell].Enthalpy_last_t = cells[i_cell].Enthalpy;
        }
    }



    // Inserting FVM face pointers to heat_solver faces
    for(size_t i = 0 ; i < faces.size() ; i++)
    {
        faces[i].left_cell_ptr = &cells[FVM->faces[i].left_cell_index];
        faces[i].right_cell_ptr = &cells[FVM->faces[i].right_cell_index];
    }

    // Inserting FVM BC_face pointers to heat_solver BC_faces
    for(size_t i = 0 ; i < BC_faces.size() ; i++)
    {
        BC_faces[i].left_cell_ptr = &cells[FVM->BC_faces[i].left_cell_index];
    }

    // Inserting boundary conditions to faces
    auto BC_face = BC_faces.begin();
    auto FVM_BC_face = FVM->BC_faces.begin();
    for(; BC_face != BC_faces.end(); BC_face++ ,FVM_BC_face++)
    {
        int i_BC = root.geometry->faces[FVM_BC_face->global_face_index]
                .boundary_marker;

        BC_face->eflow_BC_type = root.configuration->boundary_conditions[i_BC]
                .eflow_BC_type;


        if(BC_face->eflow_BC_type == enums::Eflow_BC_type::velocity_inlet)
        {
            BC_face->Velocity = root.configuration->boundary_conditions[i_BC].Velocity;
            BC_face->Pressure = BC_face->left_cell_ptr->Pressure;
            BC_face->Temperature = root.configuration->boundary_conditions[i_BC].Temperature;

            BC_face->Density = Fluid_properties->get_Density(BC_face->Pressure,BC_face->Temperature);
            BC_face->Enthalpy = Fluid_properties->get_Enthalpy(BC_face->Temperature,BC_face->Velocity);
            BC_face->Mdot = BC_face->Density * dot(BC_face->Velocity, FVM_BC_face->area_vector);
        }

        if(BC_face->eflow_BC_type == enums::Eflow_BC_type::supersonic_velocity_inlet)
        {
            BC_face->Velocity = root.configuration->boundary_conditions[i_BC].Velocity;
            BC_face->Pressure = root.configuration->boundary_conditions[i_BC].Pressure;
            BC_face->Temperature = root.configuration->boundary_conditions[i_BC].Temperature;

            BC_face->Density = Fluid_properties->get_Density(BC_face->Pressure,BC_face->Temperature);
            BC_face->Enthalpy = Fluid_properties->get_Enthalpy(BC_face->Temperature,BC_face->Velocity);
            BC_face->Mdot = BC_face->Density * dot(BC_face->Velocity, FVM_BC_face->area_vector);
        }

        if(BC_face->eflow_BC_type == enums::Eflow_BC_type::pressure_outlet)
        {
            BC_face->Temperature = BC_face->left_cell_ptr->Temperature;
            BC_face->Pressure = root.configuration->boundary_conditions[i_BC].Pressure;
            BC_face->Out_Pressure = root.configuration->boundary_conditions[i_BC].Pressure;
            BC_face->Density = Fluid_properties->get_Density(BC_face->Pressure,BC_face->Temperature);
            BC_face->Velocity = BC_face->left_cell_ptr->Velocity;
            BC_face->Mdot = BC_face->left_cell_ptr->Density * dot(BC_face->left_cell_ptr->Velocity, FVM_BC_face->area_vector);
            BC_face->Enthalpy = Fluid_properties->get_Enthalpy(BC_face->left_cell_ptr->Temperature,BC_face->left_cell_ptr->Velocity);
            if((BC_face->Mdot /BC_face->left_cell_ptr->Density / FVM_BC_face->area) >= Fluid_properties->get_speed_of_sound(BC_face->left_cell_ptr->Temperature))
            {
                BC_face->supersonic_flag = true;
                BC_face->Pressure = BC_face->left_cell_ptr->Pressure;
                BC_face->Density = Fluid_properties->get_Density(BC_face->Pressure,BC_face->Temperature);
                BC_face->Mdot = BC_face->left_cell_ptr->Density * dot(BC_face->left_cell_ptr->Velocity, FVM_BC_face->area_vector);
            }
        }

        if(BC_face->eflow_BC_type == enums::Eflow_BC_type::wall)
        {
            BC_face->ewall_BC_type = root.configuration->boundary_conditions[i_BC]
                    .ewall_BC_type;
            BC_face->eThermal_BC_type = root.configuration->boundary_conditions[i_BC]
                    .ethermal_BC_type;

            BC_face->Pressure = BC_face->left_cell_ptr->Pressure;
        }
    }

    for(auto& cell : cells)
    {
        cell.ap = 1.0;
    }

    // Calculating face mdot
    if(root.configuration->econv_flux_type == enums::Econv_flux_type::Blended_central)
    {
        auto face = faces.begin();
        auto FVM_face = FVM->faces.begin();
        for(; face != faces.end() ; face++, FVM_face++)
            calc_face_mdot<enums::Econv_flux_type::Blended_central>(*FVM_face,*face);
    }else
    {
        auto face = faces.begin();
        auto FVM_face = FVM->faces.begin();
        for(; face != faces.end() ; face++, FVM_face++)
            calc_face_mdot<enums::Econv_flux_type::First_order_upwind>(*FVM_face,*face);
    }

    BC_face = BC_faces.begin();
    FVM_BC_face = FVM->BC_faces.begin();
    for(; BC_face != BC_faces.end(); BC_face++ ,FVM_BC_face++)
        calc_face_mdot(*FVM_BC_face, *BC_face,Fluid_properties);



    // Underrelaxation ramping
    if(root.configuration->etimestepping_type == enums::Etimestepping_type::Steady)
    {
        implicit_underrelaxation = root.configuration->implicit_underrelax;
    }else
    {
        double deltaT = root.configuration->global_delta_t;
        for(auto& cell : cells)
            cell.delta_t = deltaT;
    }


}

void CFlow_solver::set_boundary_conditions(const CFVM::CFVM_BC_face &FVM_BC_face , CFlow_BC_face &BC_face,const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    if(BC_face.eflow_BC_type == enums::Eflow_BC_type::velocity_inlet)
    {
        BC_face.Pressure = BC_face.left_cell_ptr->Pressure;
        BC_face.Density = Fluid_properties->get_Density(BC_face.left_cell_ptr->Pressure,BC_face.Temperature);
        BC_face.Mdot = BC_face.Density * dot((BC_face.Velocity),FVM_BC_face.area_vector);
        BC_face.Enthalpy = Fluid_properties->get_Enthalpy(BC_face.Temperature,BC_face.Velocity);
    }

    // DO notring if supersonic velocity inlet

    if(BC_face.eflow_BC_type == enums::Eflow_BC_type::pressure_outlet)
    {
        BC_face.Temperature = BC_face.left_cell_ptr->Temperature;
        BC_face.Velocity = BC_face.left_cell_ptr->Velocity;
        BC_face.Density = Fluid_properties->get_Density(BC_face.Pressure,BC_face.left_cell_ptr->Temperature);
        BC_face.Enthalpy = Fluid_properties->get_Enthalpy(BC_face.left_cell_ptr->Temperature,BC_face.left_cell_ptr->Velocity);
        if((BC_face.Mdot /BC_face.left_cell_ptr->Density / FVM_BC_face.area) >= Fluid_properties->get_speed_of_sound(BC_face.left_cell_ptr->Temperature))
        {
            BC_face.supersonic_flag = true;
            BC_face.Pressure = BC_face.left_cell_ptr->Pressure;
        }
        else
        {
            BC_face.supersonic_flag = true;
            BC_face.Pressure = BC_face.Out_Pressure;
        }
    }

    if(BC_face.eflow_BC_type == enums::Eflow_BC_type::wall)
    {
        BC_face.Pressure = BC_face.left_cell_ptr->Pressure + dot(BC_face.left_cell_ptr->Pressure_grad,FVM_BC_face.cell_center_unit_vector)
                *FVM_BC_face.cell_center_distance;
    }
}

template<enums::Etimestepping_type etimestepping_type , enums::Eflow_type_visc eflow_type_visc , enums::Econv_flux_type econv_flux_type >
void CFlow_solver::iterate(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{

    // Set Boundary conditions
    size_t i_BC_face;
    #pragma omp parallel for
    for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
        set_boundary_conditions(FVM->BC_faces[i_BC_face], BC_faces[i_BC_face], Fluid_properties);

    // Calculate pressure gradients
    calc_press_grad(FVM,Fluid_properties);

    // Clear ap data
    size_t i_cell;
    #pragma omp parallel for
    for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
        cells[i_cell].ap = 0.0;

    //Face loop for explicit part of viscosity
    if constexpr(eflow_type_visc == enums::Eflow_type_visc::Viscous)
    {
        size_t i_face;
        #pragma omp parallel for
        for(i_face = 0; i_face < faces.size() ; i_face++)
            calc_flux_explicit_visc_part( FVM->faces[i_face], faces[i_face] , Fluid_properties);
    }

    // Solve Uncorrected Momentum equations
    for(int i_dimention = 0 ; i_dimention < 3 ; i_dimention++)
    {
        // Clear matrices
        Matrix_solver->matrix.clear();
        Matrix_solver->b_vector.clear();
        Matrix_solver->x_vector.clear();

        #pragma omp parallel
        {
            //Face loop
            size_t i_face;
            #pragma omp for
            for(i_face = 0; i_face < faces.size() ; i_face++)
                calc_flux_momentum<eflow_type_visc , econv_flux_type>(i_dimention, FVM->faces[i_face], faces[i_face] , Matrix_solver,Fluid_properties);
            #pragma omp barrier

            // BC Face Loop
            size_t i_BC_face;
            #pragma omp for
            for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
                calc_flux_momentum<eflow_type_visc>(i_dimention, FVM->BC_faces[i_BC_face], BC_faces[i_BC_face], Matrix_solver,Fluid_properties);
            #pragma omp barrier

            // Cell loop
            size_t i_cell;
            #pragma omp for
            for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
                cell_loop_momentum<etimestepping_type>(i_dimention, FVM->cells [i_cell] ,cells[i_cell], Matrix_solver, Fluid_properties);
        }

        // Solve matrix
        Matrix_solver->solve();

        // Update uncorrected velocity variables
        double dum_residual = 0.0;
        size_t i_cell;
        #pragma omp parallel for reduction (+:dum_residual)
        for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
        {
            dum_residual += abs(cells[i_cell].Velocity.val[i_dimention] - Matrix_solver->x_vector.get_element(i_cell));
            cells[i_cell].Velocity.val[i_dimention] = Matrix_solver->x_vector.get_element(i_cell);
        }

        Momentum_global_residual.val[i_dimention] = dum_residual;
    }

    // Update face Mdots
    size_t i_face;
    #pragma omp parallel for
    for(i_face = 0; i_face < faces.size() ; i_face++)
        calc_face_mdot<econv_flux_type >(FVM->faces[i_face],faces[i_face]);


    i_BC_face = 0;
    #pragma omp parallel for
    for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
        calc_face_mdot(FVM->BC_faces[i_BC_face], BC_faces[i_BC_face] , Fluid_properties);


    // Solve Pressure Equation
    // Clear matrices
    Matrix_solver->matrix.clear();
    Matrix_solver->b_vector.clear();
    Matrix_solver->x_vector.clear();

    #pragma omp parallel
    {
        //Face loop
        i_face = 0;
        #pragma omp for
        for(i_face = 0; i_face < faces.size() ; i_face++)
            calc_flux_pressure<econv_flux_type >(FVM->faces[i_face], faces[i_face], Matrix_solver,Fluid_properties);
        #pragma omp barrier

        // BC Face Loop
        i_BC_face = 0;
        #pragma omp for
        for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
            calc_flux_pressure(FVM->BC_faces[i_BC_face], BC_faces[i_BC_face], Matrix_solver,Fluid_properties);
        #pragma omp barrier

        // Cell loop
        i_cell = 0;
        #pragma omp for
        for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
            cell_loop_pressure<etimestepping_type>(FVM->cells[i_cell], cells[i_cell], Matrix_solver, Fluid_properties);
    }

    Pressure_global_residual= CMatrix_solver::ABSsum(Matrix_solver->b_vector);

    // Solve matrix
    Matrix_solver->solve();

    // Correct pressure
    i_cell = 0;
    #pragma omp parallel for
    for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
    {
        cells[i_cell].Pressure_corr = explicit_underrelaxation_P * Matrix_solver->x_vector.get_element(i_cell);
        cells[i_cell].Pressure += cells[i_cell].Pressure_corr;
    }


    i_BC_face = 0;
    #pragma omp parallel for
    for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
    {
        if(BC_faces[i_BC_face].eflow_BC_type == enums::Eflow_BC_type::velocity_inlet)
        {
            BC_faces[i_BC_face].Pressure = BC_faces[i_BC_face].left_cell_ptr->Pressure;
        }
        if(BC_faces[i_BC_face].eflow_BC_type == enums::Eflow_BC_type::pressure_outlet)
        {
            if(BC_faces[i_BC_face].supersonic_flag)
                BC_faces[i_BC_face].Pressure = BC_faces[i_BC_face].left_cell_ptr->Pressure;
        }
    }

    // Calculate pressure correction gradients
    calc_press_corr_grad(FVM,Fluid_properties);


    // Correct face Mass Flow Rates
    i_face = 0;
    #pragma omp parallel for
    for(i_face = 0; i_face < faces.size() ; i_face++)
        correct_face_mdot <econv_flux_type> (FVM->faces[i_face], faces[i_face], Fluid_properties);

    i_BC_face = 0;
    #pragma omp parallel for
    for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
        correct_face_mdot(FVM->BC_faces[i_BC_face],BC_faces[i_BC_face],Fluid_properties);



    // Correct densities
    i_cell = 0;
    #pragma omp parallel for
    for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
    {
        cells[i_cell].Density = Fluid_properties->get_Density(cells[i_cell].Pressure,cells[i_cell].Temperature);
    }

    i_BC_face = 0;
    #pragma omp parallel for
    for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
    {
        if(BC_faces[i_BC_face].eflow_BC_type == enums::Eflow_BC_type::velocity_inlet)
        {
            BC_faces[i_BC_face].Density = Fluid_properties->get_Density(BC_faces[i_BC_face].left_cell_ptr->Pressure,BC_faces[i_BC_face].Temperature);
        }
        // Do notring if supersonic velocity inlet
        if(BC_faces[i_BC_face].eflow_BC_type == enums::Eflow_BC_type::pressure_outlet)
        {
            if(BC_faces[i_BC_face].supersonic_flag)
                BC_faces[i_BC_face].Density = Fluid_properties->get_Density(BC_faces[i_BC_face].left_cell_ptr->Pressure,BC_faces[i_BC_face].left_cell_ptr->Temperature);
            else
                BC_faces[i_BC_face].Density = Fluid_properties->get_Density(BC_faces[i_BC_face].Pressure,BC_faces[i_BC_face].left_cell_ptr->Temperature);
        }
    }


    // Correct velocities
    i_cell = 0;
    #pragma omp parallel for
    for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
    {
        cells[i_cell].Velocity -= explicit_underrelaxation_u *( FVM->cells[i_cell].volume) * dot_div(cells[i_cell].Pressure_corr_grad , cells[i_cell].ap_vec);
    }


    // Solve Enthalpy Equation
    // Clear matrices
    Matrix_solver->matrix.clear();
    Matrix_solver->b_vector.clear();
    Matrix_solver->x_vector.clear();

    #pragma omp parallel
    {
        //Face loop
        i_face = 0;
        #pragma omp for
        for(i_face = 0; i_face < faces.size() ; i_face++)
            calc_flux_enthalpy <econv_flux_type> (FVM->faces[i_face], faces[i_face], Matrix_solver,Fluid_properties);

        // BC Face Loop
        i_BC_face = 0;
        #pragma omp for
        for(i_BC_face = 0; i_BC_face < BC_faces.size() ; i_BC_face++)
            calc_flux_enthalpy(FVM->BC_faces[i_BC_face], BC_faces[i_BC_face], Matrix_solver,Fluid_properties);

        // Cell loop
        i_cell = 0;
        #pragma omp for
        for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
            cell_loop_enthalpy<etimestepping_type>(FVM->cells[i_cell], cells[i_cell], Matrix_solver, Fluid_properties);
    }


    Enthalpy_global_residual= 0.0;

    // Solve matrix
    Matrix_solver->solve();

    // Update Enthalpies and calculate new temperature
    double dum_residual = 0.0;
    i_cell = 0;
    #pragma omp parallel for reduction (+:dum_residual)
    for(i_cell = 0 ; i_cell < cells.size() ; i_cell++)
    {
        dum_residual += abs(cells[i_cell].Enthalpy - Matrix_solver->x_vector.get_element(i_cell));
        cells[i_cell].Enthalpy = Matrix_solver->x_vector.get_element(i_cell);
        cells[i_cell].Temperature = Fluid_properties->get_Temperature(cells[i_cell].Temperature , cells[i_cell].Enthalpy, cells[i_cell].Velocity);
        cells[i_cell].Mach = mag(cells[i_cell].Velocity) / Fluid_properties->get_speed_of_sound(cells[i_cell].Temperature);
    }

    Enthalpy_global_residual = dum_residual;

}

template void CFlow_solver::iterate<enums::Etimestepping_type::Steady, enums::Eflow_type_visc::Viscous , enums::Econv_flux_type::First_order_upwind >
(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::iterate<enums::Etimestepping_type::Steady, enums::Eflow_type_visc::Viscous , enums::Econv_flux_type::Blended_central >
(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::iterate<enums::Etimestepping_type::Steady, enums::Eflow_type_visc::Inviscid , enums::Econv_flux_type::First_order_upwind >
(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::iterate<enums::Etimestepping_type::Steady, enums::Eflow_type_visc::Inviscid , enums::Econv_flux_type::Blended_central >
(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::iterate<enums::Etimestepping_type::Transient, enums::Eflow_type_visc::Viscous , enums::Econv_flux_type::First_order_upwind >
(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::iterate<enums::Etimestepping_type::Transient, enums::Eflow_type_visc::Viscous , enums::Econv_flux_type::Blended_central >
(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::iterate<enums::Etimestepping_type::Transient, enums::Eflow_type_visc::Inviscid , enums::Econv_flux_type::First_order_upwind >
(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::iterate<enums::Etimestepping_type::Transient, enums::Eflow_type_visc::Inviscid , enums::Econv_flux_type::Blended_central >
(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


void CFlow_solver::update_timestep()
{
    for(auto & cell : cells)
    {
        cell.Velocity_last_t = cell.Velocity;
        cell.Density_last_t = cell.Density;
        cell.Pressure_last_t = cell.Pressure;
        cell.Enthalpy_last_t = cell.Enthalpy;
    }
    // Implicit underrelaxation ramping

    if(root.configuration->etimestepping_type == enums::Etimestepping_type::Steady && this->implicit_underrelaxation<root.configuration->implicit_underrelax_max)
    {
        this->implicit_underrelaxation += (root.configuration->implicit_underrelax_max - root.configuration->implicit_underrelax)/root.configuration->underrelax_ramp_iteration;
    }

}
