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
#include "tensor.h"

template<enums::Econv_flux_type econv_flux_type >
void CFlow_solver::calc_face_mdot(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face)
{
    Flow_face.Qy = ((FVM_face.left_cell_ptr->volume / Flow_face.left_cell_ptr->ap)+(FVM_face.right_cell_ptr->volume / Flow_face.right_cell_ptr->ap))*
            dot(FVM_face.alpha,FVM_face.area_vector);

    Svector pgrad_face = solvers_common::face_average(FVM_face,Flow_face.left_cell_ptr->Pressure_grad , Flow_face.right_cell_ptr->Pressure_grad);
    double YF = Flow_face.Qy*(Flow_face.right_cell_ptr->Pressure - Flow_face.left_cell_ptr->Pressure -
                    dot(pgrad_face,(FVM_face.cell_center_distance*FVM_face.cell_center_unit_vector)));

    Svector velface = solvers_common::face_average(FVM_face,Flow_face.left_cell_ptr->Velocity, Flow_face.right_cell_ptr->Velocity);
    double mdot = dot(velface,FVM_face.area_vector) - YF;

    Flow_face.Density = solvers_common::Convective_face_interp<econv_flux_type>
            (FVM_face, mdot,Flow_face.left_cell_ptr->Density,Flow_face.right_cell_ptr->Density);

    mdot *= Flow_face.Density;

    Flow_face.Mdot = mdot;
}

template void CFlow_solver::calc_face_mdot<enums::Econv_flux_type::First_order_upwind>
(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face);

template void CFlow_solver::calc_face_mdot<enums::Econv_flux_type::Blended_central>
(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face);

void CFlow_solver::calc_face_mdot(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::wall)
    {

    }
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::velocity_inlet)
    {
        Flow_BC_face.Mdot = Fluid_properties->get_Density(Flow_BC_face.left_cell_ptr->Pressure,Flow_BC_face.Temperature)*
                dot((Flow_BC_face.Velocity),FVM_BC_face.area_vector);
    }

    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::supersonic_velocity_inlet)
    {
        Flow_BC_face.Qy = (FVM_BC_face.left_cell_ptr->volume)/(Flow_BC_face.left_cell_ptr->ap)*
                dot(FVM_BC_face.alpha,FVM_BC_face.area_vector);
    }

    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::pressure_outlet)
    {
        if(Flow_BC_face.supersonic_flag)
        {
            Flow_BC_face.Mdot = Fluid_properties->get_Density(Flow_BC_face.left_cell_ptr->Pressure,Flow_BC_face.left_cell_ptr->Temperature)*
                    dot((Flow_BC_face.left_cell_ptr->Velocity),FVM_BC_face.area_vector);
            Flow_BC_face.Qy = 1.0;
        }
        else
        {
            Flow_BC_face.Qy = (FVM_BC_face.left_cell_ptr->volume)/(Flow_BC_face.left_cell_ptr->ap)*
                    dot(FVM_BC_face.alpha,FVM_BC_face.area_vector);

            double YF = Flow_BC_face.Qy*(Flow_BC_face.Pressure - Flow_BC_face.left_cell_ptr->Pressure -
                           0.5 * dot((Flow_BC_face.left_cell_ptr->Pressure_grad),(FVM_BC_face.cell_center_distance*FVM_BC_face.cell_center_unit_vector)));

            double mdot = dot((Flow_BC_face.left_cell_ptr->Velocity),FVM_BC_face.area_vector) - YF;

            mdot *= Flow_BC_face.Density;

            Flow_BC_face.Mdot = mdot;
        }
    }
}

template<enums::Econv_flux_type econv_flux_type >
void CFlow_solver::correct_face_mdot(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    /*double T_face = solvers_common::face_average(FVM_face,Flow_face.left_cell_ptr->Temperature,Flow_face.right_cell_ptr->Temperature);
    double pcorr_face = solvers_common::face_average(FVM_face,Flow_face.left_cell_ptr->Pressure_corr , Flow_face.right_cell_ptr->Pressure_corr);*/

    double T_face = solvers_common::Convective_face_interp<econv_flux_type>
            (FVM_face, Flow_face.Mdot,Flow_face.left_cell_ptr->Temperature,Flow_face.right_cell_ptr->Temperature);

    double pcorr_face = solvers_common::Convective_face_interp<econv_flux_type>
            (FVM_face,Flow_face.Mdot, Flow_face.left_cell_ptr->Pressure_corr , Flow_face.right_cell_ptr->Pressure_corr);

    Flow_face.Mdot += Flow_face.Mdot/Flow_face.Density *
                Fluid_properties->get_Pressure_Rho_derivative(T_face)*pcorr_face;

    Flow_face.Mdot += Flow_face.Qy*Flow_face.Density * (Flow_face.left_cell_ptr->Pressure_corr - Flow_face.right_cell_ptr->Pressure_corr);

}

template void CFlow_solver::correct_face_mdot<enums::Econv_flux_type::First_order_upwind>
(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::correct_face_mdot<enums::Econv_flux_type::Blended_central>
(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

void CFlow_solver::correct_face_mdot(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &BC_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{

    // DO nothing on normal and supersonic velocity inlet
    if(BC_face.eflow_BC_type == enums::Eflow_BC_type::pressure_outlet)
    {
        if(!BC_face.supersonic_flag)
            BC_face.Mdot += BC_face.Qy*BC_face.Density * (BC_face.left_cell_ptr->Pressure_corr);
    }

}

void CFlow_solver::calc_press_grad(const std::shared_ptr<CFVM> &FVM,  const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    // Clear gradients
    for(auto& cell : cells)
        cell.Pressure_grad.clear();

    //Face loop
    auto face = faces.begin();
    auto FVM_face = FVM->faces.begin();
    for(; face != faces.end() ; face++, FVM_face++)
    {
        double face_pressure ;
        face_pressure = (face->left_cell_ptr->ap + face->right_cell_ptr->ap) ;
        face_pressure = (face->left_cell_ptr->ap * face->left_cell_ptr->Pressure +
            face->right_cell_ptr->ap * face->right_cell_ptr->Pressure) / face_pressure;

        face->left_cell_ptr->Pressure_grad = face->left_cell_ptr->Pressure_grad + face_pressure*FVM_face->area_vector;
        face->right_cell_ptr->Pressure_grad = face->right_cell_ptr->Pressure_grad - face_pressure*FVM_face->area_vector;
    }


    // BC Face Loop
    auto BC_face = BC_faces.begin();
    auto FVM_BC_face = FVM->BC_faces.begin();
    for(; BC_face != BC_faces.end() ; BC_face++, FVM_BC_face++)
    {
        if( BC_face->eflow_BC_type == enums::Eflow_BC_type::wall ||BC_face->eflow_BC_type == enums::Eflow_BC_type::symmetry
                || BC_face->eflow_BC_type == enums::Eflow_BC_type::empty || BC_face->eflow_BC_type == enums::Eflow_BC_type::velocity_inlet)
        {
            BC_face->left_cell_ptr->Pressure_grad = BC_face->left_cell_ptr->Pressure_grad + BC_face->left_cell_ptr->Pressure*FVM_BC_face->area_vector;
        }
        if( BC_face->eflow_BC_type == enums::Eflow_BC_type::pressure_outlet
                || BC_face->eflow_BC_type == enums::Eflow_BC_type::supersonic_velocity_inlet)
        {
            BC_face->left_cell_ptr->Pressure_grad = BC_face->left_cell_ptr->Pressure_grad + BC_face->Pressure*FVM_BC_face->area_vector;
        }
    }

    // Cell loop
    auto cell = cells.begin();
    auto FVM_cell = FVM->cells.begin();
    for(; cell != cells.end() ; cell++, FVM_cell++)
    {
        cell->Pressure_grad = (1.0/FVM_cell->volume)*cell->Pressure_grad;
    }
}

void CFlow_solver::calc_press_corr_grad(const std::shared_ptr<CFVM> &FVM, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    // Clear gradients
    for(auto& cell : cells)
        cell.Pressure_corr_grad.clear();

    //Face loop
    auto face = faces.begin();
    auto FVM_face = FVM->faces.begin();
    for(; face != faces.end() ; face++, FVM_face++)
    {
        double face_pressure ;
        face_pressure = solvers_common::face_average(*FVM_face, face->left_cell_ptr->Pressure_corr , face->right_cell_ptr->Pressure_corr);

        face->left_cell_ptr->Pressure_corr_grad = face->left_cell_ptr->Pressure_corr_grad + face_pressure*FVM_face->area_vector;
        face->right_cell_ptr->Pressure_corr_grad = face->right_cell_ptr->Pressure_corr_grad - face_pressure*FVM_face->area_vector;
    }


    // BC Face Loop
    auto BC_face = BC_faces.begin();
    auto FVM_BC_face = FVM->BC_faces.begin();
    for(; BC_face != BC_faces.end() ; BC_face++, FVM_BC_face++)
    {
        if(BC_face->eflow_BC_type == enums::Eflow_BC_type::wall || BC_face->eflow_BC_type == enums::Eflow_BC_type::symmetry || BC_face->eflow_BC_type == enums::Eflow_BC_type::empty
                || BC_face->eflow_BC_type == enums::Eflow_BC_type::velocity_inlet || (BC_face->supersonic_flag &&  BC_face->eflow_BC_type == enums::Eflow_BC_type::pressure_outlet) )
        {
            BC_face->left_cell_ptr->Pressure_corr_grad = BC_face->left_cell_ptr->Pressure_corr_grad + BC_face->left_cell_ptr->Pressure_corr * FVM_BC_face->area_vector;
        }


        // If boundary is pressure outlet or supersonic velocity inlet , pressure correction is zero so no operation here
    }

    // Cell loop
    auto cell = cells.begin();
    auto FVM_cell = FVM->cells.begin();
    for(; cell != cells.end() ; cell++, FVM_cell++)
    {
        cell->Pressure_corr_grad = (1.0/FVM_cell->volume)*cell->Pressure_corr_grad;
    }
}


Svector CFlow_solver::calc_wall_stress(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    Svector res;
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::wall)
    {
        if( root.configuration->eflow_type_visc == enums::Eflow_type_visc::Viscous && Flow_BC_face.ewall_BC_type == enums::Ewall_BC_type::no_slip)
        {
            Svector parvel;
            parvel = Flow_BC_face.left_cell_ptr->Velocity - (dot(Flow_BC_face.left_cell_ptr->Velocity,(1.0/FVM_BC_face.area)*FVM_BC_face.area_vector)
                                                             * ((1.0/FVM_BC_face.area)*FVM_BC_face.area_vector));

            res = Fluid_properties->get_viscosity(Flow_BC_face.left_cell_ptr->Temperature)
                    / FVM_BC_face.perp_dist * parvel;
        }
    }
    return res;
}


/*
void CFlow_solver::calc_timestep(std::shared_ptr<CFVM> &FVM, std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    // Clear timesteps
    for(auto& cell : cells)
        cell.delta_t = 0.0;

    //Face loop
    auto face = faces.begin();
    auto FVM_face = FVM->faces.begin();
    for(; face != faces.end() ; face++, FVM_face++)
    {
        Svector facevel ;
        facevel = (face->left_cell_ptr->Velocity + face->right_cell_ptr->Velocity);
        facevel = 0.5 * facevel;

        double facec;
        facec = Fluid_properties->get_speed_of_sound(face->left_cell_ptr->Temperature);
        facec = 0.5 * (facec+ Fluid_properties->get_speed_of_sound(face->right_cell_ptr->Temperature));

        facec *= FVM_face->area;
        facec = facec + abs(dot(facevel,FVM_face->area_vector));

        face->left_cell_ptr->delta_t += facec;
        face->right_cell_ptr->delta_t += facec;
    }


    // BC Face Loop
    auto BC_face = BC_faces.begin();
    auto FVM_BC_face = FVM->BC_faces.begin();
    for(; BC_face != BC_faces.end() ; BC_face++, FVM_BC_face++)
    {
        if(BC_face->eflow_BC_type == enums::Eflow_BC_type::wall || BC_face->eflow_BC_type == enums::Eflow_BC_type::empty ||  BC_face->eflow_BC_type == enums::Eflow_BC_type::pressure_outlet)
        {
            BC_face->left_cell_ptr->delta_t += abs(dot(BC_face->left_cell_ptr->Velocity,FVM_BC_face->area_vector));
            BC_face->left_cell_ptr->delta_t += Fluid_properties->get_speed_of_sound(BC_face->left_cell_ptr->Temperature) * FVM_BC_face->area;
        }

        if(BC_face->eflow_BC_type == enums::Eflow_BC_type::velocity_inlet)
        {
            BC_face->left_cell_ptr->delta_t += abs(dot(BC_face->Velocity,FVM_BC_face->area_vector));
            BC_face->left_cell_ptr->delta_t += Fluid_properties->get_speed_of_sound(BC_face->Temperature) * FVM_BC_face->area;
        }
    }

    // Cell loop
    auto cell = cells.begin();
    auto FVM_cell = FVM->cells.begin();
    for(; cell != cells.end() ; cell++, FVM_cell++)
    {
        cell->delta_t = CFL * FVM_cell->volume / cell->delta_t;
    }
}*/

