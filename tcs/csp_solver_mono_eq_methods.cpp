/**
BSD-3-Clause
Copyright 2019 Alliance for Sustainable Energy, LLC
Redistribution and use in source and binary forms, with or without modification, are permitted provided 
that the following conditions are met :
1.  Redistributions of source code must retain the above copyright notice, this list of conditions 
and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used to endorse 
or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT HOLDER, CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES 
DEPARTMENT OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, 
OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "csp_solver_core.h"
#include "numeric_solvers.h"
#include "CO2_properties.h"
#include <cmath>
#include <math.h>
#include <algorithm>

#include "lib_util.h"

int C_csp_solver::C_MEQ_cr_on__pc_q_dot_max__tes_off__defocus::operator()(double defocus /*-*/, double *diff_q_dot_pc /*MWt*/)
{
    // the last argument is just so it compiles -> need to change everything here
    C_MEQ_cr_on__pc_q_dot_max__tes_off c_eq(mpc_csp_solver, m_pc_mode, defocus);
    C_monotonic_eq_solver c_solver(c_eq);

    c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

    double T_htf_cold_guess_colder = mpc_csp_solver->m_T_htf_pc_cold_est + 30;  //[C], convert from [K]
    double T_htf_cold_guess_warmer = T_htf_cold_guess_colder + 10.0;        //[C]

    double T_htf_cold_solved, tol_solved;
    T_htf_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    int solver_code = 0;
    try
    {
        solver_code = c_solver.solve(T_htf_cold_guess_colder, T_htf_cold_guess_warmer, 0.0, T_htf_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception("C_MEQ_cr_on__pc_max__tes_off__defocus->C_MEQ_cr_on__pc_q_dot_max__tes_off received exception from mono equation solver"));
    }

    if (solver_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (solver_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) <= 0.1)
        {
            mpc_csp_solver->error_msg = util::format("At time = %lg the C_MEQ_cr_on__pc_max__tes_off__defocus->C_MEQ_cr_on__pc_q_dot_max__tes_off iteration to find the cold HTF temperature connecting the power cycle and receiver only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, mpc_csp_solver->error_msg);
        }
        else
        {
            *diff_q_dot_pc = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    *diff_q_dot_pc = (mpc_csp_solver->mc_pc_out_solver.m_q_dot_htf - m_q_dot_max) / m_q_dot_max;            //[-]

    return 0;
}

int C_csp_solver::C_MEQ_cr_on__pc_q_dot_max__tes_off::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Should not be called directly, only via C_MEQ_cr_on__pc_q_dot_max__tes_off__defocus::operator()(double defocus /*-*/, double *diff_q_dot_pc /*MWt*/)
    C_csp_tes::S_csp_tes_outputs tes_outputs;           // the aggregate of the different tes calls
    C_csp_tes::S_csp_tes_outputs tes_outputs_temp;      // output of each tes call, used to update the aggregate

    mpc_csp_solver->mc_tes.use_calc_vals(true);

    // Solve the tower model with T_htf_cold from the LT HX
    double T_htf_rec_in = T_htf_cold + 273.15;      //[K]
    double P_rec_in = mpc_csp_solver->mc_cr_htf_state_in.m_pres;    //[kPa]
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_rec_in - 273.15;      //[C]
    double P_in = mpc_csp_solver->m_P_cold_des;                 //[kPa] use the receiver design inlet pressure
    mpc_csp_solver->mc_cr_htf_state_in.m_pres = P_in;

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outputs
    double m_dot_rec_out = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double T_htf_rec_out = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15;  //[K]
    double P_rec_out = P_rec_in - mpc_csp_solver->mc_cr_out_solver.m_dP_sf * 100.;  //[kPa]
    double m_dot_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
    double T_store_in = mpc_csp_solver->mc_cr_out_solver.m_T_store_hot + 273.15;   //[K]

    // Charge storage
    // First set available charge to that coming from the tower
    mpc_csp_solver->mc_tes.set_max_charge_flow(m_dot_store);
    double T_cold_tes_K;
    mpc_csp_solver->mc_tes.charge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_store_in,
        T_cold_tes_K,
        tes_outputs);

    // First estimate available discharge in order to updated m_m_dot_tes_dc_max
    double q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est;
    mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_rec_out, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est);

    mpc_csp_solver->mc_tes.update_calc_vals(false);     // do not update calc values due to following iterations (which are within larger iterations)

    // Solve the HT HX using a steady-state media discharge (m_dot_store from the receiver)
    // This is a test call (update_calc_vals = false) using the receiver outlet temperature
    // The .calc values are not updated so discharge() is called again later to update them.
    double T_htf_hx_out, m_dot_hx_out;
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_htf_rec_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    m_dot_hx_out *= 3600.;

    double T_htf_hx_in, m_dot_hx_in, T_htf_pc_in, m_dot_pc_in, P_hx_out;
    if (m_dot_rec_out > m_dot_hx_out) {
        // Needs defocusing if this is the converged state
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    else {
        C_MEQ_cr_on_tes_dc_m_dot_tank c_eq(mpc_csp_solver, T_htf_rec_in, T_htf_rec_out, P_in, P_rec_out, m_dot_rec_out, m_dot_store);
        C_monotonic_eq_solver c_solver(c_eq);

        // Set up solver
        c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

        // Solve for cold temperature
        double T_cold_guess_low = std::min(T_htf_rec_in, T_htf_rec_out) - 273.15 - 10;  //[C]
        double T_cold_guess_high = std::max(T_htf_rec_in, T_htf_rec_out) - 273.15 + 10; //[C]

        double T_cold_solved, tol_solved;
        T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
        int iter_solved = -1;

        int T_cold_code = 0;
        try
        {
            T_cold_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
        }
        catch (C_csp_exception)
        {
            mpc_csp_solver->mc_tes.use_calc_vals(false);
            mpc_csp_solver->mc_tes.update_calc_vals(true);
            throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
        }

        if (T_cold_code != C_monotonic_eq_solver::CONVERGED)
        {
            if (T_cold_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
            {
                std::string msg = util::format("At time = %lg C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank "
                    "iteration to find the cold HTF temperature to balance energy between the CR and PC only reached a convergence "
                    "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                    mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
                mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
            }
            else
            {
                mpc_csp_solver->mc_tes.use_calc_vals(false);
                mpc_csp_solver->mc_tes.update_calc_vals(true);
                *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
                return -1;
            }
        }

        T_htf_hx_in = T_cold_solved + 273.15;
        m_dot_hx_in = m_dot_hx_out = c_eq.m_m_dot_htf_out; //[kg/hr] mass flow out of the HX on the field side
        m_dot_pc_in = m_dot_hx_out;     //[kg/hr]
        T_htf_hx_out = mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out + 273.15;     //[K]
        P_hx_out = P_rec_out * (1. - mpc_csp_solver->mc_tes_outputs.dP_perc / 100.);    //[kPa]

        // Recombine mass flow that was diverted around the CR with that after the HT HX
        double P_hx_in = P_rec_out;            //[kPa]
        double m_dot_bypassed = m_dot_hx_out - m_dot_rec_out;  //[kg/hr]

        // get enthalpy, assume sCO2 HTF
        CO2_state co2_props;
        int prop_error_code = CO2_TP(T_htf_rec_in, P_rec_in, &co2_props);
        double h_in = co2_props.enth;
        double h_out = h_in;
        prop_error_code = CO2_PH(P_hx_out, h_out, &co2_props);
        double T_htf_bypassed = co2_props.temp; //[K]

        T_htf_pc_in = (T_htf_hx_out * m_dot_hx_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_out + m_dot_bypassed);  //[K]  mix streams to get PC inlet temp
    }

    mpc_csp_solver->mc_tes.update_calc_vals(true);

    // call discharge again with calc_vals = true to update the hot and warm tank .calc values
    double T_htf_hot;  //[K] HTF temp out of the HX on the field side
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_htf_hx_in,
        T_htf_hot,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_hot_ave = tes_outputs_temp.m_T_hot_ave - 273.15;       //[C]
    tes_outputs.m_m_dot = tes_outputs_temp.m_m_dot;
    tes_outputs.m_W_dot_rhtf_pump = tes_outputs_temp.m_W_dot_rhtf_pump;
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf = tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_hot_ave = tes_outputs_temp.m_T_hot_ave;
    tes_outputs.m_T_hot_final = tes_outputs_temp.m_T_hot_final;
    tes_outputs.dP_perc = tes_outputs_temp.dP_perc;

    // Solve the PC performance at the receiver htf flow rate
    // Need to do this to get back PC T_htf_cold
    // HTF State
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_pc_in - 273.15;   //[C]
    mpc_csp_solver->mc_pc_htf_state_in.m_pres = P_hx_out;   //[kPa]
    // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_pc_in;                         //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = C_csp_power_cycle::ON;     //[-]
    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }

    // Get power cycle HTF return state
    double T_htf_pc_out = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold + 273.15;       //[K]
    double m_dot_pc_out = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;                 //[kg/hr]
    double P_pc_out = mpc_csp_solver->mc_pc_out_solver.m_P_phx_in * 1000.;              //[kPa]

    // Discharge virtual warm tank through LT HX
    //double m_dot_hx_out;    //[kg/s]
    mpc_csp_solver->mc_tes.discharge_full_lt(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_pc_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_cold_ave = tes_outputs_temp.m_T_cold_ave - 273.15;       //[C]
    m_dot_hx_out *= 3600.;      //[kg/hr]
    double P_lthx_out = P_pc_out * (1. - tes_outputs_temp.dP_perc / 100.);           //[kPa]
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf += tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_cold_ave = tes_outputs_temp.m_T_cold_ave;
    tes_outputs.m_T_cold_final = tes_outputs_temp.m_T_cold_final;
    tes_outputs.dP_perc += tes_outputs_temp.dP_perc;

    // Recombine excess mass flow from the power cycle
    double T_htf_rec_in_solved;     //[K]
    if (m_dot_pc_out > m_dot_hx_out) {
        double m_dot_bypassed = m_dot_pc_out - m_dot_hx_out;                                    //[kg/hr]

        // get enthalpy, assume sCO2 HTF
        CO2_state co2_props;
        int prop_error_code = CO2_TP(T_htf_pc_out, P_pc_out, &co2_props);
        double h_in = co2_props.enth;
        double h_out = h_in;
        prop_error_code = CO2_PH(P_lthx_out, h_out, &co2_props);
        double T_htf_bypassed = co2_props.temp; //[K]

        T_htf_rec_in_solved = (T_htf_hx_out * m_dot_hx_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_out + m_dot_bypassed);  // [K]  mix streams to get LT HX outlet temp
    }
    else {
        T_htf_rec_in_solved = T_htf_hx_out;      //[K]
    }

    //Calculate pressure difference (which is not used)
    double diff_P = (P_lthx_out - P_in) / P_in;

    // Set charging inlet/outlet temps to hot/cold ave temps, respectively
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                  //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = T_store_hot_ave;    //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_store_cold_ave;  //[C]

    // Set discharge HTF state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_hx_in;                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_hx_in - 273.15;           //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hx_out - 273.15;         //[C]

    mpc_csp_solver->mc_tes_outputs = tes_outputs;

    //Calculate diff_T_htf_cold
    *diff_T_htf_cold = (T_htf_rec_in_solved - T_htf_cold) / T_htf_cold;     //[-]

    mpc_csp_solver->mc_tes.use_calc_vals(false);
    mpc_csp_solver->mc_tes.update_calc_vals(true);
    return 0;
}

int C_csp_solver::C_mono_eq_cr_to_pc_to_cr_m_dot::operator()(double m_dot /*kg/hr*/, double *m_dot_bal /*-*/)
{
    // Converge on an HTF mass flow (m_dot) that balances the particle mass flow (m_dot_bal) between the tower outlet
    //  and the hot tank outlet, when the TES is in steady state (could be full, could be empty, could just be steady-state)
    
    C_mono_eq_cr_to_pc_to_cr c_eq(mpc_csp_solver, m_pc_mode, m_P_field_in, -1, m_field_control_in, m_dot);
    C_monotonic_eq_solver c_solver(c_eq);

    // Set up solver
    double T_cold_min = mpc_csp_solver->m_cycle_T_htf_cold_des - 273.15 - 10.;
    double T_cold_max = mpc_csp_solver->m_T_htf_cold_des - 273.15 + 40.;
    c_solver.settings(1.E-3, 50, T_cold_min, T_cold_max, false);

    // Solve for cold temperature
    double T_cold_guess_low = mpc_csp_solver->m_T_htf_cold_des - 273.15;        //[C], convert from [K]
    double T_cold_guess_high = T_cold_guess_low + 10.0;                         //[C]

    double T_cold_solved, tol_solved;
    T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    int solver_code = 0;
    try
    {
        solver_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_mono_eq_cr_to_pc_to_cr", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
    }

    if (solver_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (solver_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver::C_mono_eq_cr_to_pc_to_cr "
                "iteration to find the cold HTF temperature only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            *m_dot_bal = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    // Calculate and report mass flow rate balance
    double m_dot_rec_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
    double m_dot_tes_store = mpc_csp_solver->mc_tes_outputs.m_m_dot * 3600.;        //[kg/hr]

    *m_dot_bal = (m_dot_rec_store - m_dot_tes_store) / m_dot_rec_store;         //[-]

    return 0;
}

int C_csp_solver::C_mono_eq_cr_to_pc_to_cr::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Should not be called directly, only via C_mono_eq_cr_to_pc_to_cr_m_dot(double m_dot /*kg/hr*/, double *m_dot_bal /*-*/)
    C_csp_tes::S_csp_tes_outputs tes_outputs;           // the aggregate of the different tes calls
    C_csp_tes::S_csp_tes_outputs tes_outputs_temp;      // output of each tes call, used to update the aggregate

    // Solve the tower model with T_htf_cold from the LT HX
    double T_htf_rec_in = T_htf_cold + 273.15;      //[K]
    double P_rec_in = mpc_csp_solver->mc_cr_htf_state_in.m_pres;    //[kPa]
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_rec_in - 273.15;      //[C]
    
    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
                        mpc_csp_solver->mc_cr_htf_state_in,
                        m_field_control_in,
                        mpc_csp_solver->mc_cr_out_solver,
                        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outputs
    double m_dot_rec_out = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double T_htf_rec_out = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15;  //[K]
    double P_rec_out = mpc_csp_solver->mc_cr_htf_state_in.m_pres - mpc_csp_solver->mc_cr_out_solver.m_dP_sf * 100.;  //[kPa]
    double m_dot_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
    double T_store_in = mpc_csp_solver->mc_cr_out_solver.m_T_store_hot + 273.15;   //[K]

    // Charge storage
    // First estimate available charge
    double q_dot_tes_ch_max, m_dot_tes_ch_max, T_tes_cold_ch_max, m_dot_store_ch_max;
    q_dot_tes_ch_max = m_dot_tes_ch_max = T_tes_cold_ch_max = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.charge_avail_est(T_store_in,
        mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_tes_ch_max,
        m_dot_tes_ch_max,
        T_tes_cold_ch_max,
        m_dot_store_ch_max);

    m_dot_tes_ch_max *= 3600.0;     //[kg/hr]

    // Test if particle flow from tower is greater than tes can store, factoring in later discharge
    if (m_dot_store > m_dot_tes_ch_max + m_dot_store) {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    else {
        // update available charge
        mpc_csp_solver->mc_tes.set_max_charge_flow((m_dot_tes_ch_max + m_dot_store) / 3600.);
    }

    double T_cold_tes_K;
    bool ch_solved = mpc_csp_solver->mc_tes.charge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_store_in,
        T_cold_tes_K,
        tes_outputs);

    // Check if TES.charge method solved
    if (!ch_solved) {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -3;
    }

    mpc_csp_solver->mc_tes.use_calc_vals(true);

    // First estimate available discharge in order to updated m_m_dot_tes_dc_max
    double q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est;
    mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_rec_out, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est);

    // Solve the HT HX using a given media discharge (m_m_dot_tank from the outer MEQ)
    // This is a test call (update_calc_vals = false) using the receiver outlet temperature
    // The .calc values are not updated so discharge() is called again later to update them.
    mpc_csp_solver->mc_tes.update_calc_vals(false);     // do not update calc values due to following iterations (which are within larger iterations)
    double T_htf_hx_out, m_dot_hx_out;
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_htf_rec_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    mpc_csp_solver->mc_tes.update_calc_vals(true);
    m_dot_hx_out *= 3600.;

    double T_htf_hx_in, m_dot_hx_in, T_htf_pc_in, m_dot_pc_in, P_hx_out;
    if (m_dot_rec_out > m_dot_hx_out) {
        // Needs defocusing if this is the converged state
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    else {
        C_MEQ_cr_on_tes_dc_m_dot_tank c_eq(mpc_csp_solver, T_htf_rec_in, T_htf_rec_out, P_rec_in, P_rec_out, m_dot_rec_out, m_dot_store);
        C_monotonic_eq_solver c_solver(c_eq);

        // Set up solver
        c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

        // Solve for cold temperature
        double T_cold_guess_low = std::min(T_htf_rec_in, T_htf_rec_out) - 273.15 - 10;  //[C]
        double T_cold_guess_high = std::max(T_htf_rec_in, T_htf_rec_out) - 273.15 + 10; //[C]

        double T_cold_solved, tol_solved;
        T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
        int iter_solved = -1;

        mpc_csp_solver->mc_tes.update_calc_vals(false);
        int T_cold_code = 0;
        try
        {
            T_cold_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
        }
        catch (C_csp_exception)
        {
            mpc_csp_solver->mc_tes.use_calc_vals(false);
            mpc_csp_solver->mc_tes.update_calc_vals(true);
            throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
        }
        mpc_csp_solver->mc_tes.update_calc_vals(true);

        if (T_cold_code != C_monotonic_eq_solver::CONVERGED)
        {
            if (T_cold_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
            {
                std::string msg = util::format("At time = %lg C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank "
                    "iteration to find the cold HTF temperature to balance energy between the CR and PC only reached a convergence "
                    "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                    mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
                mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
            }
            else
            {
                mpc_csp_solver->mc_tes.use_calc_vals(false);
                mpc_csp_solver->mc_tes.update_calc_vals(true);
                *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
                return -1;
            }
        }

        T_htf_hx_in = T_cold_solved + 273.15;
        m_dot_hx_in = m_dot_hx_out = c_eq.m_m_dot_htf_out; //[kg/hr] mass flow out of the HX on the field side
        m_dot_pc_in = m_dot_hx_out;     //[kg/hr]
        T_htf_pc_in = T_htf_hx_out = mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out + 273.15;     //[K]
        P_hx_out = P_rec_out * (1. - mpc_csp_solver->mc_tes_outputs.dP_perc / 100.);    //[kPa]
    }

    // call discharge again with calc_vals = true to update the hot and warm tank .calc values
    double T_htf_hot;  //[K] HTF temp out of the HX on the field side
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_htf_hx_in,
        T_htf_hot,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_hot_ave = tes_outputs_temp.m_T_hot_ave - 273.15;       //[C]
    tes_outputs.m_m_dot = tes_outputs_temp.m_m_dot;
    tes_outputs.m_W_dot_rhtf_pump = tes_outputs_temp.m_W_dot_rhtf_pump;
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf = tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_hot_ave = tes_outputs_temp.m_T_hot_ave;
    tes_outputs.m_T_hot_final = tes_outputs_temp.m_T_hot_final;
    tes_outputs.dP_perc = tes_outputs_temp.dP_perc;

    // Solve the PC performance at the receiver htf flow rate
    // Need to do this to get back PC T_htf_cold
    // HTF State
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_pc_in - 273.15;   //[C]
    mpc_csp_solver->mc_pc_htf_state_in.m_pres = P_hx_out;   //[kPa]
    // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_pc_in;                 //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode;         //[-]
    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }

    // Get power cycle HTF return state
    double T_htf_pc_out = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold + 273.15;       //[K]
    double m_dot_pc_out = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;                 //[kg/hr]
    double P_pc_out = mpc_csp_solver->mc_pc_out_solver.m_P_phx_in * 1000.;              //[kPa]

    // Discharge virtual warm tank through LT HX
    mpc_csp_solver->mc_tes.discharge_full_lt(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_pc_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_cold_ave = tes_outputs_temp.m_T_cold_ave - 273.15;       //[C]
    m_dot_hx_out *= 3600.;      //[kg/hr]
    double P_lthx_out = P_pc_out * (1. - tes_outputs_temp.dP_perc / 100.);           //[kPa]
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf += tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_cold_ave = tes_outputs_temp.m_T_cold_ave;
    tes_outputs.m_T_cold_final = tes_outputs_temp.m_T_cold_final;
    tes_outputs.dP_perc += tes_outputs_temp.dP_perc;

    // Recombine excess mass flow from the power cycle
    double T_htf_rec_in_solved;     //[K]
    if (m_dot_pc_out > m_dot_hx_out) {
        double m_dot_bypassed = m_dot_pc_out - m_dot_hx_out;                                    //[kg/hr]

        // get enthalpy, assume sCO2 HTF
        CO2_state co2_props;
        int prop_error_code = CO2_TP(T_htf_pc_out, P_pc_out, &co2_props);
        double h_in = co2_props.enth;
        double h_out = h_in;
        prop_error_code = CO2_PH(P_lthx_out, h_out, &co2_props);
        double T_htf_bypassed = co2_props.temp; //[K]

        T_htf_rec_in_solved = (T_htf_hx_out * m_dot_hx_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_out + m_dot_bypassed);  // [K]  mix streams to get LT HX outlet temp
    }
    else {
        T_htf_rec_in_solved = T_htf_hx_out;      //[K]
    }

    //Calculate pressure difference (which is not used)
    double diff_P = (P_lthx_out - P_rec_in) / P_rec_in;

    // Set charging inlet/outlet temps to hot/cold ave temps, respectively
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                  //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = T_store_hot_ave;    //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_store_cold_ave;  //[C]

    // Set discharge HTF state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_hx_in;                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_hx_in - 273.15;           //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hx_out - 273.15;         //[C]

    mpc_csp_solver->mc_tes_outputs = tes_outputs;

    //Calculate diff_T_htf_cold
    *diff_T_htf_cold = (T_htf_rec_in_solved - 273.15 - T_htf_cold) / T_htf_cold;        //[-]

    mpc_csp_solver->mc_tes.use_calc_vals(false);
    mpc_csp_solver->mc_tes.update_calc_vals(true);
    return 0;
}

int C_csp_solver::C_mono_eq_pc_su_cont_tes_dc::operator()(double T_htf_hot /*C*/, double *diff_T_htf_hot /*-*/)
{   
    // Call the power cycle in STARTUP_CONTROLLED mode
    mpc_csp_solver->mc_pc_inputs.m_m_dot = 0.0;     //[kg/hr]
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_hot;      //[C] convert from K
    mpc_csp_solver->mc_pc_inputs.m_standby_control = C_csp_power_cycle::STARTUP_CONTROLLED;

    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
                            mpc_csp_solver->mc_pc_htf_state_in,
                            mpc_csp_solver->mc_pc_inputs,
                            mpc_csp_solver->mc_pc_out_solver,
                            mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check for new timestep, probably will find one here
    m_time_pc_su = mpc_csp_solver->mc_pc_out_solver.m_time_required_su; //[s] power cycle model returns MIN(time required to completely startup, full timestep duration)

    // Get the required mass flow rate from the power cycle outputs
    double m_dot_pc = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf / 3600.0;    //[kg/s]

    // Reset mass flow rate in 'mc_pc_htf_state'
    mpc_csp_solver->mc_pc_inputs.m_m_dot = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;    //[kg/hr]

    double T_htf_cold = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold;      //[C]

    // Estimate available discharge in order to updated m_m_dot_tes_dc_max
    double q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est;
    mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_cold + 273.15, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est);
    m_dot_field_est *= 3600.;   //[kg/hr]

    // Solve TES discharge
    double T_htf_hot_calc = std::numeric_limits<double>::quiet_NaN();
    bool is_dc_solved = mpc_csp_solver->mc_tes.discharge_both(m_time_pc_su, 
                                            mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15, 
                                            m_dot_pc,
                                            T_htf_cold + 273.15,
                                            T_htf_hot_calc,
                                            mpc_csp_solver->mc_tes_outputs);

    if (!is_dc_solved) {
        *diff_T_htf_hot = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    T_htf_hot_calc = T_htf_hot_calc - 273.15;       //[C] convert from K

    // If not actually charging (i.e. mass flow rate = 0.0), set charging inlet/outlet temps to hot/cold ave temps, respectively
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                                                          //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;        //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;      //[C]

    // Set discharge HTF state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_pc*3600.0;      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_cold;         //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hot_calc;    //[C]

    if (is_dc_solved)
    {
        *diff_T_htf_hot = (T_htf_hot_calc - T_htf_hot) / T_htf_hot;
    }
    else
    {
        *diff_T_htf_hot = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    
    return 0;
}

int C_csp_solver::C_mono_eq_pc_target_tes_dc__m_dot::operator()(double m_dot_htf /*kg/hr*/, double *q_dot_pc /*MWt*/)
{
    double T_htf_hot = std::numeric_limits<double>::quiet_NaN();
    bool is_tes_success = mpc_csp_solver->mc_tes.discharge_both(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
                                                mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
                                                m_dot_htf / 3600.0,
                                                m_T_htf_cold + 273.15,
                                                T_htf_hot,
                                                mpc_csp_solver->mc_tes_outputs);

    if (!is_tes_success)
    {
        *q_dot_pc = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    T_htf_hot -= 273.15;        //[C] convert from K

    // HTF discharging state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_htf;        //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = m_T_htf_cold;   //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hot;     //[C]

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;              //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;    //[C], convert from K
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;  //[C], convert from K

    // Solve power cycle model
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_hot;      //[C]
        // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_htf;               //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode;     //[-]
        // Performance
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
                                    mpc_csp_solver->mc_pc_htf_state_in,
                                    mpc_csp_solver->mc_pc_inputs,
                                    mpc_csp_solver->mc_pc_out_solver,
                                    mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is producing power and solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful && mpc_csp_solver->mc_pc_inputs.m_standby_control == C_csp_power_cycle::ON)
    {
        *q_dot_pc = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    *q_dot_pc = mpc_csp_solver->mc_pc_out_solver.m_q_dot_htf;   //[MWt]
    return 0;
}

int C_csp_solver::C_mono_eq_pc_target_tes_dc__T_cold::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Expect mc_pc_out_solver to be set in inner mono eq loop that converges m_dot_htf
    C_mono_eq_pc_target_tes_dc__m_dot c_eq(mpc_csp_solver, m_pc_mode, T_htf_cold);
    C_monotonic_eq_solver c_solver(c_eq);

    // Calculate the maximum mass flow rate available for discharge
    double q_dot_tes_dc_max, m_dot_tes_dc_max, T_htf_hot_dc_max, m_dot_store_dc_max;
    q_dot_tes_dc_max = m_dot_tes_dc_max = T_htf_hot_dc_max = std::numeric_limits<double>::quiet_NaN();

    mpc_csp_solver->mc_tes.discharge_avail_est_both(T_htf_cold + 273.15,
                                        mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
                                        q_dot_tes_dc_max,
                                        m_dot_tes_dc_max,
                                        T_htf_hot_dc_max,
                                        m_dot_store_dc_max);

    m_dot_tes_dc_max *= 3600.0;     //[kg/hr] convert from kg/s

    // Now take the minimum of the max possible from TES and the max mass flow rate for the PC from design info
    double m_dot_max = fmin(m_dot_tes_dc_max, mpc_csp_solver->m_m_dot_pc_max);  //[kg/hr]

    // Use this to calculate the power cycle thermal power input
    m_q_dot_calc = std::numeric_limits<double>::quiet_NaN();
    int m_dot_code = c_solver.test_member_function(m_dot_max, &m_q_dot_calc);
    if (m_dot_code != 0)
    {
        // Should be able to pass the maximum mass flow rate to the power cycle
        // So not expecting failure here
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        m_m_dot_calc = m_dot_max;   //[kg/hr]
        return -1;
    }

    if ( (m_q_dot_calc - m_q_dot_target) / m_q_dot_target < -1.E-3 )
    {   // Can't achieve target thermal power without exceeding either TES or system mass flow rate constraints
        // So calculate T_htf_cold diff and return
        *diff_T_htf_cold = (mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold - T_htf_cold) / T_htf_cold;   //[-]
        m_m_dot_calc = m_dot_max;   //[kg/hr]
        return 0;
    }

    C_monotonic_eq_solver::S_xy_pair xy_pair_1;
    xy_pair_1.x = m_dot_max;        //[kg/hr]
    xy_pair_1.y = m_q_dot_calc;     //[MWt]

    // Check against minimum mass flow rate
        // What if this is 0?
    double m_dot_min = mpc_csp_solver->m_m_dot_pc_min;  //[kg/hr]
    if (m_dot_min > 0.0)
    {
        m_dot_code = c_solver.test_member_function(m_dot_min, &m_q_dot_calc);
    }
    else
    {
        m_dot_code = -1;
    }
    
    // If heat sink solved at m_dot_min, then check solved thermal power
    if (m_dot_code == 0)
    {
        if ((m_q_dot_calc - m_q_dot_target) / m_q_dot_target > 1.E-3)
        {   // At minimum mass flow rate the thermal power is still greater than target
            // So calculate T_htf_cold diff and return
            *diff_T_htf_cold = (mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold - T_htf_cold) / T_htf_cold;   //[-]
            m_m_dot_calc = m_dot_min;   //[kg/hr]
            return 0;
        }
    }

    // Guess another realistic mass flow rate
    double m_dot_guess = std::min( 0.97*xy_pair_1.x, m_q_dot_target / xy_pair_1.y * xy_pair_1.x );  //[kg/hr]

    // And calculate a second power cycle thermal power input
    m_dot_code = c_solver.test_member_function(m_dot_guess, &m_q_dot_calc);
    if (m_dot_code != 0)
    {
        // Should be able to pass this mass flow rate estimate to the power cycle
        // So not expecting failure here
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        m_m_dot_calc = m_dot_guess; //[kg/hr]
        return -1;
    }

    C_monotonic_eq_solver::S_xy_pair xy_pair_2;
    xy_pair_2.x = m_dot_guess;      //[kg/hr]
    xy_pair_2.y = m_q_dot_calc;     //[MWt]

    // Set up solver for mass flow rate
    c_solver.settings(1.E-3, 50, mpc_csp_solver->m_m_dot_pc_min, mpc_csp_solver->m_m_dot_pc_max, true);

    // Now solve for the required mass flow rate
    double m_dot_solved, tol_solved;
    m_dot_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    m_dot_code = 0;
    try
    {
        m_dot_code = c_solver.solve(xy_pair_1, xy_pair_2, m_q_dot_target, m_dot_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_mono_eq_pc_target_tes_dc__T_cold failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time),""));
    }

    if (m_dot_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (m_dot_code > C_monotonic_eq_solver::CONVERGED && std::abs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver::C_mono_eq_pc_target_tes_dc__T_cold "
                "iteration to find a mass flow rate resulting in the target power cycle heat input only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_mono_eq_pc_target_tes_dc__T_cold failed with Eq Solver Code %d", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time, m_dot_code), ""));
        }
    }

    m_m_dot_calc = m_dot_solved;        //[kg/hr]
    m_q_dot_calc = m_q_dot_target;      //[MWt]
    *diff_T_htf_cold = (mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold - T_htf_cold) / T_htf_cold;   //[-]
    return 0;
}

int C_csp_solver::C_mono_eq_pc_match_tes_empty::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // First, get the maximum possible mass flow rate from a full TES discharge
    double T_htf_tes_hot, m_dot_tes_dc;
    T_htf_tes_hot = m_dot_tes_dc = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.discharge_full_both(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
                            mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
                            T_htf_cold + 273.15,
                            T_htf_tes_hot, 
                            m_dot_tes_dc, 
                            mpc_csp_solver->mc_tes_outputs);

    // Set TES HTF states (this needs to be less bulky...)
    // HTF discharging state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_tes_dc*3600.0;  //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_cold;         //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_tes_hot - 273.15;    //[C]

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                                  //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;    //[C] convert from K
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;//[C] convert from K

    // Solve PC model
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_tes_hot - 273.15;     //[C]
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_tes_dc*3600.0;         //[kg/hr]

    // Inputs
    mpc_csp_solver->mc_pc_inputs.m_standby_control = C_csp_power_cycle::ON;

    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    *diff_T_htf_cold = (mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold - T_htf_cold) / T_htf_cold;

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on_pc_su_tes_ch_mdot::operator()(double m_dot_tank /*kg/hr*/, double *m_dot_htf_bal /*-*/)
{
    // Converge on the power cycle mass flow that is set during controlled startup

    C_mono_eq_cr_on_pc_su_tes_ch c_eq(mpc_csp_solver, m_pc_mode, m_defocus, m_dot_tank);
    C_monotonic_eq_solver c_solver(c_eq);

    // Set up solver
    c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

    // Solve for cold temperature
    double T_cold_guess_low = mpc_csp_solver->m_T_htf_cold_des - 273.15;        //[C], convert from [K]
    double T_cold_guess_high = T_cold_guess_low + 10.0;                         //[C]

    double T_cold_solved, tol_solved;
    T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    // Solve for hot tank mass flow
    int solver_code = 0;
    try
    {
        solver_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0, T_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_mono_eq_cr_on_pc_su_tes_ch", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
    }

    if (solver_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (solver_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver::C_mono_eq_cr_on_pc_su_tes_ch "
                "iteration to find the cold HTF temperature only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            *m_dot_htf_bal = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    double m_dot_pc_in = mpc_csp_solver->mc_pc_inputs.m_m_dot;              //[kg/hr]
    double m_dot_pc_out = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;     //[kg/hr]
    *m_dot_htf_bal = (m_dot_pc_in - m_dot_pc_out) / m_dot_pc_out;           //[-]
    m_step_pc_su = c_eq.m_step_pc_su;

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on_pc_su_tes_ch::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{  
    // Should not be called directly, only via C_mono_eq_cr_on_pc_su_tes_ch_mdot::operator()(double m_dot_store /*kg/hr*/, double *m_dot_htf_bal /*-*/)
    C_csp_tes::S_csp_tes_outputs tes_outputs;           // the aggregate of the different tes calls
    C_csp_tes::S_csp_tes_outputs tes_outputs_temp;      // output of each tes call, used to update the aggregate
    
    // Solve the tower model with T_htf_cold from the LT HX
    double T_htf_rec_in = T_htf_cold + 273.15;      //[K]
    double P_rec_in = mpc_csp_solver->mc_cr_htf_state_in.m_pres;    //[kPa]
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_rec_in - 273.15;      //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outputs
    double m_dot_rec_out = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double T_htf_rec_out = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15;  //[K]
    double P_rec_out = P_rec_in - mpc_csp_solver->mc_cr_out_solver.m_dP_sf * 100.;  //[kPa]
    double m_dot_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
    double T_store_in = mpc_csp_solver->mc_cr_out_solver.m_T_store_hot + 273.15;   //[K]

    // Charge storage
    // First estimate available charge
    double q_dot_tes_ch_max, m_dot_tes_ch_max, T_tes_cold_ch_max, m_dot_store_ch_max;
    q_dot_tes_ch_max = m_dot_tes_ch_max = T_tes_cold_ch_max = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.charge_avail_est(T_store_in,
        mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_tes_ch_max,
        m_dot_tes_ch_max,
        T_tes_cold_ch_max,
        m_dot_store_ch_max);

    m_dot_tes_ch_max *= 3600.0;     //[kg/hr]

    // Test if particle flow from tower is greater than tes can store, factoring in later discharge
    if (m_dot_store > m_dot_tes_ch_max + m_m_dot_tank) {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    else {
        // update available charge
        mpc_csp_solver->mc_tes.set_max_charge_flow((m_dot_tes_ch_max + m_m_dot_tank) / 3600.);
    }

    double T_cold_tes_K;
    bool ch_solved = mpc_csp_solver->mc_tes.charge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_store_in,
        T_cold_tes_K,
        tes_outputs);

    // Check if TES.charge method solved
    if (!ch_solved) {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -3;
    }

    mpc_csp_solver->mc_tes.use_calc_vals(true);

    // First estimate available discharge in order to updated m_m_dot_tes_dc_max
    double q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est;
    mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_rec_out, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est);

    // Solve the HT HX using a given media discharge (m_m_dot_tank from the outer MEQ)
    // This is a test call (update_calc_vals = false) using the receiver outlet temperature
    // The .calc values are not updated so discharge() is called again later to update them.
    mpc_csp_solver->mc_tes.update_calc_vals(false);
    double T_htf_hx_out, m_dot_hx_out;
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_m_dot_tank / 3600.,
        T_htf_rec_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    mpc_csp_solver->mc_tes.update_calc_vals(true);
    m_dot_hx_out *= 3600.;

    double T_htf_hx_in, m_dot_hx_in, T_htf_pc_in, m_dot_pc_in, P_hx_out;
    if (m_dot_rec_out > m_dot_hx_out) {
        // Needs defocusing if this is the converged state
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    else {
        C_MEQ_cr_on_tes_dc_m_dot_tank c_eq(mpc_csp_solver, T_htf_rec_in, T_htf_rec_out, P_rec_in, P_rec_out, m_dot_rec_out, m_m_dot_tank);
        C_monotonic_eq_solver c_solver(c_eq);

        // Set up solver
        c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

        // Solve for cold temperature
        double T_cold_guess_low = std::min(T_htf_rec_in, T_htf_rec_out) - 273.15 - 10;  //[C]
        double T_cold_guess_high = std::max(T_htf_rec_in, T_htf_rec_out) - 273.15 + 10; //[C]

        double T_cold_solved, tol_solved;
        T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
        int iter_solved = -1;

        mpc_csp_solver->mc_tes.update_calc_vals(false);
        int T_cold_code = 0;
        try
        {
            T_cold_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
        }
        catch (C_csp_exception)
        {
            mpc_csp_solver->mc_tes.use_calc_vals(false);
            mpc_csp_solver->mc_tes.update_calc_vals(true);
            throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
        }
        mpc_csp_solver->mc_tes.update_calc_vals(true);

        if (T_cold_code != C_monotonic_eq_solver::CONVERGED)
        {
            if (T_cold_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
            {
                std::string msg = util::format("At time = %lg C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank "
                    "iteration to find the cold HTF temperature to balance energy between the CR and PC only reached a convergence "
                    "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                    mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
                mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
            }
            else
            {
                mpc_csp_solver->mc_tes.use_calc_vals(false);
                mpc_csp_solver->mc_tes.update_calc_vals(true);
                *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
                return -1;
            }
        }

        T_htf_hx_in = T_cold_solved + 273.15;
        m_dot_hx_in = m_dot_hx_out = c_eq.m_m_dot_htf_out; //[kg/hr] mass flow out of the HX on the field side
        m_dot_pc_in = m_dot_hx_out;     //[kg/hr]
        T_htf_pc_in = T_htf_hx_out = mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out + 273.15;     //[K]
        P_hx_out = P_rec_out * (1. - mpc_csp_solver->mc_tes_outputs.dP_perc / 100.);    //[kPa]
    }

    // call discharge again with calc_vals = true to update the hot and warm tank .calc values
    double T_htf_hot;  //[K] HTF temp out of the HX on the field side
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_m_dot_tank / 3600.,
        T_htf_hx_in,
        T_htf_hot,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_hot_ave = tes_outputs_temp.m_T_hot_ave - 273.15;       //[C]
    tes_outputs.m_m_dot = tes_outputs_temp.m_m_dot;
    tes_outputs.m_W_dot_rhtf_pump = tes_outputs_temp.m_W_dot_rhtf_pump;
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf = tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_hot_ave = tes_outputs_temp.m_T_hot_ave;
    tes_outputs.m_T_hot_final = tes_outputs_temp.m_T_hot_final;
    tes_outputs.dP_perc = tes_outputs_temp.dP_perc;

    // Solve the PC performance using the receiver or HX htf flow rate -> this may be ignored if PC is in controlled startup
    // Need to do this to get back PC T_htf_cold
    // HTF State
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_pc_in - 273.15;   //[C]
    mpc_csp_solver->mc_pc_htf_state_in.m_pres = P_hx_out;   //[kPa]
    // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_pc_in;                         //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode;     //[-]  should be C_csp_power_cycle::STARTUP_CONTROLLED
    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }

    // Get power cycle HTF return state
    // IF STARTUP_CONTROLLED, m_dot_pc_out is calculated, and != m_dot_pc_in
    double T_htf_pc_out = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold + 273.15;       //[K]
    double m_dot_pc_out = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;                 //[kg/hr]
    double P_pc_out = mpc_csp_solver->mc_pc_out_solver.m_P_phx_in * 1000.;              //[kPa]

    // Check for new PC startup timestep here
    m_step_pc_su = mpc_csp_solver->mc_pc_out_solver.m_time_required_su;     //[s] power cycle model returns MIN(time required to completely startup, full timestep duration)

    // Discharge virtual warm tank through LT HX
    //double m_dot_hx_out;    //[kg/s]
    mpc_csp_solver->mc_tes.discharge_full_lt(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_pc_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_cold_ave = tes_outputs_temp.m_T_cold_ave - 273.15;       //[C]
    m_dot_hx_out *= 3600.;      //[kg/hr]
    double P_lthx_out = P_pc_out * (1. - tes_outputs_temp.dP_perc / 100.);           //[kPa]
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf += tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_cold_ave = tes_outputs_temp.m_T_cold_ave;
    tes_outputs.m_T_cold_final = tes_outputs_temp.m_T_cold_final;
    tes_outputs.dP_perc += tes_outputs_temp.dP_perc;

    // Recombine excess mass flow from the power cycle
    double T_htf_rec_in_solved;     //[K]
    if (m_dot_pc_out > m_dot_hx_out) {
        double m_dot_bypassed = m_dot_pc_out - m_dot_hx_out;                                    //[kg/hr]

        // get enthalpy, assume sCO2 HTF
        CO2_state co2_props;
        int prop_error_code = CO2_TP(T_htf_pc_out, P_pc_out, &co2_props);
        double h_in = co2_props.enth;
        double h_out = h_in;
        prop_error_code = CO2_PH(P_lthx_out, h_out, &co2_props);
        double T_htf_bypassed = co2_props.temp; //[K]

        T_htf_rec_in_solved = (T_htf_hx_out * m_dot_hx_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_out + m_dot_bypassed);  // [K]  mix streams to get LT HX outlet temp
    }
    else {
        T_htf_rec_in_solved = T_htf_hx_out;      //[K]
    }

    //Calculate pressure difference (which is not used)
    double diff_P = (P_lthx_out - P_rec_in) / P_rec_in;

    // Set charging inlet/outlet temps to hot/cold ave temps, respectively
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = m_dot_store;                      //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = T_store_in - 273.15;            //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_store_cold_ave;              //[C]

    // Set discharge HTF state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_hx_in;                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_hx_in - 273.15;           //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hx_out - 273.15;         //[C]

    mpc_csp_solver->mc_tes_outputs = tes_outputs;

    //Calculate diff_T_htf_cold
    *diff_T_htf_cold = (T_htf_rec_in_solved - 273.15 - T_htf_cold) / T_htf_cold;        //[-]

    mpc_csp_solver->mc_tes.use_calc_vals(false);
    mpc_csp_solver->mc_tes.update_calc_vals(true);
    return 0;
}

int C_csp_solver::C_mono_eq_pc_target__m_dot::operator()(double m_dot_htf_pc /*kg/hr*/, double *q_dot_pc /*MWt*/)
{
    // Set power cycle HTF inlet state
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = m_T_htf_hot;    //[C]

    // Set power cycle inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_htf_pc;        //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode; //[-]

    // Power cycle performance call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
                                mpc_csp_solver->mc_pc_htf_state_in,
                                mpc_csp_solver->mc_pc_inputs,
                                mpc_csp_solver->mc_pc_out_solver,
                                mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is producing power or model didn't solve
    // Assumes that standby mode always solves
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful && mpc_csp_solver->mc_pc_inputs.m_standby_control == C_csp_power_cycle::ON)
    {
        *q_dot_pc = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    *q_dot_pc = mpc_csp_solver->mc_pc_out_solver.m_q_dot_htf;   //[MWt]
    return 0;
}

int C_csp_solver::C_mono_eq_cr_on_pc_target_tes_ch__T_cold::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Solve the CR
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
                                        mpc_csp_solver->mc_cr_htf_state_in,
                                        m_defocus,
                                        mpc_csp_solver->mc_cr_out_solver,
                                        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is off or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get the calculated receiver mass flow rate
    double m_dot_cr = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;    //[kg/hr]

    // Get the maximum possible mass mass flow rate to power cycle
    double m_dot_pc_max = m_dot_cr;             //[kg/hr]
    if (m_dot_pc_max > mpc_csp_solver->m_m_dot_pc_max)
    {
        m_dot_pc_max = mpc_csp_solver->m_m_dot_pc_max;  //[kg/hr]
    }

    // Try max sending max mass flow rate to power cycle and check calculated thermal power
    C_mono_eq_pc_target__m_dot c_eq(mpc_csp_solver, m_pc_mode, mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot);
    C_monotonic_eq_solver c_solver(c_eq);

    double q_dot_pc_calc = std::numeric_limits<double>::quiet_NaN();    //[MWt]
    int q_dot_pc_code = c_solver.test_member_function(m_dot_pc_max, &q_dot_pc_calc);
    if (q_dot_pc_code != 0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }

    double m_dot_pc_solved = m_dot_pc_max;      //[kg/hr]

    if ( (mpc_csp_solver->mc_pc_out_solver.m_q_dot_htf - m_q_dot_target) / m_q_dot_target > 1.E-3 )
    {   // With maximum possible mass flow rate, the power cycle is over target
        // So iterate on mass flow rate to hit the target
        C_monotonic_eq_solver::S_xy_pair xy_pair_1;
        xy_pair_1.x = m_dot_pc_max;     //[kg/hr]
        xy_pair_1.y = q_dot_pc_calc;    //[MWt]

        // Guess another mass flow rate based on target and first results
        double m_dot_pc_guess = m_q_dot_target / q_dot_pc_calc * m_dot_pc_max;  //[kg/hr]
        q_dot_pc_code = c_solver.test_member_function(m_dot_pc_guess, &q_dot_pc_calc);
        if (q_dot_pc_code != 0)
        {
            *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
            return -3;
        }

        C_monotonic_eq_solver::S_xy_pair xy_pair_2;
        xy_pair_2.x = m_dot_pc_guess;   //[kg/hr]
        xy_pair_2.y = q_dot_pc_calc;    //[MWt]

        c_solver.settings(1.E-3, 50, 0.0, m_dot_pc_max, true);

        // Solve for m_dot_pc
        m_dot_pc_solved = std::numeric_limits<double>::quiet_NaN();
        double tol_solved = std::numeric_limits<double>::quiet_NaN();
        int iter_solved = -1;

        q_dot_pc_code = 0;
        try
        {
            q_dot_pc_code = c_solver.solve(xy_pair_1, xy_pair_2, m_q_dot_target, m_dot_pc_solved, tol_solved, iter_solved);
        }
        catch (C_csp_exception)
        {
            throw(C_csp_exception("C_mono_eq_cr_on_pc_target_tes_ch__T_cold method to calculate the power cycle mass flow rate returned an unexpected exemption"));
        }
        
        if (q_dot_pc_code != C_monotonic_eq_solver::CONVERGED)
        {
            if (q_dot_pc_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) <= 0.1)
            {
                mpc_csp_solver->error_msg = util::format("At time = %lg the iteration to find the power cycle HTF mass flow rate resulting in the target thermal power only reached a convergence "
                    "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                    mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
                mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, mpc_csp_solver->error_msg);
            }
            else
            {
                *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
                return -4;
            }
        }
    }

    // Get power cycle HTF return temperature
    double T_pc_out = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold + 273.15;   //[K]

    // Get maximum charging mass flow rate
    // Knowing the receiver outlet temperature, can calculate the maximum mass flow rate available for charging
    double q_dot_tes_ch_max, m_dot_tes_ch_max, T_tes_cold_ch_max, m_dot_store_ch_max;
    q_dot_tes_ch_max = m_dot_tes_ch_max = T_tes_cold_ch_max = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.charge_avail_est(mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15, 
                                        mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step, 
                                        q_dot_tes_ch_max, 
                                        m_dot_tes_ch_max, 
                                        T_tes_cold_ch_max,
                                        m_dot_store_ch_max);

    m_dot_tes_ch_max *= 3600.0;     //[kg/hr] convert from kg/s

    // Charge storage
    double m_dot_tes = m_dot_cr - m_dot_pc_solved;      //[kg/hr]

    // If amount we want to send to storage is greater than max amount, return
    if (m_dot_tes > m_dot_tes_ch_max)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -4;
    }

    double T_tes_cold_out = std::numeric_limits<double>::quiet_NaN();   //[K]
    bool is_tes_success = mpc_csp_solver->mc_tes.charge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
                                                mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
                                                m_dot_tes / 3600.0,
                                                mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15,
                                                T_tes_cold_out,
                                                mpc_csp_solver->mc_tes_outputs);

    if (!is_tes_success)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -5;
    }

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = m_dot_tes;                                        //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;  //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_tes_cold_out - 273.15;                       //[C] convert from K

    // If not actually discharging (i.e. mass flow rate = 0.0), what should the temperatures be?
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = 0.0;                                                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;   //[C] convert from K
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;   //[C] convert from K

    // Enthalpy balancer (mixer)
    double T_htf_cold_calc = (m_dot_tes*T_tes_cold_out + m_dot_pc_solved*T_pc_out) / m_dot_cr - 273.15;     //[C]

    // Calculate diff_T_rec_in
    *diff_T_htf_cold = (T_htf_cold_calc - T_htf_cold) / T_htf_cold;     //[-]

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on_pc_match_tes_empty_mdot::operator()(double m_dot /*kg/hr*/, double *m_dot_bal /*-*/)
{
    // Converge on an HTF mass flow (m_dot) that balances the particle mass flow (m_dot_bal) between the tower outlet
    //  and the hot tank outlet, when the TES is emptying

    C_mono_eq_cr_on_pc_match_tes_empty c_eq(mpc_csp_solver, m_defocus, m_dot);
    C_monotonic_eq_solver c_solver(c_eq);

    // Set up solver
    double T_cold_min = mpc_csp_solver->m_cycle_T_htf_cold_des - 273.15 - 30.;  //[C]
    double T_cold_max = mpc_csp_solver->m_cycle_T_htf_hot_des - 273.15;         //[C]
    c_solver.settings(1.E-3, 50, T_cold_min, T_cold_max, false);

    // Solve for cold temperature
    double T_cold_guess_low = mpc_csp_solver->m_T_htf_cold_des - 273.15;        //[C], convert from [K]
    double T_cold_guess_high = T_cold_guess_low + 10.0;                         //[C]

    double T_cold_solved, tol_solved;
    T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    int solver_code = 0;
    try
    {
        solver_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_mono_eq_cr_on_pc_match_tes_empty", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
    }

    if (solver_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (solver_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver::C_mono_eq_cr_on_pc_match_tes_empty "
                "iteration to find the cold HTF temperature only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            *m_dot_bal = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    // Calculate and report mass flow rate balance
    double m_dot_rec_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
    double m_dot_tes_store = mpc_csp_solver->mc_tes_outputs.m_m_dot * 3600.;        //[kg/hr]

    *m_dot_bal = (m_dot_rec_store - m_dot_tes_store) / m_dot_rec_store;         //[-]

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on_pc_match_tes_empty::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Should not be called directly, only via C_mono_eq_cr_on_pc_match_tes_empty_mdot(double m_dot /*kg/hr*/, double *m_dot_bal /*-*/)
    C_csp_tes::S_csp_tes_outputs tes_outputs;           // the aggregate of the different tes calls
    C_csp_tes::S_csp_tes_outputs tes_outputs_temp;      // output of each tes call, used to update the aggregate

    mpc_csp_solver->mc_tes.use_calc_vals(true);

    // Solve the tower model with T_htf_cold from the LT HX
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]
    double P_in = mpc_csp_solver->m_P_cold_des;                 //[kPa] use the receiver design inlet pressure
    mpc_csp_solver->mc_cr_htf_state_in.m_pres = P_in;

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outputs
    double m_dot_rec_out = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double T_htf_rec_out = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15;  //[K]
    double P_rec_out = mpc_csp_solver->mc_cr_htf_state_in.m_pres - mpc_csp_solver->mc_cr_out_solver.m_dP_sf * 100.;  //[kPa]
    double m_dot_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
    double T_store_in = mpc_csp_solver->mc_cr_out_solver.m_T_store_hot + 273.15;   //[K]

    // Charge storage
    // First set available charge to that coming from the tower
    mpc_csp_solver->mc_tes.set_max_charge_flow(m_dot_store);
    double T_cold_tes_K;
    mpc_csp_solver->mc_tes.charge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_store_in,
        T_cold_tes_K,
        tes_outputs);

    // Recombine excess mass flow from the power cycle
    double T_htf_hx_in;     //[K]
    double m_dot_rec_in = m_m_dot;                                             //[kg/hr]  (from outer '_mdot' MEQ)
    if (m_dot_rec_in > m_dot_rec_out) {
        double T_htf_rec_in = T_htf_cold + 273.15;                              //[K]
        double P_rec_in = mpc_csp_solver->mc_cr_htf_state_in.m_pres;            //[kPa]
        double m_dot_bypassed = m_dot_rec_in - m_dot_rec_out;                   //[kg/hr]

        // get enthalpy, assume sCO2 HTF
        CO2_state co2_props;
        int prop_error_code = CO2_TP(T_htf_rec_in, P_rec_in, &co2_props);
        double h_in = co2_props.enth;
        double h_out = h_in;
        prop_error_code = CO2_PH(P_rec_out, h_out, &co2_props);
        double T_htf_bypassed = co2_props.temp; //[K]

        T_htf_hx_in = (T_htf_rec_out * m_dot_rec_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_rec_out + m_dot_bypassed);  // [K]  mix streams to get HX inlet temp
    }
    else {
        T_htf_hx_in = T_htf_rec_out;      //[K]
    }
    double m_dot_hx_in = m_dot_rec_in;      //[kg/hr]  don't let receiver create mass


    // Solve the HT HX using a full storage media discharge
    // First estimate available discharge in order to updated m_m_dot_tes_dc_max
    double q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est;
    mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_hx_in, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est);
    double T_htf_hx_out, m_dot_hx_out;
    mpc_csp_solver->mc_tes.discharge_full(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_hx_in,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_hot_ave = tes_outputs_temp.m_T_hot_ave - 273.15;       //[C]

    m_dot_hx_out *= 3600.;  //[kg/hr]
    double P_hx_out = P_rec_out * (1. - tes_outputs_temp.dP_perc / 100.);    //[kPa]

    // Recombine excess mass flow from before the HT HX (since TES is emptying, htf mass flow is constrained)
    double T_htf_pc_in;     //[K]
    if (m_dot_hx_in > m_dot_hx_out) {
        double P_hx_in = P_rec_out;            //[kPa]
        double m_dot_bypassed = m_dot_hx_in - m_dot_hx_out;                   //[kg/hr]

        // get enthalpy, assume sCO2 HTF
        CO2_state co2_props;
        int prop_error_code = CO2_TP(T_htf_hx_in, P_hx_in, &co2_props);
        double h_in = co2_props.enth;
        double h_out = h_in;
        prop_error_code = CO2_PH(P_hx_out, h_out, &co2_props);
        double T_htf_bypassed = co2_props.temp; //[K]

        T_htf_pc_in = (T_htf_hx_out * m_dot_hx_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_out + m_dot_bypassed);  //[K]  mix streams to get PC inlet temp
    }
    else {
        T_htf_pc_in = T_htf_hx_out;      //[K]
    }
    double m_dot_pc_in = m_dot_hx_in;      //[kg/hr]  don't let heat exchanger create mass


    // Solve the PC performance at the receiver htf flow rate
    // Need to do this to get back PC T_htf_cold
    // HTF State
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_pc_in - 273.15;   //[C]
    mpc_csp_solver->mc_pc_htf_state_in.m_pres = P_hx_out;   //[kPa]
    // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_pc_in;                         //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = C_csp_power_cycle::ON;     //[-]
    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }

    // Get power cycle HTF return state
    double T_htf_pc_out = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold + 273.15;       //[K]
    double m_dot_pc_out = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;                 //[kg/hr]
    double P_pc_out = mpc_csp_solver->mc_pc_out_solver.m_P_phx_in * 1000.;              //[kPa]

    // Discharge virtual warm tank through LT HX
    //double m_dot_hx_out;    //[kg/s]
    mpc_csp_solver->mc_tes.discharge_full_lt(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_pc_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_cold_ave = tes_outputs_temp.m_T_cold_ave - 273.15;       //[C]
    m_dot_hx_out *= 3600.;      //[kg/hr]
    double P_lthx_out = P_pc_out * (1. - tes_outputs_temp.dP_perc / 100.);           //[kPa]
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf += tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_cold_ave = tes_outputs_temp.m_T_cold_ave;
    tes_outputs.m_T_cold_final = tes_outputs_temp.m_T_cold_final;
    tes_outputs.dP_perc += tes_outputs_temp.dP_perc;

    // Recombine excess mass flow from the power cycle
    double T_htf_rec_in;     //[K]
    if (m_dot_pc_out > m_dot_hx_out) {
        double m_dot_bypassed = m_dot_pc_out - m_dot_hx_out;                                    //[kg/hr]

        // get enthalpy, assume sCO2 HTF
        CO2_state co2_props;
        int prop_error_code = CO2_TP(T_htf_pc_out, P_pc_out, &co2_props);
        double h_in = co2_props.enth;
        double h_out = h_in;
        prop_error_code = CO2_PH(P_lthx_out, h_out, &co2_props);
        double T_htf_bypassed = co2_props.temp; //[K]

        T_htf_rec_in = (T_htf_hx_out * m_dot_hx_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_out + m_dot_bypassed);  // [K]  mix streams to get LT HX outlet temp
    }
    else {
        T_htf_rec_in = T_htf_hx_out;      //[K]
    }

    //Calculate pressure difference (which is not used)
    double diff_P = (P_lthx_out - P_in) / P_in;

    // Set charging inlet/outlet temps to hot/cold ave temps, respectively
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                  //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = T_store_hot_ave;    //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_store_cold_ave;  //[C]

    // Set discharge HTF state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_hx_in;                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_hx_in - 273.15;           //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hx_out - 273.15;         //[C]

    mpc_csp_solver->mc_tes_outputs = tes_outputs;

    //Calculate diff_T_htf_cold
    *diff_T_htf_cold = ((T_htf_rec_in - 273.15) - T_htf_cold) / T_htf_cold;     //[-]

    mpc_csp_solver->mc_tes.use_calc_vals(false);
    mpc_csp_solver->mc_tes.update_calc_vals(true);
    return 0;
}

int C_csp_solver::C_mono_eq_pc_target__m_dot_fixed_plus_tes_dc::operator()(double m_dot_tes_dc /*kg/hr*/, double *q_dot_pc /*MWt*/)
{
    double T_htf_tes_hot = std::numeric_limits<double>::quiet_NaN();
    bool is_tes_success = mpc_csp_solver->mc_tes.discharge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
                                mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
                                m_dot_tes_dc / 3600.0,
                                m_T_htf_cold + 273.15,
                                T_htf_tes_hot,
                                mpc_csp_solver->mc_tes_outputs);

    if (!is_tes_success)
    {
        *q_dot_pc = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    T_htf_tes_hot -= 273.15;        //[C] convert from K

    // HTF discharging state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_tes_dc;     //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = m_T_htf_cold;   //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_tes_hot; //[C]

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;              //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;    //[C], convert from K
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;  //[C], convert from K

    // Enthalpy balance with 'fixed' (e.g. receiver) hot HTF
    double m_dot_htf_pc = m_dot_tes_dc + m_m_dot_htf_fixed;     //[kg/hr]
    double T_htf_pc_hot = (m_dot_tes_dc*T_htf_tes_hot + m_m_dot_htf_fixed*m_T_htf_fixed_hot) / m_dot_htf_pc;    //[C]

    // Solve power cycle model
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_pc_hot;   //[C]
    // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_htf_pc;        //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode;             //[-]
    // Performance
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is producing power and solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful && mpc_csp_solver->mc_pc_inputs.m_standby_control == C_csp_power_cycle::ON)
    {
        *q_dot_pc = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }

    *q_dot_pc = mpc_csp_solver->mc_pc_out_solver.m_q_dot_htf;   //[MWt]
    return 0;
}

int C_csp_solver::C_mono_eq_pc_target_tes_empty__x_step::operator()(double step /*s*/, double *q_dot_pc /*MWt*/)
{
    double T_htf_tes_hot, m_dot_tes_dc = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.discharge_full_both(step,
                        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
                        m_T_htf_cold + 273.15,
                        T_htf_tes_hot,
                        m_dot_tes_dc,
                        mpc_csp_solver->mc_tes_outputs);

    m_dot_tes_dc *= 3600.0;     //[kg/hr] convert from [kg/s]

    // Set TES HTF states (this needs to be less bulky...)
    // HTF discharging state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_tes_dc;             //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = m_T_htf_cold;           //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_tes_hot - 273.15;    //[C]

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                                  //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;    //[C] convert from K
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;  //[C] convert from K

    // Report thermal power to cycle
    *q_dot_pc = mpc_csp_solver->mc_tes_outputs.m_q_dot_dc_to_htf;       //[MWt]

    return 0;
}

int C_csp_solver::C_mono_eq_pc_target_tes_empty__T_cold::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Clear public member data
    m_step = std::numeric_limits<double>::quiet_NaN();

    // First, get the maximum possible mass flow rate from TES discharge 
    //   ... using the guess value for the TES cold inlet temperature
    double T_htf_tes_hot, m_dot_htf_full_ts;
    T_htf_tes_hot = m_dot_htf_full_ts = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.discharge_full_both(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
                            mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
                            T_htf_cold + 273.15,
                            T_htf_tes_hot,
                            m_dot_htf_full_ts,
                            mpc_csp_solver->mc_tes_outputs);

    // Know the mass flow rate for full discharge and the duration of the full timestep
    // so can calculate max TES discharge MASS
    double mass_tes_max = m_dot_htf_full_ts*mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step;     //[kg]
    
    C_mono_eq_pc_target_tes_empty__x_step c_eq(mpc_csp_solver, T_htf_cold);
    C_monotonic_eq_solver c_solver(c_eq);

    double time_max = mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step;       //[s]

    // Check whether the mass flow rate at the full timestep is less than the minimum PC mass flow rate
    if (m_dot_htf_full_ts*3600.0 < mpc_csp_solver->m_m_dot_pc_min)
    {
        // If it is, then calculate the time required to deplete storage at the minimum mass flow rate
        time_max = mass_tes_max / (mpc_csp_solver->m_m_dot_pc_min / 3600.0);    //[s]
    }
    // ELSE:
    // At full timestep the discharge mass flow is greater than minimum,
    //    so iterate on thermal power

    // Now use this timestep to calculate the thermal power to the power cycle  
    double q_dot_pc_m_dot_min = std::numeric_limits<double>::quiet_NaN();
    int eq_code = c_solver.test_member_function(time_max, &q_dot_pc_m_dot_min);
    if (eq_code != 0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    solve_pc(time_max);

    *diff_T_htf_cold = (mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold - T_htf_cold) / T_htf_cold;

    // Compare solved q_dot_pc to target value
    // If it is greater, than exit because we can't decrease the mass flow rate
    if ((q_dot_pc_m_dot_min - m_q_dot_pc_target) / m_q_dot_pc_target > -1.E-3)
    {
        m_step = time_max;
        return 0;
    }

    // Calculate the time required to deplete storage at the minimum mass flow rate
    // ... Be conservative with 0.75 multiplier
    double time_min = 0.75*std::max(mass_tes_max / (mpc_csp_solver->m_m_dot_pc_max / 3600.0), 0.001);   //[s]

    // Set up solver to iterate on timestep to achieve q_dot_pc_target
    c_solver.settings(1.E-3, 50, time_min, time_max, true);
    
    // Guess time required to deplete storage while delivering thermal power requirements to PC
    double time_guess_q_dot_high = mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step*
                                    (mpc_csp_solver->mc_tes_outputs.m_q_dot_dc_to_htf / m_q_dot_pc_target);     //[s]
    
    time_guess_q_dot_high = std::max(1.02*time_min, std::min(time_guess_q_dot_high, 0.98*time_max));        //[s]
    double time_guess_q_dot_low = std::max(1.01*time_min, 0.85*time_guess_q_dot_high);      //[s]

    double time_solved, tol_solved;
    time_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();    //[s]
    int iter_solved = -1;

    int solver_code = 0;

    try
    {
        solver_code = c_solver.solve(time_guess_q_dot_high, time_guess_q_dot_low, m_q_dot_pc_target, time_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception("C_mono_eq_pc_target_tes_empty__x_step method to calculate the time step required to empty TES at the target thermal power returned an exemption"));
    }

    if (solver_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (solver_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) <= 0.1)
        {
            mpc_csp_solver->error_msg = util::format("At time = %lg the iteration to find the time step resulting in emptying TES at the target thermal power only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, mpc_csp_solver->error_msg);
        }
        else
        {
            *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
            return -2;
        }
    }

    solve_pc(time_solved);

    *diff_T_htf_cold = (mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold - T_htf_cold) / T_htf_cold;
    m_step = time_solved;   //[s]

    return 0;
}

void C_csp_solver::C_mono_eq_pc_target_tes_empty__T_cold::solve_pc(double step /*s*/)
{
    // Solve PC model
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out; //[C]

    mpc_csp_solver->mc_pc_inputs.m_m_dot = mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot;     //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = C_csp_power_cycle::ON;

    // Set new local timestep
    C_csp_solver_sim_info temp_sim_info = mpc_csp_solver->mc_kernel.mc_sim_info;
    temp_sim_info.ms_ts.m_step = step;  //[s]
    temp_sim_info.ms_ts.m_time = mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time -
        mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step +
        step;   //[s]

    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        temp_sim_info);
}

int C_csp_solver::C_mono_eq_cr_on_pc_target_tes_ch_mdot::operator()(double m_dot_tank /*kg/hr*/, double *diff_pc_target /*-*/)
{
    // Converge on the hot tank mass that results in the power cycle operating at its target power

    C_mono_eq_cr_on_pc_su_tes_ch c_eq(mpc_csp_solver, m_pc_mode, m_defocus, m_dot_tank);
    C_monotonic_eq_solver c_solver(c_eq);

    // Set up solver
    c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

    // Solve for cold temperature
    double T_cold_guess_low = mpc_csp_solver->m_T_htf_cold_des - 273.15;        //[C], convert from [K]
    double T_cold_guess_high = T_cold_guess_low + 10.0;                         //[C]

    double T_cold_solved, tol_solved;
    T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    // Solve for hot tank mass flow
    int solver_code = 0;
    try
    {
        solver_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0, T_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_mono_eq_cr_on_pc_target_tes_ch_mdot", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
    }

    if (solver_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (solver_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver::C_mono_eq_cr_on_pc_target_tes_ch_mdot "
                "iteration to find the cold HTF temperature only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            *diff_pc_target = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    if (m_pc_mode == C_csp_power_cycle::ON) {
        *diff_pc_target = (mpc_csp_solver->mc_pc_out_solver.m_P_cycle - m_W_dot_pc_target) / m_W_dot_pc_target;         //[-]
    }
    else {
        *diff_pc_target = (mpc_csp_solver->mc_pc_out_solver.m_q_dot_htf - m_q_dot_pc_target) / m_q_dot_pc_target;           //[-]
    }

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on_pc_target_tes_dc_mdot::operator()(double m_dot_tank /*kg/hr*/, double *diff_pc_target /*-*/)
{
    // Converge on the hot tank mass that results in the power cycle operating at its target power

    C_mono_eq_cr_on_pc_target_tes_dc c_eq(mpc_csp_solver, m_pc_mode, m_defocus, m_dot_tank);
    C_monotonic_eq_solver c_solver(c_eq);
         
    // Set up solver
    double T_cold_min = mpc_csp_solver->m_cycle_T_htf_cold_des - 273.15 - 10.;
    double T_cold_max = mpc_csp_solver->m_T_htf_cold_des - 273.15 + 40.;
    c_solver.settings(2.E-3, 50, T_cold_min, T_cold_max, false);

    // Solve for cold temperature
    double T_cold_guess_low = mpc_csp_solver->m_T_htf_cold_des - 273.15;        //[C], convert from [K]
    double T_cold_guess_high = T_cold_guess_low + 10.0;                         //[C]

    double T_cold_solved, tol_solved;
    T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    // Solve for hot tank mass flow
    int solver_code = 0;
    try
    {
        solver_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0, T_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_mono_eq_cr_on_pc_target_tes_dc_mdot", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
    }

    if (solver_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (solver_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver::C_mono_eq_cr_on_pc_target_tes_dc_mdot "
                "iteration to find the cold HTF temperature only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            *diff_pc_target = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }


    *diff_pc_target = (mpc_csp_solver->mc_pc_out_solver.m_P_cycle - m_W_dot_pc_target) / m_W_dot_pc_target;         //[-]

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on_pc_target_tes_dc::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Should not be called directly, only via C_mono_eq_cr_on_pc_target_tes_dc_mdot(double m_dot_tank /*kg/hr*/, double *q_dot_pc /*MWt*/)
    C_csp_tes::S_csp_tes_outputs tes_outputs;           // the aggregate of the different tes calls
    C_csp_tes::S_csp_tes_outputs tes_outputs_temp;      // output of each tes call, used to update the aggregate

    // Solve the tower model with T_htf_cold from the LT HX
    double T_htf_rec_in = T_htf_cold + 273.15;      //[K]
    double P_rec_in = mpc_csp_solver->mc_cr_htf_state_in.m_pres;    //[kPa]
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_rec_in - 273.15;      //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outputs
    double m_dot_rec_out = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double T_htf_rec_out = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15;  //[K]
    double P_rec_out = P_rec_in - mpc_csp_solver->mc_cr_out_solver.m_dP_sf * 100.;  //[kPa]
    double m_dot_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
    double T_store_in = mpc_csp_solver->mc_cr_out_solver.m_T_store_hot + 273.15;   //[K]

    // Charge storage
    // First estimate available charge
    double q_dot_tes_ch_max, m_dot_tes_ch_max, T_tes_cold_ch_max, m_dot_store_ch_max;
    q_dot_tes_ch_max = m_dot_tes_ch_max = T_tes_cold_ch_max = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.charge_avail_est(T_store_in,
        mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_tes_ch_max,
        m_dot_tes_ch_max,
        T_tes_cold_ch_max,
        m_dot_store_ch_max);

    m_dot_tes_ch_max *= 3600.0;     //[kg/hr]

    // Test if particle flow from tower is greater than tes can store, factoring in later discharge
    if (m_dot_store > m_dot_tes_ch_max + m_m_dot_tank) {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    else {
        // update available charge
        mpc_csp_solver->mc_tes.set_max_charge_flow((m_dot_tes_ch_max + m_m_dot_tank) / 3600.);
    }

    double T_cold_tes_K;
    bool ch_solved = mpc_csp_solver->mc_tes.charge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_dot_store / 3600.,
        T_store_in,
        T_cold_tes_K,
        tes_outputs);

    // Check if TES.charge method solved
    if (!ch_solved) {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -3;
    }

    mpc_csp_solver->mc_tes.use_calc_vals(true);

    // First estimate available discharge in order to updated m_m_dot_tes_dc_max
    double q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est;
    mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_rec_out, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est);

    // Solve the HT HX using a given media discharge (m_m_dot_tank from the outer MEQ)
    // This is a test call (update_calc_vals = false) using the receiver outlet temperature
    // The .calc values are not updated so discharge() is called again later to update them.
    mpc_csp_solver->mc_tes.update_calc_vals(false);     // do not update calc values due to following iterations (which are within larger iterations)
    double T_htf_hx_out, m_dot_hx_out;
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_m_dot_tank / 3600.,
        T_htf_rec_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    mpc_csp_solver->mc_tes.update_calc_vals(true);
    m_dot_hx_out *= 3600.;

    double T_htf_hx_in, m_dot_hx_in, T_htf_pc_in, m_dot_pc_in, P_hx_out;
    if (m_dot_rec_out > m_dot_hx_out) {
        // Needs defocusing if this is the converged state
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    else {
        C_MEQ_cr_on_tes_dc_m_dot_tank c_eq(mpc_csp_solver, T_htf_rec_in, T_htf_rec_out, P_rec_in, P_rec_out, m_dot_rec_out, m_m_dot_tank);
        C_monotonic_eq_solver c_solver(c_eq);

        // Set up solver
        c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

        // Solve for cold temperature
        double T_cold_guess_low = std::min(T_htf_rec_in, T_htf_rec_out) - 273.15 - 10;  //[C]
        double T_cold_guess_high = std::max(T_htf_rec_in, T_htf_rec_out) - 273.15 + 10; //[C]

        double T_cold_solved, tol_solved;
        T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
        int iter_solved = -1;

        mpc_csp_solver->mc_tes.update_calc_vals(false);
        int T_cold_code = 0;
        try
        {
            T_cold_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
        }
        catch (C_csp_exception)
        {
            mpc_csp_solver->mc_tes.use_calc_vals(false);
            mpc_csp_solver->mc_tes.update_calc_vals(true);
            throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
        }
        mpc_csp_solver->mc_tes.update_calc_vals(true);

        if (T_cold_code != C_monotonic_eq_solver::CONVERGED)
        {
            if (T_cold_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
            {
                std::string msg = util::format("At time = %lg C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank "
                    "iteration to find the cold HTF temperature to balance energy between the CR and PC only reached a convergence "
                    "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                    mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
                mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
            }
            else
            {
                mpc_csp_solver->mc_tes.use_calc_vals(false);
                mpc_csp_solver->mc_tes.update_calc_vals(true);
                *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
                return -1;
            }
        }

        T_htf_hx_in = T_cold_solved + 273.15;
        m_dot_hx_in = m_dot_hx_out = c_eq.m_m_dot_htf_out; //[kg/hr] mass flow out of the HX on the field side
        m_dot_pc_in = m_dot_hx_out;     //[kg/hr]
        T_htf_pc_in = T_htf_hx_out = mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out + 273.15;     //[K]
        P_hx_out = P_rec_out * (1. - mpc_csp_solver->mc_tes_outputs.dP_perc / 100.);    //[kPa]
    }

    // call discharge again with calc_vals = true to update the hot and warm tank .calc values
    double T_htf_hot;  //[K] HTF temp out of the HX on the field side
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_m_dot_tank / 3600.,
        T_htf_hx_in,
        T_htf_hot,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_hot_ave = tes_outputs_temp.m_T_hot_ave - 273.15;       //[C]
    tes_outputs.m_m_dot = tes_outputs_temp.m_m_dot;
    tes_outputs.m_W_dot_rhtf_pump = tes_outputs_temp.m_W_dot_rhtf_pump;
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf = tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_hot_ave = tes_outputs_temp.m_T_hot_ave;
    tes_outputs.m_T_hot_final = tes_outputs_temp.m_T_hot_final;
    tes_outputs.dP_perc = tes_outputs_temp.dP_perc;

    // Solve the PC performance at the receiver htf flow rate
    // Need to do this to get back PC T_htf_cold
    // HTF State
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_pc_in - 273.15;   //[C]
    mpc_csp_solver->mc_pc_htf_state_in.m_pres = P_hx_out;   //[kPa]
    // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_pc_in;                         //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = C_csp_power_cycle::ON;     //[-]
    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful)
    {
        mpc_csp_solver->mc_tes.use_calc_vals(false);
        mpc_csp_solver->mc_tes.update_calc_vals(true);
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }

    // Get power cycle HTF return state
    double T_htf_pc_out = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold + 273.15;       //[K]
    double m_dot_pc_out = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;                 //[kg/hr]
    double P_pc_out = mpc_csp_solver->mc_pc_out_solver.m_P_phx_in * 1000.;              //[kPa]

    // Discharge virtual warm tank through LT HX
    //double m_dot_hx_out;    //[kg/s]
    mpc_csp_solver->mc_tes.discharge_full_lt(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_pc_out,
        T_htf_hx_out,
        m_dot_hx_out,
        tes_outputs_temp);
    double T_store_cold_ave = tes_outputs_temp.m_T_cold_ave - 273.15;       //[C]
    m_dot_hx_out *= 3600.;      //[kg/hr]
    double P_lthx_out = P_pc_out * (1. - tes_outputs_temp.dP_perc / 100.);           //[kPa]
    // don't double count heater power and thermal losses, already accounted for during charging
    tes_outputs.m_q_dot_dc_to_htf += tes_outputs_temp.m_q_dot_dc_to_htf;
    tes_outputs.m_T_cold_ave = tes_outputs_temp.m_T_cold_ave;
    tes_outputs.m_T_cold_final = tes_outputs_temp.m_T_cold_final;
    tes_outputs.dP_perc += tes_outputs_temp.dP_perc;

    // Recombine excess mass flow from the power cycle
    double T_htf_rec_in_solved;     //[K]
    if (m_dot_pc_out > m_dot_hx_out) {
        double m_dot_bypassed = m_dot_pc_out - m_dot_hx_out;                                    //[kg/hr]

        // get enthalpy, assume sCO2 HTF
        CO2_state co2_props;
        int prop_error_code = CO2_TP(T_htf_pc_out, P_pc_out, &co2_props);
        double h_in = co2_props.enth;
        double h_out = h_in;
        prop_error_code = CO2_PH(P_lthx_out, h_out, &co2_props);
        double T_htf_bypassed = co2_props.temp; //[K]

        T_htf_rec_in_solved = (T_htf_hx_out * m_dot_hx_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_out + m_dot_bypassed);  // [K]  mix streams to get LT HX outlet temp
    }
    else {
        T_htf_rec_in_solved = T_htf_hx_out;      //[K]
    }

    //Calculate pressure difference (which is not used)
    double diff_P = (P_lthx_out - P_rec_in) / P_rec_in;

    // Set charging inlet/outlet temps to hot/cold ave temps, respectively
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                  //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = T_store_hot_ave;    //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_store_cold_ave;  //[C]

    // Set discharge HTF state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_hx_in;                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_hx_in - 273.15;           //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hx_out - 273.15;         //[C]

    mpc_csp_solver->mc_tes_outputs = tes_outputs;

    //Calculate diff_T_htf_cold
    *diff_T_htf_cold = (T_htf_rec_in_solved - 273.15 - T_htf_cold) / T_htf_cold;        //[-]

    mpc_csp_solver->mc_tes.use_calc_vals(false);
    mpc_csp_solver->mc_tes.update_calc_vals(true);
    return 0;
}

int C_csp_solver::C_mono_eq_cr_on__pc_match__tes_full::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Solve the receiver model with T_htf_cold
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
                                        mpc_csp_solver->mc_cr_htf_state_in,
                                        m_defocus,
                                        mpc_csp_solver->mc_cr_out_solver,
                                        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outputs
    double m_dot_receiver = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;  //[kg/hr]
    double T_htf_rec_hot = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;       //[C]

    // Solve TES for *full* charge
    double T_htf_tes_cold, m_dot_tes;
    T_htf_tes_cold = m_dot_tes = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.charge_full(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
                            mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
                            T_htf_rec_hot + 273.15,
                            T_htf_tes_cold,
                            m_dot_tes,
                            mpc_csp_solver->mc_tes_outputs);

    T_htf_tes_cold -= 273.15;   //[C] convert from K
    m_dot_tes *= 3600.0;        //[kg/hr] convert from kg/s

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = m_dot_tes;                //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;      //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_htf_tes_cold;        //[C]

    // If not actually discharging (i.e. mass flow rate = 0.0), what should the temperatures be?
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = 0.0;                                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;   //[C] convert from K
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;   //[C] convert from K

    double m_dot_pc = m_dot_receiver - m_dot_tes;   //[kg/hr]

    if (m_dot_pc < 0.0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }
    // Do we want to check PC HTF mass flow rate here?
    //if (m_dot_pc > mpc_csp_solver->m_m_dot_pc_max)
    //{
    //  *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
    //  return -3;
    //}

    // Solve the PC performance
        // HTF State
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;  //[C]
        // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_pc;                //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode;     //[-]
        // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
                            mpc_csp_solver->mc_pc_htf_state_in,
                            mpc_csp_solver->mc_pc_inputs,
                            mpc_csp_solver->mc_pc_out_solver,
                            mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is producing power and solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful && mpc_csp_solver->mc_pc_inputs.m_standby_control == C_csp_power_cycle::ON)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -3;
    }

    // Get power cycle HTF return temperature...
    double T_htf_pc_cold = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold;       //[C]

    // Enthalpy balance to get T_htf_cold_calc
    double T_htf_cold_calc = (m_dot_tes*T_htf_tes_cold + m_dot_pc*T_htf_pc_cold) / m_dot_receiver;  //[C]

    //Calculate diff_T_htf_cold
    *diff_T_htf_cold = (T_htf_cold_calc - T_htf_cold) / T_htf_cold;     //[-]

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on__pc_max_m_dot__tes_full::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Solve the receiver model with T_htf_cold
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outputs
    double m_dot_receiver = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;  //[kg/hr]
    double T_htf_rec_hot = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;       //[C]

    // First, call power cycle, because if it's in startup mode, we need the new timestep
    // Solve the PC performance at MAX PC HTF FLOW RATE
    // Need to do this to get back PC T_htf_cold
    // HTF State
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;  //[C]
    // Inputs
    if (m_pc_mode == C_csp_power_cycle::ON)
    {
        mpc_csp_solver->mc_pc_inputs.m_m_dot = mpc_csp_solver->m_m_dot_pc_max;      //[kg/hr]
    }
    else
    {
        mpc_csp_solver->mc_pc_inputs.m_m_dot = 0.0;             //[kg/hr]
    }
    mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode;     //[-]
    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is producing power and solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful && mpc_csp_solver->mc_pc_inputs.m_standby_control == C_csp_power_cycle::ON)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -2;
    }

    // Get power cycle HTF return temperature...
    double T_htf_pc_cold = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold;       //[C]
    double m_dot_pc = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;             //[kg/hr]

    // Check for new timestep
    double step_calc = mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step;      //[s]
    if (m_pc_mode == C_csp_power_cycle::STARTUP_CONTROLLED)
    {
        if (mpc_csp_solver->mc_pc_out_solver.m_time_required_su < step_calc)
        {
            step_calc = mpc_csp_solver->mc_pc_out_solver.m_time_required_su;
        }
    }

    // Solve TES for *full* charge
    double T_htf_tes_cold, m_dot_tes;
    T_htf_tes_cold = m_dot_tes = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.charge_full(step_calc,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_rec_hot + 273.15,
        T_htf_tes_cold,
        m_dot_tes,
        mpc_csp_solver->mc_tes_outputs);

    T_htf_tes_cold -= 273.15;   //[C] convert from K
    m_dot_tes *= 3600.0;        //[kg/hr] convert from kg/s

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = m_dot_tes;                //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;      //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_htf_tes_cold;        //[C]

    // If not actually discharging (i.e. mass flow rate = 0.0), what should the temperatures be?
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = 0.0;                                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;   //[C] convert from K
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;   

    // Enthalpy balance to get T_htf_cold_calc
        // So mass flow rate is probably *not* balanced here, because we're fixing both TES and PC
    double T_htf_cold_calc = (m_dot_tes*T_htf_tes_cold + m_dot_pc*T_htf_pc_cold) / (m_dot_tes + m_dot_pc);  //[C]

    //Calculate diff_T_htf_cold
    *diff_T_htf_cold = (T_htf_cold_calc - T_htf_cold) / T_htf_cold;     //[-]

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on__pc_target__tes_full__defocus::operator()(double defocus /*-*/, double *q_dot_pc /*MWt*/)
{
    int T_htf_cold_code = mpc_csp_solver->solver_cr_on__pc_match__tes_full(m_pc_mode, defocus);

    if (T_htf_cold_code != 0)
    {
        *q_dot_pc = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    *q_dot_pc = mpc_csp_solver->mc_pc_out_solver.m_q_dot_htf;   //[MWt]

    return 0;
}

int C_csp_solver::C_mono_eq_cr_on__pc_m_dot_max__tes_full_defocus::operator()(double defocus /*-*/, double *m_dot_bal /*-*/)
{
    C_mono_eq_cr_on__pc_max_m_dot__tes_full c_eq(mpc_csp_solver, m_pc_mode, defocus);
    C_monotonic_eq_solver c_solver(c_eq);

    // Set up solver
    c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

    // Solve for cold temperature
    double T_cold_guess_low = mpc_csp_solver->m_T_htf_pc_cold_est;  //[C]
    double T_cold_guess_high = T_cold_guess_low + 10.0;     //[C]

    double T_cold_solved, tol_solved;
    T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    int T_cold_code = 0;
    try
    {
        T_cold_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver:::solver_pc_fixed__tes_dc failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
    }

    if (T_cold_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (T_cold_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver:::solver_pc_fixed__tes_dc failed "
                "iteration to find the cold HTF temperature to balance energy between the TES and PC only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            *m_dot_bal = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    // Calculate and report mass flow rate balance
    double m_dot_rec = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double m_dot_pc = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;         //[kg/hr]
    double m_dot_tes = mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot;         //[kg/hr]

    *m_dot_bal = (m_dot_rec - (m_dot_pc + m_dot_tes)) / m_dot_rec;          //[-]
    return 0;
}

int C_csp_solver::C_mono_eq_cr_on__pc_match_m_dot_ceil__tes_full::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Solve the receiver model with T_htf_cold
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outputs
    double m_dot_receiver = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;  //[kg/hr]
    double T_htf_rec_hot = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;       //[C]

    // Solve TES for *full* charge
    double T_htf_tes_cold, m_dot_tes;
    T_htf_tes_cold = m_dot_tes = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.charge_full(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_rec_hot + 273.15, 
        T_htf_tes_cold,
        m_dot_tes,
        mpc_csp_solver->mc_tes_outputs);

    T_htf_tes_cold -= 273.15;   //[C] convert from K
    m_dot_tes *= 3600.0;        //[kg/hr] convert from kg/s

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = m_dot_tes;                //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;      //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_htf_tes_cold;        //[C]

    // If not actually discharging (i.e. mass flow rate = 0.0), what should the temperatures be?
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = 0.0;                                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;   //[C] convert from K
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;

    // Solve the PC performance at MIN( calculated pc mass flow rate, max allowable pc mass flow rate)
    double m_dot_pc = std::min(m_dot_receiver - m_dot_tes, mpc_csp_solver->m_m_dot_pc_max);     //[kg/hr]
    // Need to do this to get back PC T_htf_cold
    // HTF State
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;  //[C]
    // Inputs
    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_pc;                //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode;     //[-]
    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check that power cycle is producing power and solving without errors
    if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful && mpc_csp_solver->mc_pc_inputs.m_standby_control == C_csp_power_cycle::ON)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -4;
    }

    // Get power cycle HTF return temperature...
    double T_htf_pc_cold = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold;       //[C]

    // Enthalpy balance to get T_htf_cold_calc
    // So mass flow rate is probably *not* balanced here, because we're fixing both TES and PC
    double T_htf_cold_calc = (m_dot_tes*T_htf_tes_cold + m_dot_pc*T_htf_pc_cold) / (m_dot_tes + m_dot_pc);  //[C]

    //Calculate diff_T_htf_cold
    *diff_T_htf_cold = (T_htf_cold_calc - T_htf_cold) / T_htf_cold;     //[-]

    return 0;
}

int C_csp_solver::C_MEQ_cr_on__pc_m_dot_max__tes_off__defocus::operator()(double defocus /*-*/, double *m_dot_bal /*-*/)
{
    C_MEQ_cr_on__pc_max_m_dot__tes_off__T_htf_cold c_eq(mpc_csp_solver, m_pc_mode, defocus);
    C_monotonic_eq_solver c_solver(c_eq);
    double T_cold_guess_low, T_cold_guess_high;
    
    double T_cold_guess = mpc_csp_solver->m_cycle_T_htf_cold_des - 273.15 - 10.;            //[C]

    double diff_T_htf_cold_calc = std::numeric_limits<double>::quiet_NaN();
    int T_htf_code = c_solver.test_member_function(T_cold_guess, &diff_T_htf_cold_calc);
    if (T_htf_code != 0) {
        *m_dot_bal = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }
    
    if (diff_T_htf_cold_calc < 0) { // if the first guess is not cold enough
        return -1;  // don't anticipate this happening
    }
    else { // it's cold enough
        T_cold_guess_low = T_cold_guess;    //[C]
        T_cold_guess_high = T_cold_guess + 40.;     //[C]
    }


    // Set up solver
    c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

    double T_cold_solved, tol_solved;
    T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    int T_cold_code = 0;
    try
    {
        T_cold_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::cr_on__pc_m_dot_max__tes_off__defocus failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
    }

    if (T_cold_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (T_cold_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver::cr_on__pc_m_dot_max__tes_off__defocus "
                "iteration to find the cold HTF temperature to balance energy between the CR and PC only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            *m_dot_bal = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    // Calculate and report mass flow rate balance
    double m_dot_rec_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
    double m_dot_tes_store = mpc_csp_solver->mc_tes_outputs.m_m_dot * 3600.;        //[kg/hr]

    *m_dot_bal = (m_dot_rec_store - m_dot_tes_store) / m_dot_rec_store;         //[-]

    return 0;
}

int C_csp_solver::C_MEQ_cr_on__pc_max_m_dot__tes_off__T_htf_cold::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // This MEQ doesn't seem valid. You can't constrain the PC mass flow and temperature and TES discharge. Constraining the mass
    //  flow and PC inlet temperature also constrains the power. And if you're not defocusing, you can't adjust the power because
    //  you can't adjust the TES discharge.

    // More practically, if you constrain the PC mass flow you are also contraining the HT HX sco2 mass flow. And you're constraining
    //  the particle mass flow (because it's steady state), so you're overconstrained unless you let the HX operate non-optimally,
    //  but then you're not hitting the target PC inlet temperature.
    

    *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
    return -1;
    
 //   mpc_csp_solver->mc_tes.use_calc_vals(true);

 //   double m_dot_pc = mpc_csp_solver->m_m_dot_pc_max;   //[kg/hr]

    //// Solve the tower model with T_htf_cold from the LT HX
 //   double T_htf_rec_in = T_htf_cold + 273.15;      //[K]
 //   double P_rec_in = mpc_csp_solver->mc_cr_htf_state_in.m_pres;    //[kPa]
 //   mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_rec_in - 273.15;        //[C]
    ////mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]
 ////   double P_in = mpc_csp_solver->m_P_cold_des;   //[kPa] use the receiver design inlet pressure
 ////   mpc_csp_solver->mc_cr_htf_state_in.m_pres = P_in;

    //mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
    //  mpc_csp_solver->mc_cr_htf_state_in,
    //  m_defocus,
    //  mpc_csp_solver->mc_cr_out_solver,
    //  mpc_csp_solver->mc_kernel.mc_sim_info);

    //// Check if receiver is OFF or didn't solve
    //if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    //{
    //  mpc_csp_solver->mc_tes.use_calc_vals(false);
    //  mpc_csp_solver->mc_tes.update_calc_vals(true);
    //  *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
    //  return -1;
    //}

    //// Get receiver HTF outputs
    //double m_dot_rec_out = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot; //[kg/hr]
    //double T_htf_rec_out = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15;    //[K]
 //   double P_rec_out = mpc_csp_solver->mc_cr_htf_state_in.m_pres - mpc_csp_solver->mc_cr_out_solver.m_dP_sf * 100.;  //[kPa]
 //   double m_dot_store = mpc_csp_solver->mc_cr_out_solver.m_m_dot_store_tot;    //[kg/hr]
 //   double T_store_in = mpc_csp_solver->mc_cr_out_solver.m_T_store_hot + 273.15;   //[K]

 //   if (m_dot_rec_out > m_dot_pc) {
 //       // Might be able to take some flow from the tower exit instead
 //       mpc_csp_solver->mc_tes.use_calc_vals(false);
 //       mpc_csp_solver->mc_tes.update_calc_vals(true);
 //       *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
 //       return -1;
 //   }

 //   // Charge storage
 //   // First estimate available charge
 //   double q_dot_tes_ch_max, m_dot_tes_ch_max, T_tes_cold_ch_max;
 //   q_dot_tes_ch_max = m_dot_tes_ch_max = T_tes_cold_ch_max = std::numeric_limits<double>::quiet_NaN();
 //   mpc_csp_solver->mc_tes.charge_avail_est(T_store_in,
 //       mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
 //       q_dot_tes_ch_max,
 //       m_dot_tes_ch_max,
 //       T_tes_cold_ch_max);

 //   m_dot_tes_ch_max *= 3600.0;       //[kg/hr]

 //   // Test if particle flow from tower is greater than tes can store, factoring in later discharge
 //   if (m_dot_store > m_dot_tes_ch_max + m_dot_store) {
 //       mpc_csp_solver->mc_tes.use_calc_vals(false);
 //       mpc_csp_solver->mc_tes.update_calc_vals(true);
 //       *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
 //       return -1;
 //   }
 //   else {
 //       // update available charge
 //       mpc_csp_solver->mc_tes.set_max_charge_flow((m_dot_tes_ch_max + m_dot_store) / 3600.);
 //   }

 //   double T_cold_tes_K;
 //   bool ch_solved = mpc_csp_solver->mc_tes.charge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
 //                                   mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
 //                                   m_dot_store / 3600.,
 //                                   T_store_in,
 //                                   T_cold_tes_K,
 //                                   mpc_csp_solver->mc_tes_outputs);

 //   // Check if TES.charge method solved
 //   if (!ch_solved) {
 //       mpc_csp_solver->mc_tes.use_calc_vals(false);
 //       mpc_csp_solver->mc_tes.update_calc_vals(true);
 //       *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
 //       return -3;
 //   }

 //   // First estimate available discharge in order to updated m_m_dot_tes_dc_max
 //   double q_dot_dc_est, m_dot_field_est, T_hot_field_est;
 //   mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_rec_out, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
 //       q_dot_dc_est, m_dot_field_est, T_hot_field_est);

 //   // Solve the HT HX using a given media discharge (m_dot_store from the receiver as TES is off (SS))
 //   // This is a test call (update_calc_vals = false) using the receiver outlet temperature
 //   // The .calc values are not updated so discharge() is called again later to update them.
 //   mpc_csp_solver->mc_tes.update_calc_vals(false);     // do not update calc values due to following iterations (which are within larger iterations)
 //   double T_htf_hx_out, m_dot_hx_out;
 //   mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
 //       mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
 //       m_dot_store / 3600.,
 //       T_htf_rec_out,
 //       T_htf_hx_out,
 //       mpc_csp_solver->mc_tes_outputs);
 //   mpc_csp_solver->mc_tes.update_calc_vals(true);
 //   m_dot_hx_out = mpc_csp_solver->mc_tes_outputs.m_m_dot * 3600.;

 //   double T_htf_hx_in, m_dot_hx_in, T_htf_pc_in, m_dot_pc_in, P_hx_out;
 //   if (m_dot_rec_out > m_dot_hx_out) {
 //       //Needs defocusing if this is the converged state
 //       mpc_csp_solver->mc_tes.use_calc_vals(false);
 //       mpc_csp_solver->mc_tes.update_calc_vals(true);
 //       *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
 //       return -1;
 //   }
 //   else {
 //       C_MEQ_cr_on_tes_dc_m_dot_tank c_eq(mpc_csp_solver, T_htf_rec_in, T_htf_rec_out, P_in, P_rec_out, m_dot_rec_out, m_m_dot_tank);
 //       C_monotonic_eq_solver c_solver(c_eq);

 //       // Set up solver
 //       c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

 //       // Solve for cold temperature
 //       double T_cold_guess_low = std::min(T_htf_rec_in, T_htf_rec_out) - 273.15 - 10;    //[C]
 //       double T_cold_guess_high = std::max(T_htf_rec_in, T_htf_rec_out) - 273.15 + 10;   //[C]

 //       double T_cold_solved, tol_solved;
 //       T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
 //       int iter_solved = -1;

 //       mpc_csp_solver->mc_tes.update_calc_vals(false);
 //       int T_cold_code = 0;
 //       try
 //       {
 //           T_cold_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
 //       }
 //       catch (C_csp_exception)
 //       {
 //           mpc_csp_solver->mc_tes.use_calc_vals(false);
 //           mpc_csp_solver->mc_tes.update_calc_vals(true);
 //           throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
 //       }
 //       mpc_csp_solver->mc_tes.update_calc_vals(true);

 //       if (T_cold_code != C_monotonic_eq_solver::CONVERGED)
 //       {
 //           if (T_cold_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
 //           {
 //               std::string msg = util::format("At time = %lg C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank "
 //                   "iteration to find the cold HTF temperature to balance energy between the CR and PC only reached a convergence "
 //                   "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
 //                   mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
 //               mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
 //           }
 //           else
 //           {
 //               mpc_csp_solver->mc_tes.use_calc_vals(false);
 //               mpc_csp_solver->mc_tes.update_calc_vals(true);
 //               *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
 //               return -1;
 //           }
 //       }

 //       T_htf_hx_in = T_cold_solved + 273.15;
 //       m_dot_hx_in = m_dot_hx_out = mpc_csp_solver->mc_tes_outputs.m_m_dot * 3600.; //[kg/hr] mass flow out of the HX on the field side
 //       m_dot_pc_in = m_dot_hx_out;     //[kg/hr]
 //       T_htf_pc_in = T_htf_hx_out = mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out + 273.15;     //[K]
 //       P_hx_out = P_rec_out * (1. - mpc_csp_solver->mc_tes_outputs.dP_perc / 100.);    //[kPa]
 //   }











 //   // Recombine excess mass flow from the power cycle
 //   double T_htf_hx_in;     //[K]
 //   double m_dot_rec_in = mpc_csp_solver->m_m_dot_pc_max;                       //[kg/hr]
 //   if (m_dot_rec_in > m_dot_rec_out) {
 //       mpc_csp_solver->mc_tes.use_calc_vals(false);
 //       mpc_csp_solver->mc_tes.update_calc_vals(true);
 //       *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
 //       return -1;
//   }
 //   else {
 //       T_htf_hx_in = T_htf_rec_out;      //[K]
 //   }
 //   double m_dot_hx_in = m_dot_rec_in;  //[kg/hr]  don't let receiver create mass

 //   // Estimate available discharge in order to updated m_m_dot_tes_dc_max
 //   double q_dot_dc_est, m_dot_field_est, T_hot_field_est;
 //   mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_hx_in, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
 //       q_dot_dc_est, m_dot_field_est, T_hot_field_est);
 //   m_dot_field_est *= 3600.;   //[kg/hr]

 //   // Divert HTF around HT HX if it requires more than what hot tank contains
 //   double m_dot_bypassed = 0.;
 //   if (m_dot_hx_in > m_dot_field_est * 0.99) {
 //       m_dot_bypassed = m_dot_hx_in - m_dot_field_est * 0.99;      //[kg/hr]
 //       m_dot_hx_in = m_dot_field_est * 0.99;                       //[kg/hr]
 //   }
 //   
 //   // Solve the HT HX using a steady-state storage media mass flow, equal to the receiver storage media mass flow
 //   double T_htf_hx_out;
 //   mpc_csp_solver->mc_tes.discharge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
 //                                   mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
 //                                   m_dot_hx_in / 3600.,
 //                                   T_htf_hx_in,
 //                                   T_htf_hx_out,
 //                                   mpc_csp_solver->mc_tes_outputs);
 //   double T_store_hot_ave = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;       //[C]
 //   double P_hx_out = P_rec_out * (1. - mpc_csp_solver->mc_tes_outputs.dP_perc / 100.); //[kPa]

 //   // If any HTF was diverted around the HT HX, mix with HTF after HT HX
 //   double T_htf_pc_in;     //[K]
 //   if (m_dot_bypassed > 0) {
 //       // get enthalpy, assume sCO2 HTF
 //       CO2_state co2_props;
 //       int prop_error_code = CO2_TP(T_htf_hx_in, P_rec_out, &co2_props);
 //       double h_in = co2_props.enth;
 //       double h_out = h_in;
 //       prop_error_code = CO2_PH(P_hx_out, h_out, &co2_props);
 //       double T_htf_bypassed = co2_props.temp; //[K]

 //       T_htf_pc_in = (T_htf_hx_out * m_dot_hx_in + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_in + m_dot_bypassed);  // [K]  mix streams to get HX inlet temp
 //   }
 //   else {
 //       T_htf_pc_in = T_htf_hx_out;     //[K]
 //   }
 //   double m_dot_pc_in = m_dot_hx_in + m_dot_bypassed;

    //// Solve the PC performance at MAX PC HTF FLOW RATE
    //// Need to do this to get back PC T_htf_cold
    //// HTF State
    //mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_hx_out - 273.15;        //[C]
 //   mpc_csp_solver->mc_pc_htf_state_in.m_pres = P_hx_out;                   //[kPa]
    //// Inputs
    //mpc_csp_solver->mc_pc_inputs.m_m_dot = mpc_csp_solver->m_m_dot_pc_max;    //[kg/hr]
    //mpc_csp_solver->mc_pc_inputs.m_standby_control = m_pc_mode;               //[-]
    //// Performance Call
    //mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
    //  mpc_csp_solver->mc_pc_htf_state_in,
    //  mpc_csp_solver->mc_pc_inputs,
    //  mpc_csp_solver->mc_pc_out_solver,
    //  mpc_csp_solver->mc_kernel.mc_sim_info);

    //// Check that power cycle is producing power and solving without errors
    //if (!mpc_csp_solver->mc_pc_out_solver.m_was_method_successful && mpc_csp_solver->mc_pc_inputs.m_standby_control == C_csp_power_cycle::ON)
    //{
    //  mpc_csp_solver->mc_tes.use_calc_vals(false);
    //  mpc_csp_solver->mc_tes.update_calc_vals(true);
    //  *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
    //  return -2;
    //}

    //// Get power cycle HTF return state
    //double T_htf_pc_out = mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold + 273.15;     //[K]
 //   double m_dot_pc_out = mpc_csp_solver->mc_pc_out_solver.m_m_dot_htf;                 //[kg/hr]
 //   double P_pc_out = mpc_csp_solver->mc_pc_out_solver.m_P_phx_in * 1000.;              //[kPa]

 //   // Discharge virtual warm tank through LT HX
 //   double m_dot_hx_out;    //[kg/s]
 //   mpc_csp_solver->mc_tes.discharge_full_lt(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
 //       mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
 //       T_htf_pc_out,
 //       T_htf_hx_out,
 //       m_dot_hx_out,
 //       mpc_csp_solver->mc_tes_outputs);
 //   double T_store_cold_ave = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;       //[C]

 //   m_dot_hx_out *= 3600.;      //[kg/hr]
 //   double P_lthx_out = P_pc_out * (1. - mpc_csp_solver->mc_tes_outputs.dP_perc / 100.);           //[kPa]

 //   // Recombine excess mass flow from the power cycle
 //   double T_htf_rec_in;     //[K]
 //   if (m_dot_pc_out > m_dot_hx_out) {
 //       double m_dot_bypassed = m_dot_pc_out - m_dot_hx_out;                                    //[kg/hr]

 //       // get enthalpy, assume sCO2 HTF
 //       CO2_state co2_props;
 //       int prop_error_code = CO2_TP(T_htf_pc_out, P_pc_out, &co2_props);
 //       double h_in = co2_props.enth;
 //       double h_out = h_in;
 //       prop_error_code = CO2_PH(P_lthx_out, h_out, &co2_props);
 //       double T_htf_bypassed = co2_props.temp; //[K]

 //       T_htf_rec_in = (T_htf_hx_out * m_dot_hx_out + T_htf_bypassed * m_dot_bypassed) / (m_dot_hx_out + m_dot_bypassed);  // [K]  mix streams to get LT HX outlet temp
 //   }
 //   else {
 //       T_htf_rec_in = T_htf_hx_out;      //[K]
 //   }

 //   //Calculate pressure difference (which is not used)
 //   double diff_P = (P_lthx_out - P_in) / P_in;

 //   // Set charging inlet/outlet temps to hot/cold ave temps, respectively
 //   mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                    //[kg/hr]
 //   mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = T_store_hot_ave;  //[C]
 //   mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_store_cold_ave;  //[C]

 //   // Set discharge HTF state
 //   mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = mpc_csp_solver->m_m_dot_pc_max; //[kg/hr]
 //   mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_hx_in - 273.15;         //[C]
 //   mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hx_out - 273.15;           //[C]

    ////Calculate diff_T_htf_cold
    //*diff_T_htf_cold = ((T_htf_rec_in - 273.15) - T_htf_cold) / T_htf_cold;       //[-]

 //   mpc_csp_solver->mc_tes.use_calc_vals(false);
 //   mpc_csp_solver->mc_tes.update_calc_vals(true);
    //return 0;
}

int C_csp_solver::C_MEQ_cr_on__pc_off__tes_ch__T_htf_cold::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Solve the collector-receiver
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Now, solved TES charge with CR outputs
    double T_htf_tes_cold_out = std::numeric_limits<double>::quiet_NaN();
    bool tes_charge_success = mpc_csp_solver->mc_tes.charge(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step, 
                                mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15, 
                                mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot / 3600.0, 
                                mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot + 273.15,
                                T_htf_tes_cold_out, 
                                mpc_csp_solver->mc_tes_outputs);
    T_htf_tes_cold_out -= 273.15;       //[C] convert back from K

    // Set charge htf state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;    //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;      //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_htf_tes_cold_out;            //[C]

    // Set discharge htf state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = 0.0;                                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;   //[C] convert from K
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;   //[C] convert from K

    if (!tes_charge_success)
    {   // If receiver output overcharges storage during iteration, then assume we need some defocus and break loop
        // Receiver thermal output is *roughly* constant for varying receiver inlet temperatures,
        // ... and we don't want to try to throttle thermal power output by controlling this value
        return -2;
    }

    *diff_T_htf_cold = (T_htf_tes_cold_out - T_htf_cold) / T_htf_cold;

    return 0;
}

int C_csp_solver::C_MEQ_cr_on__pc_target__tes_empty__T_htf_cold::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // Clear public member data  
    m_step = std::numeric_limits<double>::quiet_NaN();

    // Solve CR at full timestep
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or model didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get the receiver mass flow rate
    double m_dot_rec_full_ts = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double T_htf_rec_hot = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;   //[C]

    // Get the maximum possible mass flow rate from TES discharge
    // ... using the guess value for the TES cold inlet temperature
    double T_htf_tes_hot, m_dot_htf_full_ts;
    T_htf_tes_hot = m_dot_htf_full_ts = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.discharge_full(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_cold + 273.15,
        T_htf_tes_hot,
        m_dot_htf_full_ts,
        mpc_csp_solver->mc_tes_outputs);

    // Know the mass flow rate for full discharge and the duration of the full timestep
    // so can calculate max TES discharge MASS
    double mass_tes_max = m_dot_htf_full_ts*mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step;     //[kg]

    C_MEQ_cr_on__pc_target__tes_empty__step c_eq(mpc_csp_solver, m_defocus, T_htf_cold);
    C_monotonic_eq_solver c_solver(c_eq);

    double time_max = mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step;       //[s]

    // Check whether the mass flow rate at the full timestep is less than the minimum PC mass flow rate
    if (m_dot_htf_full_ts*3600.0 + m_dot_rec_full_ts < mpc_csp_solver->m_m_dot_pc_min)
    {
        // If it is, then calculate the time required to deplete storage at the minimum mass flow rate
        time_max = mass_tes_max / ((mpc_csp_solver->m_m_dot_pc_min - m_dot_rec_full_ts) / 3600.0);  //[s]
    }
    // ELSE:
    // At full timestep the discharge mass flow is greater than minimum,
    //    so iterate on thermal power

    // Now use this timestep to calculate the thermal power to the power cycle
    double q_dot_pc_m_dot_min = std::numeric_limits<double>::quiet_NaN();
    int eq_code = c_solver.test_member_function(time_max, &q_dot_pc_m_dot_min);
    if (eq_code != 0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Solve PC model with calculated inlet and timestep values
    solve_pc(time_max, c_eq.m_T_htf_pc_hot, c_eq.m_m_dot_pc);

    *diff_T_htf_cold = (mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold - T_htf_cold) / T_htf_cold;

    // Compare solved q_dot_pc to target value
    // If it is greater, than exit because we can't decrease the mass flow rate
    if ((q_dot_pc_m_dot_min - m_q_dot_pc_target) / m_q_dot_pc_target > -1.E-3)
    {
        m_step = time_max;
        return 0;
    }

    // Calculate the time required to deplete storage at the maximum PC HTF mass flow rate
    // ... Be conservative with 0.75 multiplier
    double time_min = 0.75*std::max(mass_tes_max / ( (mpc_csp_solver->m_m_dot_pc_max - m_dot_rec_full_ts) / 3600.0), 0.001);    //[s]

    //// Now use this timestep to calculate the thermal power to the power cycle
    //double q_dot_pc_m_dot_max = std::numeric_limits<double>::quiet_NaN();
    //eq_code = c_solver.test_member_function(time_min, &q_dot_pc_m_dot_max);
    //if (eq_code != 0)
    //{
    //  *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
    //  return -2;
    //}
    //
    //// Compare solved q_dot_pc to target value
    //// If it is less than target, than exit because we can't increase the mass flow rate
    //if ((q_dot_pc_m_dot_max - m_q_dot_pc_target) / m_q_dot_pc_target < -1.E-3)
    //{
    //  *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
    //  return -3;
    //}

    //time_min = 0.001;
    
    // Set up solver to iterate on timestep to achieve q_dot_pc_target
    c_solver.settings(1.E-3, 50, time_min, time_max, true);

    // Guess time required to deplete storage while delivering thermal power requirements to PC

    double time_guess_q_dot_high = mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step*
        (mpc_csp_solver->mc_tes_outputs.m_q_dot_dc_to_htf / (m_q_dot_pc_target - mpc_csp_solver->mc_cr_out_solver.m_q_thermal));    //[s]   

    time_guess_q_dot_high = std::max(1.02*time_min, std::min(time_guess_q_dot_high, 0.98*time_max));        //[s]
    double time_guess_q_dot_low = std::max(1.01*time_min, 0.85*time_guess_q_dot_high);      //[s]

    double time_solved, tol_solved;
    time_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();    //[s]
    int iter_solved = -1;

    int solver_code = 0;

    try
    {
        solver_code = c_solver.solve(time_guess_q_dot_high, time_guess_q_dot_low, m_q_dot_pc_target, time_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception("C_mono_eq_cr_on__pc_target__tes_empty__step method to calculate the time step required to empty TES at the target thermal power returned an exemption"));
    }

    if (solver_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (solver_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) <= 0.1)
        {
            mpc_csp_solver->error_msg = util::format("At time = %lg the iteration to find the time step resulting in emptying TES at the target thermal power only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, mpc_csp_solver->error_msg);
        }
        else
        {
            *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
            return -2;
        }
    }

    // Solve PC model with calculated inlet and timestep values
    solve_pc(time_solved, c_eq.m_T_htf_pc_hot, c_eq.m_m_dot_pc);

    *diff_T_htf_cold = (mpc_csp_solver->mc_pc_out_solver.m_T_htf_cold - T_htf_cold) / T_htf_cold;
    m_step = time_solved;   //[s]

    return 0;
}

int C_csp_solver::C_MEQ_cr_on__pc_target__tes_empty__step::operator()(double step /*s*/, double *q_dot_pc /*MWt*/)
{
    m_m_dot_pc = std::numeric_limits<double>::quiet_NaN();      //[kg/hr]
    m_T_htf_pc_hot = std::numeric_limits<double>::quiet_NaN();  //[MWt]

    // Set new local timestep
    C_csp_solver_sim_info temp_sim_info = mpc_csp_solver->mc_kernel.mc_sim_info;
    temp_sim_info.ms_ts.m_step = step;  //[s]
    temp_sim_info.ms_ts.m_time = mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time -
        mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step +
        step;   //[s]

    mpc_csp_solver->mc_cr_htf_state_in.m_temp = m_T_htf_cold;       //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        temp_sim_info);

    // Check if receiver is OFF or model didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        *q_dot_pc = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get the receiver mass flow rate
    double m_dot_rec = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double T_htf_rec_hot = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;   //[C]
    double q_dot_rec = mpc_csp_solver->mc_cr_out_solver.m_q_thermal;        //[MWt]

    double T_htf_tes_hot, m_dot_tes_dc = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.discharge_full(step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_T_htf_cold + 273.15,
        T_htf_tes_hot,
        m_dot_tes_dc,
        mpc_csp_solver->mc_tes_outputs);

    m_dot_tes_dc *= 3600.0;     //[kg/hr] convert from [kg/s]
    T_htf_tes_hot -= 273.15;    //[C] convert from [K]
    double q_dot_tes = mpc_csp_solver->mc_tes_outputs.m_q_dot_dc_to_htf;    //[MWt]

    // Set TES HTF states (this needs to be less bulky...)
    // HTF discharging state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_dot_tes_dc;         //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = m_T_htf_cold;       //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_tes_hot;     //[C]

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                                  //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;    //[C] convert from K
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;  //[C] convert from K

    // Enthalpy balance on HTF mix from CR and TES
    m_m_dot_pc = m_dot_rec + m_dot_tes_dc;      //[kg/hr]
    m_T_htf_pc_hot = (T_htf_rec_hot*m_dot_rec + T_htf_tes_hot*m_dot_tes_dc) / m_m_dot_pc;   //[C]
    *q_dot_pc = q_dot_rec + q_dot_tes;      //[MWt]

    return 0;
}

void C_csp_solver::C_MEQ_cr_on__pc_target__tes_empty__T_htf_cold::solve_pc(double step /*s*/, double T_htf_pc_hot /*C*/, double m_dot_htf_pc /*kg/hr*/)
{
    // Solve PC model
    mpc_csp_solver->mc_pc_htf_state_in.m_temp = T_htf_pc_hot;   //[C]

    mpc_csp_solver->mc_pc_inputs.m_m_dot = m_dot_htf_pc;        //[kg/hr]
    mpc_csp_solver->mc_pc_inputs.m_standby_control = C_csp_power_cycle::ON;

    // Set new local timestep
    C_csp_solver_sim_info temp_sim_info = mpc_csp_solver->mc_kernel.mc_sim_info;
    temp_sim_info.ms_ts.m_step = step;  //[s]
    temp_sim_info.ms_ts.m_time = mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time -
        mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step +
        step;   //[s]

    // Performance Call
    mpc_csp_solver->mc_power_cycle.call(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_pc_htf_state_in,
        mpc_csp_solver->mc_pc_inputs,
        mpc_csp_solver->mc_pc_out_solver,
        temp_sim_info);
}

int C_csp_solver::C_MEQ_cr_df__pc_off__tes_full__T_cold::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    mpc_csp_solver->mc_cr_htf_state_in.m_temp = T_htf_cold;     //[C]

    mpc_csp_solver->mc_collector_receiver.on(mpc_csp_solver->mc_weather.ms_outputs,
        mpc_csp_solver->mc_cr_htf_state_in,
        m_defocus,
        mpc_csp_solver->mc_cr_out_solver,
        mpc_csp_solver->mc_kernel.mc_sim_info);

    // Check if receiver is OFF or model didn't solve
    if (mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot == 0.0 || mpc_csp_solver->mc_cr_out_solver.m_q_thermal == 0.0)
    {
        *diff_T_htf_cold = std::numeric_limits<double>::quiet_NaN();
        return -1;
    }

    // Get receiver HTF outlet temperature
    double T_htf_rec_hot = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;       //[C]

    // Solve TES for *full* charge
    double T_htf_tes_cold, m_dot_tes;
    T_htf_tes_cold = m_dot_tes = std::numeric_limits<double>::quiet_NaN();
    mpc_csp_solver->mc_tes.charge_full(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        T_htf_rec_hot + 273.15,
        T_htf_tes_cold,
        m_dot_tes,
        mpc_csp_solver->mc_tes_outputs);

    T_htf_tes_cold -= 273.15;   //[C] convert from K
    m_dot_tes *= 3600.0;        //[kg/hr] convert from kg/s

    // HTF charging state
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = m_dot_tes;                //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_cr_out_solver.m_T_salt_hot;      //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = T_htf_tes_cold;        //[C]

    // If not actually discharging (i.e. mass flow rate = 0.0), what should the temperatures be?
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = 0.0;                                      //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;   //[C] convert from K
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;   //[C] convert from K

    *diff_T_htf_cold = (T_htf_tes_cold - T_htf_cold) / T_htf_cold;  //[-]
    
    return 0;
}

int C_csp_solver::C_MEQ_cr_df__pc_off__tes_full__defocus::operator()(double defocus /*-*/, double *diff_m_dot /*-*/)
{
    C_MEQ_cr_df__pc_off__tes_full__T_cold c_eq(mpc_csp_solver, defocus);
    C_monotonic_eq_solver c_solver(c_eq);

    // Set up solver
    c_solver.settings(1.E-3, 50, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(), false);

    // Solve for cold temperature
    double T_cold_guess_low = mpc_csp_solver->m_T_htf_pc_cold_est - 10.0;   //[C]
    double T_cold_guess_high = T_cold_guess_low + 10.0;     //[C]

    double T_cold_solved, tol_solved;
    T_cold_solved = tol_solved = std::numeric_limits<double>::quiet_NaN();
    int iter_solved = -1;

    int T_cold_code = 0;
    try
    {
        T_cold_code = c_solver.solve(T_cold_guess_low, T_cold_guess_high, 0.0, T_cold_solved, tol_solved, iter_solved);
    }
    catch (C_csp_exception)
    {
        throw(C_csp_exception(util::format("At time = %lg, C_csp_solver::C_MEQ_cr_df__pc_off__tes_full__defocus failed", mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time), ""));
    }

    if (T_cold_code != C_monotonic_eq_solver::CONVERGED)
    {
        if (T_cold_code > C_monotonic_eq_solver::CONVERGED && fabs(tol_solved) < 0.1)
        {
            std::string msg = util::format("At time = %lg C_csp_solver::C_MEQ_cr_df__pc_off__tes_full__defocus "
                "iteration to find the cold HTF temperature to balance energy between the CR and PC only reached a convergence "
                "= %lg. Check that results at this timestep are not unreasonably biasing total simulation results",
                mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_time / 3600.0, tol_solved);
            mpc_csp_solver->mc_csp_messages.add_message(C_csp_messages::NOTICE, msg);
        }
        else
        {
            *diff_m_dot = std::numeric_limits<double>::quiet_NaN();
            return -1;
        }
    }

    // Calculate and report mass flow rate balance
    double m_dot_rec = mpc_csp_solver->mc_cr_out_solver.m_m_dot_salt_tot;   //[kg/hr]
    double m_dot_tes = mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot;         //[kg/hr]

    *diff_m_dot = (m_dot_rec - m_dot_tes) / m_dot_rec;

    return 0;
}

int C_csp_solver::C_MEQ_cr_on_tes_dc_m_dot_tank::operator()(double T_htf_cold /*C*/, double *diff_T_htf_cold /*-*/)
{
    // T_htf_cold is the guessed temperature into the high-temp (HT) HX
    
    double q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est;
    mpc_csp_solver->mc_tes.discharge_avail_est(T_htf_cold + 273.15, mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        q_dot_dc_est, m_dot_field_est, T_hot_field_est, m_dot_store_est);

    double T_htf_hot;  //[K] HTF temp out of the HX on the field side
    mpc_csp_solver->mc_tes.discharge_tes_side(mpc_csp_solver->mc_kernel.mc_sim_info.ms_ts.m_step,
        mpc_csp_solver->mc_weather.ms_outputs.m_tdry + 273.15,
        m_m_dot_store / 3600.,
        T_htf_cold + 273.15,
        T_htf_hot,
        m_m_dot_htf_out,
        mpc_csp_solver->mc_tes_outputs);
    m_m_dot_htf_out *= 3600.;

    if (m_m_dot_htf_out < m_m_dot_rec_out) {
        return -1; // in this MEQ we're assuming m_dot_HX > m_dot_CR
    }
    double m_dot_pc = m_m_dot_htf_out;     //[kg/hr]  
    double m_dot_bypassed = std::max(0., m_dot_pc - m_m_dot_rec_out);   //[kg/hr] amount of HTF from before the CR that is mixed after the CR before the HX

    // Recombine excess mass flow from the power cycle
    // get enthalpy, assume sCO2 HTF
    CO2_state co2_props;
    int prop_error_code = CO2_TP(m_T_htf_rec_in, m_P_htf_rec_in, &co2_props);
    double h_in = co2_props.enth;
    double h_out = h_in;
    prop_error_code = CO2_PH(m_P_htf_rec_out, h_out, &co2_props);
    double T_htf_bypassed = co2_props.temp; //[K]

    double T_htf_hx_in = (m_T_htf_rec_out * m_m_dot_rec_out + T_htf_bypassed * m_dot_bypassed) / (m_m_dot_rec_out + m_dot_bypassed);  // [K]  mix streams to get LT HX outlet temp

    // Set charging inlet/outlet temps to hot/cold ave temps, respectively
    mpc_csp_solver->mc_tes_ch_htf_state.m_m_dot = 0.0;                  //[kg/hr]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_in = mpc_csp_solver->mc_tes_outputs.m_T_hot_ave - 273.15;    //[C]
    mpc_csp_solver->mc_tes_ch_htf_state.m_temp_out = mpc_csp_solver->mc_tes_outputs.m_T_cold_ave - 273.15;  //[C]

    // Set discharge HTF state
    mpc_csp_solver->mc_tes_dc_htf_state.m_m_dot = m_m_dot_store;    //[kg/hr]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_in = T_htf_hx_in - 273.15;           //[C]
    mpc_csp_solver->mc_tes_dc_htf_state.m_temp_out = T_htf_hot - 273.15;            //[C]


    *diff_T_htf_cold = (T_htf_hx_in - 273.15 - T_htf_cold) / T_htf_cold;    //[-]

    return 0;
}
