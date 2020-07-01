/**
BSD-3-Clause
Copyright 2019 Alliance for Sustainable Energy, LLC
Redistribution and use in source and binary forms, with or without modification, are permitted provided 
that the following conditions are met :
1.	Redistributions of source code must retain the above copyright notice, this list of conditions 
and the following disclaimer.
2.	Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
and the following disclaimer in the documentation and/or other materials provided with the distribution.
3.	Neither the name of the copyright holder nor the names of its contributors may be used to endorse 
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

#ifndef __csp_solver_tower_collector_receiver_
#define __csp_solver_tower_collector_receiver_

#include "csp_solver_core.h"
#include "csp_solver_mspt_collector_receiver.h"
#include "csp_solver_two_tank_tes.h"
#include "CO2_properties.h"
//#include "csp_solver_pt_sf_perf_interp.h"
//#include "csp_solver_pt_receiver.h"



class C_csp_tower_collector_receiver : public C_csp_collector_receiver
{

private:
    //const double dP_rec_perc = 3.5;         // [%] HTF pressure drop in individual receiver as percent of inlet pressure
    //const double dP_recHX_perc = 0.15;       // [%] HTF pressure drop in individual receiver HX as percent of inlet pressure

    std::vector<C_csp_mspt_collector_receiver> collector_receivers;
    std::vector<C_heat_exchanger> hxs;
    C_csp_two_tank_two_hx_tes *tes;

    sco2Properties mc_field_htfProps;		// Instance of HTFProperties class for field HTF
    HTFProperties mc_store_htfProps;		// Instance of HTFProperties class for storage HTF

    std::string error_msg;                  // member string for exception messages

    double m_P_rec_in_des;        // [kPa] Receiver (after riser) inlet pressure at design

public:
	
	enum
    {
        E_FIELD_AREA_TOT,            //[m2] Total heliostat field area
        E_FIELD_AREA1,              //[m2] Heliostat field area 1
        E_FIELD_AREA2,              //[m2] Heliostat field area 1
        E_FIELD_AREA3,              //[m2] Heliostat field area 1
        E_FIELD_Q_DOT_INC,          //[MWt] Field incident thermal power
        E_FIELD_ETA_TOT,            //[-] Average optical efficiency including tower reflectivity
        E_FIELD_ETA_OPT1,            //[-] Optical efficiency including Tower refl
        E_FIELD_ETA_OPT2,            //[-] Optical efficiency including Tower refl
        E_FIELD_ETA_OPT3,            //[-] Optical efficiency including Tower refl
        E_FIELD_ADJUST,             //[-] Field adjustment factor
        
        E_Q_DOT_INC,                //[MWt] Tower incident thermal power
        E_ETA_THERMAL,              //[-] Tower thermal efficiency
        E_Q_DOT_THERMAL,            //[MWt] Field incident thermal power
        E_M_DOT_HTF,                //[kg/hr] Tower mass flow rate
        E_Q_DOT_STARTUP,            //[MWt] Tower startup thermal power consumed
        E_T_HTF_IN,                 //[C] Tower HTF inlet temperature
        E_T_HTF_OUT,                //[C] Tower HTF outlet temperature
        E_Q_DOT_PIPE_LOSS,          //[MWt] Tower piping losses
        E_Q_DOT_LOSS,               //[MWt] Tower convection and radiation losses
        E_P_HEATTRACE,              //[MWe] Tower heat trace parasitic
        E_T_HTF_OUT_END,            //[C] Instantaneous tower HTF outlet temperature at the end of the time step
        E_T_HTF_OUT_MAX,            //[C] Tower maximum HTF outlet temperature at any point during time step
        E_T_HTF_PANEL_OUT_MAX,      //[C] Tower panel maximum HTF outlet temperature at any point during time step
        E_T_WALL_INLET,             //[C] Tower inlet wall temperature at end of time step
        E_T_WALL_OUTLET,            //[C] Tower inlet wall temperature at end of time step
        E_T_RISER,                  //[C] Riser temperature at the end of the time step
        E_T_DOWNC,                  //[C] Downcomer temperature at the end of the time step

        E_Q_DOT_INC1,               //[MWt] Receiver 1 incident thermal power
        E_Q_DOT_INC2,               //[MWt] Receiver 2 incident thermal power
        E_Q_DOT_INC3,               //[MWt] Receiver 3 ncident thermal power
        E_M_DOT_HTF1,               //[kg/hr] Receiver 1 mass flow rate
        E_M_DOT_HTF2,               //[kg/hr] Receiver 2 mass flow rate
        E_M_DOT_HTF3,               //[kg/hr] Receiver 3 mass flow rate
        E_T_HTF_IN1,                //[C] Receiver 1 HTF inlet temperature
        E_T_HTF_IN2,                //[C] Receiver 2 HTF inlet temperature
        E_T_HTF_IN3,                //[C] Receiver 3 HTF inlet temperature
        E_T_HTF_OUT1,               //[C] Receiver 1 HTF outlet temperature
        E_T_HTF_OUT2,               //[C] Receiver 2 HTF outlet temperature
        E_T_HTF_OUT3,               //[C] Receiver 3 HTF outlet temperature
        E_Q_DOT_HX1,                //[MWt] HX 1 transferred power
        E_M_DOT_HX1,                //[kg/hr] HX 1 TES mass flow rate
        E_T_HX_OUT1,                //[C] HX 1 TES outlet temperature
        E_Q_DOT_HX2,                //[MWt] HX 1 transferred power
        E_M_DOT_HX2,                //[kg/hr] HX 1 TES mass flow rate
        E_T_HX_OUT2,                //[C] HX 1 TES outlet temperature
        E_Q_DOT_HX3,                //[MWt] HX 1 transferred power
        E_M_DOT_HX3,                //[kg/hr] HX 1 TES mass flow rate
        E_T_HX_OUT3,                //[C] HX 1 TES outlet temperature
        E_ETA_THERM1,               //[-] Receiver 1 thermal efficiency
        E_ETA_THERM2,               //[-] Receiver 2 thermal efficiency
        E_ETA_THERM3,               //[-] Receiver 3 thermal efficiency
        E_DP_REC1,                  //[kPa] Receiver 1 pressure drop
        E_DP_REC2,                  //[kPa] Receiver 2 pressure drop
        E_DP_REC3,                  //[kPa] Receiver 3 pressure drop
        E_DP_RISER,                 //[kPa] Riser pressure drop
        E_DP_DOWNCOMER,             //[kPa] Downcomer pressure drop
    };
	
    int m_field_fl;
    util::matrix_t<double> m_field_fl_props;
    int m_tes_fl;
    util::matrix_t<double> m_tes_fl_props;
    double hx_duty;
    double m_dt_hot;
    double T_rec_hot_des;
    double T_hx_cold_des;       // design cold inlet temperature on storage side
    double h_lift;              // [m] Effective lift height
    double riser_diam;          // [m] Riser inner diameter
    double downcomer_diam;          // [m] Downcomer inner diameter
    double riser_length;        // [m] Total riser length for pressure drop calc
    double dP_recHX_perc;       // [%] HTF pressure drop in individual receiver HX as percent of inlet pressure

	C_csp_reported_outputs mc_reported_outputs;
	
	C_csp_tower_collector_receiver(std::vector<C_csp_mspt_collector_receiver> & collector_receivers);

	~C_csp_tower_collector_receiver();

	virtual void init(const C_csp_collector_receiver::S_csp_cr_init_inputs init_inputs, 
			C_csp_collector_receiver::S_csp_cr_solved_params & solved_params);

	virtual int get_operating_state();

    virtual double get_startup_time();
    virtual double get_startup_energy(); //MWh
    virtual double get_pumping_parasitic_coef();  //MWe/MWt
    virtual double get_min_power_delivery();    //MWt
	virtual double get_tracking_power();		//MWe
	virtual double get_col_startup_power();		//MWe-hr
    void set_tes(C_csp_two_tank_two_hx_tes *tes);

    virtual void off(const C_csp_weatherreader::S_outputs &weather,
		const C_csp_solver_htf_1state &htf_state_in,
		C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
		//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
		const C_csp_solver_sim_info &sim_info);

	virtual void startup(const C_csp_weatherreader::S_outputs &weather,
		const C_csp_solver_htf_1state &htf_state_in,
		C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
		//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
		const C_csp_solver_sim_info &sim_info);

	virtual void on(const C_csp_weatherreader::S_outputs &weather,
		const C_csp_solver_htf_1state &htf_state_in,
		double field_control,
		C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
		//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
		const C_csp_solver_sim_info &sim_info);

	virtual void estimates(const C_csp_weatherreader::S_outputs &weather,
		const C_csp_solver_htf_1state &htf_state_in,
		C_csp_collector_receiver::S_csp_cr_est_out &est_out,
		const C_csp_solver_sim_info &sim_info);

	virtual void converged();

	virtual void write_output_intervals(double report_time_start,
		const std::vector<double> & v_temp_ts_time_end, double report_time_end);

    virtual double calculate_optical_efficiency( const C_csp_weatherreader::S_outputs &weather, const C_csp_solver_sim_info &sim );
  
    virtual double calculate_thermal_efficiency_approx( const C_csp_weatherreader::S_outputs &weather, double q_incident /*MW*/ );

    virtual double get_collector_area();

    double conveyor_power(double m_dot_particle);


	void call(const C_csp_weatherreader::S_outputs &weather,
		const C_csp_solver_htf_1state &htf_state_in,
		const C_csp_collector_receiver::S_csp_cr_inputs &inputs,
		C_csp_collector_receiver::S_csp_cr_out_solver &cr_out_solver,
		//C_csp_collector_receiver::S_csp_cr_out_report &cr_out_report,
		const C_csp_solver_sim_info &sim_info);

};

#endif //__csp_solver_tower_collector_receiver_
