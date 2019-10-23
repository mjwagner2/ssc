#include "csp_solver_pt_receiver.h"
#include "csp_solver_mspt_receiver.h"
#include "csp_solver_core.h"
#include "sam_csp_util.h"
#include <set>

C_mspt_receiver::C_mspt_receiver()
{    
	m_w_rec = std::numeric_limits<double>::quiet_NaN();
	m_h_rec = std::numeric_limits<double>::quiet_NaN();

	m_pipe_loss_per_m = std::numeric_limits<double>::quiet_NaN();
	m_pipe_length_add = std::numeric_limits<double>::quiet_NaN();
	m_pipe_length_mult = std::numeric_limits<double>::quiet_NaN();

	m_T_salt_hot_target = std::numeric_limits<double>::quiet_NaN();
	m_eta_pump = std::numeric_limits<double>::quiet_NaN();
	m_night_recirc = -1;
	m_hel_stow_deploy = std::numeric_limits<double>::quiet_NaN();

	m_field_fl = -1;

	m_Q_dot_piping_loss = std::numeric_limits<double>::quiet_NaN();
	m_m_dot_htf_max = std::numeric_limits<double>::quiet_NaN();

	m_itermode = -1;
	m_od_control = std::numeric_limits<double>::quiet_NaN();
	m_eta_field_iter_prev = std::numeric_limits<double>::quiet_NaN();
	m_tol_od = std::numeric_limits<double>::quiet_NaN();
	m_q_dot_inc_min = std::numeric_limits<double>::quiet_NaN();

	m_E_su = std::numeric_limits<double>::quiet_NaN();
	m_E_su_prev = std::numeric_limits<double>::quiet_NaN();
	m_t_su = std::numeric_limits<double>::quiet_NaN();
	m_t_su_prev = std::numeric_limits<double>::quiet_NaN();

	m_ncall = -1;
}

C_mspt_receiver::~C_mspt_receiver()
{}

void C_mspt_receiver::init()
{
	ambient_air.SetFluid(ambient_air.Air);

	// Declare instance of fluid class for FIELD fluid
	if( m_field_fl != HTFProperties::User_defined && m_field_fl < HTFProperties::End_Library_Fluids )
	{
		if( !field_htfProps.SetFluid( m_field_fl ) )
		{
			throw(C_csp_exception("Receiver HTF code is not recognized", "sCO2PT receiver"));
		}
	}
	else if( m_field_fl == HTFProperties::User_defined )
	{
		// Check that 'm_field_fl_props' is allocated and correct dimensions
		int n_rows = (int)m_field_fl_props.nrows();
		int n_cols = (int)m_field_fl_props.ncols();
		if( n_rows > 2 && n_cols == 7 )
		{
			if( !field_htfProps.SetUserDefinedFluid(m_field_fl_props) )
			{
				error_msg = util::format(field_htfProps.UserFluidErrMessage(), n_rows, n_cols);
				throw(C_csp_exception(error_msg, "sCO2PT receiver"));
			}
		}
		else
		{
			error_msg = util::format("The user defined field HTF table must contain at least 3 rows and exactly 7 columns. The current table contains %d row(s) and %d column(s)", n_rows, n_cols);
			throw(C_csp_exception(error_msg, "sCO2PT receiver"));
		}
	}
	else
	{
		throw(C_csp_exception("Receiver HTF code is not recognized", "sCO2PT receiver"));
	}

	
	// Unit Conversions
	m_T_htf_hot_des += 273.15;	//[K] Convert from input in [C]
	m_T_htf_cold_des += 273.15;	//[K] Convert from input in [C]
	m_q_rec_des *= 1.E6;		//[W] Convert from input in [MW]

	m_mode = C_csp_collector_receiver::OFF;					//[-] 0 = requires startup, 1 = starting up, 2 = running
	m_itermode = 1;			//[-] 1: Solve for design temp, 2: solve to match mass flow restriction
	m_od_control = 1.0;			//[-] Additional defocusing for over-design conditions
	m_tol_od = 0.001;		//[-] Tolerance for over-design iteration

	double c_htf_des = field_htfProps.Cp((m_T_htf_hot_des + m_T_htf_cold_des) / 2.0)*1000.0;		//[J/kg-K] Specific heat at design conditions
	m_m_dot_htf_des = m_q_rec_des / (c_htf_des*(m_T_htf_hot_des - m_T_htf_cold_des));					//[kg/s]
	double eta_therm_des = 0.9;
	m_q_dot_inc_min = m_q_rec_des * m_f_rec_min / eta_therm_des;	//[W] Minimum receiver thermal power

	if (m_m_dot_htf_max_frac != m_m_dot_htf_max_frac)
	{
		// if max frac not set, then max mass flow (absolute) needs to be defined
		if (m_m_dot_htf_max != m_m_dot_htf_max)
		{
			throw(C_csp_exception("maximum rec htf mass flow rate not defined", "sCO2PT receiver"));
		}
		m_m_dot_htf_max /= 3600.0;	//[kg/s] Convert from input in [kg/hr]
	}
	m_m_dot_htf_max = m_m_dot_htf_max_frac * m_m_dot_htf_des;	//[kg/s]

	m_mode_prev = m_mode;
	m_E_su_prev = m_q_rec_des * m_rec_qf_delay;	//[W-hr] Startup energy
	m_t_su_prev = m_rec_su_delay;				//[hr] Startup time requirement
	m_eta_field_iter_prev = 1.0;				//[-] Set to largest possible value

	m_T_salt_hot_target += 273.15;			//[K] convert from C
	

    //set up the lookup tables
    m_efficiency_lookup.Set_2D_Lookup_Table(m_efficiency_lookup_input);
    m_pressure_lookup.Set_2D_Lookup_Table(m_pressure_lookup_input);

    //calculate reference receiver efficiency
    m_eta_rec_des = m_efficiency_lookup.bilinear_2D_interp(1., m_T_htf_cold_des - 273.15); 


	param_inputs.T_amb = std::numeric_limits<double>::quiet_NaN();
    param_inputs.T_sky = std::numeric_limits<double>::quiet_NaN();
    param_inputs.c_htf = std::numeric_limits<double>::quiet_NaN();
    param_inputs.rho_htf = std::numeric_limits<double>::quiet_NaN();
    param_inputs.mu_htf = std::numeric_limits<double>::quiet_NaN();
    param_inputs.k_htf = std::numeric_limits<double>::quiet_NaN();
    param_inputs.Pr_htf = std::numeric_limits<double>::quiet_NaN();
	
    //force paramters
    m_night_recirc = 0;

	return;
}

void C_mspt_receiver::call(const C_csp_weatherreader::S_outputs &weather, 
	const C_csp_solver_htf_1state &htf_state_in,
	const C_mspt_receiver::S_inputs &inputs,
	const C_csp_solver_sim_info &sim_info)
{
	// Increase call-per-timestep counter
	// Converge() sets it to -1, so on first call this line will adjust it = 0
	m_ncall++;
	
	// Get inputs
	double field_eff = inputs.m_field_eff;					//[-]
	//const util::matrix_t<double> *flux_map_input = inputs.m_flux_map_input;
		// When this function is called from TCS solver, input_operation_mode should always be == 2
	int input_operation_mode = inputs.m_input_operation_mode;

	if(input_operation_mode < C_csp_collector_receiver::OFF || input_operation_mode > C_csp_collector_receiver::STEADY_STATE)
	{
		error_msg = util::format("Input operation mode must be either [0,1,2], but value is %d", input_operation_mode);
		throw(C_csp_exception(error_msg, "MSPT receiver timestep performance call"));
	}

	// Get sim info 
	double step = sim_info.ms_ts.m_step;			//[s]
	double time = sim_info.ms_ts.m_time;	//[s]

	// Get applicable htf state info
	double T_salt_cold_in = htf_state_in.m_temp;		//[C]
    double P_in = htf_state_in.m_pres;

	// Complete necessary conversions/calculations of input variables
	T_salt_cold_in += 273.15;				//[K] Cold salt inlet temp, convert from C
	double P_amb = weather.m_pres*100.0;	//[Pa] Ambient pressure, convert from mbar
	double hour = time / 3600.0;			//[hr] Hour of the year
	double T_dp = weather.m_tdew + 273.15;	//[K] Dewpoint temperature, convert from C
	double T_amb = weather.m_tdry + 273.15;	//[K] Dry bulb temperature, convert from C
	// **************************************************************************************

	// Read in remaining weather inputs from weather output structure
	double zenith = weather.m_solzen;
	double azimuth = weather.m_solazi;
	double v_wind_10 = weather.m_wspd;
	double I_bn = weather.m_beam;


	double T_sky = CSP::skytemp(T_amb, T_dp, hour);

	// Set current timestep stored values to NaN so we know that code solved for them
	m_mode = -1;
	m_E_su = std::numeric_limits<double>::quiet_NaN();
	m_t_su = std::numeric_limits<double>::quiet_NaN();

	m_itermode = 1;

	double v_wind = log((m_h_tower + m_h_rec / 2) / 0.003) / log(10.0 / 0.003)*v_wind_10;

    bool rec_is_off = false;
	bool rec_is_defocusing = false;
	double field_eff_adj = 0.0;


	// ************* Outputs for ISCC model ****************
	//double q_thermal_ss = 0.0;
	double f_rec_timestep = 1.0;
	// *****************************************************

	// Do an initial check to make sure the solar position called is valid
	// If it's not, return the output equal to zeros. Also check to make sure
	// the solar flux is at a certain level, otherwise the correlations aren't valid
	if( input_operation_mode == C_csp_collector_receiver::OFF )
	{
		rec_is_off = true;
	}

	if( zenith>(90.0 - m_hel_stow_deploy) || I_bn <= 1.E-6 || (zenith == 0.0 && azimuth == 180.0) )
	{
		m_mode = C_csp_collector_receiver::OFF;
		rec_is_off = true;
	}

	double T_coolant_prop = (m_T_salt_hot_target + T_salt_cold_in) / 2.0;		//[K] The temperature at which the coolant properties are evaluated. Validated as constant (mjw)
	double c_p_coolant = field_htfProps.Cp(T_coolant_prop)*1000.0;						//[J/kg-K] Specific heat of the coolant

	double m_dot_htf_max = m_m_dot_htf_max;
    
	//double q_abs_sum = 0.0;
	double err_od = 999.0;	// Reset error before iteration

    // Suggests controller applied defocus, so reset *controller* defocus
	if (field_eff < m_eta_field_iter_prev && m_od_control < 1.0)
		m_od_control = fmin(m_od_control + (1.0 - field_eff / m_eta_field_iter_prev), 1.0);

    double m_dot_salt, q_dot_inc_sum, eta_therm, rho_coolant, T_salt_hot_rec;

	do
	{
		if( rec_is_off )
			break;

		field_eff_adj = field_eff*m_od_control;

        //calculate total power on the reciever
        q_dot_inc_sum = field_eff_adj * inputs.m_A_sf * I_bn;      //W

        //calculate fractional load
        m_q_dot_inc = q_dot_inc_sum;
        double load = m_q_dot_inc / m_q_rec_des;

        //look up receiver efficiency
        eta_therm = m_efficiency_lookup.bilinear_2D_interp(load, T_salt_cold_in - 273.15); 
        if (eta_therm < 0. || eta_therm != eta_therm)
        {
            rec_is_off = true;
            break;
        }

        m_q_dot_loss = (1. - eta_therm) * m_q_dot_inc; 		//[W] Total reciever losses
		m_q_dot_abs = m_q_dot_inc - m_q_dot_loss;	//[W] Absorbed flux at each node
				
		rho_coolant = field_htfProps.dens(T_coolant_prop, 1.0);			//[kg/m^3] Density of the coolant

        m_dot_salt = m_q_dot_abs / (c_p_coolant*(m_T_salt_hot_target - T_salt_cold_in));			//[kg/s]

        T_salt_hot_rec = m_T_salt_hot_target;   // T_salt_cold_in + m_q_dot_abs / (m_dot_salt*c_p_coolant);	//[K] Energy balance for each node																																																

		if( rec_is_off )
			break;

        if(m_dot_salt < 1.e-5)
		{
			m_mode = C_csp_collector_receiver::OFF;				//[-] Set the startup mode
			rec_is_off = true;
		}

		if( rec_is_off )
			break;

		// Limit the HTF mass flow rate to the maximum, if needed
		if( (m_dot_salt > m_dot_htf_max) || m_itermode == 2 )
		{
			err_od = (m_dot_salt - m_dot_htf_max) / m_dot_htf_max;
			if( err_od < m_tol_od )
			{
				m_itermode = 1;
				m_od_control = 1.0;
				rec_is_defocusing = false;
			}
			else
			{
                if (m_od_control < 0.001)
                {
                    rec_is_off = true;
                    break;
                }

				m_od_control = m_od_control*pow((m_dot_htf_max / m_dot_salt), 0.8);	//[-] Adjust the over-design defocus control by modifying the current value
				m_itermode = 2;
				rec_is_defocusing = true;
				// GOTO 15
			}
		}
	} while( rec_is_defocusing );

	double Pres_D, q_thermal, q_startup;
	Pres_D = q_thermal = q_startup = std::numeric_limits<double>::quiet_NaN();

	double q_startup_energy = 0.0;		//[J]
	q_startup = 0.0;

	double time_required_su = step/3600.0;

	if( !rec_is_off )
	{
		double m_dot_rec_des = m_q_rec_des / (c_p_coolant*(m_T_htf_hot_des - m_T_htf_cold_des)); // Design point receiver mass flow rate (kg/s)

		switch( input_operation_mode )
		{
		    case C_csp_collector_receiver::STARTUP:
		    {

			    // Startup model based on fixed time and energy requirements
				double time_require_su_energy = m_E_su_prev / (m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in));	//[hr]
				double time_require_su_ramping = m_t_su_prev;

				double time_required_max = fmax(time_require_su_energy, time_require_su_ramping);	//[hr]

				double time_step_hrs = step / 3600.0;		//[hr]

				if( time_required_max  > time_step_hrs )		// Can't completely startup receiver in maximum allowable timestep
				{											// Need to advance timestep and try again
					time_required_su = time_step_hrs;		
					m_mode = C_csp_collector_receiver::STARTUP;
					q_startup = m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in)*step / 3600.0;
				}
				else
				{
					time_required_su = time_required_max;		//[hr]
					m_mode = C_csp_collector_receiver::ON;

					double q_startup_energy_req = m_E_su_prev;	//[W-hr]
					double q_startup_ramping_req = m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in)*m_t_su;	//[W-hr]
					q_startup = fmax(q_startup_energy_req, q_startup_ramping_req);
				}

				m_E_su = fmax(0.0, m_E_su_prev - m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in)*step / 3600.0);
				m_t_su = fmax(0.0, m_t_su_prev - step / 3600.0);

				rec_is_off = true;

                //look up receiver pressure drop
                Pres_D = m_pressure_lookup.bilinear_2D_interp(m_dot_salt / m_dot_rec_des, P_in / m_P_cold_des);

			    break;
			}
		case C_csp_collector_receiver::ON:
			{

			    // Steady state receiver model
				if (m_E_su_prev > 0.0 || m_t_su_prev > 0.0)
				{
					m_E_su = fmax(0.0, m_E_su_prev - m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in)*step / 3600.0);	//[W-hr]
					m_t_su = fmax(0.0, m_t_su_prev - step / 3600.0);	//[hr]

					if (m_E_su + m_t_su > 0.0)
					{
						m_mode = C_csp_collector_receiver::STARTUP;		// If either are greater than 0, we're staring up but not finished

						// 4.28.15 twn: Startup energy also needs to consider energy consumed during time requirement, if that is greater than energy requirement
						//q_startup = (m_E_su_prev - m_E_su) / (step / 3600.0)*1.E-6;
						q_startup = m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in)*step / 3600.0;
						rec_is_off = true;
						f_rec_timestep = 0.0;
					}
					else
					{
						m_mode = C_csp_collector_receiver::ON;

						double q_startup_energy_req = m_E_su_prev;	//[W-hr]
						double q_startup_ramping_req = m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in)*m_t_su;	//[W-hr]
						q_startup = fmax(q_startup_energy_req, q_startup_ramping_req);

						// Adjust the available mass flow to reflect startup
						m_dot_salt = fmin((1.0 - m_t_su_prev / (step / 3600.0))*m_dot_salt, m_dot_salt - m_E_su_prev / ((step / 3600.0)*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in)));
						f_rec_timestep = fmax(0.0, fmin(1.0 - m_t_su_prev / (step / 3600.0), 1.0 - m_E_su_prev / (m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in))));
					}
					//4.28.15 twn: Startup energy needs to consider
					//q_startup = (m_E_su_prev - m_E_su) / (step / 3600.0)*1.E-6;
				}
				else
				{
					m_E_su = m_E_su_prev;
					m_t_su = m_t_su_prev;
					m_mode = C_csp_collector_receiver::ON;
					q_startup = 0.0;

					if (q_dot_inc_sum*1.E3 < m_q_dot_inc_min)
					{
						// If output here is less than specified allowed minimum, then need to shut off receiver
						m_mode = C_csp_collector_receiver::OFF;

						// Include here outputs that are ONLY set to zero if receiver completely off, and not attempting to start-up
                        Pres_D = 0.0; 
					}
				}
				q_thermal = m_dot_salt*c_p_coolant*(T_salt_hot_rec - T_salt_cold_in);
                Pres_D = m_pressure_lookup.bilinear_2D_interp(m_dot_salt / m_dot_rec_des, P_in / m_P_cold_des);

			    if (q_dot_inc_sum*1.E3 < m_q_dot_inc_min)
				    rec_is_off = true;
			    
                break;
			}

		case C_csp_collector_receiver::STEADY_STATE:
            {
                m_mode = C_csp_collector_receiver::STEADY_STATE;
                f_rec_timestep = 1.0;
                q_thermal = m_dot_salt * c_p_coolant*(T_salt_hot_rec - T_salt_cold_in);

                if (q_dot_inc_sum*1.E3 < m_q_dot_inc_min && m_mode_prev == C_csp_collector_receiver::ON)
                    rec_is_off = true;

                break;
            }
		}	
	}
	else
	{	// If receiver was off BEFORE startup deductions
		m_mode = C_csp_collector_receiver::OFF;

		// Include here outputs that are ONLY set to zero if receiver completely off, and not attempting to start-up
        Pres_D = 0.0; 
	}

	if( rec_is_off )
	{
        // Receiver isn't producing usable energy
        m_dot_salt = 0.0; eta_therm = 0.0; q_thermal = 0.0;
        // Set the receiver outlet temperature equal to the inlet design temperature
        T_salt_hot_rec = m_T_htf_cold_des;
        q_dot_inc_sum = 0.0;
        f_rec_timestep = 0.0; 

        // Reset m_od_control
        m_od_control = 1.0;		//[-]

	}

	// Steady state outputs
	outputs.m_m_dot_salt_tot = m_dot_salt*3600.0;		//[kg/hr] convert from kg/s
	outputs.m_eta_therm = eta_therm;							//[-] RECEIVER thermal efficiency (includes radiation and convective losses. reflection losses are contained in receiver flux model)
    outputs.m_W_dot_pump = 0.; // W_dot_pump / 1.E6;				//[MW] convert from W
	//outputs.m_q_conv_sum = q_conv_sum / 1.E6;				//[MW] convert from W
	//outputs.m_q_rad_sum = q_rad_sum / 1.E6;					//[MW] convert from W
	outputs.m_Q_thermal = q_thermal / 1.E6;					//[MW] convert from W
	outputs.m_T_salt_hot = T_salt_hot_rec - 273.15;		//[C] convert from K
	outputs.m_field_eff_adj = field_eff_adj;					//[-]
	outputs.m_component_defocus = m_od_control;				//[-]
	outputs.m_q_dot_rec_inc = q_dot_inc_sum / 1.E3;			//[MW] convert from kW
	outputs.m_q_startup = q_startup/1.E6;					//[MW-hr] convert from W-hr
	//outputs.m_dP_receiver = DELTAP*m_n_panels / m_n_lines / 1.E5;	//[bar] receiver pressure drop, convert from Pa
	outputs.m_dP_total = Pres_D*10.0;						//[bar] total pressure drop, convert from MPa
	//outputs.m_vel_htf = u_coolant;							//[m/s]
	outputs.m_T_salt_cold = T_salt_cold_in - 273.15;			//[C] convert from K
	//outputs.m_m_dot_ss = m_dot_salt_tot_ss*3600.0;			//[kg/hr] convert from kg/s
	//outputs.m_q_dot_ss = q_thermal_ss / 1.E6;				//[MW] convert from W
	outputs.m_f_timestep = f_rec_timestep;					//[-]
	outputs.m_time_required_su = time_required_su*3600.0;	//[s], convert from hr in code
	if(q_thermal > 0.0)
		outputs.m_q_dot_piping_loss = m_Q_dot_piping_loss/1.E6;	//[MWt]
	else
		outputs.m_q_dot_piping_loss = 0.0;		//[MWt]


	outputs.m_inst_T_salt_hot = 
	    outputs.m_max_T_salt_hot = 
	    outputs.m_min_T_salt_hot = 
	    outputs.m_max_rec_tout = T_salt_hot_rec - 273.15;
	outputs.m_q_heattrace = 0.0;
	outputs.m_Twall_inlet = 0.0;
	outputs.m_Twall_outlet = 0.0;
	outputs.m_Triser = 0.0;
	outputs.m_Tdownc = 0.0;

    ms_outputs = outputs;

	m_eta_field_iter_prev = field_eff;	//[-]
}

void C_mspt_receiver::off(const C_csp_weatherreader::S_outputs &weather,
	const C_csp_solver_htf_1state &htf_state_in,
	const C_csp_solver_sim_info &sim_info)
{
	// Don't currently need *any* of these inputs, but if we add recirculation or thermal capacitance it would be helpful to have in place
	m_mode = C_csp_collector_receiver::OFF;

	// Assuming no night recirculation, so... these should be zero
	outputs.m_m_dot_salt_tot = 0.0;		//[kg/hr] convert from kg/s
	outputs.m_eta_therm = 0.0;			//[-] RECEIVER thermal efficiency (includes radiation and convective losses. reflection losses are contained in receiver flux model)
	outputs.m_W_dot_pump = 0.0;			//[MW] convert from W
	outputs.m_q_conv_sum = 0.0;			//[MW] convert from W
	outputs.m_q_rad_sum = 0.0;			//[MW] convert from W
	outputs.m_Q_thermal = 0.0;			//[MW] convert from W
	outputs.m_T_salt_hot = 0.0;			//[C] convert from K
	outputs.m_field_eff_adj = 0.0;		//[-]
	outputs.m_component_defocus = 1.0;	//[-]
	outputs.m_q_dot_rec_inc = 0.0;		//[MW] convert from kW
	outputs.m_q_startup = 0.0;			//[MW-hr] convert from W-hr
	outputs.m_dP_receiver = 0.0;			//[bar] receiver pressure drop, convert from Pa
	outputs.m_dP_total = 0.0;			//[bar] total pressure drop, convert from MPa
	outputs.m_vel_htf = 0.0;				//[m/s]
	outputs.m_T_salt_cold = 0.0;			//[C] convert from K
	outputs.m_m_dot_ss = 0.0;			//[kg/hr] convert from kg/s
	outputs.m_q_dot_ss = 0.0;			//[MW] convert from W
	outputs.m_f_timestep = 0.0;			//[-]
	outputs.m_time_required_su = sim_info.ms_ts.m_step;	//[s], convert from hr in code
	outputs.m_q_dot_piping_loss = 0.0;	//[MWt]
	
	outputs.m_inst_T_salt_hot = 0.0;
	outputs.m_max_T_salt_hot = 0.0;
	outputs.m_min_T_salt_hot = 0.0;
	outputs.m_max_rec_tout = 0.0;
	outputs.m_q_heattrace = 0.0;
	outputs.m_Twall_inlet = 0.0;
	outputs.m_Twall_outlet = 0.0;
	outputs.m_Triser = 0.0;
	outputs.m_Tdownc = 0.0;

    ms_outputs = outputs;

	return;
}

void C_mspt_receiver::converged()
{
	// Check HTF props?
	//!MJW 9.8.2010 :: Call the property range check subroutine with the inlet and outlet HTF temps to make sure they're in the valid range
	//call check_htf(Coolant,T_salt_hot)
	//call check_htf(Coolant,T_salt_cold)

	if( m_mode == C_csp_collector_receiver::STEADY_STATE )
	{
		throw(C_csp_exception("Receiver should only be run at STEADY STATE mode for estimating output. It must be run at a different mode before exiting a timestep",
			"MSPT receiver converged method"));
	}

	if( m_mode == C_csp_collector_receiver::OFF )
	{
		m_E_su = m_q_rec_des * m_rec_qf_delay;
		m_t_su = m_rec_su_delay;
	}

	m_mode_prev = m_mode;
	m_E_su_prev = m_E_su;
	m_t_su_prev = m_t_su;

	m_itermode = 1;
	m_od_control = 1.0;
	m_eta_field_iter_prev = 1.0;		//[-]

	m_ncall = -1;

    ms_outputs = outputs;
}

void C_mspt_receiver::calc_pump_performance(double rho_f, double mdot, double ffact, double &PresDrop_calc, double &WdotPump_calc)
{

    // Pressure drop calculations
    PresDrop_calc = WdotPump_calc = std::numeric_limits<double>::quiet_NaN();

}

double C_mspt_receiver::get_pumping_parasitic_coef()
{
    return std::numeric_limits<double>::quiet_NaN();
}

double C_mspt_receiver::area_proj()
{
    return m_w_rec * m_h_rec; //[m^2] projected or aperture area of the receiver
}

double C_mspt_receiver::calc_external_convection_coeff(const parameter_eval_inputs &pinputs, double Twall)
{
    return std::numeric_limits<double>::quiet_NaN();
}

void C_mspt_receiver::calc_header_size(double pdrop, double mdot, double rhof, double muf, double Lh, double &id_calc, double &th_calc, double &od_calc)
{
    th_calc = od_calc = std::numeric_limits<double>::quiet_NaN();
}

double C_mspt_receiver::interpolate(double x, const std::vector<double> &xarray, const std::vector<double> &yarray, int klow, int khigh)
{
	// Linear interpolation between discrete tabulated points
	// x = value of independent variable
	// xarray = independent variable points (assumed to be increasing between klow and khigh), 
	// yarray = dependent variable points
	// klow, khigh = lowest,highest indicies of interest 

	int jl = klow, ju = khigh;
	int jm;
	while (ju - jl > 1)
	{
		jm = (ju + jl) / 2;	//middle index of the range
		if (x < xarray.at(jm)) ju = jm;
		else jl = jm;
	}
	double yinterp = yarray.at(jl) + (yarray.at(ju) - yarray.at(jl)) / (xarray.at(ju) - xarray.at(jl)) * (x - xarray.at(jl));
	return yinterp;
}

double C_mspt_receiver::integrate(double xlow, double xhigh, const std::vector<double> &xarray, const std::vector<double> &yarray, int klow, int khigh)
{

	// Numerical integral between upper and lower bounds xlow and xhigh
	// xarray = independent variable points (assumed to be increasing between klow and khigh), 
	// yarray = dependent variable points
	// klow, khigh = lowest,highest indicies of interest 

	int i = klow; int j = khigh - 1;
	while (i < khigh && xarray.at(i) < xlow)		// i = first point > lower integration bound
		i++;
	while (j >= klow && xarray.at(i) > xhigh)		// j = last point < upper integration bound
		j--;

	// Interpolate to find values at lower and upper integration bounds
	double y1 = yarray.at(i);
	if (i>klow)   y1 = yarray.at(i) + (yarray.at(i) - yarray.at(i - 1)) / (xarray.at(i) - xarray.at(i - 1)) * (xlow - xarray.at(i));
	double y2 = yarray.at(j);
	if (j<khigh)   y2 = yarray.at(j) + (yarray.at(j) - yarray.at(j + 1)) / (xarray.at(j) - xarray.at(j + 1)) * (xhigh - xarray.at(j));

	double inteval = 0.0;
	for (int k = i; k < j; k++)		// Intergral between tabulated points entirely included in the integration range
		inteval = inteval + (xarray.at(k + 1) - xarray.at(k)) * 0.5 * (yarray.at(k) + yarray.at(k + 1));
	inteval = inteval + (xarray.at(i) - xlow) * 0.5 * (y1 + yarray.at(i));
	if (j >= i)
		inteval = inteval + (xhigh - xarray.at(j)) * 0.5 * (yarray.at(j) + y2);

	return inteval;
}

void C_mspt_receiver::cubic_splines(const std::vector<double> &xarray, const std::vector<double> &yarray, util::matrix_t<double> &splines)
{
	// Fit cubic splines to data points in xarray, yarray
	int n = xarray.size()-1;
	splines.resize_fill(n, 5, 0.0);

	vector<double> a(n + 1, 0.0);
	vector<double> b(n, 0.0);
	vector<double> d(n, 0.0);
	vector<double> h(n, 0.0);
	vector<double> alpha(n, 0.0);
	vector<double> c(n + 1, 0.0);
	vector<double> l(n + 1, 0.0);
	vector<double> mu(n + 1, 0.0);
	vector<double> z(n + 1, 0.0);

	a = yarray;
	l.at(0) = 1.0;
	mu.at(0) = 0.0;
	z.at(0) = 0.0;
	for (int i = 0; i < n; i++)
	{
		h.at(i) = xarray.at(i + 1) - xarray.at(i);
		if (i > 0)
		{
			alpha.at(i) = (3.0 / h.at(i)) * (a.at(i + 1) - a.at(i)) - (3.0 / h.at(i - 1))*(a.at(i) - a.at(i - 1));
			l.at(i) = 2.0*(xarray.at(i + 1) - xarray.at(i-1)) - h.at(i - 1)*mu.at(i - 1);
			mu.at(i) = h.at(i) / l.at(i);
			z.at(i) = (alpha.at(i) - h.at(i - 1)*z.at(i - 1)) / l.at(i);
		}
	}

	l.at(n) = 1.0;
	z.at(n) = 0.0;
	c.at(n) = 0.0;
	for (int i = n - 1; i >= 0; i--)
	{
		c.at(i) = z.at(i) - mu.at(i)*c.at(i + 1);
		b.at(i) = (a.at(i + 1) - a.at(i)) / h.at(i) - h.at(i)*(c.at(i + 1) + 2.0*c.at(i)) / 3.0;
		d.at(i) = (c.at(i + 1) - c.at(i)) / 3.0 / h.at(i);
	}

	for (int i = 0; i < n; i++)
	{
		splines.at(i, 0) = a.at(i);
		splines.at(i, 1) = b.at(i);
		splines.at(i, 2) = c.at(i);
		splines.at(i, 3) = d.at(i);
		splines.at(i, 4) = xarray.at(i);
	}
}

double C_mspt_receiver::get_startup_time()
{
    return m_rec_su_delay;
}

double C_mspt_receiver::get_startup_energy()
{
    return m_q_rec_des * m_rec_qf_delay;
}