#ifndef __csp_solver_mspt_receiver_
#define __csp_solver_mspt_receiver_

#include "ngcc_powerblock.h"
#include "csp_solver_pt_receiver.h"
#include "csp_solver_util.h"

class C_mspt_receiver : public C_pt_receiver
{
// The transient receiver, including legacy steady-state receiver code for either non-transient startup or non-transient operation

private:

	int m_itermode;
	double m_od_control;
	double m_eta_field_iter_prev;	//[-] Efficiency from heliostat on last iteration. Maybe change if CR gets defocus signal from controller
	double m_tol_od;

	/* declare storage variables here */
	double m_E_su;
	double m_E_su_prev;
	double m_t_su;
	double m_t_su_prev;

    double m_q_dot_inc;
    double m_q_dot_loss;
    double m_q_dot_abs;
    double m_eta_rec_des;
    
    //structures for performance lookup tables
    Bilinear_Interp m_efficiency_lookup;
    Bilinear_Interp m_pressure_lookup;

	// track number of calls per timestep, reset = -1 in converged() call
	int m_ncall;

	struct parameter_eval_inputs
	{
		double T_amb, T_sky, pres, wspd, c_htf, rho_htf, mu_htf, k_htf, Pr_htf, mflow_tot, finitial, ffinal, ramptime;
		std::vector<double> tm;
		util::matrix_t<double> Tfeval, Tseval, qinc;
	} param_inputs;

	double calc_external_convection_coeff(const parameter_eval_inputs &pinputs, double Twall);
	void calc_header_size(double pdrop, double mdot, double rhof, double muf, double Lh, double &id_calc, double &th_calc, double &od_calc);
	double interpolate(double x, const std::vector<double> &xarray, const std::vector<double> &yarray, int klow, int khigh);
	double integrate(double xlow, double xhigh, const std::vector<double> &xarray, const std::vector<double> &yarray, int klow, int khigh);
	void cubic_splines(const std::vector<double> &xarray, const std::vector<double> &yarray, util::matrix_t<double> &splines);
	
    enum startup_modes
	{
		HEAT_TRACE = 0,		// No flux on receiver, riser/downcomer heated with heat tracing
		PREHEAT,			// Low flux on receiver, no HTF flow
		CIRCULATE,			// Full available power on receiver, HTF mass flow rate selected to hit target hot at SS
		HOLD				// Models predict that startup has been completed, but minimum startup time has not yet been reached.  Fluid continues to circulate through the receiver.  
	};

public:
	// Class to save messages for up stream classes
	C_csp_messages csp_messages;

	// Data
	double m_w_rec;					//[m]
	double m_h_rec;					//[m]
	//double m_A_sf;					//[m2]

	// 8.10.2015 twn: add tower piping thermal losses to receiver performance
	double m_pipe_loss_per_m;		//[Wt/m]
	double m_pipe_length_add;		//[m]
	double m_pipe_length_mult;		//[-]

	// 7.13.17 twn: keep this public for now so iscc can calculate
	double m_m_dot_htf_max;			//[kg/s]

    util::matrix_t<double> m_efficiency_lookup_input; 
    util::matrix_t<double> m_pressure_lookup_input;
    
		// 4.17.15 twn: former TCS inputs, moved to member data because are constant throughout simulation
	double m_T_salt_hot_target;			//[C], convert to K in init() call
	double m_hel_stow_deploy;			//[-]

		// Added for csp_solver/tcs wrappers:
	int m_field_fl;
	util::matrix_t<double> m_field_fl_props;	
	
	S_outputs outputs;

	// Methods
	C_mspt_receiver();

	~C_mspt_receiver();

	virtual void init();

	virtual void call(const C_csp_weatherreader::S_outputs &weather, 
		const C_csp_solver_htf_1state &htf_state_in, 
		const C_pt_receiver::S_inputs &inputs,
		const C_csp_solver_sim_info &sim_info);

	virtual void off(const C_csp_weatherreader::S_outputs &weather,
		const C_csp_solver_htf_1state &htf_state_in,
		const C_csp_solver_sim_info &sim_info);

	virtual void converged();

    void calc_pump_performance(double rho_f, double mdot, double ffact, double &PresDrop_calc,
        double &WdotPump_calc);

    virtual double get_pumping_parasitic_coef();

	//void est_startup_time_energy(double fract, double &est_time, double &est_energy);

	//double est_heattrace_energy();

    virtual double get_startup_time();

    virtual double get_startup_energy();

    virtual double area_proj();
};

#endif 
