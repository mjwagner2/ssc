from PySSC import PySSC
import numpy as np
import os
import pandas as pd

#modules with cost/performance functions
import piping
import receiver
import tes
import cycle

"""
####################################################
Optimization variables
####################################################
"""

class Variables:
    def __init__(self):
        self.initialize()

    def initialize(self):
        #initialize variable values
        self.cycle_design_power = 100.        # MWe
        # self.receiver_design_power = 100/.43*3   # MWt
        self.solar_multiple = 3.
        self.dni_design_point = 976.          # W/m2
        self.receiver_height = 5.3           # m
        self.riser_inner_diam = 0.490      # m
        self.downcomer_inner_diam = 0.490      # m
        self.hours_tes = 13                   # hr        
        self.dT_approach_charge_hx = 15       # C  charge hx approach temp
        self.dT_approach_disch_hx = 15        # C  discharge hx total approach temp

class Settings:
    def __init__(self):
        self.print_ssc_messages = False
        self.print_summary_output = False
        self.save_hourly_results = False
        self.dni_des_ref = 976.
        self.cycle_temperature_drop = 700 - 562
        self.lift_technology = 'skip'    #or 'bucket'
        self.is_north = False                  # is field north or surround
        self.cycle_efficiency_nominal = 0.43  #must correspond to the nominal efficiency used to develop the cycle lookup tables
        self.scale_hx_cost = 1.

class Gen3opt:
    def __init__(self):
        self.clear()

    def clear(self):
        self.variables = Variables()
        self.settings = Settings()
        self.summary_results = []
        self.hourly_data = None
    
    #----------------------------------------------------------------
    def load_ssc_constants(self, ssc, data):

        """
        ####################################################
            Configurations and simulation options
        ####################################################
        """

        ssc.data_set_number( data, b'ppa_multiplier_model', 0 );
        ssc.data_set_number( data, b'is_udpc_co2', 1);
        ssc.data_set_number( data, b'time_start', 0);    #24*3600*79 );
        ssc.data_set_number( data, b'time_stop', 31536000 );
        ssc.data_set_number( data, b'is_dispatch', 0 );
        ssc.data_set_number( data, b'is_wlim_series', 0 );
        ssc.data_set_number( data, b'is_dispatch_series', 0 );
        ssc.data_set_number( data, b'store_htf', 34 );

        """
        ####################################################
                    files
        ####################################################
        """
        ssc.data_set_array_from_csv( data, b'wlim_series', b'resource/wlim_series.csv');
        ssc.data_set_array_from_csv( data, b'dispatch_factors_ts', b'resource/dispatch_factors_ts.csv');
        
        """
        ####################################################
                    parameters
        ####################################################
        """
        #design
        ssc.data_set_number( data, b'gross_net_conversion_factor', 0.9 );
        ssc.data_set_number( data, b'P_ref', 100 );
        ssc.data_set_number( data, b'design_eff', .43);     # 0.43; 
        ssc.data_set_number( data, b'tshours', 13 );     # 10.3
        ssc.data_set_number( data, b'solarm',  3);

        #heliostat field
        ssc.data_set_number( data, b'helio_width', 8.66 );
        ssc.data_set_number( data, b'helio_height', 8.66 );
        ssc.data_set_number( data, b'helio_active_fraction', 0.99 );
        ssc.data_set_number( data, b'dens_mirror', 0.97 );
        ssc.data_set_number( data, b'helio_reflectance', 999 );     #0.9 not used in user field mode
        ssc.data_set_number( data, b'rec_absorptance', 999 );        #0.94 not used in user field mode
        ssc.data_set_number( data, b'rec_hl_perm2', 0 );
        ssc.data_set_number( data, b'land_max', 9.5 );
        ssc.data_set_number( data, b'land_min', 0.75 );
        ssc.data_set_number( data, b'dni_des', 950 );
        ssc.data_set_number( data, b'p_start', 0.025 );
        ssc.data_set_number( data, b'p_track', 0.055 );
        ssc.data_set_number( data, b'hel_stow_deploy', 8 );
        ssc.data_set_number( data, b'v_wind_max', 15 );

        #total height and width of all recievers (cost calculation)
        ssc.data_set_number( data, b'rec_height', 22 );     #524.67 m^2
        ssc.data_set_number( data, b'D_rec', 15 );
        ssc.data_set_number( data, b'h_tower', 230 );

        ssc.data_set_number( data, b'water_usage_per_wash', 0.7 );
        ssc.data_set_number( data, b'washing_frequency', 63 );

        ssc.data_set_number( data, b'tower_fixed_cost', 2.3602 * 0.78232e6 );
        ssc.data_set_number( data, b'tower_exp', 0.0113 );
        ssc.data_set_number( data, b'foundation_fixed_cost', 6684590 );
        ssc.data_set_number( data, b'foundation_cost_scaling_quadratic', 154.343 );
        ssc.data_set_number( data, b'foundation_cost_scaling_linear', 115727. );
        ssc.data_set_number( data, b'particle_lift_cost', 0. )  #  60e6 );
        ssc.data_set_number( data, b'riser_and_downcomer_cost', 0. );

        ssc.data_set_number( data, b'rec_ref_cost', 1e7 );
        ssc.data_set_number( data, b'rec_ref_area', 1110 );
        ssc.data_set_number( data, b'rec_cost_exp', 0.7 );

        #field costs
        ssc.data_set_number( data, b'site_spec_cost', 10. );
        ssc.data_set_number( data, b'heliostat_spec_cost', 75. );

        #Plant and BOP
        ssc.data_set_number( data, b'plant_spec_cost', 600 );
        ssc.data_set_number( data, b'bop_spec_cost', 0 );        #<<<< hx included in tes_spec_cost, lift is separate line item

        #TES
        ssc.data_set_number( data, b'tes_spec_cost', 80.)  #152895955./(100000/0.372*15.5) );  #$/kwht

        #land
        ssc.data_set_number( data, b'land_spec_cost', 0); #10000 );
        ssc.data_set_number( data, b'csp.pt.sf.fixed_land_area', 0 );
        ssc.data_set_number( data, b'csp.pt.sf.land_overhead_factor', 1 );
        ssc.data_set_number( data, b'land_area_base', 2822 );         #from spreadsheet

        ssc.data_set_number( data, b'contingency_rate', 7 );
        ssc.data_set_number( data, b'sales_tax_rate', 5 );
        ssc.data_set_number( data, b'sales_tax_frac', 0 );
        ssc.data_set_number( data, b'cost_sf_fixed', 0 );
        ssc.data_set_number( data, b'fossil_spec_cost', 0 );
        ssc.data_set_number( data, b'csp.pt.cost.epc.per_acre', 0 );
        ssc.data_set_number( data, b'csp.pt.cost.epc.percent', 17.6 );
        ssc.data_set_number( data, b'csp.pt.cost.epc.per_watt', 0 );
        ssc.data_set_number( data, b'csp.pt.cost.epc.fixed', 0 );
        ssc.data_set_number( data, b'csp.pt.cost.plm.percent', 0 );
        ssc.data_set_number( data, b'csp.pt.cost.plm.per_watt', 0 );
        ssc.data_set_number( data, b'csp.pt.cost.plm.fixed', 0 );


        #receiver parameters
        ssc.data_set_number( data, b'rec_htf', 5 );
        field_fl_props = [ [ 1, 7, 0, 0, 0, 0, 0, 0, 0 ] ]
        ssc.data_set_matrix( data, b'field_fl_props', field_fl_props );
        ssc.data_set_number( data, b'f_rec_min', 0.05 );
        ssc.data_set_number( data, b'rec_su_delay', 0.1); #0.2 );
        ssc.data_set_number( data, b'rec_qf_delay', 0.1); #0.25 );
        ssc.data_set_number( data, b'csp.pt.rec.max_oper_frac', 2 );
        ssc.data_set_number( data, b'piping_loss', 10200 );
        ssc.data_set_number( data, b'piping_length_mult', 2.6 );
        ssc.data_set_number( data, b'piping_length_const', 0 );
        ssc.data_set_number( data, b'eta_pump', 0.85 );


        ssc.data_set_number( data, b'T_rec_hot_des', 730 );
        ssc.data_set_number( data, b'T_rec_cold_des', 590 );
        store_fl_props = [ [ 1, 7, 0, 0, 0, 0, 0, 0, 0 ] ]
        ssc.data_set_matrix( data, b'store_fl_props', store_fl_props );
        ssc.data_set_number( data, b'tes_pump_coef', 0.15 );
        ssc.data_set_number( data, b'T_tes_hot_des', 715 );
        ssc.data_set_number( data, b'T_tes_warm_des', 615 );
        ssc.data_set_number( data, b'T_tes_cold_des', 565 );
        ssc.data_set_number( data, b'csp.pt.tes.init_hot_htf_percent', 0 );
        ssc.data_set_number( data, b'h_tank', 20 );
        ssc.data_set_number( data, b'cold_tank_max_heat', 0 );
        ssc.data_set_number( data, b'dt_charging', 15 );
        ssc.data_set_number( data, b'dt_ht_discharging', 12 );
        ssc.data_set_number( data, b'dt_lt_discharging', 3 );
        ssc.data_set_number( data, b'dP_LTHX_perc', 0.5 );
        ssc.data_set_number( data, b'dP_HTHX_perc', 1.5 );

        ssc.data_set_number( data, b'u_tank', 0.4 );
        ssc.data_set_number( data, b'tank_pairs', 1 );
        ssc.data_set_number( data, b'cold_tank_Thtr', 280 );
        ssc.data_set_number( data, b'h_tank_min', 1 );
        ssc.data_set_number( data, b'hot_tank_Thtr', 500 );
        ssc.data_set_number( data, b'hot_tank_max_heat', 0 );

        ssc.data_set_number( data, b'T_pc_hot_des', 700 );
        ssc.data_set_number( data, b'T_pc_cold_des', 535 );
        ssc.data_set_number( data, b'pb_pump_coef', 0.55 );
        ssc.data_set_number( data, b'startup_time', 0.5 );
        ssc.data_set_number( data, b'startup_frac', 0.5 );
        ssc.data_set_number( data, b'cycle_max_frac', 1.2); #1.05 );
        ssc.data_set_number( data, b'cycle_cutoff_frac', 0.2 );
        ssc.data_set_number( data, b'q_sby_frac', 0.2 );

        ssc.data_set_number( data, b'ud_T_amb_des', 35 );
        ssc.data_set_number( data, b'ud_f_W_dot_cool_des', 2 );    # was set to 0. Guess for now
        ssc.data_set_number( data, b'ud_m_dot_water_cool_des', 0 );
        ssc.data_set_number( data, b'ud_T_htf_low', 680 );
        ssc.data_set_number( data, b'ud_T_htf_high', 715 );
        ssc.data_set_number( data, b'ud_T_amb_low', 0 );
        ssc.data_set_number( data, b'ud_T_amb_high', 45 );
        ssc.data_set_number( data, b'ud_m_dot_htf_low', 0.5 );
        ssc.data_set_number( data, b'ud_m_dot_htf_high', 1.05 );
        ud_T_htf_ind_od = [[0]]
        ssc.data_set_matrix( data, b'ud_T_htf_ind_od', ud_T_htf_ind_od );
        ud_T_amb_ind_od = [[0]]
        ssc.data_set_matrix( data, b'ud_T_amb_ind_od', ud_T_amb_ind_od );
        ud_m_dot_htf_ind_od = [[0]]
        ssc.data_set_matrix( data, b'ud_m_dot_htf_ind_od', ud_m_dot_htf_ind_od );
        ssc.data_set_number( data, b'P_phx_in_co2_des', 24750.625 );
        ssc.data_set_number( data, b'P_turb_in_co2_des', 20790.525 );
        ssc.data_set_number( data, b'P_turb_in_co2_off_sun_des', 24504.);

        ssc.data_set_number( data, b'pb_fixed_par', 0); #0.0055 );
        ssc.data_set_number( data, b'bop_par', 0 );
        ssc.data_set_number( data, b'bop_par_f', 1 );
        ssc.data_set_number( data, b'bop_par_0', 0 );
        ssc.data_set_number( data, b'bop_par_1', 0.483 );
        ssc.data_set_number( data, b'bop_par_2', 0 );
        f_turb_tou_periods = [ 0.99, 1, 1, 1, 1, 1, 1, 1, 1 ]
        ssc.data_set_array( data, b'f_turb_tou_periods', f_turb_tou_periods );
        weekday_schedule = \
            [ [ 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5 ], 
            [ 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3 ], 
            [ 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3 ], 
            [ 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3 ], 
            [ 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5 ] ] 
        ssc.data_set_matrix( data, b'weekday_schedule', weekday_schedule );
        weekend_schedule = \
            [ [ 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ], 
            [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ], 
            [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ], 
            [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ], 
            [ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ], 
            [ 6, 6, 6, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ] ] 
        ssc.data_set_matrix( data, b'weekend_schedule', weekend_schedule);
        ssc.data_set_number( data, b'disp_horizon', 48 );
        ssc.data_set_number( data, b'disp_frequency', 24 );
        ssc.data_set_number( data, b'disp_max_iter', 35000 );
        ssc.data_set_number( data, b'disp_timeout', 5 );
        ssc.data_set_number( data, b'disp_mip_gap', 0.001 );
        ssc.data_set_number( data, b'disp_time_weighting', 0.99 );
        ssc.data_set_number( data, b'disp_rsu_cost', 950 );
        ssc.data_set_number( data, b'disp_csu_cost', 10000 );
        ssc.data_set_number( data, b'disp_pen_delta_w', 0.1 );
        ssc.data_set_matrix( data, b'dispatch_sched_weekday', weekday_schedule);
        ssc.data_set_matrix( data, b'dispatch_sched_weekend', weekend_schedule);
        ssc.data_set_number( data, b'dispatch_factor1', 1 );
        ssc.data_set_number( data, b'dispatch_factor2', 1 );
        ssc.data_set_number( data, b'dispatch_factor3', 1 );
        ssc.data_set_number( data, b'dispatch_factor4', 1 );
        ssc.data_set_number( data, b'dispatch_factor5', 1 );
        ssc.data_set_number( data, b'dispatch_factor6', 1 );
        ssc.data_set_number( data, b'dispatch_factor7', 1 );
        ssc.data_set_number( data, b'dispatch_factor8', 1 );
        ssc.data_set_number( data, b'dispatch_factor9', 1 );
        dispatch_series = [ 0 ]
        ssc.data_set_array( data, b'dispatch_series', dispatch_series );
        ssc.data_set_number( data, b'const_per_interest_rate1', 4 );
        ssc.data_set_number( data, b'const_per_interest_rate2', 0 );
        ssc.data_set_number( data, b'const_per_interest_rate3', 0 );
        ssc.data_set_number( data, b'const_per_interest_rate4', 0 );
        ssc.data_set_number( data, b'const_per_interest_rate5', 0 );
        ssc.data_set_number( data, b'const_per_months1', 24 );
        ssc.data_set_number( data, b'const_per_months2', 0 );
        ssc.data_set_number( data, b'const_per_months3', 0 );
        ssc.data_set_number( data, b'const_per_months4', 0 );
        ssc.data_set_number( data, b'const_per_months5', 0 );
        ssc.data_set_number( data, b'const_per_percent1', 100 );
        ssc.data_set_number( data, b'const_per_percent2', 0 );
        ssc.data_set_number( data, b'const_per_percent3', 0 );
        ssc.data_set_number( data, b'const_per_percent4', 0 );
        ssc.data_set_number( data, b'const_per_percent5', 0 );
        ssc.data_set_number( data, b'const_per_upfront_rate1', 1 );
        ssc.data_set_number( data, b'const_per_upfront_rate2', 0 );
        ssc.data_set_number( data, b'const_per_upfront_rate3', 0 );
        ssc.data_set_number( data, b'const_per_upfront_rate4', 0 );
        ssc.data_set_number( data, b'const_per_upfront_rate5', 0 );
        ssc.data_set_number( data, b'adjust:constant', 4 );
        ssc.data_set_number( data, b'sf_adjust:constant', 0 );

        #------------------------------------------------------------------------------

        ssc.data_set_number( data, b'analysis_period', 30 ); 
        federal_tax_rate = [ 35 ]
        ssc.data_set_array( data, b'federal_tax_rate', federal_tax_rate );
        state_tax_rate = [ 5 ]
        ssc.data_set_array( data, b'state_tax_rate', state_tax_rate );
        ssc.data_set_number( data, b'property_tax_rate', 0 );
        ssc.data_set_number( data, b'prop_tax_cost_assessed_percent', 100 );
        ssc.data_set_number( data, b'prop_tax_assessed_decline', 0 );
        ssc.data_set_number( data, b'real_discount_rate', 4.4 );
        ssc.data_set_number( data, b'inflation_rate', 2.5 );
        ssc.data_set_number( data, b'insurance_rate', 0.5 );
        ssc.data_set_number( data, b'system_capacity', 1e5 );
        om_fixed = [ 0 ]
        ssc.data_set_array( data, b'om_fixed', om_fixed );
        ssc.data_set_number( data, b'om_fixed_escal', 0 );
        om_production = [4]
        ssc.data_set_array( data, b'om_production', om_production); 
        ssc.data_set_number( data, b'om_production_escal', 0 );
        om_capacity = [ 40 ]
        ssc.data_set_array( data, b'om_capacity', om_capacity );    
        ssc.data_set_number( data, b'om_capacity_escal', 0 );
        om_fuel_cost = [ 0 ]
        ssc.data_set_array( data, b'om_fuel_cost', om_fuel_cost );
        ssc.data_set_number( data, b'om_fuel_cost_escal', 0 );
        om_replacement_cost1 = [ 0 ]
        ssc.data_set_array( data, b'om_replacement_cost1', om_replacement_cost1 );
        ssc.data_set_number( data, b'om_replacement_cost_escal', 0 );

        ssc.data_set_number( data, b'itc_fed_amount', 0 );
        ssc.data_set_number( data, b'itc_fed_amount_deprbas_fed', 1 );
        ssc.data_set_number( data, b'itc_fed_amount_deprbas_sta', 1 );
        ssc.data_set_number( data, b'itc_sta_amount', 0 );
        ssc.data_set_number( data, b'itc_sta_amount_deprbas_fed', 0 );
        ssc.data_set_number( data, b'itc_sta_amount_deprbas_sta', 0 );
        ssc.data_set_number( data, b'itc_fed_percent', 0 );
        ssc.data_set_number( data, b'itc_fed_percent_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'itc_fed_percent_deprbas_fed', 1 );
        ssc.data_set_number( data, b'itc_fed_percent_deprbas_sta', 1 );
        ssc.data_set_number( data, b'itc_sta_percent', 0 );
        ssc.data_set_number( data, b'itc_sta_percent_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'itc_sta_percent_deprbas_fed', 0 );
        ssc.data_set_number( data, b'itc_sta_percent_deprbas_sta', 0 );
        ptc_fed_amount = [ 0 ]
        ssc.data_set_array( data, b'ptc_fed_amount', ptc_fed_amount );
        ssc.data_set_number( data, b'ptc_fed_term', 10 );
        ssc.data_set_number( data, b'ptc_fed_escal', 0 );
        ptc_sta_amount = [ 0 ]
        ssc.data_set_array( data, b'ptc_sta_amount', ptc_sta_amount );
        ssc.data_set_number( data, b'ptc_sta_term', 10 );
        ssc.data_set_number( data, b'ptc_sta_escal', 0 );
        ssc.data_set_number( data, b'ibi_fed_amount', 0 );
        ssc.data_set_number( data, b'ibi_fed_amount_tax_fed', 1 );
        ssc.data_set_number( data, b'ibi_fed_amount_tax_sta', 1 );
        ssc.data_set_number( data, b'ibi_fed_amount_deprbas_fed', 0 );
        ssc.data_set_number( data, b'ibi_fed_amount_deprbas_sta', 0 );
        ssc.data_set_number( data, b'ibi_sta_amount', 0 );
        ssc.data_set_number( data, b'ibi_sta_amount_tax_fed', 1 );
        ssc.data_set_number( data, b'ibi_sta_amount_tax_sta', 1 );
        ssc.data_set_number( data, b'ibi_sta_amount_deprbas_fed', 0 );
        ssc.data_set_number( data, b'ibi_sta_amount_deprbas_sta', 0 );
        ssc.data_set_number( data, b'ibi_uti_amount', 0 );
        ssc.data_set_number( data, b'ibi_uti_amount_tax_fed', 1 );
        ssc.data_set_number( data, b'ibi_uti_amount_tax_sta', 1 );
        ssc.data_set_number( data, b'ibi_uti_amount_deprbas_fed', 0 );
        ssc.data_set_number( data, b'ibi_uti_amount_deprbas_sta', 0 );
        ssc.data_set_number( data, b'ibi_oth_amount', 0 );
        ssc.data_set_number( data, b'ibi_oth_amount_tax_fed', 1 );
        ssc.data_set_number( data, b'ibi_oth_amount_tax_sta', 1 );
        ssc.data_set_number( data, b'ibi_oth_amount_deprbas_fed', 0 );
        ssc.data_set_number( data, b'ibi_oth_amount_deprbas_sta', 0 );
        ssc.data_set_number( data, b'ibi_fed_percent', 0 );
        ssc.data_set_number( data, b'ibi_fed_percent_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'ibi_fed_percent_tax_fed', 1 );
        ssc.data_set_number( data, b'ibi_fed_percent_tax_sta', 1 );
        ssc.data_set_number( data, b'ibi_fed_percent_deprbas_fed', 0 );
        ssc.data_set_number( data, b'ibi_fed_percent_deprbas_sta', 0 );
        ssc.data_set_number( data, b'ibi_sta_percent', 0 );
        ssc.data_set_number( data, b'ibi_sta_percent_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'ibi_sta_percent_tax_fed', 1 );
        ssc.data_set_number( data, b'ibi_sta_percent_tax_sta', 1 );
        ssc.data_set_number( data, b'ibi_sta_percent_deprbas_fed', 0 );
        ssc.data_set_number( data, b'ibi_sta_percent_deprbas_sta', 0 );
        ssc.data_set_number( data, b'ibi_uti_percent', 0 );
        ssc.data_set_number( data, b'ibi_uti_percent_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'ibi_uti_percent_tax_fed', 1 );
        ssc.data_set_number( data, b'ibi_uti_percent_tax_sta', 1 );
        ssc.data_set_number( data, b'ibi_uti_percent_deprbas_fed', 0 );
        ssc.data_set_number( data, b'ibi_uti_percent_deprbas_sta', 0 );
        ssc.data_set_number( data, b'ibi_oth_percent', 0 );
        ssc.data_set_number( data, b'ibi_oth_percent_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'ibi_oth_percent_tax_fed', 1 );
        ssc.data_set_number( data, b'ibi_oth_percent_tax_sta', 1 );
        ssc.data_set_number( data, b'ibi_oth_percent_deprbas_fed', 0 );
        ssc.data_set_number( data, b'ibi_oth_percent_deprbas_sta', 0 );
        ssc.data_set_number( data, b'cbi_fed_amount', 0 );
        ssc.data_set_number( data, b'cbi_fed_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'cbi_fed_tax_fed', 1 );
        ssc.data_set_number( data, b'cbi_fed_tax_sta', 1 );
        ssc.data_set_number( data, b'cbi_fed_deprbas_fed', 0 );
        ssc.data_set_number( data, b'cbi_fed_deprbas_sta', 0 );
        ssc.data_set_number( data, b'cbi_sta_amount', 0 );
        ssc.data_set_number( data, b'cbi_sta_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'cbi_sta_tax_fed', 1 );
        ssc.data_set_number( data, b'cbi_sta_tax_sta', 1 );
        ssc.data_set_number( data, b'cbi_sta_deprbas_fed', 0 );
        ssc.data_set_number( data, b'cbi_sta_deprbas_sta', 0 );
        ssc.data_set_number( data, b'cbi_uti_amount', 0 );
        ssc.data_set_number( data, b'cbi_uti_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'cbi_uti_tax_fed', 1 );
        ssc.data_set_number( data, b'cbi_uti_tax_sta', 1 );
        ssc.data_set_number( data, b'cbi_uti_deprbas_fed', 0 );
        ssc.data_set_number( data, b'cbi_uti_deprbas_sta', 0 );
        ssc.data_set_number( data, b'cbi_oth_amount', 0 );
        ssc.data_set_number( data, b'cbi_oth_maxvalue', 1.e38 );
        ssc.data_set_number( data, b'cbi_oth_tax_fed', 1 );
        ssc.data_set_number( data, b'cbi_oth_tax_sta', 1 );
        ssc.data_set_number( data, b'cbi_oth_deprbas_fed', 0 );
        ssc.data_set_number( data, b'cbi_oth_deprbas_sta', 0 );
        pbi_fed_amount = [ 0 ]
        ssc.data_set_array( data, b'pbi_fed_amount', pbi_fed_amount );
        ssc.data_set_number( data, b'pbi_fed_term', 0 );
        ssc.data_set_number( data, b'pbi_fed_escal', 0 );
        ssc.data_set_number( data, b'pbi_fed_tax_fed', 1 );
        ssc.data_set_number( data, b'pbi_fed_tax_sta', 1 );
        pbi_sta_amount = [ 0 ]
        ssc.data_set_array( data, b'pbi_sta_amount', pbi_sta_amount );
        ssc.data_set_number( data, b'pbi_sta_term', 0 );
        ssc.data_set_number( data, b'pbi_sta_escal', 0 );
        ssc.data_set_number( data, b'pbi_sta_tax_fed', 1 );
        ssc.data_set_number( data, b'pbi_sta_tax_sta', 1 );
        pbi_uti_amount = [ 0 ]
        ssc.data_set_array( data, b'pbi_uti_amount', pbi_uti_amount );
        ssc.data_set_number( data, b'pbi_uti_term', 0 );
        ssc.data_set_number( data, b'pbi_uti_escal', 0 );
        ssc.data_set_number( data, b'pbi_uti_tax_fed', 1 );
        ssc.data_set_number( data, b'pbi_uti_tax_sta', 1 );
        pbi_oth_amount = [ 0 ]
        ssc.data_set_array( data, b'pbi_oth_amount', pbi_oth_amount );
        ssc.data_set_number( data, b'pbi_oth_term', 0 );
        ssc.data_set_number( data, b'pbi_oth_escal', 0 );
        ssc.data_set_number( data, b'pbi_oth_tax_fed', 1 );
        ssc.data_set_number( data, b'pbi_oth_tax_sta', 1 );
        degradation = [ 0 ]
        ssc.data_set_array( data, b'degradation', degradation );
        roe_input = [ 0 ]
        ssc.data_set_array( data, b'roe_input', roe_input );
        ssc.data_set_number( data, b'loan_moratorium', 0 );
        ssc.data_set_number( data, b'system_use_recapitalization', 0 );
        ssc.data_set_number( data, b'system_use_lifetime_output', 0 );
        #ssc.data_set_number( data, b'total_installed_cost', 673465472 );  #calculated by performance module
        ssc.data_set_number( data, b'reserves_interest', 0 );
        ssc.data_set_number( data, b'equip1_reserve_cost', 0 );
        ssc.data_set_number( data, b'equip1_reserve_freq', 12 );
        ssc.data_set_number( data, b'equip2_reserve_cost', 0 );
        ssc.data_set_number( data, b'equip2_reserve_freq', 15 );
        ssc.data_set_number( data, b'equip3_reserve_cost', 0 );
        ssc.data_set_number( data, b'equip3_reserve_freq', 3 );
        ssc.data_set_number( data, b'equip_reserve_depr_sta', 0 );
        ssc.data_set_number( data, b'equip_reserve_depr_fed', 0 );
        ssc.data_set_number( data, b'salvage_percentage', 0 );
        ssc.data_set_number( data, b'ppa_soln_mode', 0 );
        ppa_price_input = [0.13]
        ssc.data_set_array( data, b'ppa_price_input', ppa_price_input );
        ssc.data_set_number( data, b'ppa_escalation', 1 );
        #ssc.data_set_number( data, b'construction_financing_cost', 33673272 ); #calculated by performance module
        ssc.data_set_number( data, b'term_tenor', 18 );
        ssc.data_set_number( data, b'term_int_rate', 7 );
        ssc.data_set_number( data, b'dscr', 1.3 );
        ssc.data_set_number( data, b'dscr_reserve_months', 6 );
        ssc.data_set_number( data, b'debt_percent', 50 );
        ssc.data_set_number( data, b'debt_option', 1 );
        ssc.data_set_number( data, b'payment_option', 0 );
        ssc.data_set_number( data, b'cost_debt_closing', 0 );
        ssc.data_set_number( data, b'cost_debt_fee', 0 );
        ssc.data_set_number( data, b'months_working_reserve', 6 );
        ssc.data_set_number( data, b'months_receivables_reserve', 0 );
        ssc.data_set_number( data, b'cost_other_financing', 0 );
        ssc.data_set_number( data, b'flip_target_percent', 7 );
        ssc.data_set_number( data, b'flip_target_year', 30 );
        ssc.data_set_number( data, b'depr_alloc_macrs_5_percent', 90 );
        ssc.data_set_number( data, b'depr_alloc_macrs_15_percent', 1.5 );
        ssc.data_set_number( data, b'depr_alloc_sl_5_percent', 0 );
        ssc.data_set_number( data, b'depr_alloc_sl_15_percent', 2.5 );
        ssc.data_set_number( data, b'depr_alloc_sl_20_percent', 3 );
        ssc.data_set_number( data, b'depr_alloc_sl_39_percent', 0 );
        ssc.data_set_number( data, b'depr_alloc_custom_percent', 0 );
        depr_custom_schedule = [ 0 ]
        ssc.data_set_array( data, b'depr_custom_schedule', depr_custom_schedule );
        ssc.data_set_number( data, b'depr_bonus_sta', 0 );
        ssc.data_set_number( data, b'depr_bonus_sta_macrs_5', 1 );
        ssc.data_set_number( data, b'depr_bonus_sta_macrs_15', 1 );
        ssc.data_set_number( data, b'depr_bonus_sta_sl_5', 0 );
        ssc.data_set_number( data, b'depr_bonus_sta_sl_15', 0 );
        ssc.data_set_number( data, b'depr_bonus_sta_sl_20', 0 );
        ssc.data_set_number( data, b'depr_bonus_sta_sl_39', 0 );
        ssc.data_set_number( data, b'depr_bonus_sta_custom', 0 );
        ssc.data_set_number( data, b'depr_bonus_fed', 0 );
        ssc.data_set_number( data, b'depr_bonus_fed_macrs_5', 1 );
        ssc.data_set_number( data, b'depr_bonus_fed_macrs_15', 1 );
        ssc.data_set_number( data, b'depr_bonus_fed_sl_5', 0 );
        ssc.data_set_number( data, b'depr_bonus_fed_sl_15', 0 );
        ssc.data_set_number( data, b'depr_bonus_fed_sl_20', 0 );
        ssc.data_set_number( data, b'depr_bonus_fed_sl_39', 0 );
        ssc.data_set_number( data, b'depr_bonus_fed_custom', 0 );
        ssc.data_set_number( data, b'depr_itc_sta_macrs_5', 1 );
        ssc.data_set_number( data, b'depr_itc_sta_macrs_15', 0 );
        ssc.data_set_number( data, b'depr_itc_sta_sl_5', 0 );
        ssc.data_set_number( data, b'depr_itc_sta_sl_15', 0 );
        ssc.data_set_number( data, b'depr_itc_sta_sl_20', 0 );
        ssc.data_set_number( data, b'depr_itc_sta_sl_39', 0 );
        ssc.data_set_number( data, b'depr_itc_sta_custom', 0 );
        ssc.data_set_number( data, b'depr_itc_fed_macrs_5', 1 );
        ssc.data_set_number( data, b'depr_itc_fed_macrs_15', 0 );
        ssc.data_set_number( data, b'depr_itc_fed_sl_5', 0 );
        ssc.data_set_number( data, b'depr_itc_fed_sl_15', 0 );
        ssc.data_set_number( data, b'depr_itc_fed_sl_20', 0 );
        ssc.data_set_number( data, b'depr_itc_fed_sl_39', 0 );
        ssc.data_set_number( data, b'depr_itc_fed_custom', 0 );
        ssc.data_set_number( data, b'pbi_fed_for_ds', 0 );
        ssc.data_set_number( data, b'pbi_sta_for_ds', 0 );
        ssc.data_set_number( data, b'pbi_uti_for_ds', 0 );
        ssc.data_set_number( data, b'pbi_oth_for_ds', 0 );
        ssc.data_set_number( data, b'depr_stabas_method', 1 );
        ssc.data_set_number( data, b'depr_fedbas_method', 1 );
        ssc.data_set_number( data, b'cp_capacity_payment_esc', 0 );
        ssc.data_set_number( data, b'cp_capacity_payment_type', 0 );
        cp_capacity_payment_amount = [ 0 ]
        ssc.data_set_array( data, b'cp_capacity_payment_amount', cp_capacity_payment_amount );
        cp_capacity_credit_percent = [ 0 ]
        ssc.data_set_array( data, b'cp_capacity_credit_percent', cp_capacity_credit_percent );
        ssc.data_set_number( data, b'cp_system_nameplate', 90 );
        ssc.data_set_number( data, b'cp_battery_nameplate', 0 );
        grid_curtailment_price = [ 0 ]
        ssc.data_set_array( data, b'grid_curtailment_price', grid_curtailment_price );
        ssc.data_set_number( data, b'grid_curtailment_price_esc', 0 );

        return

    #----------------------------------------------------------------
    def exec(self):

        ssc = PySSC()
        if self.settings.print_ssc_messages:
            print('Process ID ', os.getpid())
            print ('Current folder = ' + os.getcwd() )
            print ('SSC Version = ', ssc.version())
            print ('SSC Build Information = ', ssc.build_info().decode("utf - 8"))
            ssc.module_exec_set_print(1)
        else:
            ssc.module_exec_set_print(0)

        data = ssc.data_create()

        #load in the default parameters
        self.load_ssc_constants(ssc, data)

        """
        ####################################################
                    files
        ####################################################
        """
        ssc.data_set_string( data, b'solar_resource_file', b'resource/daggett_ca_34.865371_-116.783023_psmv3_60_tmy.csv' );
        
        #calculate system temperatures
        T_rec_hot_des = 730  #C this is fixed
        T_pc_hot_des = T_rec_hot_des - (self.variables.dT_approach_charge_hx + self.variables.dT_approach_disch_hx)
        T_rec_cold_des = T_rec_hot_des - self.settings.cycle_temperature_drop
        T_tes_hot_des = T_rec_hot_des - self.variables.dT_approach_charge_hx
        T_tes_cold_des = T_tes_hot_des - self.settings.cycle_temperature_drop 
        T_tes_warm_des = T_tes_cold_des + self.variables.dT_approach_disch_hx
        T_pc_cold_des = T_pc_hot_des - self.settings.cycle_temperature_drop

        #cycle efficiency
        cycle_efficiency = self.settings.cycle_efficiency_nominal/0.489 * cycle.calculate_nominal_efficiency(T_pc_hot_des, self.variables.cycle_design_power*1000.)
        ssc.data_set_matrix( data, b'ud_ind_od', cycle.create_updc_lookup('resource/ud_ind_od.csv', T_pc_hot_des) );
        ssc.data_set_matrix( data, b'ud_ind_od_off_sun', cycle.create_updc_lookup('resource/ud_ind_od_off_sun.csv', T_pc_hot_des) );

        # Do initial calculations for parameters used in cost/performance models
        helio_area = 8.66**2*.97
        q_pb_des = self.variables.cycle_design_power/cycle_efficiency  #MWt
        receiver_design_power = q_pb_des * self.variables.solar_multiple
        receiver_eff_des = receiver.calculate_efficiency(self.variables.receiver_height)
        q_sf_des = receiver_design_power / receiver_eff_des * self.settings.dni_des_ref / self.variables.dni_design_point 

        #get the heliostat field for the rated solar field power
        eta_map = receiver.create_heliostat_field_lookup('resource/eta_lookup_{:s}.csv'.format('north' if self.settings.is_north else 'surround'), 
                                                q_sf_des*1000, helio_area)
        ssc.data_set_matrix( data, b'eta_map', eta_map);

        #tower height
        tht = receiver.calculate_tower_height(q_sf_des*1000, self.settings.is_north)

        #riser cost
        L_riser = tht * ssc.data_get_number(data, b'piping_length_mult') + ssc.data_get_number(data, b'piping_length_const')
        riser_cost = piping.solve(self.variables.riser_inner_diam, L_riser, ssc.data_get_number( data, b'P_phx_in_co2_des'))['cost']
        downcomer_cost = piping.solve(self.variables.downcomer_inner_diam, L_riser, ssc.data_get_number( data, b'P_phx_in_co2_des'))['cost']

        #receiver
        ntd = receiver.calculate_n_tubes(receiver_design_power*1000, T_rec_cold_des, T_rec_hot_des, self.variables.receiver_height)
        n_tubes = ntd['n_tubes']
        recd = receiver.calculate_cost(self.variables.receiver_height, n_tubes)
        rec_total_cost = recd['total_cost'] 
        rec_area = recd['A_rec']
        D_rec = recd['W_rec']
        ssc.data_set_matrix(data, b'rec_efficiency_lookup', receiver.create_receiver_efficiency_lookup("resource/rec_efficiency.csv", self.variables.receiver_height) )
        ssc.data_set_matrix(data, b'rec_pressure_lookup', receiver.create_receiver_pressure_lookup("resource/rec_pressure.csv", self.variables.receiver_height) )

        #lift power and cost
        m_dot_p = receiver_design_power*1e3 / (tes.cp_particle() * (T_tes_hot_des - T_tes_cold_des))  #kg/s
        lift_cost = tes.calculate_lift_cost(receiver_design_power*1000, self.settings.lift_technology)
        lift_eff = tes.calculate_lift_efficiency(q_sf_des*1000, receiver_design_power*1000, m_dot_p, self.settings.lift_technology)

        #TES costs
        e_tes = self.variables.hours_tes * q_pb_des * 1000  #kWh
        dtes = tes.calculate_silo_cost( q_pb_des*1000, self.variables.hours_tes, self.settings.cycle_temperature_drop)
        dhx = tes.calculate_hx_cost(q_pb_des*1000, self.variables.dT_approach_charge_hx, self.variables.dT_approach_disch_hx, T_rec_hot_des, T_rec_cold_des, self.settings.scale_hx_cost)
        hx_cost = dhx['total_cost']
        # tes_spec_cost = (hx_cost + dtes['silo_cost'] + dtes['media_cost'])/e_tes + 5  #need to assume some overage (5) for balance of equipment. 
        tes_spec_cost = (hx_cost + dtes['media_cost'])/e_tes + 5  #need to assume some overage (5) for balance of equipment. 

        """
        ####################################################
                    parameters
        ####################################################
        """
        #design
        ssc.data_set_number( data, b'P_ref', self.variables.cycle_design_power );
        ssc.data_set_number( data, b'design_eff', cycle_efficiency);     # 0.43; 
        ssc.data_set_number( data, b'tshours', self.variables.hours_tes );     # 10.3
        ssc.data_set_number( data, b'solarm',  self.variables.solar_multiple);

        #heliostat field
        ssc.data_set_number( data, b'helio_width', 8.66 );
        ssc.data_set_number( data, b'helio_height', 8.66 );
        ssc.data_set_number( data, b'dens_mirror', 0.97 );
        ssc.data_set_number( data, b'dni_des', self.variables.dni_design_point );

        #total height and width of all recievers (cost calculation)
        ssc.data_set_number( data, b'rec_height', self.variables.receiver_height );     #524.67 m^2
        ssc.data_set_number( data, b'D_rec', D_rec );
        ssc.data_set_number( data, b'h_tower', tht );

        ssc.data_set_number( data, b'tower_fixed_cost', 2.3602 * 0.78232e6 );
        ssc.data_set_number( data, b'tower_exp', 0.0113 );
        ssc.data_set_number( data, b'foundation_fixed_cost', 6684590 );
        ssc.data_set_number( data, b'foundation_cost_scaling_quadratic', 154.343 );
        ssc.data_set_number( data, b'foundation_cost_scaling_linear', 115727. );
        ssc.data_set_number( data, b'particle_lift_cost', lift_cost )  #  60e6 );
        ssc.data_set_number( data, b'riser_and_downcomer_cost',  riser_cost + downcomer_cost );

        ssc.data_set_number( data, b'rec_ref_cost', rec_total_cost );
        ssc.data_set_number( data, b'rec_ref_area', rec_area );

        #field costs
        ssc.data_set_number( data, b'site_spec_cost', 10. );
        ssc.data_set_number( data, b'heliostat_spec_cost', 75. );

        #Plant and BOP
        ssc.data_set_number( data, b'plant_spec_cost', 600 );

        #TES
        ssc.data_set_number( data, b'tes_spec_cost', tes_spec_cost)  #$/kwht

        #land
        ssc.data_set_number( data, b'contingency_rate', 7 );
        ssc.data_set_number( data, b'csp.pt.cost.epc.percent', 17.6 );


        #receiver parameters
        ssc.data_set_number( data, b'f_rec_min', 0.05 );
        ssc.data_set_number( data, b'rec_su_delay', 0.1); #0.2 );
        ssc.data_set_number( data, b'rec_qf_delay', 0.1); #0.25 );
        ssc.data_set_number( data, b'csp.pt.rec.max_oper_frac', 1.1 );
        ssc.data_set_number( data, b'piping_loss', 10200 );
        ssc.data_set_number( data, b'piping_length_mult', 1.5 );
        ssc.data_set_number( data, b'piping_length_const', 50 );
        ssc.data_set_number( data, b'piping_riser_diam', self.variables.riser_inner_diam );
        ssc.data_set_number( data, b'piping_downcomer_diam', self.variables.downcomer_inner_diam );
        ssc.data_set_number( data, b'eta_pump', lift_eff );


        ssc.data_set_number( data, b'T_rec_hot_des', T_rec_hot_des );
        ssc.data_set_number( data, b'T_rec_cold_des', T_rec_cold_des );
        ssc.data_set_number( data, b'T_tes_hot_des', T_tes_hot_des );
        ssc.data_set_number( data, b'T_tes_warm_des', T_tes_warm_des );
        ssc.data_set_number( data, b'T_tes_cold_des', T_tes_cold_des );
        ssc.data_set_number( data, b'dt_charging', self.variables.dT_approach_charge_hx );
        ssc.data_set_number( data, b'dt_ht_discharging', 0.8 * self.variables.dT_approach_disch_hx );
        ssc.data_set_number( data, b'dt_lt_discharging', 0.2 * self.variables.dT_approach_disch_hx );
        ssc.data_set_number( data, b'dP_LTHX_perc', 0.5 );
        ssc.data_set_number( data, b'dP_HTHX_perc', 1.5 );

        ssc.data_set_number( data, b'T_pc_hot_des', T_pc_hot_des );
        ssc.data_set_number( data, b'T_pc_cold_des', T_pc_cold_des );

        #------------------------------------------------------------------------------

        module = ssc.module_create(b'tcsmolten_salt') 
        ssc.module_exec_set_print( 0 );
        if ssc.module_exec(module, data) == 0:
            print ('tcsmolten_salt simulation error')
            idx = 1
            msg = ssc.module_log(module, 0)
            while (msg != None):
                print ('    : ' + msg.decode("utf - 8"))
                msg = ssc.module_log(module, idx)
                idx = idx + 1
            SystemExit( "Simulation Error" );
        ssc.module_free(module)

        #------------------------------------------------------------------------------
        #------------------------------------------------------------------------------

        p_ref = ssc.data_get_number( data, b'P_ref' );
        ssc.data_set_number( data, b'system_capacity', p_ref*1000. );

        #ssc.data_set_number( data, b'construction_financing_cost', 33673272 ); #calculated by performance module
        gross_to_net = ssc.data_get_number( data, b'gross_net_conversion_factor' );
        ssc.data_set_number( data, b'cp_system_nameplate', p_ref * gross_to_net );

        #------------------------------------------------------------------------------
        module = ssc.module_create(b'singleowner')
        ssc.module_exec_set_print( 0 );
        if ssc.module_exec(module, data) == 0:
            print ('singleowner simulation error')
            idx = 1
            msg = ssc.module_log(module, 0)
            while (msg != None):
                print ('    : ' + msg.decode("utf - 8"))
                msg = ssc.module_log(module, idx)
                idx = idx + 1
            SystemExit( "Simulation Error" );
        ssc.module_free(module)

        #------------------------------------------------------------------------------
        #------------------------------------------------------------------------------

        if self.settings.save_hourly_results:
        
            #outputs
            self.hourly_output_columns=[
                ["a_sf1", "m2", 1],
                ["a_sf2", "m2", 1],
                ["a_sf3", "m2", 1],
                ["beam", "W/m2", 1],
                ["defocus", "-", 1],
                ["dp_rec1", "kPa", 1],
                ["dp_rec2", "kPa", 1],
                ["dp_rec3", "kPa", 1],
                ["dp_riser", "kPa", 1],
                ["dp_downcomer", "kPa", 1],
                ["e_ch_tes", "MWh", 1],
                ["eta", "-", 1],
                ["eta_field_tot", "-", 1],
                ["eta_field1", "-", 1],
                ["eta_field2", "-", 1],
                ["eta_field3", "-", 1],
                ["eta_therm", "-", 1],
                ["eta_rec_therm1", "-", 1],
                ["eta_rec_therm2", "-", 1],
                ["eta_rec_therm3", "-", 1],
                ["htf_pump_power", "MW", 1],
                ["m_dot_rec", "kg/s", 1],
                ["m_dot_rec1", "kg/s", 0.000277778],
                ["m_dot_rec2", "kg/s", 0.000277778],
                ["m_dot_rec3", "kg/s", 0.000277778],
                ["m_dot_pc", "kg/s", 0.000277778],
                ["m_dot_tes_dc", "kg/s", 0.000277778],
                ["m_dot_tes_ch", "kg/s", 0.000277778],
                ["op_mode_1", "", 1],
                ["op_mode_2", "", 1],
                ["op_mode_3", "", 1],
                ["p_cooling_tower_tot", "MW", 1],
                ["p_cycle", "MW", 1],
                ["p_out_net", "MW", 1],
                ["p_phx_in", "kPa", 1],
                ["p_phx_out", "kPa", 1],
                ["p_tower_pump", "MW", 1],
                ["pparasi", "MW", 1],
                ["q_ch_tes", "MW", 1],
                ["q_dc_tes", "MW", 1],
                ["q_dot_hx1", "MW", 1],
                ["q_dot_hx2", "MW", 1],
                ["q_dot_hx3", "MW", 1],
                ["q_dot_rec_inc", "MW", 0.001],
                ["q_dot_rec_inc1", "MW", 0.001],
                ["q_dot_rec_inc2", "MW", 0.001],
                ["q_dot_rec_inc3", "MW", 0.001],
                ["q_pb", "MW", 1],
                ["q_thermal", "MW", 1],
                ["q_startup", "MW", 1],
                ["tank_losses", "MW", 1],
                ["t_hx_tes_out1", "C", 1],
                ["t_hx_tes_out2", "C", 1],
                ["t_hx_tes_out3", "C", 1],
                ["t_pc_in", "C", 1],
                ["t_pc_out", "C", 1],
                ["t_rec_in1", "C", 1],
                ["t_rec_in2", "C", 1],
                ["t_rec_in3", "C", 1],
                ["t_rec_out1", "C", 1],
                ["t_rec_out2", "C", 1],
                ["t_rec_out3", "C", 1],
                ["t_tes_cold", "C", 1],
                ["t_tes_hot", "C", 1],
                ["tdry", "C", 1]
            ];

            alldata = []
            colnames = []
            for name, units, mult in self.hourly_output_columns:
                newcol = np.array(ssc.data_get_array(data, name.encode()))*mult
                alldata.append( newcol )
                colnames.append(name)
            alldata = np.array(alldata).T
            self.hourly_data = df = pd.DataFrame(alldata, columns=colnames, index=pd.date_range('1/1/2019', periods=8760, freq='h'))

            # print(df[df.index.month == 3].describe())

            dfsun = df[df.q_dot_rec_inc > 0]
            dfnight = df[df.q_pb > 0][df.q_dot_rec_inc == 0]

        self.summary_results = []

        printouts = [
            ['Annual energy', 'annual_energy'],
            ['Capacity factor', 'capacity_factor'],
            ['LCOE (real)', 'lcoe_real'],
            ['Site improvement', 'csp.pt.cost.site_improvements'],
            ['Heliostats', 'csp.pt.cost.heliostats'],
            ['Tower', 'csp.pt.cost.tower'],
            ['Receiver', 'csp.pt.cost.receiver'],
            ['Storage', 'csp.pt.cost.storage'],
            ['Power block', 'csp.pt.cost.power_block'],
            ['Contingency', 'csp.pt.cost.contingency'],
            ['Direct costs subtotal', 'total_direct_cost'],
            # ['EPC', 'csp.pt.cost.epc.total'],
            # ['Total land cost', 'csp.pt.cost.plm.total'],
            # ['Sales tax', 'csp.pt.cost.sales_tax.total'],
            ['Indirect costs subtotal', 'total_indirect_cost'],
            ['Net capital cost', 'cost_installed'],
            ['Cost per capacity', 'csp.pt.cost.installed_per_capacity'],
            ['Tower height', 'h_tower'],
        ]

        for lab,var in printouts:
            val = ssc.data_get_number(data, var.encode())
            self.summary_results.append([lab, val])

        if self.settings.save_hourly_results:
            #dni-weighted annual field efficiency
            dni_sum = df.beam.sum()
            self.summary_results.append(['Field efficiency - annual', (df.beam * df.eta_field_tot).sum() / dni_sum*100])
            #... and by subfield
            self.summary_results.append(['Field 1 efficiency - annual', (df.beam * df.eta_field1).sum() / dni_sum*100])
            self.summary_results.append(['Field 2 efficiency - annual', (df.beam * df.eta_field2).sum() / dni_sum*100])
            self.summary_results.append(['Field 3 efficiency - annual', (df.beam * df.eta_field3).sum() / dni_sum*100])
            #annual receiver efficiency
            self.summary_results.append(['Receiver efficiency - annual', (df.eta_therm * df.q_dot_rec_inc).sum() / df.q_dot_rec_inc.sum()*100])
            #... and by subfield
            self.summary_results.append(['Receiver 1 efficiency - annual', (df.eta_rec_therm1 * df.q_dot_rec_inc1).sum() / df.q_dot_rec_inc1.sum()*100])
            self.summary_results.append(['Receiver 2 efficiency - annual', (df.eta_rec_therm2 * df.q_dot_rec_inc2).sum() / df.q_dot_rec_inc2.sum()*100])
            self.summary_results.append(['Receiver 3 efficiency - annual', (df.eta_rec_therm3 * df.q_dot_rec_inc3).sum() / df.q_dot_rec_inc3.sum()*100])
            #annual cycle efficiency
            self.summary_results.append(['Cycle efficiency - annual', (df.eta * df.q_pb).sum() / df.q_pb.sum()*100])
            #... on sun
            self.summary_results.append(['Cycle on-sun efficiency - annual', (dfsun.eta * dfsun.q_pb).sum() / dfsun.q_pb.sum()*100])
            #... off sun
            self.summary_results.append(['Cycle off-sun efficiency annual', (dfnight.eta * dfnight.q_pb).sum() / dfnight.q_pb.sum()*100])
            #

        #if printing summary results...
        if self.settings.print_summary_output:
            for lab,val in self.summary_results:
                if val > 1000:
                    vals = "{:,d}".format(int(val))
                else:
                    vals = "{:.3f}".format(val)
                print("{:35s}\t{:>15s}".format(lab, vals))

        return

    #----------------------------------------------------------------
    def get_result_value(self, name):
        for row in self.summary_results:
            if row[0] == name:
                return row[1]
        print("Variable not found: " + name)
    
    #----------------------------------------------------------------
    def write_hourly_results_to_file(self, file_path="output_dview.csv"):

        if not self.settings.save_hourly_results:
            print("Hourly results were not saved for this simulation and connot be written to file.")
            return

        allout = []
        for name,units,_ in self.hourly_output_columns:
            allout.append( np.concatenate( ([name, units], self.hourly_data[name].values) ) )

        allout = np.array(allout).T.astype(np.str)

        with open(file_path, 'w') as fout:
            txt = "\n".join([",".join(line) for line in allout])
            fout.write(txt)

        return

#------------------------------------------------------------------------------

if __name__ == "__main__":

    # 100% HX Cost
    cases = [
        ['base', 'surround', 'skip', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
        ['optimal', 'surround', 'skip', 78.712, 2.724, 766.321, 5.082, 0.633, 0.597, 17.555, 42.446, 40.663],
        ['base', 'surround', 'bucket', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
        ['optimal', 'surround', 'bucket', 24.713, 2.535, 721.602, 5.759, 0.262, 0.23, 15.185, 47.676, 31.107],
        ['base', 'north', 'skip', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
        ['optimal', 'north', 'skip', 48.626, 2.47, 701.404, 5.768, 0.339, 0.264, 24.214, 42.807, 30.383],
        ['base', 'north', 'bucket', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
        ['optimal', 'north', 'bucket', 44.906, 2.435, 784.226, 5.931, 0.496, 0.281, 13.389, 41.356, 40.221],
    ]

    # 75% HX Cost
    # cases = [
    #     ['base', 'surround', 'skip', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
    #     ['optimal', 'surround', 'skip', 98.5, 2.18, 704, 7.94, 0.392, 0.674, 14.4, 26.7, 36.2],
    #     ['base', 'surround', 'bucket', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
    #     ['optimal', 'surround', 'bucket', 52.2, 2.46, 730, 8.47, 0.394, 0.367, 15.9, 33.8, 30.5],
    #     ['base', 'north', 'skip', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
    #     ['optimal', 'north', 'skip', 68.6, 2.38, 777, 7.32, 0.302, 0.388, 14.2, 43.2, 38.9],
    #     ['base', 'north', 'bucket', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
    #     ['optimal', 'north', 'bucket', 30.7, 2.59, 737, 5.31, 0.27, 0.294, 18.1, 47.1, 25.1],
    # ]

    # 50% HX Cost
    # cases = [
        # ['base', 'surround', 'skip', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
        # ['optimal', 'surround', 'skip', 89.6, 2.51, 769, 6.08, 0.363, 0.446, 16.5, 41.4, 22.4],
        # ['base', 'surround', 'bucket', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
        # ['optimal', 'surround', 'bucket', 29.145,     2.450,   740.121,     6.284,     0.234,     0.192,    18.375,    33.216,    39.647],
        # ['base', 'north', 'skip', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
        # ['optimal', 'north', 'skip', 82.6, 2.38, 746, 6.42, 0.354, 0.382, 12.9, 32.4, 34.0],
        # ['base', 'north', 'bucket', 100, 3, 976, 5.3, 0.45, 0.45, 13.3, 15, 15],
        # ['optimal', 'north', 'bucket', 50.5, 2.43, 756, 6.18, 0.315, 0.395, 16.9, 42.8, 29.9],
    # ]

    all_sum_results = {}
    casenames = []

    for case in cases:
    
        g = Gen3opt()

        g.settings.print_summary_output = True
        g.settings.save_hourly_results = True
        # g.settings.print_ssc_messages = True

        g.settings.scale_hx_cost = 0.5
        
        evaltype, \
        northstr, \
        g.settings.lift_technology, \
        g.variables.cycle_design_power, \
        g.variables.solar_multiple, \
        g.variables.dni_design_point, \
        g.variables.receiver_height, \
        g.variables.riser_inner_diam, \
        g.variables.downcomer_inner_diam, \
        g.variables.hours_tes, \
        g.variables.dT_approach_charge_hx, \
        g.variables.dT_approach_disch_hx = case

        g.settings.is_north = 'north' in northstr

        g.exec()

        if case == cases[0]:
            keyord = []
            for res,v in g.summary_results:
                keyord.append(res)
                all_sum_results[res] = [v]
        else:
            for res,v in g.summary_results:
                all_sum_results[res].append(v)

        casename = northstr + '-' + g.settings.lift_technology + '-' + evaltype
        casenames.append(casename)

        g.write_hourly_results_to_file( casename + '.csv')

        

    fsum = open('optimal-summary-results.csv', 'w')
    fsum.write("," + ",".join(casenames) + '\n')
    
    for key in keyord:
        fsum.write( ','.join([key] + [str(v) for v in all_sum_results[key]]) + '\n')

    fsum.close()


