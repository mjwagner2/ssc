from math import exp, log
from numpy import array
from scipy.interpolate import interp2d
from receiver import calculate_tower_height, specheat_co2

__particle_specheat = 1.3 #[kJ/kg-K]
__particle_density = 1600 #[kg/m3]
__media_cost = 50. #$/ton

__raw_data_dphx = {
    "m_dot" : [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.07, 0.1, 0.13, 0.15],
    "T_avg" : [500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800],
    #dp (fractional) vs temperature (vertical axis) and mass flow (horizontal axis)
    "dp" : [ \
                [3.40916e-05, 6.81832e-05, 1.02275e-04, 1.36366e-04, 1.70458e-04, 2.75829e-04, 3.53229e-04, 4.40094e-04, 5.36399e-04, 1.15833e-03, 2.20751e-03, 3.57548e-03, 4.65771e-03],
                [3.56127e-05, 7.12254e-05, 1.06838e-04, 1.42451e-04, 1.78063e-04, 2.13676e-04, 3.65512e-04, 4.55140e-04, 5.54484e-04, 1.19577e-03, 2.27730e-03, 3.68752e-03, 4.80337e-03],
                [3.71558e-05, 7.43116e-05, 1.11467e-04, 1.48623e-04, 1.85779e-04, 2.22935e-04, 3.77834e-04, 4.70219e-04, 5.72594e-04, 1.23315e-03, 2.34689e-03, 3.79916e-03, 4.94846e-03],
                [3.87206e-05, 7.74411e-05, 1.16162e-04, 1.54882e-04, 1.93603e-04, 2.32323e-04, 3.90200e-04, 4.85336e-04, 5.90735e-04, 1.27051e-03, 2.41630e-03, 3.91044e-03, 5.09303e-03],
                [4.03065e-05, 8.06130e-05, 1.20919e-04, 1.61226e-04, 2.01532e-04, 2.41839e-04, 4.02610e-04, 5.00495e-04, 6.08912e-04, 1.30784e-03, 2.48558e-03, 4.02140e-03, 5.23716e-03],
                [4.19132e-05, 8.38264e-05, 1.25740e-04, 1.67653e-04, 2.09566e-04, 2.51479e-04, 4.15069e-04, 5.15698e-04, 6.27130e-04, 1.34517e-03, 2.55475e-03, 4.13210e-03, 5.38088e-03],
                [4.35403e-05, 8.70806e-05, 1.30621e-04, 1.74161e-04, 2.17701e-04, 2.61242e-04, 4.27577e-04, 5.30949e-04, 6.45390e-04, 1.38249e-03, 2.62382e-03, 4.24255e-03, 5.52426e-03],
                [4.51874e-05, 9.03749e-05, 1.35562e-04, 1.80750e-04, 2.25937e-04, 2.71125e-04, 4.40137e-04, 5.46250e-04, 6.63699e-04, 1.41984e-03, 2.69281e-03, 4.35281e-03, 5.66732e-03],
                [4.68543e-05, 9.37087e-05, 1.40563e-04, 1.87417e-04, 2.34272e-04, 2.81126e-04, 4.52749e-04, 5.61603e-04, 6.82056e-04, 1.45720e-03, 2.76175e-03, 4.46289e-03, 5.81012e-03],
                [4.85407e-05, 9.70813e-05, 1.45622e-04, 1.94163e-04, 2.42703e-04, 2.91244e-04, 4.65416e-04, 5.77010e-04, 7.00467e-04, 1.49458e-03, 2.83065e-03, 4.57282e-03, 5.95268e-03],
                [5.02461e-05, 1.00492e-04, 1.50738e-04, 2.00985e-04, 2.51231e-04, 3.01477e-04, 4.78138e-04, 5.92473e-04, 7.18932e-04, 1.53201e-03, 2.89951e-03, 4.68263e-03, 6.09503e-03],
                [5.19705e-05, 1.03941e-04, 1.55912e-04, 2.07882e-04, 2.59853e-04, 3.11823e-04, 4.90917e-04, 6.07992e-04, 7.37455e-04, 1.56947e-03, 2.96837e-03, 4.79234e-03, 6.23720e-03],
                [5.37135e-05, 1.07427e-04, 1.61141e-04, 2.14854e-04, 2.68568e-04, 3.22281e-04, 3.75995e-04, 6.23571e-04, 7.56036e-04, 1.60697e-03, 3.03722e-03, 4.90196e-03, 6.37922e-03],
                [5.54750e-05, 1.10950e-04, 1.66425e-04, 2.21900e-04, 2.77375e-04, 3.32850e-04, 3.88325e-04, 6.39209e-04, 7.74677e-04, 1.64453e-03, 3.10607e-03, 5.01151e-03, 6.52111e-03],
                [5.72546e-05, 1.14509e-04, 1.71764e-04, 2.29018e-04, 2.86273e-04, 3.43528e-04, 4.00782e-04, 6.54908e-04, 7.93380e-04, 1.68213e-03, 3.17494e-03, 5.12102e-03, 6.66289e-03],
                [5.90523e-05, 1.18105e-04, 1.77157e-04, 2.36209e-04, 2.95261e-04, 3.54314e-04, 4.13366e-04, 6.70668e-04, 8.12146e-04, 1.71980e-03, 3.24384e-03, 5.23048e-03, 6.80458e-03],
            ]
}

__dp_interp_f = interp2d(array(__raw_data_dphx["m_dot"]), array(__raw_data_dphx["T_avg"]), array(__raw_data_dphx["dp"]))

def cp_particle():
    """
    kJ/kg-K
    """
    return __particle_specheat

#----------------------------------------------------------------------
def calculate_hx_cost(q_cycle_in_kw, dT_approach_chg, dT_approach_dis, T_rec_out_C, T_rec_in_C, scale_cost=1.):
    """
    Inputs: 
    q_cycle_in_kw - cycle design thermal rating / discharge rating at design (kWt)
    dT_approach_chg - nominal charge heat exchanger approach !! Either 20C or 15C
    dt_approach_dis - nominal discharge heat exchangers approach (C)
    T_rec_out_C - nominal receiver hot outlet temperature (C)
    T_rec_in_C - nominal receiver cold inlet temperature (C)

    Returns:
    total_cost ($)
    cost_charge ($)
    cost_hot_discharge ($)
    cost_cold_discharge ($)
    eta_charge
    eta_hot_discharge
    eta_cold_discharge
    UA_charge
    UA_hot_discharge
    UA_cold_discharge
    """
    
    #calculate state point temperatures
    dT_cycle = dT_receiver = T_rec_out_C - T_rec_in_C
    T_charge_particle_out = T_rec_out_C - dT_approach_chg       #hot particles to storage
    T_charge_particle_in = T_charge_particle_out - dT_cycle     #cold particles to the receiver charge HX. assume CR=1
    T_hot_disch_co2_in = T_rec_in_C
    T_hot_disch_co2_out = T_cycle_in = T_rec_out_C - dT_approach_chg - dT_approach_dis
    T_cycle_out = T_cycle_in - dT_cycle
    T_cold_disch_co2_in = T_cycle_out
    T_cold_disch_co2_out = T_rec_in_C
    T_hot_disch_particle_in = T_charge_particle_out
    T_hot_disch_particle_out = T_hot_disch_co2_in + dT_approach_dis
    
    #heat exchanger mass flows
    cp_co2_cycle = specheat_co2(T_cycle_in + dT_cycle/2)    #kJ/kg-K
    m_dot_co2 = q_cycle_in_kw / (cp_co2_cycle * dT_cycle)
    m_dot_particle = q_cycle_in_kw / (__particle_specheat * dT_cycle)

    #duty (kW)
    q_charge_duty = q_cycle_in_kw

    cp_co2_hot_disch = specheat_co2( (T_hot_disch_co2_out + T_hot_disch_co2_in)/2 )
    q_hot_disch_duty = m_dot_co2*cp_co2_hot_disch*(T_hot_disch_co2_out - T_hot_disch_co2_in)
    
    cp_co2_cold_disch = specheat_co2( (T_cold_disch_co2_out + T_cold_disch_co2_in)/2 )
    q_cold_disch_duty = m_dot_co2*cp_co2_cold_disch*(T_cold_disch_co2_out - T_cold_disch_co2_in)


    #effectiveness
    eta_charge = q_charge_duty / (m_dot_co2 * cp_co2_cycle * (T_rec_out_C - T_charge_particle_in))
    
    cr_min_hot_disch = min([m_dot_co2*cp_co2_hot_disch, m_dot_particle*__particle_specheat])
    eta_hot_disch = q_hot_disch_duty / (cr_min_hot_disch*(T_hot_disch_particle_in - T_rec_in_C))
    
    cr_min_cold_disch = min([m_dot_co2*cp_co2_cold_disch, m_dot_particle*__particle_specheat])
    eta_cold_disch = q_cold_disch_duty / (cr_min_cold_disch*(T_hot_disch_particle_out - T_cycle_out))

    #UA
    UA_charge = m_dot_co2 * cp_co2_cycle * (eta_charge / (1.-eta_charge))  
    
    cr_hot_disch = cr_min_hot_disch / max([m_dot_co2*cp_co2_cold_disch, m_dot_particle*__particle_specheat])
    UA_hot_disch = cr_min_hot_disch * log( (1 - eta_hot_disch*cr_hot_disch) / (1.-eta_hot_disch)) / (1 - cr_hot_disch)
    
    cr_cold_disch = cr_min_cold_disch / max([m_dot_co2*cp_co2_cold_disch, m_dot_particle*__particle_specheat])
    UA_cold_disch = cr_min_cold_disch * log( (1 - eta_cold_disch*cr_cold_disch) / (1.-eta_cold_disch)) / (1 - cr_cold_disch)

    cost_charge = 3400 * scale_cost * UA_charge * 3         # $3.40/UA
    cost_hot_disch = 3400 * scale_cost * UA_hot_disch
    cost_cold_disch = 9800 * scale_cost * UA_cold_disch  # $9.80/UA --- costs from Brayton study
    
    ### Specify parameters and factors

    #Performance confidence derate (1 -> full confidence)
    x_safety = 1.0
    #assumed future cost reduction multiplier
    x_future = 0.7      
    #reduced cost factor at low temperature
    x_lt_cost = 0.75
    #max particle flow rate per cell
    m_dot_particle_c_max = 0.056    #kg/s       --> rho_b * sqrt(g) * (w_cell / tan(theta^prime)) * (w_outlet - 1.5 d_p)^(3/2)
    #Max cell length
    L_cell_max = 1.65       #m
    # cell width
    W_cell = 0.2032         #m
    
    #Calculate cell surface area based on assumed conductance
    U = 0.450               #kW/m2-K
    A_charge = UA_charge / U
    A_hot_disch = UA_hot_disch / U
    A_cold_disch = UA_cold_disch / U

    #Determine the minimum number of cells
    N_cells_min = m_dot_particle / (x_safety * m_dot_particle_c_max)

    #Calculate the length of each cell, limiting to the maximum cell length
    L_cell_charge =     min([A_charge     / (2 * W_cell * N_cells_min), L_cell_max])
    L_cell_hot_disch =  min([A_hot_disch  / (2 * W_cell * N_cells_min), L_cell_max])
    L_cell_cold_disch = min([A_cold_disch / (2 * W_cell * N_cells_min), L_cell_max])

    #Ensure cell count matches if the max length constraint is active
    N_cells_charge     = A_charge     / (2 * W_cell * L_cell_charge)
    N_cells_hot_disch  = A_hot_disch  / (2 * W_cell * L_cell_hot_disch)
    N_cells_cold_disch = A_cold_disch / (2 * W_cell * L_cell_cold_disch)

    #Fraction of total costs associated with low cost materials
    x_m_charge     = (0.1466 * L_cell_charge     + 0.2681) * x_lt_cost \
                        if max([ T_rec_out_C, T_charge_particle_out]) < 600 else 0.
    x_m_hot_disch  = (0.1466 * L_cell_hot_disch  + 0.2681) * x_lt_cost \
                        if max([ T_hot_disch_co2_out, T_hot_disch_particle_in]) < 600 else 0.
    x_m_cold_disch = (0.1466 * L_cell_cold_disch + 0.2681) * x_lt_cost \
                        if max([ T_cold_disch_co2_out, T_hot_disch_particle_out]) < 600 else 0.

    #Cost of the cells ($/cell)
    c_cell_charge     = 271.1 * L_cell_charge     + 603
    c_cell_hot_disch  = 271.1 * L_cell_hot_disch  + 603
    c_cell_cold_disch = 271.1 * L_cell_cold_disch + 603

    #Total heat exchanger cost
    cost_charge     = c_cell_charge     * N_cells_charge     * x_future * (1. - x_m_charge) * 3
    cost_hot_disch  = c_cell_hot_disch  * N_cells_hot_disch  * x_future * (1. - x_m_hot_disch)
    cost_cold_disch = c_cell_cold_disch * N_cells_cold_disch * x_future * (1. - x_m_cold_disch)

    #Fractional pressure loss
    dp_charge     = __dp_interp_f(m_dot_co2 / N_cells_charge,     (T_rec_out_C          + T_hot_disch_co2_in )/2.)[0]*L_cell_charge
    dp_hot_disch  = __dp_interp_f(m_dot_co2 / N_cells_hot_disch,  (T_hot_disch_co2_out  + T_hot_disch_co2_in )/2.)[0]*L_cell_hot_disch
    dp_cold_disch = __dp_interp_f(m_dot_co2 / N_cells_cold_disch, (T_cold_disch_co2_out + T_cold_disch_co2_in)/2.)[0]*L_cell_cold_disch

    return { 
        'total_cost':(cost_charge + cost_hot_disch + cost_cold_disch), 
        'cost_charge':cost_charge,
        'cost_hot_discharge':cost_hot_disch,
        'cost_cold_discharge':cost_cold_disch,
        'eta_charge':eta_charge,
        'eta_hot_discharge':eta_hot_disch,
        'eta_cold_discharge':eta_cold_disch,
        'UA_charge':UA_charge,
        'UA_hot_discharge':UA_hot_disch,
        'UA_cold_discharge':UA_cold_disch,
        'duty_charge':q_charge_duty,
        'duty_discharge_hot':q_hot_disch_duty,
        'duty_discharge_cold':q_cold_disch_duty,
        'eta_charge':eta_charge,
        'eta_hot_disch':eta_hot_disch,
        'eta_cold_disch':eta_cold_disch,
        'L_cell_charge':L_cell_charge,
        'L_cell_hot_disch':L_cell_hot_disch,
        'L_cell_cold_disch':L_cell_cold_disch,
        'N_cells_charge':N_cells_charge,
        'N_cells_hot_disch':N_cells_hot_disch,
        'N_cells_cold_disch':N_cells_cold_disch,
        'dp_charge':dp_charge,
        'dp_hot_disch':dp_hot_disch,
        'dp_cold_disch':dp_cold_disch,
    }

#----------------------------------------------------------------------

def calculate_hx_base_dp(m_dot_cell, T_avg_C, L):
    """
    Inputs:
        m_dot_cell - (kg/s) mass flow rate of CO2 through a single cell
        T_avg_C - (C) average temperature of the CO2 through the HX
        L - (m) cell length

    Returns
        fractional pressure drop
    """

    return __dp_interp_f(m_dot_cell, T_avg_C)[0]*L

#----------------------------------------------------------------------

def calculate_balance_tes_cost(q_cycle_in_kw):
    """
    Inputs:
        q_cycle_in_kw - cycle design thermal rating / discharge rating at design (kWt)

    Returns:
        Specific cost ($/kWht) of thermal storage balance of system costs
    """

    #the correlation 
    cycle_power_mw = 0.001 * q_cycle_in_kw * 0.43       #convert to equivalent cycle power using original assumed efficiency

    return -2.27525e-5 * cycle_power_mw**3 + 6.42615e-3 * cycle_power_mw**2 - 6.76546e-1 * cycle_power_mw + 4.61456e1


#----------------------------------------------------------------------

def calculate_silo_cost(q_cycle_in_kw, hours_tes, dt_cycle):
    """
    Inputs:
        q_cycle_in_kw - cycle design thermal rating / discharge rating at design (kWt)
        hours_tes - hours of full load thermal storage (hr)
        dt_cycle - nominal temperature drop across the cycle (C)

    Returns:
        silo_cost - Cost of a pair of silos, each able to handle full particle capacity ($)
        media_cost - Bulk media purchase cost ($)
    """

    #assume material properties from UW
    cp = __particle_specheat #[kJ/kg-K]
    rho = __particle_density #[kg/m3]

    #total energy stored (kJ)
    E = q_cycle_in_kw * hours_tes * 3600

    #mass and volume of particles
    m = E / (cp * dt_cycle)  #kg
    V = m/rho
    m_ton = m * 0.001102311  

    #cost curve from Megan/Jack's old spreadsheet cost-model-jh_ROM BOM Baseload100MW.xlsx
    cost = 1.551e6 * exp(4.2662e-5 * m_ton) * 2

    #media cost
    media = __media_cost * (m_ton + m_ton*.02*30 + m_ton*0.1)   #assume 2% annual replacement, 10% extra inventory

    return {'silo_cost':cost, 'media_cost':media}
#----------------------------------------------------------------------

def calculate_lift_efficiency(q_solarfield_in_kw, q_solarfield_out_kw, m_dot_p, lift_type):
    """
    Calculate lift cost as a function of solar field rating and lift type. 

    Inputs
        q_solarfield_out_kw - total absorbed power at design from all receivers (kWt)
        q_solarfield_in_kw - total incident power at design from all receivers (kWt)
        m_dot_p - lift particle mass flow rate (kg/s)
        lift_type - one of 'bucket' or 'skip'
    Returns
        Lift cost ($)
    """

    #look up power from WP curve
    x = q_solarfield_out_kw/1000
    if lift_type == 'bucket':
        lift_power = 1.487E-02*x**2 + 1.025E+01*x - 2.901E+02
    elif lift_type == 'skip':
        lift_power = 6.931E-03*x**2 + 1.373E+01*x - 6.519E+02
    else:
        raise Exception("Invalid lift_type. Must be one of 'bucket' or 'skip'")
    
    #assumed tower height is WP model
    tht = calculate_tower_height(q_solarfield_in_kw, wp_data=True)
    bulk_power = m_dot_p * 9.81 * tht*1.1 / 1e3 #kW
    
    return bulk_power / lift_power

def calculate_lift_cost(q_solarfield_out_kw, lift_type):
    """
    Calculate lift cost as a function of solar field rating and lift type. 

    Inputs
        q_solarfield_out_kw - total absorbed power at design from all receivers (kWt)
        lift_type - one of 'bucket' or 'skip'
    Returns
        Lift cost ($)
    """
    x = q_solarfield_out_kw /1000. * 0.43 / 3       #Convert to nominal cycle power assuming SM=3 and 43% cycle eff.

    if lift_type == 'bucket':
        return 8.65829e3 * x**2 + 5.41087e5 * x - 5.38127e5
    elif lift_type == 'skip':
        return 1.84708e3 * x**2 + 3.91687e5 * x + 1.27537e6
    else:
        raise Exception("Invalid lift_type. Must be one of 'bucket' or 'skip'")

#----------------------------------------------------------------------

def calculate_lift_availability(q_cycle_in_kw, lift_type):
    """
    Inputs
        q_cycle_in_kw - cycle design thermal rating / discharge rating at design (kWt)
        lift_type - one of 'bucket' or 'skip'

    Returns
        Availability (-) of the lift system
    """

    x = q_cycle_in_kw /1000. * 0.43      #Convert to nominal cycle power assuming 43% cycle eff.

    if lift_type == 'bucket':
        return 2.21083E-10 * x**5 - 5.90973E-08 * x**4 + 5.33490E-06 * x**3 - 1.79084E-04 * x**2 + 3.24384E-04 * x + 9.79899E-01
    elif lift_type == 'skip':
        return -3.18735E-10 * x**4 + 2.88137E-08 * x**3 - 7.44566E-07 * x**2 + 5.07427E-06 * x + 9.79998E-01
    else:
        raise Exception("Invalid lift_type. Must be one of 'bucket' or 'skip'")


#----------------------------------------------------------------------
if __name__ == "__main__":

    qc=100000/.43
    
    # print(calculate_hx_cost(qc*3, qc, 15))

    # print(calculate_silo_cost(qc, 13, 715-560))

    # print( calculate_hx_cost(qc*3, 20, 15, 730, 592) )

    print( calculate_hx_base_dp(0.13, 660, 1.5) )

    