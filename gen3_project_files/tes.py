from math import exp, log
from receiver import calculate_tower_height, specheat_co2

__particle_specheat = 1.3 #[kJ/kg-K]
__particle_density = 1600 #[kg/m3]
__media_cost = 50. #$/ton

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
    }

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
    x = q_solarfield_out_kw /1000.

    if lift_type == 'bucket':
        return -7.640E-02*x**3 + 1.624E+02*x**2 + 1.806E+04*x + 2.923E+05
    elif lift_type == 'skip':
        return 5.536E-02*x**3 - 1.516E+01*x**2 + 6.626E+04*x + 1.354E+06
    else:
        raise Exception("Invalid lift_type. Must be one of 'bucket' or 'skip'")


#----------------------------------------------------------------------
if __name__ == "__main__":

    qc=100000/.43
    
    # print(calculate_hx_cost(qc*3, qc, 15))

    # print(calculate_silo_cost(qc, 13, 715-560))

    print( calculate_hx_cost(qc*3, 20, 15, 730, 592) )
    