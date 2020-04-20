from math import exp
from receiver import calculate_tower_height, specheat_co2

__particle_specheat = 1.3 #[kJ/kg-K]
__particle_density = 1600 #[kg/m3]
__media_cost = 50. #$/ton

#----------------------------------------------------------------------
def calculate_hx_cost(q_cycle_in_kw, dT_approach_chg, dT_approach_dis, T_rec_out_C, T_rec_in_C):
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
    UA_charge = m_dot_co2 * cp_co2_cycle * (eta_charge / (1.-eta_charge)) * 3 
    UA_hot_disch = cr_min_hot_disch * (eta_hot_disch / (1.-eta_hot_disch))
    UA_cold_disch = cr_min_cold_disch * (eta_cold_disch / (1.-eta_cold_disch))

    cost_charge = 3400 * UA_charge          # $3.40/UA
    cost_hot_disch = 3400 * UA_hot_disch
    cost_cold_disch = 9800 * UA_cold_disch  # $9.80/UA --- costs from Brayton study

    #the fraction of discharge HX duty on the high-temp unit
    # f_ht = (176405)/(176405 + 44956)

    # if dT_approach_chg == 20:
    #     cost_charge = 182. * q_solarfield_out_kw
    #     cost_discharge = 368.*f_ht*q_cycle_in_kw + 988*(1-f_ht)*q_cycle_in_kw

    # elif dT_approach_chg == 15:
    #     cost_charge = 245. * q_solarfield_out_kw
    #     cost_discharge = 232.*f_ht*q_cycle_in_kw + 637*(1-f_ht)*q_cycle_in_kw
    # else:
    #     raise Exception("Invalid approach temperature passed to heat exchanger cost calculation. Requires 15 or 20C, got " + str(dT_approach_chg) + " C.")

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

def calculate_lift_efficiency(q_solarfield_in_kw, q_solarfield_out_kw, lift_type):
    """
    Calculate lift cost as a function of solar field rating and lift type. 

    Inputs
        q_solarfield_out_kw - total absorbed power at design from all receivers (kWt)
        q_solarfield_in_kw - total incident power at design from all receivers (kWt)
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
    bulk_power = __particle_density * 9.81 * tht*1.1 / 1000 #kW
    
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

    print( calculate_hx_cost(qc*3, qc, 15, 15, 730, 592) )
    