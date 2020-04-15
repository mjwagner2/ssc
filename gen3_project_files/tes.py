from math import exp

#----------------------------------------------------------------------
def calculate_hx_cost(q_solarfield_out_kw, q_cycle_in_kw, dT_approach_chg):
    """
    Inputs: 
    q_solarfield_out_kw - total absorbed power at design from all receivers (kWt)
    q_cycle_in_kw - cycle design thermal rating / discharge rating at design (kWt)
    dT_approach_chg - nominal charge heat exchanger approach !! Either 20C or 15C

    Returns:
    total_cost ($)
    cost_charge ($)
    cost_discharge ($)
    """
    
    #the fraction of discharge HX duty on the high-temp unit
    f_ht = (176405)/(176405 + 44956)

    if dT_approach_chg == 20:
        cost_charge = 182. * q_solarfield_out_kw
        cost_discharge = 368.*f_ht*q_cycle_in_kw + 988*(1-f_ht)*q_cycle_in_kw

    elif dT_approach_chg == 15:
        cost_charge = 245. * q_solarfield_out_kw
        cost_discharge = 232.*f_ht*q_cycle_in_kw + 637*(1-f_ht)*q_cycle_in_kw

    else:
        raise Exception("Invalid approach temperature passed to heat exchanger cost calculation. Requires 15 or 20C, got " + str(dT_approach_chg) + " C.")

    return {'total_cost':(cost_charge + cost_discharge), 'cost_charge':cost_charge, 'cost_discharge':cost_discharge}

#----------------------------------------------------------------------

def calculate_silo_cost(q_cycle_in_kw, hours_tes, dt_cycle):
    """
    Inputs:
    q_cycle_in_kw - cycle design thermal rating / discharge rating at design (kWt)
    hours_tes - hours of full load thermal storage (hr)
    dt_cycle - nominal temperature drop across the cycle (C)

    Returns:
    Cost of a pair of silos, each able to handle full particle capacity ($)
    """

    #assume material properties from UW
    cp = 1.3 #[kJ/kg-K]
    rho = 1600 #[kg/m3]

    #total energy stored (kJ)
    E = q_cycle_in_kw * hours_tes * 3600

    #mass and volume of particles
    m = E / (cp * dt_cycle)  #kg
    V = m/rho
    m_ton = m * 0.001102311  

    #cost curve from Megan/Jack's old spreadsheet cost-model-jh_ROM BOM Baseload100MW.xlsx
    cost = 1.551e6 * exp(4.2662e-5 * m_ton) * 2

    return {'cost':cost}


#----------------------------------------------------------------------


if __name__ == "__main__":

    qc=100000/.43
    
    print(calculate_hx_cost(qc*3, qc, 15))

    print(calculate_silo_cost(qc, 13, 715-560))

