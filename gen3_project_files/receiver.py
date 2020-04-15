from numpy import interp, pi, array, argsort
from math import ceil


#----------------------------------------------------------------------------
#lookup for receiver performance versus height (length) in m
__raw_data = {
    'length':[3, 3.556, 4.111, 4.2, 4.5, 4.667, 5.222, 5.3, 5.778, 6.3, 6.333, 6.889, 7.444, 8],
    'efficiency':[0.6047, 0.6821, 0.7351, 0.742, 0.7627, 0.7728, 0.8002, 0.8035, 0.8195, 0.824, 0.8243, 0.8281, 0.8312, 0.8338],  # -
    'pressure_drop':[0.01905, 0.05688, 0.141, 0.1609, 0.2459, 0.3071, 0.6068, 0.6636, 1.086, 1.478, 1.507, 2.027, 2.666, 3.435],  # MPa
    'm_dot_tube':[0.01627, 0.02775, 0.04279, 0.04555, 0.05554, 0.06157, 0.08409, 0.08757, 0.109, 0.1229, 0.1238, 0.1386, 0.1535, 0.1682], #kg/s
}

def specheat_co2(T_C):
    """
    Inputs
        T_C - temperature in deg C
    Returns
        specific heat (kJ/kg-K)
    """
    return 0.412972 * T_C**0.165755

def density_co2(T_C, P_MPa):
    """
    Inputs
        T_C - temperature of the gas (C)
        P_MPa - pressure of the gas (MPa)

    Returns
        Density - kg/m3
    """
    return P_MPa*1e6 / ((T_C+273.15) * 188.92)

#----------------------------------------------------------------------------
def calculate_efficiency(L):
    return interp(L, __raw_data['length'], __raw_data['efficiency'])

#----------------------------------------------------------------------------
def calculate_pressure_drop(L):
    return interp(L, __raw_data['length'], __raw_data['pressure_drop'])

#----------------------------------------------------------------------------
def calculate_m_dot_tube(L):
    return interp(L, __raw_data['length'], __raw_data['m_dot_tube'])

#----------------------------------------------------------------------------
def calculate_n_tubes(q_solarfield_out_kw, T_rec_in_des_C, T_rec_out_des_C, L):
    """
    Calculate the receiver number of tubes.

    Inputs
        q_solarfield_out_kw - Design-point receiver power output (kW)
        T_rec_in_des_C - Design-point receiver CO2 inlet temperature (C)
        T_rec_out_des_C - Design-point receiver CO2 outlet temperature (C)
        L - Specified receiver tube height (m)

    Outputs
        n_tubes - number of parallel tubes in the receiver (-)
    """

    m_dot_rec_tot = q_solarfield_out_kw / (specheat_co2((T_rec_in_des_C + T_rec_out_des_C)/2) * (T_rec_out_des_C - T_rec_in_des_C))

    m_dot_tube = calculate_m_dot_tube(L)

    n_tubes = ceil(m_dot_rec_tot/m_dot_tube)

    return {'n_tubes':n_tubes, 'm_dot_rec_tot':m_dot_rec_tot, 'm_dot_tube':m_dot_tube}

#----------------------------------------------------------------------------
def calculate_cost(L, N_tubes):
    """
    Calculate the receiver cost as a function of height (L) and number of tubes (N_tubes)

    The defaults are:
    L = 5.3 [m]
    N_tubes = 14124 

    This model provided by Brayton
    """
    
    rho = 0.291 * 27679.905 #  [lbm/in3] * convert(lbm/in3,kg/m3)
    D_out = 0.375 * 0.0254 #  [in] * convert(in,m)
    th = 0.07 * 0.0254 #  [in] * convert(in,m)
    D_in = D_out - 2 *th
    
    L_extra = 1.3 #  [m]	"extra to account for expansion loop"
    
    a_cs = pi/4*(D_out**2-D_in**2)
    V_tube = a_cs * (L + L_extra) * N_tubes 
    m_tube = V_tube * rho 
    
    #"header"
    f_extra = 1.25
    L_header = N_tubes * D_out * f_extra
    
    
    D_out_h = 2.875 * 0.0254 #  [in] * convert(in,m)
    th_h_in = 0.343 * 0.0254 #  [in]* convert(in,m)
    th_h_out = 0.688 * 0.0254 #  [in]* convert(in,m)
    D_in_h_in = D_out_h - 2 * th_h_in
    D_in_h_out = D_out_h - 2 * th_h_out
    
    V_h_in = pi/4*(D_out_h**2-D_in_h_in**2) * L_header
    V_h_out = pi/4*(D_out_h**2-D_in_h_out**2) * L_header
    
    m_header = rho * (V_h_in + V_h_out)
    
    # "cap"
    N_caps = L_header * 2 #  [1/m]
    m_cap_in_spec = 2.4 * 0.4536 #  [lbm] * convert(lbm,kg)
    m_cap_out_spec = 0.77 * 0.4536 #  [lbm] * convert(lbm,kg)
    
    m_cap = m_cap_in_spec * N_caps + m_cap_out_spec * N_caps
    
    # "tube_stub"
    m_tube_stub_spec = 0.18 * 0.4536 #  [lbm] * convert(lbm,kg)
    
    N_stub = N_tubes * 2 
    
    m_tube_stub = N_stub * m_tube_stub_spec
    
    # "plate"
    m_plate_spec = 0.02 * 0.4536 #  [lbm] * convert(lbm,kg)
    N_plates = N_tubes * 24
    
    m_plate = N_plates * m_plate_spec
    
    # "sum" 
    m_total = m_tube + m_header + m_cap + m_tube_stub + m_plate
    spec_cost = 40 * 2.2046226 #  [$/lbm] * convert(1/lbm,1/kg)
    mat_cost = m_total * spec_cost
    
    # "Toll Processing"
    stub_cost_spec = 8
    stub_toll_cost = N_stub * stub_cost_spec
    
    spec_cost_header = 4000 #  [$]
    header_toll_cost = L_header * spec_cost_header * 2 #  [1/m]
    
    spec_cost_cap = 2000 #  [$]
    cap_toll_cost = N_caps * spec_cost_cap
    
    plate_toll_spec = 0.1 #  [$]
    plate_toll_cost = plate_toll_spec * N_plates
    
    tube_bend_toll_spec = 5 #  [$]
    tube_bend_toll_cost = N_tubes * tube_bend_toll_spec
    
    toll_cost_total = stub_toll_cost + header_toll_cost + cap_toll_cost + plate_toll_cost + tube_bend_toll_cost
    
    # "Welding"
    cap_weld_spec = 109.76 * 3.28084 #  [$/ft] * convert(1/ft,1/m)
    cap_weld_cost = pi * D_out_h * cap_weld_spec * N_caps
    
    stub_weld_spec = 33.47 * 3.28084 #  [$/ft] * convert(1/ft,1/m)
    stub_weld_cost = N_stub * pi * 0.75 * 0.0254 * stub_weld_spec #  [in] *  convert(in,m) 
    
    stay_plate_spec = 8.76 * 3.28084#  [$/ft] * convert(1/ft,1/m)
    stay_plate_weld_cost = stay_plate_spec * 0.316 * 0.0254 * N_plates #  [in]  * convert(in,m)
    
    weld_cost = cap_weld_cost + stub_weld_cost + stay_plate_weld_cost
    
    # "frame"
    #                    [$]      [m]    [in] 
    frame_spec_cost = 46039.68 / (5.3 * 0.375 * 0.0254 * 14124)
    frame_cost = frame_spec_cost * L * N_tubes * D_out
    
    # "insulation"
    #                [$]   [ft][ft] convert(ft2,m2)
    ins_spec_cost = 1000 / (3 * 4  * 0.09290304)
    A_rec = L * N_tubes * D_out
    
    # A_rec = W_rec * L_rec
    L_rec = (A_rec/2.)**(0.5)
    W_rec = 2 * L_rec	#"assume worst case apsect ratio of 2"

    per = 2*L_rec + 2*W_rec
    L_ins_sur = 1.5 #  [m]	"1.5 m for spillage"
    
    A_ins = per * L_ins_sur
    ins_cost = ins_spec_cost * A_ins
    
    # "total cost"
    total_cost = mat_cost + toll_cost_total + weld_cost + frame_cost +  ins_cost
    return {'total_cost':total_cost, 'A_rec':A_rec}


#----------------------------------------------------------------------------
def __poly_coefs(X,Y):
    """
    Solve for polynomial coefficients that give an exact fit on a 3-point array.

    Input
        X - dim 3 array of x-axis values
        Y - dim 3 array of y-axis values
    
    Returns
        array[3] : polynomial coefficients [const, linear, quad]
    """

    if len(X) != len(Y) or len(X) != 3:
        raise Exception("__poly_coefs(X,Y) requires X,Y arrays of length 3")

    x1,x2,x3 = X
    y1,y2,y3 = Y

    #solve explicitly for coefficients based on x,y data
    c2 = (y3-y1 - (x3-x1)/(x2-x1)*(y2-y1))/(x3**2 - x1**2 + (x3-x1)/(x2-x1)*(x1**2 - x2**2))
    c1 = (y2-y1 + c2*(x1**2 - x2**2))/(x2-x1)
    c0 = y1 - c2*x1**2 - c1*x1

    return [c0, c1, c2]

def __interp_poly(X,Y,x):
    c0,c1,c2 = __poly_coefs(X,Y)
    return (lambda v: c0 + c1*v + c2*v**2)(x)
#----------------------------------------------------------------------------
    

def create_heliostat_field_lookup(efficiency_file_path, q_solarfield_in_kw, heliostat_area_m2):
    """
    Read the efficiency file and interpolate between power levels to produce a single, power-appropriate
    efficiency lookup table for use in the performance simulation.\n\n

    Inputs\n
        efficiency_file_path\n
            Relative path to efficiency file. The first row of the file lists included power levels.
            The second and following rows list efficiency and active area of each subfield as a function
            of sun position. Columns are:\n
                                 | Power level 1 .........................................| Power level 2.......\n
            SolAz(deg)|SolEl(deg)|Sf1Eff(-)|Sf2Eff(-)|Sf3Eff(-)|Sf1A(m2)|Sf2A(m2)|Sf3A(m3)|Sf1Eff(-)|Sf2Eff(-)|...\n
        q_solarfield_in_kw\n
            Nominal power level for the solar field at design (kw)\n
        heliostat_area_m2\n
            Active reflective area of a *single* heliostat (m2)\n\n

    Returns:\n
        Lookup table expressing efficiency and number of heliostats active for each receiver. The table
        contains the following for an array of solar angles:\n
        SolAz(deg)|SolEl(deg)|Sf1Eff(-)|Sf2Eff(-)|Sf3Eff(-)|Sf1A(m2)|Sf2A(m2)|Sf3A(m3)
    """

    #load the file
    fdat = open(efficiency_file_path, 'r').readlines()
    
    #collect the power levels
    power_levels = fdat[0].replace('\n','').split(",")
    while ("" in power_levels):
        power_levels.remove("")
    power_levels = [float(v) for v in power_levels]

    #read raw data
    raw_data = [[float(v) for v in line.split(",")] for line in fdat[1:]]

    #there should be 2 + 6*len(power_levels) columns
    if len(raw_data[0]) != (2 + 6*len(power_levels)):
        raise Exception("The number of columns in the heliostat efficiency data file does not match "
                        "the number of power levels. Expected {:d} columns for {:d} power levels but "
                        "read in {:d}.".format(2+6*len(power_levels), len(power_levels), len(raw_data[0])))
    
    raw_data_T = array(raw_data).T

    #structure data
    all_data = {}
    
    interp_data = {
        "azimuth" : raw_data_T[0],
        "zenith"  : raw_data_T[1],
    }

    cols = ['azimuth','zenith','eta_1','eta_2','eta_3','nh_1','nh_2','nh_3']

    #find the 3-point region closest to the given power level
    power_level_indices = argsort( array([(q-q_solarfield_in_kw)**2 for q in power_levels]) )[:3]
    power_levels_used = [power_levels[i] for i in power_level_indices]

    #collate all available data by sun position, with a list of values corresponding to each power level
    for j,lab in enumerate(cols[2:]):
        zcols = []
        #convert heliostat total area to heliostat count, if applicable
        scale = heliostat_area_m2 if "nh" in str(lab) else 1.
        #get all of the related columns by power level
        for i in power_level_indices:
            zcols.append(raw_data_T[2+6*i+j]/scale)

        #zip columns into lists
        all_data[lab] = list(zip(*zcols))

        interp_data[lab] = []

        for datset in all_data[lab]:
            interp_data[lab].append( __interp_poly(power_levels_used, datset, q_solarfield_in_kw) )
        
    return array([interp_data[col] for col in cols]).T.tolist()


#----------------------------------------------------------------------------
if __name__ == "__main__":
    #test
    # print(calculate_cost(5.3, 14124))
    
    # import matplotlib.pyplot as plt
    # import numpy

    # X = [3, 12, 33]
    # Y = [15, 5, 40]    
    # c0,c1,c2 = __poly_coefs(X,Y)
    # Xa = numpy.arange(0, 50, .2)
    # Ya = [(lambda x: c0 + c1*x + c2*x**2)(x) for x in Xa]
    # plt.plot(X,Y, 'ks')
    # plt.plot(Xa,Ya, 'b-')
    # plt.show()

    # tab = create_heliostat_field_lookup("resource/eta_lookup_north.csv", 100, 8.6**2)

    # print(specheat_co2(550))
    # print(density_co2(550, 25))