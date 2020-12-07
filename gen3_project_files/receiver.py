from numpy import interp, pi, array, argsort, zeros
from math import ceil
from scipy.interpolate import SmoothBivariateSpline, interp1d, interp2d, Rbf, griddata
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pandas
from globalspline import GlobalSpline2D

#----------------------------------------------------------------------------
def specheat_co2(T_C):
    """
    Inputs
        T_C - temperature in deg C
    Returns
        specific heat (kJ/kg-K)
    """
    return 0.412972 * T_C**0.165755

def density_co2(T_C, P_kPa):
    """
    Inputs
        T_C - temperature of the gas (C)
        P_kPa - pressure of the gas (kPa)

    Returns
        Density - kg/m3
    """
    return P_kPa*1e3 / ((T_C+273.15) * 188.92)

#----------------------------------------------------------------------------
def calculate_tower_height(q_solarfield_in_kw, is_north = True, wp_data = False):
    """
    Inputs
        q_solarfield_in_kw - design point incident power on all receivers (kW)
        is_north - use north-field configuration
        wp_data - use Worley-Parsons curve
    
    Returns
        Tower height (m)
    """
    q = q_solarfield_in_kw/1000.

    if wp_data:
        tht = 3.3627e-07*q**3 - 6.4330E-04*q**2 + 5.3928E-01*q + 3.7892E+01
    else:
        if is_north:
            tht = 11.376 * q**0.4554
        else:
            tht = 13.231 * q**0.4031

    return tht

def ReceiverHeightRange(D_receiver_tube):
    """
    Returns the nominal minimum and maximum receiver height as function of the receiver tube size

    D_receiver_tube     [in] Receiver single tube outer (nominal) diameter, valid from 1/4 to 3/8
    """

    x = D_receiver_tube
    H_rec_min = -118.14*x**2 + 95.992*x - 14.883
    H_rec_max = -74.88*x**2 + 77.816*x - 12.351

    return [H_rec_min, H_rec_max]

def ReceiverTubeDiameterRange(H_rec):
    """
    Calculate the allowable tube diameter range given a receiver height (inverse of the function ReceiverHeightRange() )

    D_tube_max is undefined above ~H>4.6m,
    D_tube_min is undefined above ~H>7.8m

    Returns tube diameter [in]
    """

    if H_rec < 4.6:
        D_tube_max = 0.4062637549 - 0.00001692906721*(1.36334599*10**8 - 2.9535000*10**7*H_rec)**0.5
    else:
        D_tube_max = float('inf')

    if H_rec < 7.8:
        D_tube_min = 0.5196047009 - 0.00005341880342*(3.6811849*10**7 - 4.680000*10**6*H_rec)**0.5
    else:
        D_tube_min = float('inf')

    return [D_tube_min, D_tube_max]


def ReceiverMinimumTubeLength(Q_rec):
    """
    Calculate the minimum tube length as a function of single receiver incident thermal power using the data from Brayton

    Q_rec       [kWt] Single receiver power rating
    """

    # q_dot_rec     L_min
    #  MWt           m
    # --------------------
    # 22.3          1.7
    # 111           3.5
    # 276           5.3

    Q_rec_MWt = Q_rec * 1.e-3
    L_min = 0.41786 * Q_rec_MWt ** 0.45174

    return L_min

#----------------------------------------------------------------------------
def ReceiverTubeDesignMassFlow(D_tube, L_tube):
    """
    D_tube      [in] Outer diameter of tube
    L_tube      [m] Length of tube

    Returns the standard mass flow of a single receiver tube [kg/s]
    """

    kTubeOuterDiameters = (0.25, 0.375, 0.5);

    # These are simple linear regressions done in Excel of the 'OD Look-up Tables' worksheet
    # in 'optimization worksheetV5.xlsx'
    def ReceiverTubeDesignMassFlow_0p25_0p055(L_tube):
        return 0.0338*L_tube - 0.0278;       # [kg/s]

    def ReceiverTubeDesignMassFlow_0p375_0p070(L_tube):
        return 0.0252*L_tube - 0.04;      # [kg/s]

    def ReceiverTubeDesignMassFlow_0p5_0p100(L_tube):
        return 0.0139*L_tube - 0.0219;      # [kg/s]

    m_dot_tube_des_0p25_0p055 = ReceiverTubeDesignMassFlow_0p25_0p055(L_tube)
    m_dot_tube_des_0p375_0p070 = ReceiverTubeDesignMassFlow_0p375_0p070(L_tube)
    m_dot_tube_des_0p5_0p100 = ReceiverTubeDesignMassFlow_0p5_0p100(L_tube)

    m_dot_tubes_des = (m_dot_tube_des_0p25_0p055, m_dot_tube_des_0p375_0p070,
                       m_dot_tube_des_0p5_0p100)

    f_m_dot_tube = interp1d(kTubeOuterDiameters, m_dot_tubes_des, kind='quadratic', fill_value='extrapolate')

    return f_m_dot_tube(D_tube)[()]

def ReceiverTubeThickness(D_tube):
    """
    D_tube        [in] Tube diameter

    Returns the tube wall thickness in [in]
    """

    kTubeOuterDiameters = [1/4, 3/8, 1/2]
    kTubeWallThicknesses = [0.055, 0.070, 0.100]
    f_wall_thickness = interp1d(kTubeOuterDiameters, kTubeWallThicknesses, kind='quadratic', fill_value='extrapolate')
    return f_wall_thickness(D_tube)

#----------------------------------------------------------------------------
def calculate_cost(D, L, N_tubes):
    """
    Calculate the receiver cost 
    
    Inputs
        D - tube outer (nominal) diameter (in)
        L - Receiver tube length / aka, receiver height (m)
        N_tubes - Number of parallel tubes in the receiver (-)

    Returns
        'total_cost' - Total cost for all receivers ($)
        'A_rec' - Receiver aperture area (m2)
        'W_rec' - Receiver width (m)

    The defaults are:
    D = 0.375 [in]
    L = 5.3 [m]
    N_tubes = 14124 
    This model provided by Brayton
    """
    th_in = ReceiverTubeThickness(D)
    
    rho = 0.291 * 27679.905 #  [lbm/in3] * convert(lbm/in3,kg/m3)
    D_out = D * 0.0254 #  [in] * convert(in,m)
    th = th_in * 0.0254 #  [in] * convert(in,m)
    D_in = D_out - 2 *th
    
    L_extra = 31.44 * D_out + 0.7312 #  [m]	"extra to account for expansion loop"
    
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
    total_cost = (mat_cost + toll_cost_total + weld_cost + frame_cost +  ins_cost)*3  #3 receivers
    return {'total_cost':total_cost, 'A_rec':A_rec, 'W_rec':W_rec}


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

#----------------------------------------------------------------------------
def __interp_poly(X,Y,x):
    c0,c1,c2 = __poly_coefs(X,Y)
    return (lambda v: c0 + c1*v + c2*v**2)(x)

#----------------------------------------------------------------------------
def load_receiver_interpolator_provider(receiver_file_path, mdot_adj_factor_tube_to_rec):
    """
    receiver_file_path:
    Path to the CSV file of receiver efficiency and pressure drop 
    as a function of tube diameter, tube length, mass flow fraction, and
    inlet temperature (for efficiency) and inlet pressure (for pressure drop).
    Needs columns:
        'D_tube_inch'
        'L_tube_inch'
        'm_dot_frac_eta'
        'T_in_C'
        'mdot_Tin_id'       # ID for each unique combination of m_dot_frac_eta and T_in_C
        'eta'               # dependent variable
        'm_dot_frac_dP'
        'P_in_kPa'
        'mdot_Pin_id'       # ID for each unique combination of m_dot_frac_dP and P_in_kPa
        'dP_kPa'            # dependent variable

    mdot_adj_factor_tube_to_rec:
    The correction factor to convert the design tube mass flow used to create the receiver tables to
    the total receiver mass flow at design (which will have lower individual tube mass flows because
    there will be *more* tubes than optimal, due to rounding up to the next whole tube). This factor
    should be equal to N_tubes/N_tubes_frac, where N_tubes_frac = m_dot_rec_des / m_dot_tube_des

    Outputs:
    Two dictionaries with keys:
        'm_dot_frac_eta','T_in_C','eta' (for rec_eta_lookup)
        'm_dot_frac_dP','P_in_kPa','dP_kPa' (for rec_dP_lookup)

    'm_dot_frac_', 'T_in_C' and 'P_in_kPa' are lists of type float
    and correspond to the fluid inlet condition, in [-] and [C] (or [kPa]), respectively

    'eta and 'dP_kPa' are lists of functions:
        eta[0] = f(D_tube, L_tube) for mdot/Tin inlet condition 0
        eta[1] = f(D_tube, L_tube) for mdot/Tin inlet condition 1
        ...
        eta[120] = f(D_tube, L_tube) for mdot/Tin inlet condition 120

        dP_kPa[0] = f(D_tube, L_tube) for mdot/Pin inlet condition 0
        dP_kPa[1] = f(D_tube, L_tube) for mdot/Pin inlet condition 1
        ...
        dP_kPa[120] = f(D_tube, L_tube) for mdot/Pin inlet condition 120

    Usage example:
    >> rec_eta_lookup, rec_dP_lookup = load_receiver_interpolator_provider(<data file path> <mdot_adj_factor_tube_to_rec>)
    >> f_eta_inlet_condition_120 = rec_eta_lookup['eta'][120]
    >> eta_inlet_condition_120 = f_eta_inlet_condition_120( <tube diameter>, <tube length>)[0][0]
    """

    df = ReadAndFilterCsv(receiver_file_path)

    # rec_eta_lookup ----------------------------------------------------------------------------
    cols_eta = ['m_dot_frac_eta','T_in_C','eta']
    interp_funcs = {}
    for c in cols_eta:
        interp_funcs[c] = []

    mdot_Tin_id_values = list(set(df.mdot_Tin_id.values))
    mdot_Tin_id_values.sort()

    for id in mdot_Tin_id_values:           # for each unique m_dot_frac\T_in combination
        df_group = df[df.mdot_Tin_id == id]

        interp_funcs['m_dot_frac_eta'].append(df_group.m_dot_frac_eta.values[0] * mdot_adj_factor_tube_to_rec)
        interp_funcs['T_in_C'].append(df_group.T_in_C.values[0])

        diameters = df_group.D_tube_inch.values
        lengths = df_group.L_tube_m.values
        
        # there will only be one loop as there's only one dependent value column (eta)
        for col in cols_eta[2:]:
            interp_funcs[col].append( 
                    # SmoothBivariateSpline(diameters, lengths, df_group[col].values, kx=2, ky=2)     # works great when kx=ky=2 for just data, but not
                                                                                                    #  between D's and L's. Doesn't work well when kx=ky=1
                    # interp2d(diameters, lengths, df_group[col].values, kind='linear')   # works a lot better than SmoothBivariateSpline
                    GlobalSpline2D(diameters, lengths, df_group[col].values, kind='linear')   # adds extrapolation to interp2d
                    # Rbf(diameters, lengths, df_group[col].values, function='linear', smooth=1)      # doesn't interpolate lengths well, using any 'function'
                    # Note: griddata() returns points, not a function
                )

    rec_eta_lookup = interp_funcs


    # rec_dP_lookup ----------------------------------------------------------------------------
    cols_dP = ['m_dot_frac_dP','P_in_kPa','dP_kPa']
    interp_funcs = {}
    for c in cols_dP:
        interp_funcs[c] = []

    mdot_Pin_id_values = list(set(df.mdot_Pin_id.values))
    mdot_Pin_id_values.sort()

    for id in mdot_Pin_id_values:
        df_group = df[df.mdot_Pin_id == id]

        interp_funcs['m_dot_frac_dP'].append(df_group.m_dot_frac_dP.values[0] * mdot_adj_factor_tube_to_rec)
        interp_funcs['P_in_kPa'].append(df_group.P_in_kPa.values[0])

        diameters = df_group.D_tube_inch.values
        lengths = df_group.L_tube_m.values
        
        # there will only be one loop as there's only one dependent value column (dP_kPa)
        for col in cols_dP[2:]:
            interp_funcs[col].append( 
                    # SmoothBivariateSpline(diameters, lengths, df_group[col].values, kx=2, ky=2) 
                    # interp2d(diameters, lengths, df_group[col].values, kind='linear')
                    GlobalSpline2D(diameters, lengths, df_group[col].values, kind='linear')
                    # Rbf(diameters, lengths, df_group[col].values, function='linear', smooth=1)
                )
    rec_dP_lookup = interp_funcs


    return (rec_eta_lookup, rec_dP_lookup)

#----------------------------------------------------------------------------
def create_receiver_eta_lookup(receiver_eta_interp_provider, D_tube, L_tube):
    """
    Load the efficiency lookup provider and interpolate between tube sizes to produce a single
    efficiency lookup table for use in the performance simulation.

    Inputs
        receiver_eta_interp_provider    returned from load_receiver_interpolator_provider)
        D_tube:                         diameter of a receiver tube
        L_tube:                         length of a receiver tube

    Returns:
        Lookup table expressing receiver efficiency
            m_dot_frac_eta  |   (-) Fractional receiver mass flow
            T_in_C          |   (C) Receiver inlet temperature
            eta             |   (-) Receiver efficiency
    """
    
    interp_data = [
        receiver_eta_interp_provider['m_dot_frac_eta'],
        receiver_eta_interp_provider['T_in_C'],
    ]

    cols = ['eta']
    for i,col in enumerate(cols):
        interp_data.append( [ fun(D_tube, L_tube)[0] for fun in receiver_eta_interp_provider[col] ] )
        # interp_data.append( [ fun(D_tube, L_tube)[()] for fun in receiver_eta_interp_provider[col] ] )  # for rbf
        
    interp_data_list = array(interp_data).T.tolist()

    # Sort list, needed for proper use by SSC
    df_interp = pandas.DataFrame(interp_data_list, columns=['m_dot_frac_eta', 'T_in_C', 'eta'])
    df_interp.sort_values(by=['T_in_C', 'm_dot_frac_eta'], inplace=True)
    list_interp = df_interp.values.tolist()        # convert back to list of lists

    # Convert single element numpy arrays to value, if given by interpolator
    if isinstance(list_interp[0][2], np.ndarray):
        for i, val in enumerate(list_interp):
            list_interp[i][2] = list_interp[i][2][0]

    return list_interp

#----------------------------------------------------------------------------
def create_receiver_dP_lookup(receiver_dP_interp_provider, D_tube, L_tube):
    """
    Load the pressure drop lookup provider and interpolate between tube sizes to produce a single
    pressure lookup table for use in the performance simulation.

    Inputs
        receiver_dP_interp_provider     returned from load_receiver_interpolator_provider)
        D_tube:                         diameter of a receiver tube
        L_tube:                         length of a receiver tube

    Returns:
        Lookup table expressing receiver pressure drop
            m_dot_frac_dP   |   (-) Fractional receiver mass flow
            P_in_kPa        |   (kPa) Receiver inlet pressure
            dP_kPa          |   (kPa) Receiver pressure drop
    """
    
    interp_data = [
        receiver_dP_interp_provider['m_dot_frac_dP'],
        receiver_dP_interp_provider['P_in_kPa'],
    ]

    cols = ['dP_kPa']
    for i,col in enumerate(cols):
        interp_data.append( [ fun(D_tube, L_tube)[0] for fun in receiver_dP_interp_provider[col] ] )
        # interp_data.append( [ fun(D_tube, L_tube)[()] for fun in receiver_dP_interp_provider[col] ] )  # for rbf
        
    interp_data_list = array(interp_data).T.tolist()

    # Sort list, needed for proper use by SSC
    df_interp = pandas.DataFrame(interp_data_list, columns=['m_dot_frac_dP', 'P_in_kPa', 'eta'])
    df_interp.sort_values(by=['P_in_kPa', 'm_dot_frac_dP'], inplace=True)
    list_interp = df_interp.values.tolist()        # convert back to list of lists

    # Convert single element numpy arrays to value, if given by interpolator
    if isinstance(list_interp[0][2], np.ndarray):
        for i, val in enumerate(list_interp):
            list_interp[i][2] = list_interp[i][2][0]

    return list_interp


def ReadAndFilterCsv(receiver_file_path):
    df = pandas.read_csv(receiver_file_path)

    cond1 = np.isclose(df.D_tube_inch, 0.25) & np.isclose(df.L_tube_m, 1.6)
    cond2 = np.isclose(df.D_tube_inch, 0.25) & np.isclose(df.L_tube_m, 2.42)
    cond3 = np.isclose(df.D_tube_inch, 0.375) & np.isclose(df.L_tube_m, 4.11)
    cond4 = np.isclose(df.D_tube_inch, 0.375) & np.isclose(df.L_tube_m, 6.89)
    cond5 = np.isclose(df.D_tube_inch, 0.5) & np.isclose(df.L_tube_m, 4)
    cond6 = np.isclose(df.D_tube_inch, 0.5) & np.isclose(df.L_tube_m, 14)

    allcond = cond1 | cond2 | cond3 | cond4 | cond5 | cond6

    return df[allcond].reset_index(drop=True)

#----------------------------------------------------------------------------
def PlotReceiverVariousTubes(receiver_file_path, m_dot_frac_eta, T_in, m_dot_frac_dP, P_in):
    """
    receiver_file_path:     path to CSV file of receiver efficiency and pressure drop

    """

    # Inputs
    D_min = 1/4
    L_maxAtDmin = 3.35
    D_max = 1/2
    L_maxAtDmax = 14
    add_modeled_points_to_surface = False           # for testing surface fit of *modeled* points
    add_data_points_to_surface = True               # for testing overall interpolation to known data points

    rec_eta_lookup, rec_dP_lookup = load_receiver_interpolator_provider(receiver_file_path, 1)
    df_out = pandas.DataFrame(columns=['D_tube', 'L_tube', 'eta', 'dP_kPa'])

    TubeDiameters = np.linspace(D_min, D_max, num=30, endpoint=True)      # [in] outer diameter
    TubeLengths = np.linspace(1, L_maxAtDmax, num=30, endpoint=True)      # [m]

    for D_tube in TubeDiameters:
        for L_tube in TubeLengths:
            # Eta
            rec_eta_table_modld = create_receiver_eta_lookup(rec_eta_lookup, D_tube, L_tube)
            df_eta_modld = pandas.DataFrame(rec_eta_table_modld, columns=['m_dot_frac_eta', 'T_in_C', 'eta'])
            df_eta_modld_f = df_eta_modld[np.isclose(df_eta_modld.m_dot_frac_eta, m_dot_frac_eta) & np.isclose(df_eta_modld.T_in_C, T_in)]

            # dP
            rec_dP_table_modld = create_receiver_dP_lookup(rec_dP_lookup, D_tube, L_tube)
            df_dP_modld = pandas.DataFrame(rec_dP_table_modld, columns=['m_dot_frac_dP', 'P_in_kPa', 'dP_kPa'])
            df_dP_modld_f = df_dP_modld[np.isclose(df_dP_modld.m_dot_frac_dP, m_dot_frac_dP) & np.isclose(df_dP_modld.P_in_kPa, P_in)]

            df_out = df_out.append({'D_tube': D_tube, 'L_tube': L_tube, 'eta': df_eta_modld_f['eta'].iloc[0], 'dP_kPa': df_dP_modld_f['dP_kPa'].iloc[0]},
                ignore_index=True)

    # Filter out lengths that are greater than the assumed max per the diameter
    #  so the plot is more easily read.
    #  I.e., filter out lengths above the line that connections the max length
    #  at the lowest diameter and the max length at the highest diameter
    def filter_fn(row):
        line_slope = (L_maxAtDmax - L_maxAtDmin) / (D_max - D_min)
        line_intercept = L_maxAtDmin - line_slope * D_min
        D = row['D_tube']
        L = row['L_tube']
        L_max = line_slope * row['D_tube'] + line_intercept
        if row['L_tube'] > L_max:
            return False         # False means don't keep it
        else:
            return True

    to_filter = df_out.apply(filter_fn, axis=1)
    df_out = df_out[to_filter]
    df_out = df_out.reset_index(drop=True)

    # Overall Figure
    fig = plt.figure(figsize=(12,6))    # width, height in inches

    if add_data_points_to_surface:
        df_data = ReadAndFilterCsv(receiver_file_path)

    # Subplot 1, Eta modeled
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_trisurf(df_out['D_tube'], df_out['L_tube'], df_out['eta'], cmap=plt.cm.viridis, linewidth=0.2)
    if add_modeled_points_to_surface:
        modld_pts = ax.scatter(df_out['D_tube'], df_out['L_tube'], df_out['eta'], c='black', s=15)
    if add_data_points_to_surface:
        df_data_eta = df_data[np.isclose(df_data.m_dot_frac_eta, m_dot_frac_eta) & np.isclose(df_data.T_in_C, T_in)].reset_index(drop=True)
        data_pts = ax.scatter(df_data_eta['D_tube_inch'], df_data_eta['L_tube_m'], df_data_eta['eta'], c='red', s=10)
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlabel('D_tube [in]')
    ax.set_ylabel('L_tube [m]')
    ax.set_zlabel('eta')
    ax.set_title('Eta\n\
        m_dot_frac_eta = {m_dot_frac_eta:.3f} [-]\n\
        T_in = {T_in:.1f} [C]'\
        .format(m_dot_frac_eta=m_dot_frac_eta, T_in=T_in))

    # Subplot 2, dP modeled
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    surf = ax2.plot_trisurf(df_out['D_tube'], df_out['L_tube'], df_out['dP_kPa'], cmap=plt.cm.viridis, linewidth=0.2)
    if add_modeled_points_to_surface:
        modld_pts = ax2.scatter(df_out['D_tube'], df_out['L_tube'], df_out['dP_kPa'], c='black', s=15)
    if add_data_points_to_surface:
        df_data_dP = df_data[np.isclose(df_data.m_dot_frac_dP, m_dot_frac_dP) & np.isclose(df_data.P_in_kPa, P_in)].reset_index(drop=True)
        data_pts = ax2.scatter(df_data_dP['D_tube_inch'], df_data_dP['L_tube_m'], df_data_dP['dP_kPa'], c='red', s=10)
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    ax2.set_xlabel('D_tube [in]')
    ax2.set_ylabel('L_tube [m]')
    ax2.set_zlabel('dP [kPa]')
    ax2.set_title('dP\n\
        m_dot_frac_dP = {m_dot_frac_dP:.3f} [-]\n\
        P_in = {P_in:.0f} [kPa]'\
        .format(m_dot_frac_dP=m_dot_frac_dP, P_in=P_in))

    # Link rotation
    def on_move(event):
        if event.inaxes == ax:
            ax2.view_init(elev=ax.elev, azim=ax.azim)
        elif event.inaxes == ax2:
            ax.view_init(elev=ax2.elev, azim=ax2.azim)
        else:
            return
        fig.canvas.draw_idle()
    c1 = fig.canvas.mpl_connect('motion_notify_event', on_move)

    plt.show()

#----------------------------------------------------------------------------
def PlotReceiverTables(receiver_file_path, tube_config):
    """
    receiver_file_path:     path to CSV file of receiver efficiency and pressure drop
    tube_config:            can be either:
                                '1/4-0.055"-1.6m')  
                                '1/4-0.055"-2.42m')
                                '1/4-0.055"-3.35m')
                                '3/8-0.070"-4.11m')
                                '3/8-0.070"-5.78m')
                                '3/8-0.070"-6.89m')
                                '1/2-0.100"-4m')
                                '1/2-0.100"-9.385m')
                                '1/2-0.100"-14m')                                
    """

    if tube_config == '1/4-0.055"-1.6m':
        kTubeConfig = 0
    elif tube_config == '1/4-0.055"-2.42m':
        kTubeConfig = 1
    elif tube_config == '1/4-0.055"-3.35m':
        kTubeConfig = 2
    elif tube_config == '3/8-0.070"-4.11m':
        kTubeConfig = 3
    elif tube_config == '3/8-0.070"-5.78m':
        kTubeConfig = 4
    elif tube_config == '3/8-0.070"-6.89m':
        kTubeConfig = 5
    elif tube_config == '1/2-0.100"-4m':
        kTubeConfig = 6
    elif tube_config == '1/2-0.100"-9.385m':
        kTubeConfig = 7
    elif tube_config == '1/2-0.100"-14m':
        kTubeConfig = 8
    else:
        raise Exception('Tube configuration not supported')

    D_tube_frac = ['1/4', '1/4', '1/4', '3/8', '3/8', '3/8', '1/2', '1/2', '1/2']
    L_tube = [1.6, 2.42, 3.35, 4.11, 5.78, 6.89, 4, 9.385, 14]
    thick_tube = [0.055, 0.055, 0.055, 0.070, 0.070, 0.070, 0.100, 0.100, 0.100]

    # Reading actual values and filtering for just single tube config
    D_tube = [eval(D) for D in D_tube_frac]
    df_meas = ReadAndFilterCsv(receiver_file_path)
    df_meas_f = df_meas[np.isclose(df_meas.D_tube_inch, D_tube[kTubeConfig]) & np.isclose(df_meas.L_tube_m, L_tube[kTubeConfig])]

    # Modeled values
    rec_eta_lookup, rec_dP_lookup = load_receiver_interpolator_provider(receiver_file_path, 1)
    rec_eta_table_modld = create_receiver_eta_lookup(rec_eta_lookup, D_tube[kTubeConfig], L_tube[kTubeConfig])
    rec_dP_table_modld = create_receiver_dP_lookup(rec_dP_lookup, D_tube[kTubeConfig], L_tube[kTubeConfig])
    df_eta_modld = pandas.DataFrame(rec_eta_table_modld, columns=['m_dot_frac_eta', 'T_in_C', 'eta'])
    df_meas_f.sort_values(by=['T_in_C', 'm_dot_frac_eta'], inplace=True)
    df_eta_modld['eta_diff'] = df_eta_modld['eta'].values - df_meas_f['eta'].values       # add modeled minus measured
    df_dP_modld = pandas.DataFrame(rec_dP_table_modld, columns=['m_dot_frac_dP', 'P_in_kPa', 'dP_kPa'])
    df_meas_f.sort_values(by=['P_in_kPa', 'm_dot_frac_dP'], inplace=True)
    df_dP_modld['dP_diff'] = df_dP_modld['dP_kPa'].values - df_meas_f['dP_kPa'].values    # add modeled minus measured
    df_modld = pandas.concat([df_eta_modld, df_dP_modld], axis=1)

    # Actual overlaid on modeled
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_trisurf(df_eta_modld['m_dot_frac_eta'], df_eta_modld['T_in_C'], df_eta_modld['eta'], cmap=plt.cm.viridis, linewidth=0.2)
    # pts = ax.scatter(df_meas_f['m_dot_frac_eta'], df_meas_f['T_in_C'], df_meas_f['eta'], c='black', s=15)
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    # ax.set_xlabel('m_dot_frac')
    # ax.set_ylabel('T_in [C]')
    # ax.set_zlabel('eta')
    # ax.set_title('{diameter}-{thickness:.3f}\"-{length:.3f}m'.format(diameter=D_tube_frac[kTubeConfig],
    # thickness=thick_tube[kTubeConfig], length=L_tube[kTubeConfig]))
    # plt.show()

    # Overall Figure
    fig = plt.figure(figsize=(9,7))    # width, height in inches
    fig.suptitle('{diameter}-{thickness:.3f}\"-{length:.3f}m'.format(diameter=D_tube_frac[kTubeConfig],
        thickness=thick_tube[kTubeConfig], length=L_tube[kTubeConfig]))

    # Subplot 1, Eta actual overlaid on modeled
    ax = fig.add_subplot(2, 2, 1, projection='3d')
    surf = ax.plot_trisurf(df_modld['m_dot_frac_eta'], df_modld['T_in_C'], df_modld['eta'], cmap=plt.cm.viridis, linewidth=0.2)
    pts = ax.scatter(df_meas_f['m_dot_frac_eta'], df_meas_f['T_in_C'], df_meas_f['eta'], c='black', s=15)
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlabel('m_dot_frac')
    ax.set_ylabel('T_in [C]')
    ax.set_zlabel('eta')
    ax.set_title('Modeled and Measured')

    # Subplot 2, Eta error (absolute difference in the fractions)
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    pts = ax2.scatter(df_modld['m_dot_frac_eta'], df_modld['T_in_C'], df_modld['eta_diff'], c='red', s=15)
    ax2.set_xlabel('m_dot_frac')
    ax2.set_ylabel('T_in [C]')
    ax2.set_zlabel('eta_diff')
    ax2.set_title('Modeled - Measured [-]\n(max={max:.3f})'.format(max=max(df_modld['eta_diff'], key=abs)))     # maximum absolute value

    # Subplot 3, dP actual overlaid on modeled
    ax3 = fig.add_subplot(2, 2, 3, projection='3d')
    surf = ax3.plot_trisurf(df_modld['m_dot_frac_dP'], df_modld['P_in_kPa'], df_modld['dP_kPa'], cmap=plt.cm.viridis, linewidth=0.2)
    pts = ax3.scatter(df_meas_f['m_dot_frac_dP'], df_meas_f['P_in_kPa'], df_meas_f['dP_kPa'], c='black', s=15)
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    ax3.set_xlabel('m_dot_frac')
    ax3.set_ylabel('P_in [kPa]')
    ax3.set_zlabel('dP [kPa]')
    ax3.set_title('Modeled and Measured')

    # Subplot 4, dP error (absolute difference in kPa)
    ax4 = fig.add_subplot(2, 2, 4, projection='3d')
    pts = ax4.scatter(df_modld['m_dot_frac_dP'], df_modld['P_in_kPa'], df_modld['dP_diff'], c='red', s=15)
    ax4.set_xlabel('m_dot_frac')
    ax4.set_ylabel('P_in [kPa]')
    ax4.set_zlabel('dP_diff [kPa]')
    ax4.set_title('Modeled - Measured [kPa]\n(max={max:.3f})'.format(max=max(df_modld['dP_diff'], key=abs)))     # maximum absolute value

    # Link rotation
    def on_move(event):
        if event.inaxes == ax:
            ax2.view_init(elev=ax.elev, azim=ax.azim)
        elif event.inaxes == ax2:
            ax.view_init(elev=ax2.elev, azim=ax2.azim)
        elif event.inaxes == ax3:
            ax4.view_init(elev=ax3.elev, azim=ax3.azim)
        elif event.inaxes == ax4:
            ax3.view_init(elev=ax4.elev, azim=ax4.azim)
        else:
            return
        fig.canvas.draw_idle()
    c1 = fig.canvas.mpl_connect('motion_notify_event', on_move)

    plt.show()


#----------------------------------------------------------------------------
if __name__ == "__main__":

    # print(calculate_cost(0.375, 5.3, 14124))
    
    # X = [3, 12, 33]
    # Y = [15, 5, 40]    
    # c0,c1,c2 = __poly_coefs(X,Y)
    # Xa = numpy.arange(0, 50, .2)
    # Ya = [(lambda x: c0 + c1*x + c2*x**2)(x) for x in Xa]
    # plt.plot(X,Y, 'ks')
    # plt.plot(Xa,Ya, 'b-')
    # plt.show()

    # print(specheat_co2(550))
    # print(density_co2(550, 25))
    # create_receiver_pressure_lookup("resource/rec_pressure.csv", 4)
    # create_receiver_efficiency_lookup(5.3)

    # rec_eta_lookup, rec_dP_lookup = load_receiver_interpolator_provider('resource/rec_lookup_all.csv', 1)
    # create_receiver_eta_lookup(rec_eta_lookup, D_tube=0.25, L_tube=1.7, N_tubes=20)
    # create_receiver_dP_lookup(rec_dP_lookup, D_tube=0.25, L_tube=1.7, N_tubes=20)


    #---------------------------------------------------------------------------------------------------------------------
    #---Testing tube wall thickness function------------------------------------------------------------------------------
    # kTubeOuterDiameters = [1/4, 3/8, 1/2]
    # tube_wall_thicknesses = [ReceiverTubeThickness(D_tube) for D_tube in kTubeOuterDiameters]

    # tube_diameters = np.arange(3/16, 9/16, 1/64)                     # [in] outer diameter, TODO: replace with linspace
    # tube_wall_thicknesses_mod = [ReceiverTubeThickness(D_tube) for D_tube in tube_diameters]

    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # ax.plot(kTubeOuterDiameters, tube_wall_thicknesses, 'o')
    # ax.plot(tube_diameters, tube_wall_thicknesses_mod, '-')
    # ax.set_xlabel('D_tube [in]')
    # ax.set_ylabel('Thickness_tube [in]')
    # ax.set_title('Tube Wall Thickness [in]')
    # plt.show()

    #---------------------------------------------------------------------------------------------------------------------
    #---Testing receiver table generation---------------------------------------------------------------------------------
    # PlotReceiverTables('resource/rec_lookup_all.csv', '1/4-0.055"-1.6m')
    # PlotReceiverTables('resource/rec_lookup_all.csv', '1/4-0.055"-2.42m') 
    # # PlotReceiverTables('resource/rec_lookup_all.csv', '1/4-0.055"-3.35m') 
    # PlotReceiverTables('resource/rec_lookup_all.csv', '3/8-0.070"-4.11m')
    # # PlotReceiverTables('resource/rec_lookup_all.csv', '3/8-0.070"-5.78m') 
    # PlotReceiverTables('resource/rec_lookup_all.csv', '3/8-0.070"-6.89m') 
    # PlotReceiverTables('resource/rec_lookup_all.csv', '1/2-0.100"-4m')
    # # PlotReceiverTables('resource/rec_lookup_all.csv', '1/2-0.100"-9.385m')
    # PlotReceiverTables('resource/rec_lookup_all.csv', '1/2-0.100"-14m')   

    #---------------------------------------------------------------------------------------------------------------------
    #---Testing receiver table generation for different diameters and lengths---------------------------------------------
    # receiver_file_path = 'resource/rec_lookup_all.csv'
    # df = ReadAndFilterCsv(receiver_file_path)

    # m_dot_frac_eta_max = df['m_dot_frac_eta'].max()
    # m_dot_frac_eta_min = df['m_dot_frac_eta'].min()
    # T_in_max = df['T_in_C'].max()
    # T_in_min = df['T_in_C'].min()
    # m_dot_frac_dP_max = df['m_dot_frac_dP'].max()
    # m_dot_frac_dP_min = df['m_dot_frac_dP'].min()
    # P_in_max = df['P_in_kPa'].max()
    # P_in_min = df['P_in_kPa'].min()

    # PlotReceiverVariousTubes(receiver_file_path,
    #     m_dot_frac_eta_max, T_in_min, m_dot_frac_dP_max, P_in_max)
    # PlotReceiverVariousTubes(receiver_file_path,
    #     m_dot_frac_eta_max, T_in_max, m_dot_frac_dP_max, P_in_min)
    # PlotReceiverVariousTubes(receiver_file_path,
    #     m_dot_frac_eta_min, T_in_min, m_dot_frac_dP_min, P_in_max)
    # PlotReceiverVariousTubes(receiver_file_path,
    #     m_dot_frac_eta_min, T_in_max, m_dot_frac_dP_min, P_in_min)

    #---------------------------------------------------------------------------------------------------------------------
    #---Testing ReceiverTubeDesignMassFlow()------------------------------------------------------------------------------
    # # 1/4 - 0.055
    # L1 = [1.6, 1.7, 1.962, 2.192, 2.423, 2.654, 2.885, 3.115, 3.346, 3.577, 3.808]
    # M1 = [0.02674, 0.0298, 0.03818, 0.04592, 0.05386, 0.06192, 0.07003, 0.0781, 0.086, 0.09352, 0.1002]
    # # 3/8 - 0.070
    # L2 = [3, 3.556, 4.111, 4.2, 4.5, 4.667, 5.222, 5.3, 5.778, 6.3, 6.333, 6.889, 7.444, 8]
    # M2 = [0.03779, 0.04998, 0.06295, 0.06508, 0.07242, 0.07655, 0.09056, 0.09256, 0.1049, 0.1183, 0.1192, 0.1336, 0.1479, 0.1622]
    # # 1/2 - 0.100
    # L3 = [4, 4.769, 5.538, 6.308, 7.077, 7.846, 8.615, 9.385, 10.15, 10.92, 11.69, 12.46, 13.23, 14]
    # M3 = [0.035, 0.04497, 0.05523, 0.06567, 0.0762, 0.08688, 0.09752, 0.1083, 0.1191, 0.13, 0.141, 0.1519, 0.1628, 0.1738]

    # L_new = np.arange(1.731, 14, 0.1), TODO: replace with linspace
    # Ma_new = [ReceiverTubeDesignMassFlow(7/32, L) for L in L_new]
    # M1_new = [ReceiverTubeDesignMassFlow(1/4, L) for L in L_new]        # has data
    # Mb_new = [ReceiverTubeDesignMassFlow(9/32, L) for L in L_new]
    # Mc_new = [ReceiverTubeDesignMassFlow(5/16, L) for L in L_new]
    # Md_new = [ReceiverTubeDesignMassFlow(11/32, L) for L in L_new]
    # M2_new = [ReceiverTubeDesignMassFlow(3/8, L) for L in L_new]        # has data
    # Me_new = [ReceiverTubeDesignMassFlow(13/32, L) for L in L_new]
    # Mf_new = [ReceiverTubeDesignMassFlow(7/16, L) for L in L_new]
    # Mg_new = [ReceiverTubeDesignMassFlow(15/32, L) for L in L_new]
    # M3_new = [ReceiverTubeDesignMassFlow(1/2, L) for L in L_new]        # has data

    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # plt.plot(L_new, Ma_new, '-', label='7/32')

    # plt.plot(L1, M1, 'ro', label='1/4 - 0.055')
    # plt.plot(L_new, M1_new, 'r-', label='1/4')

    # plt.plot(L_new, Mb_new, '-', label='9/32')

    # plt.plot(L_new, Mc_new, 'g-', label='5/16')

    # plt.plot(L_new, Md_new, '-', label='11/32')

    # plt.plot(L2, M2, 'bo', label='3/8 - 0.070')
    # plt.plot(L_new, M2_new, 'b-', label='3/8')

    # plt.plot(L_new, Me_new, '-', label='13/32')

    # plt.plot(L_new, Mf_new, '-', label='7/16')

    # plt.plot(L_new, Mg_new, '-', label='15/32')

    # plt.plot(L3, M3, 'ko', label='1/2 - 0.100')
    # plt.plot(L_new, M3_new, 'k-', label='1/2')

    # ax.set_xlim(1.02, 14)
    # ax.set_ylim(-0.02, .55)
    # ax.set_xlabel('L_tube [m]')
    # ax.set_ylabel('m_dot_std [kg/s]')
    # ax.set_title('Standard Tube Mass Flow [kg/s]')
    # plt.legend()
    # plt.show()
    #---------------------------------------------------------------------------------------------------------------------
    x=None
