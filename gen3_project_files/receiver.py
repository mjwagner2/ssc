from numpy import interp, pi


#lookup for receiver performance versus height (length) in m
__raw_data = {
    'length':[3, 3.556, 4.111, 4.2, 4.5, 4.667, 5.222, 5.3, 5.778, 6.3, 6.333, 6.889, 7.444, 8],
    'efficiency':[0.6047, 0.6821, 0.7351, 0.742, 0.7627, 0.7728, 0.8002, 0.8035, 0.8195, 0.824, 0.8243, 0.8281, 0.8312, 0.8338],  # -
    'pressure_drop':[0.01905, 0.05688, 0.141, 0.1609, 0.2459, 0.3071, 0.6068, 0.6636, 1.086, 1.478, 1.507, 2.027, 2.666, 3.435],  # MPa
    'm_dot_tube':[0.01627, 0.02775, 0.04279, 0.04555, 0.05554, 0.06157, 0.08409, 0.08757, 0.109, 0.1229, 0.1238, 0.1386, 0.1535, 0.1682], #kg/s
}

def calculate_efficiency(L):
    return interp(L, __raw_data['length'], __raw_data['efficiency'])

def calculate_pressure_drop(L):
    return interp(L, __raw_data['length'], __raw_data['pressure_drop'])

def calculate_m_dot_tube(L):
    return interp(L, __raw_data['length'], __raw_data['m_dot_tube'])

def calculate_area():
    return

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


if __name__ == "__main__":
    #test
    print(calculate_cost(5.3, 14124))