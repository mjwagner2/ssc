from numpy import interp, pi

#-------------------------------------------------------------------------------------

__raw_data = {
    'T':[37.8, 93.3, 148.9, 204.4, 260.0, 315.6, 343.3, 371.1, 398.9, 426.7, 454.4, 482.2, 510.0, 537.8, 565.6, 593.3, 621.1, 648.9],  #C
    'S_allow':[171.6855, 170.3065, 151.69, 137.2105, 127.5575, 122.0415, 119.973, 118.594, 117.215, 115.836, 114.457, 113.078, 111.699, 102.7355, 79.982, 62.055, 47.5755, 35.854],  #MPa
    'y':[0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.5, 0.7, 0.7],
    'W':[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.95, 0.91, 0.86, 0.82, 0.77],
}

#-------------------------------------------------------------------------------------

def solve(D_in, L_LT, P):
    """
    material is 347 seam welded pipe with filler material, UNS code S34700, Class 1 & 3, deformation allowed
    
    Defaults:
    D_in = 0.490 [m]
    L_LT = 550 [m]
    P = 27.5 [MPa]	"needs to be 110% of maximum operating pressure"
    """
    T = 580 #[C]	"Should be the maximum expected operating temperature, however code allows some excursions"
    # T_f = converttemp(C,F,T)
    
    # sigma_allow = interpolate1('347', 'T', 'S_allow', 'T' = T_f) * convert(ksi,MPa)
    sigma_allow = interp(T, __raw_data['T'], __raw_data['S_allow'])
    
    # y =interpolate1('347', 'T', 'y', 'T' = T_f)
    y = interp(T, __raw_data['T'], __raw_data['y'])
    # W =interpolate1('347', 'T', 'W', 'T' = T_f)
    W = interp(T, __raw_data['T'], __raw_data['W'])
    
    A = 0 #[m]
    th = (P * D_in + 2 * sigma_allow * W* A + 2 * y * P * A)/(2*(sigma_allow * W+ P*y - P))
    
    D_out = D_in + 2 * th
    
    A_cs = pi/4*(D_out**2-D_in**2)
    V = A_cs * L_LT
    # rho = density('Stainless_AISI347', T = 20 [C])
    rho=7999 #[kg/m^3]
    mass = V*rho
    cost_rate = 7*2.2046226 #[$/lbm] * convert(1/lbm,1/kg)
        
    cost = cost_rate * mass
    
    cost_length = cost/L_LT

    return {'th':th, 'D_out':D_out, 'cost':cost, 'cost_length':cost_length}

#-------------------------------------------------------------------------------------
if __name__ == "__main__":

    print(solve(0.490, 550, 27.5))
