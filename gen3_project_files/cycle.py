from numpy import interp

__power_lookup = \
{
    "power":[10,20,30,40,50,60,70,80,90,100],
    "correction":[-0.054,-0.034,-0.0235,-0.017,-0.012,-0.008,-0.006,-0.004,-0.002,0]    #absolute reduction in efficiency
}


#----------------------------------------------------------------------------
def calculate_nominal_efficiency(T_cycle_in_C, W_cycle_des_kW):
    #maximum efficiency based on temperature
    eta_max = 3.857e-4 * T_cycle_in_C + 0.2269  #from Echogen chart provided to Brayton

    #correction for reduced power block size
    eta_correction = interp(W_cycle_des_kW/1000., __power_lookup["power"], __power_lookup["correction"])

    return eta_max - eta_correction


#----------------------------------------------------------------------------
def create_updc_lookup(file_path, T_cycle_in_C):

    cycle_adj =  T_cycle_in_C - 700

    #load table
    raw_data = [[float(v) for v in line.split(",")] for line in open(file_path,'r').readlines()]
    for row in raw_data:
        row[0] += cycle_adj

    return raw_data

def create_indirect_updc_lookup(file_path, T_cycle_HTF_in_C):

    cycle_adj =  T_cycle_HTF_in_C - 715

    #load table
    raw_data = [[float(v) for v in line.split(",")] for line in open(file_path,'r').readlines()]
    for row in raw_data:
        row[0] += cycle_adj

    return raw_data
