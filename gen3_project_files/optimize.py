import gen3gas as G3
import scipy.optimize 
import multiprocessing


#-------------------------------
def f_eval(x, data):

    # z = sum([v*v for v in x])

    # return z
    x_unscaled = [x[i]*data.x_initial[i] for i in range(len(x))]

    data.variables.cycle_design_power,\
        data.variables.solar_multiple,\
        data.variables.dni_design_point,\
        data.variables.receiver_height,\
        data.variables.riser_inner_diam, \
        data.variables.downcomer_inner_diam, \
        data.variables.hours_tes,\
        data.variables.dT_approach_charge_hx,\
        data.variables.dT_approach_disch_hx = x_unscaled
    
    data.exec()
    lcoe = data.get_result_value('LCOE (real)')
    if lcoe < 0.1 or lcoe > 35:
        lcoe = float('nan')
    
    
    logline = "{:20s}Iter: {:04d}".format(data.casename, data.current_iteration) + ("{:>10s}"*(len(x)+1)).format(*["{:.3f}".format(v) for v in [lcoe]+x_unscaled])
    data.optimization_log += "\n" + logline

    print(logline)

    data.current_iteration += 1
    return lcoe

def f_callback(xk): #, data):
    # print(".... Iteration complete")
    xyz=None
    return

def optimize(casenum):

    casename = ["sm-north-bucket", "sm-surround-bucket", "sm-north-skip", "sm-surround-skip"]
    case = casename[casenum]
    # for case in casename:

    print(">>>>> Running case " + case)
    g = G3.Gen3opt()

    g.settings.is_north = 'north' in case
    g.optimization_log = ""
    g.settings.lift_technology = 'bucket' if 'bucket' in case else 'skip'
    g.casename = case 
    # g.variables.cycle_design_power = 50

    g.current_iteration = 0
    
    x0 = [
        g.variables.cycle_design_power,
        g.variables.solar_multiple,
        g.variables.dni_design_point,
        g.variables.receiver_height,
        g.variables.riser_inner_diam, 
        g.variables.downcomer_inner_diam, 
        g.variables.hours_tes,
        g.variables.dT_approach_charge_hx,
        g.variables.dT_approach_disch_hx,
    ]

    g.x_initial = [v for v in x0]

    x0 = [1. for v in x0]

    res = scipy.optimize.fmin(f_eval, x0, args = ((g,)), ftol=0.01, xtol=0.01, maxfun=250, callback=f_callback) #, callback=f_update)

    fout = open('optimization-log-'+case+'.txt', 'w')
    fout.write(g.optimization_log)
    fout.close()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    nthreads = 4
    all_args = [[i] for i in range(4)]

    pool = multiprocessing.Pool(processes=nthreads)
    pool.starmap(optimize, all_args)
    
    
    # optimize(0)
