import gen3gas as G3
import scipy.optimize 
import multiprocessing
import random


#-------------------------------
def log_entry(x, z, iternum, label):
    return "{:20s}Iter: {:04d}".format(label, iternum) + ("{:>10s}"*(len(x)+1)).format(*["{:.3f}".format(v) for v in [z]+x])

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
    
    try:
        data.exec()
    except:
        return float('nan')
        
    lcoe = data.get_result_value('LCOE (real)')
    if lcoe < 0.1 or lcoe > 35:
        lcoe = float('nan')

    if lcoe < data.z_best['z']:
        data.z_best['z'] = lcoe
        data.z_best['xk'] = [v for v in x_unscaled]
        data.z_best['iter'] = data.current_iteration
    
    # logline = "{:20s}Iter: {:04d}".format(data.casename, data.current_iteration) + ("{:>10s}"*(len(x)+1)).format(*["{:.3f}".format(v) for v in [lcoe]+x_unscaled])
    logline = log_entry(x_unscaled, lcoe, data.current_iteration, data.casename)
    data.optimization_log += "\n" + logline

    print(logline)

    data.current_iteration += 1
    return lcoe #+ max([ (data.variables.dT_approach_charge_hx + data.variables.dT_approach_disch_hx - 30), 0 ])*2.

def f_callback(xk): #, data):
    # print(".... Iteration complete")
    xyz=None
    return

def optimize(thread_id):

    casename = ["north-bucket", "surround-bucket", "north-skip", "surround-skip"][thread_id % 4]
    case = "{:03d}_{:s}".format(thread_id, casename)

    # print(">>>> " + case)

    g = G3.Gen3opt()

    g.settings.is_north = 'north' in case
    g.optimization_log = ""
    g.settings.lift_technology = 'bucket' if 'bucket' in case else 'skip'
    g.casename = case 
    # g.variables.cycle_design_power = 50

    g.current_iteration = 0
    
    #set variable bounds
    xb = [
        [   15  ,  150  ],   # cycle_design_power
        [   2.5 ,  3.5  ],   # solar_multiple
        [   650 ,  1200 ],   # dni_design_point
        [   3   ,  8    ],   # receiver_height
        [   .25 ,  .75  ],   # riser_inner_diam
        [   .25 ,  .75  ],   # downcomer_inner_diam
        [   4   ,  20   ],   # hours_tes
        [   10  ,  40   ],   # dT_approach_charge_hx
        [   10  ,  40   ],   # dT_approach_disch_hx
    ]
    
    #initial guess variable values
    x0 = [random.uniform(x[0], x[1]) for x in xb]
    
    for i in range(len(x0)):
        for j in range(2):
            xb[i][j] /= x0[i]

    #save 
    g.x_initial = [v for v in x0]
    #initialize best point tracker
    g.z_best = {'z':float('inf'), 'xk':[v for v in x0], 'iter':-1}

    #variables will be normalized
    x0 = [1. for v in x0]

    #call optimize
    scipy.optimize.minimize(f_eval, x0, args = ((g,)), method='SLSQP', tol=0.001, 
                            options={'maxiter':200, 'eps':0.1, 'ftol':0.001}, 
                            bounds=xb, callback=f_callback)

    logline = log_entry(g.z_best['xk'], g.z_best['z'], g.z_best['iter'], "***Best point:")
    g.optimization_log += "\n\n" + logline

    fout = open('runs/optimization-log-'+case+'.txt', 'w')
    fout.write(g.optimization_log)
    fout.close()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    nthreads = 14
    all_args = [[i] for i in range(100)]

    pool = multiprocessing.Pool(processes=nthreads)
    pool.starmap(optimize, all_args)
    
    
    # optimize(0)
