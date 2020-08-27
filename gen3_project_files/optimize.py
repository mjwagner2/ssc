import gen3gas as G3
import receiver as G3rec
import field as G3field
import scipy.optimize 
import multiprocessing
import random
import time


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
        data.variables.h_tower,\
        data.variables.dni_design_point,\
        data.variables.receiver_height,\
        data.variables.riser_inner_diam, \
        data.variables.downcomer_inner_diam, \
        data.variables.hours_tes,\
        data.variables.dT_approach_charge_hx,\
        data.variables.dT_approach_ht_disch_hx,\
        data.variables.dT_approach_lt_disch_hx = x_unscaled
    
    try:
        simok = data.exec()
    except Exception as E:
        print("{:s}: {0}".format(data.casename, E))
        simok = False

    if not simok:
        data.current_iteration += 1
        return float('nan')
        
    lcoe = data.get_result_value('LCOE (real)')
    if lcoe < 0.1 or lcoe > 50:
        lcoe = float('nan')

    if lcoe < data.z_best['z']:
        data.z_best['z'] = lcoe
        data.z_best['xk'] = [v for v in x_unscaled]
        data.z_best['iter'] = data.current_iteration
    

    logline = log_entry(x_unscaled, lcoe, data.current_iteration, data.casename)
    data.optimization_log += "\n" + logline

    time_elapsed = time.time() - data.clock_time_start
    timestamp = "{:03d}:{:02d}   ".format( int(time_elapsed/60), int(time_elapsed % 60) )

    print(timestamp + logline)

    data.current_iteration += 1
    return lcoe #+ max([ (data.variables.dT_approach_charge_hx + data.variables.dT_approach_disch_hx - 30), 0 ])*2.

def f_callback(xk): #, data):
    # print(".... Iteration complete")
    return

def optimize(thread_id, sf_interp_provider):

    #choose the case based on the thread_id integer
    casename = ["north-bucket", "surround-bucket", "north-skip", "surround-skip"][thread_id % 4]
    case = "{:03d}_{:s}".format(thread_id, casename)

    #instantiate the case
    g = G3.Gen3opt(sf_interp_provider = sf_interp_provider)

    g.settings.is_north = 'north' in case
    
    #force attributes for optimization
    g.optimization_log = ""
    g.settings.lift_technology = 'bucket' if 'bucket' in case else 'skip'
    g.casename = case 
    g.current_iteration = 0
    g.clock_time_start = time.time()
    
    #set variable bounds
    xb = [
        [   15  ,  150  ],   # cycle_design_power
        [   2.5 ,  3.5  ],   # solar_multiple
        [   50  ,  250  ],   # h_tower
        [   650 ,  1200 ],   # dni_design_point
        [   3   ,  8    ],   # receiver_height
        [   .25 ,  .75  ],   # riser_inner_diam
        [   .25 ,  .75  ],   # downcomer_inner_diam
        [   4   ,  20   ],   # hours_tes
        [   10  ,  40   ],   # dT_approach_charge_hx
        [   10  ,  40   ],   # dT_approach_ht_disch_hx
        [   10  ,  40   ],   # dT_approach_lt_disch_hx
    ]
    
    #initial guess variable values
    x0 = [random.uniform(x[0], x[1]) for x in xb]
    
    #tie tower height guess to power
    x0[2] = g.variables.guess_h_tower(cycle_design_power = x0[0], solar_multiple = x0[1], is_north = g.settings.is_north)

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

    #report best point
    logline = log_entry(g.z_best['xk'], g.z_best['z'], g.z_best['iter'], "***Best point:")
    g.optimization_log += "\n\n" + logline

    #write a summary log
    fout = open('runs/optimization-log-'+case+'.txt', 'w')
    fout.write(g.optimization_log)
    fout.close()


if __name__ == "__main__":


    north_interp_provider = G3field.load_heliostat_interpolator_provider('resource/eta_lookup_all.csv', 'north')
    surr_interp_provider = G3field.load_heliostat_interpolator_provider('resource/eta_lookup_all.csv', 'surround')
    
    # nthreads = 4
    # nreplicates = 1

    # all_args = []
    # for i in range(nreplicates*4):
    #     all_args.append([i, north_interp_provider if i % 2 == 0 else surr_interp_provider]) 

    # pool = multiprocessing.Pool(processes=nthreads)
    # pool.starmap(optimize, all_args)
    

    id=3
    optimize(id, north_interp_provider if id%2 == 0 else surr_interp_provider)
