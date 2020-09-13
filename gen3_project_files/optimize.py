import gen3gas as G3
import receiver as G3rec
import field as G3field
import scipy.optimize 
import multiprocessing
import random
import time
import os


#-------------------------------
def log_entry(x, z, iternum, label):
    return "{:20s}Iter: {:04d}".format(label, iternum) + ("{:>10s}"*(len(x)+1)).format(*["{:.3f}".format(v) for v in [z]+x])

#-------------------------------
# Constraint function evaluation
def fconst_eval(x, data):

    # Unscale scaled (normalized) parameters from optimizer
    x_unscaled = [x[i]*data.x_scalers[i] for i in range(len(x))]
    # unscale select initially scaled values:
    # OR IS THIS THE PROBLEM WITH THE GRADIENTS MOVING UNDER THE OPTIMIZER'S FEET??
    [receiver_height_min, receiver_height_max] = G3rec.ReceiverHeightRange(x_unscaled[5])       # [5] = receiver_tube_diam
    # data.variables.receiver_height = data.variables.receiver_height * (receiver_height_max - receiver_height_min) + receiver_height_min
    x_unscaled[4] = x_unscaled[4] * (receiver_height_max - receiver_height_min) + receiver_height_min   # [4] = receiver_height

    # Assign variables before calling exec
    data.variables.cycle_design_power,\
        data.variables.solar_multiple,\
        data.variables.h_tower,\
        data.variables.dni_design_point,\
        data.variables.receiver_height,\
        data.variables.receiver_tube_diam,\
        data.variables.riser_inner_diam, \
        data.variables.downcomer_inner_diam, \
        data.variables.hours_tes,\
        data.variables.dT_approach_charge_hx,\
        data.variables.dT_approach_ht_disch_hx,\
        data.variables.dT_approach_lt_disch_hx = x_unscaled

    # data.variables.dT_approach_ht_disch_hx = data.variables.dT_approach_lt_disch_hx
    
    q_sf_des = data.exec(sf_des_only=True)
    
    L_min = G3rec.ReceiverMinimumTubeLength(q_sf_des * 1.e3 / 3)

    #The inequality constraints are feasible if positive
    return data.variables.receiver_height - L_min


#-------------------------------
# Objective function evaluation
def f_eval(x, data):

    # Unscale scaled (normalized) parameters from optimizer
    x_unscaled = [x[i]*data.x_scalers[i] for i in range(len(x))]
    # unscale select initially scaled values:
    # OR IS THIS THE PROBLEM WITH THE GRADIENTS MOVING UNDER THE OPTIMIZER'S FEET??
    [receiver_height_min, receiver_height_max] = G3rec.ReceiverHeightRange(x_unscaled[5])       # [5] = receiver_tube_diam
    # data.variables.receiver_height = data.variables.receiver_height * (receiver_height_max - receiver_height_min) + receiver_height_min
    x_unscaled[4] = x_unscaled[4] * (receiver_height_max - receiver_height_min) + receiver_height_min   # [4] = receiver_height

    data.variables.cycle_design_power,\
        data.variables.solar_multiple,\
        data.variables.h_tower,\
        data.variables.dni_design_point,\
        data.variables.receiver_height,\
        data.variables.receiver_tube_diam,\
        data.variables.riser_inner_diam, \
        data.variables.downcomer_inner_diam, \
        data.variables.hours_tes,\
        data.variables.dT_approach_charge_hx,\
        data.variables.dT_approach_ht_disch_hx,\
        data.variables.dT_approach_lt_disch_hx = x_unscaled

    # data.variables.dT_approach_ht_disch_hx = data.variables.dT_approach_lt_disch_hx
    
    try:
        simok = data.exec()
    except Exception as E:
        print("{:s}: {0}".format(data.casename, E))
        simok = False

    if not simok:
        data.current_iteration += 1
        print("function evaluation NAN", x_unscaled)
        return float('nan')
        
    lcoe = data.get_result_value('LCOE (real)')
    if lcoe < 0.1 or lcoe > 50:
        lcoe = float('nan')

    # Track parameters of current best LCOE
    if lcoe < data.z_best['z']:
        data.z_best['z'] = lcoe
        data.z_best['xk'] = [v for v in x_unscaled]
        data.z_best['iter'] = data.current_iteration
    

    logline = log_entry(x_unscaled, lcoe, data.current_iteration, data.casename)
    data.optimization_log += "\n" + logline

    time_elapsed = time.time() - data.clock_time_start
    timestamp = "{:03d}:{:02d}   ".format( int(time_elapsed/60), int(time_elapsed % 60) )

    L_min = G3rec.ReceiverMinimumTubeLength(data.variables.q_sf_des * 1.e3 / 3)
    q_sf_des_single_and_lmin = " {q_sf_des_single:.1f} -> {L_min:.2f}   ".format(q_sf_des_single=data.variables.q_sf_des/3, L_min=L_min)
    print(timestamp + q_sf_des_single_and_lmin + logline)

    data.current_iteration += 1
    return lcoe #+ max([ (data.variables.dT_approach_charge_hx + data.variables.dT_approach_disch_hx - 30), 0 ])*2.

def f_callback(xk): #, data):
    # print(".... Iteration complete")
    return

def optimize(thread_id, sf_interp_provider):

    # choose the case based on the thread_id integer
    casename = ["north-bucket", "surround-bucket", "north-skip", "surround-skip"][thread_id % 4]
    case = "{:03d}_{:s}".format(thread_id, casename)

    # instantiate the case
    g = G3.Gen3opt(sf_interp_provider = sf_interp_provider)

    g.settings.is_north = 'north' in case
    
    # force attributes for optimization
    g.optimization_log = ""
    g.settings.lift_technology = 'bucket' if 'bucket' in case else 'skip'
    g.casename = case 
    g.current_iteration = 0
    g.clock_time_start = time.time()
    
    # set variable bounds
    xb = [
        [   15  ,  150  ],   # [0] cycle_design_power
        [   2.5 ,  3.5  ],   # [1] solar_multiple
        [   50  ,  200  ],   # [2] h_tower
        [   650 ,  1200 ],   # [3] dni_design_point
        # [   3   ,  8    ],   # [4] receiver_height
        [   0   ,  1    ],   # [4] receiver_height, normalized
        [   .25 ,  .375 ],   # [5] receiver tube outside diameter
        [   .3 ,   .75  ],   # [6] riser_inner_diam
        [   .3 ,   .75  ],   # [7] downcomer_inner_diam
        [   4   ,  20   ],   # [8] hours_tes
        [   10  ,  40   ],   # [9] dT_approach_charge_hx
        [   10  ,  40   ],   # [10] dT_approach_ht_disch_hx
        [   10  ,  40   ],   # [11] dT_approach_lt_disch_hx
    ]
    
    # Randomly choose initial guess values for most variables, but correlate some
    x0 = [random.uniform(x[0], x[1]) for x in xb]
    x0[2] = g.variables.guess_h_tower(cycle_design_power = x0[0], solar_multiple = x0[1], is_north = g.settings.is_north) # correlate tower height guess to power
    x0[6] = xb[6][0] + (x0[2] - xb[2][0]) * (xb[6][1] - xb[6][0]) / (xb[2][1] - xb[2][0])    # set riser_inner_diam guess at the same fraction across its range as 
                                                                                             #  the receiver height is across its range b/c they're positively correlated
    x0[7] = xb[7][0] + (x0[2] - xb[2][0]) * (xb[7][1] - xb[7][0]) / (xb[2][1] - xb[2][0])    # set downcomer_inner_diam guess at the same fraction across its range as 
                                                                                             #  the receiver height is across its range b/c they're positively correlated

    # normalize bounds using respective guess values
    for i in range(len(x0)):
        for j in range(2):
            xb[i][j] /= x0[i]

    g.x_scalers = [v for v in x0]   # used in the objective function for unscaling the normalized values sent by the optimizer
    x0 = [1. for v in x0]           # normalized initial guesses

    g.z_best = {'z':float('inf'), 'xk':[v for v in x0], 'iter':-1}  # initialize best point tracker

    # print column labels for terminal output
    print("time      Q_rec -> L_min  ###_field-lift      Iter: ####    LCOE     P_ref       solarm  H_tower  dni_des      H_rec     D_rec_tb  D_riser   D_dcomr  tshours   dT_chrg   dt_ht_dchrg   dt_lt_dchrg")

    # call optimize
    # scipy.optimize.minimize(f_eval, x0, args = ((g,)), method='SLSQP', tol=0.001, 
    #                         options={'maxiter':10, 'eps':0.1, 'ftol':0.001}, 
    #                         bounds=xb, callback=f_callback)

    scipy.optimize.fmin_slsqp(f_eval, x0, bounds=xb, args=((g,)), ieqcons=[fconst_eval], epsilon=0.1, iter=20, acc=0.001)

    # report best point
    logline = log_entry(g.z_best['xk'], g.z_best['z'], g.z_best['iter'], "***Best point:")
    g.optimization_log += "\n\n" + logline

    # write a summary log
    if not os.path.exists('runs'):
        os.makedirs('runs')
    fout = open('runs/optimization-log-'+case+'.txt', 'w')
    fout.write(g.optimization_log)
    fout.close()


if __name__ == "__main__":

    north_interp_provider = G3field.load_heliostat_interpolator_provider('resource/eta_lookup_all.csv', 'north')
    surr_interp_provider = G3field.load_heliostat_interpolator_provider('resource/eta_lookup_all.csv', 'surround')
    
    # Run different field-lift combinations on different threads
    # nthreads = 4
    # nreplicates = 1
    # -------
    # all_args = []
    # for i in range(nreplicates*4):
    #     all_args.append([i, north_interp_provider if i % 2 == 0 else surr_interp_provider]) 
    # -------
    # pool = multiprocessing.Pool(processes=nthreads)
    # pool.starmap(optimize, all_args)
    
    id=3
    optimize(id, north_interp_provider if id%2 == 0 else surr_interp_provider)
