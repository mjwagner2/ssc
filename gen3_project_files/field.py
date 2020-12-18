from numpy import interp, pi, array, argsort, zeros
from math import ceil
from scipy.interpolate import SmoothBivariateSpline, interp1d, interp2d, Rbf, griddata
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import pandas
from globalspline import GlobalSpline2D

def load_heliostat_interpolator_provider(efficiency_file_path, field_type):
    """
    Path to the database (csv) file for heliostat field efficiency and active area 
    as a function of sun position, power, tower height, and field type.

    Specify the data provider to generate by:
        field_type      'surround' OR 'north'
        NOTE!  Brayton's convention is to use North = 0 and Surround = 1 while this code and file is using the opposite

    Produces a dictionary with keys:
        'az','zen','eta1','eta2','eta3','a1','a2','a3'
    and value lists with an entry for each solar position. 

    The data in lists 'az' and 'zen' are of type float and correspond to sun position 
    angles in degrees. Entries in the remaining lists are function providers that
    can be called directly with arguments for power and tower height. For example,

    >> p = load_heliostat_interpolator_provider(<data file path>, <field type>)
    >> p['eta2'] ( <field power>, <tower height> )
        <return type float>

    Interpolation is provided using order=2 smooth bivariate splining from the scipy.interpolate library.
    """

    df = pandas.read_csv(efficiency_file_path)

    #get only the relevant field configuration rows
    # ============
    ## WARNING!!! Brayton's convention is to use North = 0 and Surround = 1 while this code and file is using the opposite
    # ============
    df = df[df.type == (0 if field_type == 'surround' else 1)]

    cols = ['az','zen','eta1','eta2','eta3','a1','a2','a3']
    interp_funcs = {}
    for c in cols:
        interp_funcs[c] = []


    sunid_values = list(set(df.sunid.values))
    sunid_values.sort()

    for id in sunid_values:
        df_group = df[df.sunid == id]

        interp_funcs['az'].append(df_group.az.values[0])
        interp_funcs['zen'].append(df_group.zen.values[0])

        power_levels = df_group.power.values
        tower_heights = df_group.tht.values
        

        for col in cols[2:]:
            interp_funcs[col].append( 
                    #works better than scipy.interpolate.interp2d in this case
                    # SmoothBivariateSpline(power_levels, tower_heights, df_group[col].values, kx=2, ky=2)
                    # GlobalSpline2D(power_levels, tower_heights, df_group[col].values, kind='linear')   # adds extrapolation to interp2d
                    Rbf(power_levels, tower_heights, df_group[col].values, function='thin_plate', smooth=0.0)
                    # Rbf(power_levels, tower_heights, df_group[col].values, function='linear', smooth=0.0)
                    # interp2d(power_levels, tower_heights, df_group[col].values, kind='linear')
                )    

    return interp_funcs

#----------------------------------------------------------------------------
def create_heliostat_field_lookup(field_interp_provider, q_solarfield_in_kw, h_tower, heliostat_area_m2):
    """
    Read the efficiency file and interpolate between power levels to produce a single, power-appropriate
    efficiency lookup table for use in the performance simulation.

    Inputs
        efficiency_file_path
            Relative path to efficiency file. The first row of the file lists included power levels.
            The second and following rows list efficiency and active area of each subfield as a function
            of sun position. 
            
            Columns are:
                type    |   0=surround, 1=north
                power   |   (MWt) Receiver nominal power output
                tht     |   (m) Tower optical height
                az      |   (deg) Solar azimuth angle
                zen     |   (deg) Solar zenith angle
                sunid   |   (-) order index of az/zen combinations
                eta1    |   (-) Optical efficiency of subfield 1
                eta2    |   (-) Optical efficiency of subfield 2
                eta3    |   (-) Optical efficiency of subfield 3
                a1      |   (m2) Active heliostat area of subfield 1
                a2      |   (m2) Active heliostat area of subfield 2
                a3      |   (m2) Active heliostat area of subfield 3

        q_solarfield_in_kw
            Nominal power level for the solar field at design (kw)
        h_tower
            Tower optical height (m)
        field_type
            'surround' or 'north'
        heliostat_area_m2
            Active reflective area of a *single* heliostat (m2)

    Returns:
        Lookup table expressing efficiency and number of heliostats active for each receiver. The table
        contains the following columns for an array of solar az-zen angle combinations:
            az      |   (deg) Solar azimuth angle
            zen     |   (deg) Solar zenith angle
            eta1    |   (-) Optical efficiency of subfield 1
            eta2    |   (-) Optical efficiency of subfield 2
            eta3    |   (-) Optical efficiency of subfield 3
            nh1     |   (-) Active heliostat count for subfield 1
            nh2     |   (-) Active heliostat count for subfield 2
            nh3     |   (-) Active heliostat count for subfield 3
    """

    #structure data
    all_data = {}
    
    interp_data = [
        field_interp_provider['az'],
        field_interp_provider['zen'],
    ]

    cols = ['eta1','eta2','eta3','a1','a2','a3']

    q_solarfield_in = q_solarfield_in_kw/1000.


    #collate all available data by sun position, with a list of values corresponding to each power level
    for i,col in enumerate(cols):
        
        #convert heliostat total area to heliostat count, if applicable
        scale = heliostat_area_m2 if i > 2 else 1.

        if type(field_interp_provider['eta1'][0]).__module__ == 'scipy.interpolate.rbf':
            interp_data.append( [ f(q_solarfield_in, h_tower)[()]/scale for f in field_interp_provider[col] ] )  # for rbf
        else:
            interp_data.append( [ f(q_solarfield_in, h_tower)[0]/scale for f in field_interp_provider[col] ] )
    
    list_interp = array(interp_data).T.tolist()

    # Convert single element numpy arrays to value, if given by interpolator
    if isinstance(list_interp[0][2], np.ndarray):
        for i, val in enumerate(list_interp):
            for j in range(2, len(val)):
                list_interp[i][j] = list_interp[i][j][0]

    return list_interp


def ReadAndFilterCsv(field_file_path):
    df = pandas.read_csv(field_file_path)

    # cond1 = np.isclose(df.D_tube_inch, 0.25) & np.isclose(df.L_tube_m, 1.6)
    # cond2 = np.isclose(df.D_tube_inch, 0.25) & np.isclose(df.L_tube_m, 2.42)

    # allcond = cond1 | cond2

    # return df[allcond].reset_index(drop=True)

    return df

#----------------------------------------------------------------------------
def PlotFieldTables(field_file_path, field_config):
    """
    field_file_path:        path to CSV file of field efficiency and areas
    tube_config:            can be either:
                                '66.8_MWt'
                                '333_MWt' 
                                '828_MWt' 
    """

    if field_config == '66.8_MWt':
        kFieldConfig = 0
    elif field_config == '333_MWt':
        kFieldConfig = 1
    elif field_config == '828_MWt':
        kFieldConfig = 2
    else:
        raise Exception('Field configuration not supported')

    power = [66.8, 333, 828]

    # Reading actual values and filtering for just power
    df_meas = ReadAndFilterCsv(field_file_path)
    df_meas = df_meas[np.isclose(df_meas.type, 0) & np.isclose(df_meas.power, power[kFieldConfig])].reset_index(drop=True)
    df_meas.sort_values(by=['tht', 'sunid'], inplace=True)

    # Modeled values
    interp_provider = load_heliostat_interpolator_provider(field_file_path, 'surround')
    helio_area = 8.66**2*.97
    #
    # Create eta_maps for the different heights
    column_names = ['az', 'zen', 'eta1', 'eta2', 'eta3', 'n_hel1', 'n_hel2', 'n_hel3', 'power', 'tht', 'sunid']
    df_modld = pandas.DataFrame(columns=column_names)
    tower_heights = list(set(df_meas.tht))
    for tht in tower_heights:
        eta_map = create_heliostat_field_lookup(interp_provider, power[kFieldConfig]*1000, tht, helio_area)
        df_modld_x = pandas.DataFrame(eta_map, columns=['az', 'zen', 'eta1', 'eta2', 'eta3', 'n_hel1', 'n_hel2', 'n_hel3'])
        df_modld_x['power'] = power[kFieldConfig]
        df_modld_x['tht'] = tht
        n_rows = len(df_modld_x.index)
        df_modld_x['sunid'] = np.linspace(0, n_rows, num=n_rows, endpoint=False, dtype=int)
        df_modld = df_modld.append(df_modld_x, ignore_index=True)
    df_modld.sort_values(by=['tht', 'sunid'], inplace=True)

    # add modeled minus measured
    df_modld['eta1_diff'] = df_modld['eta1'].values - df_meas['eta1'].values
    df_modld['eta2_diff'] = df_modld['eta2'].values - df_meas['eta2'].values
    df_modld['eta3_diff'] = df_modld['eta3'].values - df_meas['eta3'].values
    df_modld['a1_diff'] = (df_modld['n_hel1'].values * helio_area - df_meas['a1'].values) / df_meas['a1'].values * 100.
    df_modld['a2_diff'] = (df_modld['n_hel2'].values * helio_area - df_meas['a2'].values) / df_meas['a2'].values * 100.
    df_modld['a3_diff'] = (df_modld['n_hel3'].values * helio_area - df_meas['a3'].values) / df_meas['a3'].values * 100.

    def plot_subfield_single_height(subfield, tht):
        # Overall Figure
        fig = plt.figure(figsize=(9,7))    # width, height in inches
        fig.suptitle('q_sf_des = {power}, tower_height = {tht} [m], subfield = {subfield}'.format(power=field_config, tht=tht, subfield=subfield))

        eta_col = 'eta' + str(subfield)
        eta_diff_col = 'eta' + str(subfield) + '_diff'
        a_col_modld = 'n_hel' + str(subfield)
        a_col_meas = 'a' + str(subfield)
        a_diff_col = 'a' + str(subfield) + '_diff'

        df_modld_f = df_modld[np.isclose(list(df_modld.tht), tht)].reset_index(drop=True)
        df_meas_f = df_meas[np.isclose(list(df_meas.tht), tht)].reset_index(drop=True)

        # Subplot 1, Eta actual overlaid on modeled
        ax = fig.add_subplot(2, 2, 1, projection='3d')
        surf = ax.plot_trisurf(df_modld_f['az'], df_modld_f['zen'], df_modld_f[eta_col], cmap=plt.cm.viridis, linewidth=0.2)
        pts = ax.scatter(df_meas_f['az'], df_meas_f['zen'], df_meas_f[eta_col], c='black', s=15)
        #fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set_xlabel('azimuth')
        ax.set_ylabel('zenith')
        ax.set_zlabel('eta')
        ax.set_title('Modeled and Measured')

        # Subplot 2, Eta error (absolute difference in the fractions)
        ax2 = fig.add_subplot(2, 2, 2, projection='3d')
        pts = ax2.scatter(df_modld_f['az'], df_modld_f['zen'], df_modld_f[eta_diff_col], c='red', s=15)
        ax2.set_xlabel('azimuth')
        ax2.set_ylabel('zenith')
        ax2.set_zlabel('eta_diff')
        ax2.set_title('Modeled - Measured [-]\n(max={max:.3f})'.format(max=max(df_modld_f[eta_diff_col], key=abs)))     # maximum absolute value

        # Subplot 3, area actual overlaid on modeled
        ax3 = fig.add_subplot(2, 2, 3, projection='3d')
        surf = ax3.plot_trisurf(df_modld_f['az'], df_modld_f['zen'], df_modld_f[a_col_modld] * helio_area, cmap=plt.cm.viridis, linewidth=0.2)
        pts = ax3.scatter(df_meas_f['az'], df_meas_f['zen'], df_meas_f[a_col_meas], c='black', s=15)
        #fig.colorbar(surf, shrink=0.5, aspect=5)
        ax3.set_xlabel('azimuth')
        ax3.set_ylabel('zenith')
        ax3.set_zlabel('area_diff')
        ax3.set_title('Modeled and Measured')

        # Subplot 4, area error (absolute difference in area)
        ax4 = fig.add_subplot(2, 2, 4, projection='3d')
        pts = ax4.scatter(df_modld_f['az'], df_modld_f['zen'], df_modld_f[a_diff_col], c='red', s=15)
        ax4.set_xlabel('azimuth')
        ax4.set_ylabel('zenith')
        ax4.set_zlabel('area_diff')
        ax4.set_title('(Modeled - Measured)/Measured [%]\n(max={max:.1f})'.format(max=max(df_modld_f[a_diff_col], key=abs)))     # maximum absolute value of the percent

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

    tht = tower_heights[0]
    plot_subfield_single_height(subfield=1, tht=tht)
    plot_subfield_single_height(subfield=2, tht=tht)
    plot_subfield_single_height(subfield=3, tht=tht)

    tht = tower_heights[1]
    plot_subfield_single_height(subfield=1, tht=tht)
    plot_subfield_single_height(subfield=2, tht=tht)
    plot_subfield_single_height(subfield=3, tht=tht)

    tht = tower_heights[2]
    plot_subfield_single_height(subfield=1, tht=tht)
    plot_subfield_single_height(subfield=2, tht=tht)
    plot_subfield_single_height(subfield=3, tht=tht)

#----------------------------------------------------------------------------
def PlotFieldVariousPowersHeights(field_file_path, subfield, sunid):
    """
    field_file_path:     path to CSV file of subfield efficiencies and areas
    """

    # Inputs
    #   P_min like D_tube
    #   tht like L_tube
    P_min = 66.8
    tht_maxAtPmin = 175 # 85
    P_max = 828
    tht_maxAtPmax = 235
    helio_area = 8.66**2*.97
    add_modeled_points_to_surface = False           # for testing surface fit of *modeled* points
    add_data_points_to_surface = True               # for testing overall interpolation to known data points

    interp_provider = load_heliostat_interpolator_provider(field_file_path, 'surround')
    df_out = pandas.DataFrame(columns=['az', 'zen', 'eta1', 'eta2', 'eta3', 'n_hel1', 'n_hel2', 'n_hel3', 'power', 'tht', 'sunid'])

    powers = np.linspace(P_min, P_max, num=30, endpoint=True)         # [in] outer diameter
    thts = np.linspace(57, tht_maxAtPmax, num=30, endpoint=True)      # [m]

    for power in powers:
        for tht in thts:
            # Eta and area
            field_eta_table_modld = create_heliostat_field_lookup(interp_provider, power*1000, tht, helio_area)
            df_modld = pandas.DataFrame(field_eta_table_modld, columns=['az', 'zen', 'eta1', 'eta2', 'eta3', 'n_hel1', 'n_hel2', 'n_hel3'])
            df_modld['power'] = power
            df_modld['tht'] = tht
            n_rows = len(df_modld.index)
            df_modld['sunid'] = np.linspace(0, n_rows, num=n_rows, endpoint=False, dtype=int)
            df_modld_f = df_modld[np.isclose(df_modld.sunid, sunid)]       # keep just the sunid specified
            df_out = df_out.append(df_modld_f, ignore_index=True)

    # Filter out tower heights that are greater than the assumed max per the power
    #  so the plot is more easily read.
    #  I.e., filter out heights above the line that connections the max height
    #  at the lowest power and the max height at the highest power
    def filter_fn(row):
        line_slope = (tht_maxAtPmax - tht_maxAtPmin) / (P_max - P_min)
        line_intercept = tht_maxAtPmin - line_slope * P_min
        P = row['power']
        H = row['tht']
        tht_max = line_slope * row['power'] + line_intercept
        if row['tht'] > tht_max:
            return False         # False means don't keep it
        else:
            return True

    to_filter = df_out.apply(filter_fn, axis=1)
    df_out = df_out[to_filter]
    df_out = df_out.reset_index(drop=True)

    # Overall Figure
    fig = plt.figure(figsize=(12,6))    # width, height in inches

    if add_data_points_to_surface:
        df_data = ReadAndFilterCsv(field_file_path)

    eta_col = 'eta' + str(subfield)
    n_hel_col_modld = 'n_hel' + str(subfield)
    n_hel_col_meas = 'a' + str(subfield)

    # Subplot 1, Eta modeled
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    surf = ax.plot_trisurf(df_out['power'], df_out['tht'], df_out[eta_col], cmap=plt.cm.viridis, linewidth=0.2)
    if add_modeled_points_to_surface:
        modld_pts = ax.scatter(df_out['power'], df_out['tht'], df_out[eta_col], c='black', s=15)
    if add_data_points_to_surface:
        df_data_eta = df_data[np.isclose(df_data.type, 0) & np.isclose(df_data.sunid, sunid)].reset_index(drop=True)
        data_pts = ax.scatter(df_data_eta['power'], df_data_eta['tht'], df_data_eta[eta_col], c='red', s=10)

        # #2-D plot of just points, eta vs. tht for different powers
        # fig_2D = plt.figure()
        # ax = fig_2D.add_subplot(1, 1, 1)
        # df_data_eta.sort_values(by=['power', 'tht'], inplace=True)
        # powers = list(set(df_data_eta.power))
        # for power in powers:
        #     ax.plot(df_data_eta['tht'][np.isclose(df_data_eta.power, power)], df_data_eta[eta_col][np.isclose(df_data_eta.power, power)],\
        #         'o-', label='{power} MWt'.format(power=power))
        # ax.set_xlabel('Tower Height')
        # ax.set_ylabel('Eta')
        # plt.legend()
        # ax.set_title('Subfield {field_num} Efficiency at Sunid {sunid}'.format(field_num=subfield, sunid=sunid))
        # plt.show()

    # fig.colorbar(surf, shrink=0.5, aspect=5)
    ax.set_xlabel('power [MWt]')
    ax.set_ylabel('tht [m]')
    ax.set_zlabel('eta')
    ax.set_title('Eta\n\
        subfield = {subfield} [-]\n\
        sunid = {sunid} [-]'\
        .format(subfield=subfield, sunid=sunid))

    # Subplot 2, area modeled
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    surf = ax2.plot_trisurf(df_out['power'], df_out['tht'], df_out[n_hel_col_modld] * helio_area, cmap=plt.cm.viridis, linewidth=0.2)
    if add_modeled_points_to_surface:
        modld_pts = ax2.scatter(df_out['power'], df_out['tht'], df_out[n_hel_col_modld] * helio_area, c='black', s=15)
    if add_data_points_to_surface:
        df_data_area = df_data[np.isclose(df_data.type, 0) & np.isclose(df_data.sunid, sunid)].reset_index(drop=True)
        data_pts = ax2.scatter(df_data_area['power'], df_data_area['tht'], df_data_area[n_hel_col_meas], c='red', s=10)
    # fig.colorbar(surf, shrink=0.5, aspect=5)
    ax2.set_xlabel('power [MWt]')
    ax2.set_ylabel('tht [m]')
    ax2.set_zlabel('area [m2]')
    ax2.set_title('area\n\
        subfield = {subfield} [-]\n\
        sunid = {sunid} [-]'\
        .format(subfield=subfield, sunid=sunid))

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
if __name__ == "__main__":
    intp = load_heliostat_interpolator_provider('resource/eta_lookup_all.csv', 'surround')
    create_heliostat_field_lookup(intp, 660000, 215, 88)

    #---------------------------------------------------------------------------------------------------------------------
    #---Testing field table generation---------------------------------------------------------------------------------
    # PlotFieldTables('resource/eta_lookup_all.csv', '66.8_MWt')
    # PlotFieldTables('resource/eta_lookup_all.csv', '333_MWt')
    # PlotFieldTables('resource/eta_lookup_all.csv', '828_MWt')

    #---------------------------------------------------------------------------------------------------------------------
    #---Testing field table generation for different diameters and lengths------------------------------------------------
    field_file_path = 'resource/eta_lookup_all.csv'
    sunid_1 = 0     # az =  70, zen = 77
    sunid_2 = 3     # az = 180, zen = 11.4
    sunid_3 = 5     # az = 275, zen = 53

    # PlotFieldVariousPowersHeights(field_file_path, subfield=1, sunid=sunid_1)
    # PlotFieldVariousPowersHeights(field_file_path, subfield=1, sunid=sunid_2)
    PlotFieldVariousPowersHeights(field_file_path, subfield=1, sunid=sunid_3)

    # PlotFieldVariousPowersHeights(field_file_path, subfield=2, sunid=sunid_1)
    # PlotFieldVariousPowersHeights(field_file_path, subfield=2, sunid=sunid_2)
    # PlotFieldVariousPowersHeights(field_file_path, subfield=2, sunid=sunid_3)

    # PlotFieldVariousPowersHeights(field_file_path, subfield=3, sunid=sunid_1)
    # PlotFieldVariousPowersHeights(field_file_path, subfield=3, sunid=sunid_2)
    # PlotFieldVariousPowersHeights(field_file_path, subfield=3, sunid=sunid_3)

    x=None