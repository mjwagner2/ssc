import pandas
from numpy import array
from scipy.interpolate import SmoothBivariateSpline

def load_heliostat_interpolator_provider(efficiency_file_path, field_type):
    """
    Path to the database (csv) file for heliostat field efficiency and active area 
    as a function of sun position, power, tower height, and field type.

    Specify the data provider to generate by:
        field_type      'surround' OR 'north'

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
                    SmoothBivariateSpline(power_levels, tower_heights, df_group[col].values, kx=2, ky=2 ) 
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

        interp_data.append( [ f(q_solarfield_in, h_tower)[0][0]/scale for f in field_interp_provider[col] ] )
        
    return array(interp_data).T.tolist()