# import xlrd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def PlotScatter(file_loc):
    """Open xlsx file and plot data from all sheets, separate scatter plots for each column
    see: https://gist.github.com/jiffyclub/9ab668f63c3d0f9adf3e730dc37cd419"""

    # Gather data from all sheets into a single dataframe
    xls_file = pd.ExcelFile(file_loc)
    sheet_names = xls_file.sheet_names
    df = pd.DataFrame()
    for sheet in sheet_names:
        df = df.append(xls_file.parse(sheet), ignore_index=True)

    # Overall Figure
    fig = plt.figure(figsize=(17,10))
    fig.suptitle('surround-skip -- LCOE, normalized\n\
    black = within all constraints\n\
    blue = just exceeds receiver height range per diameter\n\
    red = below receiver minimum height per power rating')

    # Subplots
    def PlotSubplot(df, fig, grid_dims, grid_index, x_col_name, x_axis_label):
        ax = fig.add_subplot(grid_dims[0], grid_dims[1], grid_index)
        # pts = ax.scatter(df[x_col_name], df['LCOE'], c='black', s=1)
        
        H_rec_below_min = df['H_rec_above_min'] == 0
        H_rec_above_min_OOB = (df['H_rec_above_min'] == 1) & (df['H_rec_in_bounds'] == 0)
        H_rec_above_min_in_bounds =  (df['H_rec_above_min'] == 1) & (df['H_rec_in_bounds'] == 1)

        LCOE_min = df['LCOE'].min()

        pts = ax.scatter(df.loc[H_rec_below_min, [x_col_name]], df.loc[H_rec_below_min, ['LCOE']] / LCOE_min, c='red', s=3)
        pts = ax.scatter(df.loc[H_rec_above_min_OOB, [x_col_name]], df.loc[H_rec_above_min_OOB, ['LCOE']] / LCOE_min, c='blue', s=3)
        pts = ax.scatter(df.loc[H_rec_above_min_in_bounds, [x_col_name]], df.loc[H_rec_above_min_in_bounds, ['LCOE']] / LCOE_min, c='black', s=3)
        # ax.set_ylim(7, 16)
        ax.set_ylim(0.9, 2)
        ax.set_xlabel(x_axis_label)
        return

    grid_dims = [3, 3]
    PlotSubplot(df, fig, grid_dims, 1, 'P_ref',     'cycle_design_power')
    PlotSubplot(df, fig, grid_dims, 2, 'solarm',    'solar_multiple')
    PlotSubplot(df, fig, grid_dims, 3, 'dni_des',   'dni_design_point')
    PlotSubplot(df, fig, grid_dims, 4, 'H_rec',     'receiver_height')
    PlotSubplot(df, fig, grid_dims, 5, 'D_riser',   'riser_inner_diam')
    PlotSubplot(df, fig, grid_dims, 6, 'tshours',   'hours_tes')
    PlotSubplot(df, fig, grid_dims, 7, 'dT_chrg',   'dT_approach_charge_hx')
    PlotSubplot(df, fig, grid_dims, 8, 'dT_htHX',   'dT_approach_HT_disch_hx')
    PlotSubplot(df, fig, grid_dims, 9, 'dT_ltHX',   'dT_approach_LT_disch_hx')

    plt.show()

    pass

def PlotBubble(file_loc):
    """Open xlsx file and plot data from all sheets, separate bubble plots for each column"""

    # Gather data from all sheets into a single dataframe
    xls_file = pd.ExcelFile(file_loc)
    sheet_names = xls_file.sheet_names
    df = pd.DataFrame()
    for sheet in sheet_names:
        df = df.append(xls_file.parse(sheet), ignore_index=True)

    # Overall Figure
    fig = plt.figure(figsize=(12,8))
    fig.suptitle('surround-skip\n\
        Smaller marker = Lower LCOE')
    # black = within all constraints\n\
    # blue = above receiver minimum height, but exceeds receiver height upper 'limit' per diameter\n\

    # Subplots
    def PlotSubplot(df, fig, grid_dims, grid_index, x_col_name, x_axis_label, y_col_name, y_axis_label):
        ax = fig.add_subplot(grid_dims[0], grid_dims[1], grid_index)
        
        # Sort ascending by LCOE
        df.sort_values(by=['LCOE'], ascending=True, inplace=True)

        # Plot only the first 50 points, excluding those with the receiver height below the minimum      
        H_rec_below_min = df['H_rec_above_min'] == 0
        df = df.loc[df['H_rec_above_min'] == 1].head(n=200)
        df.reset_index(inplace=True)


        # Use 'try' in case there are no values that fit the filters
        try:
            H_rec_above_min_OOB = (df['H_rec_above_min'] == 1) & (df['H_rec_in_bounds'] == 0)
            LCOEs = df.loc[H_rec_above_min_OOB, ['LCOE']]
            LCOE_min = LCOEs.min()
            LCOE_max = LCOEs.max()
            sizes = 1 + (LCOEs - LCOE_min) * (10 - 1) / (LCOE_max - LCOE_min)       # max size = 10, min size = 1
            pts = ax.scatter(df.loc[H_rec_above_min_OOB, [x_col_name]], df.loc[H_rec_above_min_OOB, [y_col_name]],\
                c='magenta', s=sizes * 5)      # Larger marker = Larger LCOE
        except:
            pass

        try:
            H_rec_above_min_in_bounds =  (df['H_rec_above_min'] == 1) & (df['H_rec_in_bounds'] == 1)
            LCOEs = df.loc[H_rec_above_min_in_bounds, ['LCOE']]
            LCOE_min = LCOEs.min()
            LCOE_max = LCOEs.max()
            sizes = 1 + (LCOEs - LCOE_min) * (10 - 1) / (LCOE_max - LCOE_min)       # max size = 10, min size = 1
            pts = ax.scatter(df.loc[H_rec_above_min_in_bounds, [x_col_name]], df.loc[H_rec_above_min_in_bounds, [y_col_name]],\
                c='#1f77b4', s=sizes * 5)
        except:
            pass

        ax.set_xlim(5, 45)
        ax.set_ylim(5, 45)
        ax.set_xlabel(x_axis_label)
        ax.set_ylabel(y_axis_label)
        return

    grid_dims = [2, 2]
    PlotSubplot(df, fig, grid_dims, 1, 'dT_chrg',   'dT_approach_charge_hx',    'dT_htHX',  'dT_approach_ht_disch_hx')
    PlotSubplot(df, fig, grid_dims, 2, 'dT_chrg',   'dT_approach_charge_hx',    'dT_ltHX',  'dT_approach_lt_disch_hx')
    PlotSubplot(df, fig, grid_dims, 3, 'dT_ltHX',   'dT_approach_lt_disch_hx',  'dT_htHX',  'dT_approach_ht_disch_hx')

    plt.show()

    pass


if __name__ == '__main__':
    terminal_output_file = 'C:/Users/mboyd/Documents/Project Docs/Gen3_Gas/Brayton Model/optimum_results_rev5_indep_dTs.xlsx'
    # PlotScatter(terminal_output_file)
    PlotBubble(terminal_output_file)
