import os
import sys
import gpflow
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.patches import Circle, RegularPolygon
from matplotlib.path import Path
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
from matplotlib.spines import Spine
from matplotlib.transforms import Affine2D
import matplotlib.patches as mpatch

from fffit.pareto import find_pareto_set, is_pareto_efficient
sys.path.append("../")

from utils.ap import AP
from utils.prepare_samples import prepare_df

############################# QUANTITIES TO EDIT #############################
##############################################################################

iternum = 4

##############################################################################
##############################################################################

temperatures = [10, 78, 298]

csv_path = "/scratch365/bbefort/ap-fffit/ap-fffit/analysis/csv/"
in_csv_names = [
    "uc-lattice-iter" + str(i) + "-results.csv" for i in range(1, iternum+1)
]
out_csv_name = "uc-lattice-final-params.csv"

# Read files
df_csvs = [
    pd.read_csv(csv_path + in_csv_name, index_col=0)
    for in_csv_name in in_csv_names
]
df_csv = pd.concat(df_csvs)
df_all = prepare_df(df_csv, AP)


def main():

    ###########################################################
    ####################   Find Best Points    ###################
    ###########################################################
    
    #Compare all results from each iteration to hand tuned results
    df1 = df_all.loc[(df_all["temperature"] == 10) & (df_all["uc_mean_distance"] <= 0.156)].filter(list(AP.param_names))
    df2 = df_all.loc[(df_all["temperature"] == 10) & (df_all["lattice_ape"] <= 0.94706632)].filter(list(AP.param_names))
    df3 = df_all.loc[(df_all["temperature"] == 78) & (df_all["lattice_ape"] <= 1.62047081)].filter(list(AP.param_names))
    df4 = df_all.loc[(df_all["temperature"] == 298) & (df_all["lattice_ape"] <= 1.69615265)].filter(list(AP.param_names))
    df5 = df_all.loc[(df_all["temperature"] == 10) & (df_all["abs(HB3-HB4)"] <= 0.001)].filter(list(AP.param_names))
    df6 = df_all.loc[(df_all["temperature"] == 10) & (df_all["abs(HA3-HA4)"] <= 0.3)].filter(list(AP.param_names))
    overall_best = pd.merge(df1, pd.merge(df2, pd.merge(df3,pd.merge(df4,pd.merge(df5,df6)))))
    print(f"Intersection is {len(overall_best)} for uc_mean_distance, lattice_ape at all three temperatures, and HBond & HAngle Symmetry that is less than or equal to GT's best hand-tuned parameter set results")
    overall_best = pd.merge(df_all, overall_best[list(AP.param_names)], left_on=list(AP.param_names), right_on=list(AP.param_names))
    # Save the final pareto parameters
    overall_best.to_csv(csv_path + out_csv_name)
    
    #Calculate Lattice MAPE
    overall_best_group = overall_best.groupby(by=['sigma_Cl']).mean()
    overall_best_MAPE = overall_best_group['lattice_ape']
    overall_best_MAPE = np.asarray(overall_best_MAPE)
    
    #Drop NAN
    overall_best = overall_best.dropna()
    
    #Sort by sigma Cl - this is needed because the groupby function sorts
    overall_best = overall_best.sort_values(by='sigma_Cl')
    
    #Add Lattice MAPE to dataframe
    overall_best['Lattice_MAPE'] = overall_best_MAPE

    ###########################################################
    ##################   Find Pareto Points    ##################
    ###########################################################

    final_results = overall_best
    final_results = final_results.reset_index()
    
    #Pareto Metrics
    costs1 = final_results[['uc_mean_distance','Lattice_MAPE']]
    costs1 = costs1.to_numpy()
    
    result1 = is_pareto_efficient(costs1)
    costs1_list = costs1.tolist()
    paretoPoints=[]
    dominatedPoints=[]
    for i in range(len(costs1_list)):
        if result1[i] == True:
            paretoPoints.append(costs1_list[i])
        else:
            dominatedPoints.append(costs1_list[i])
    paretoPoints_efficient1 = np.array(paretoPoints)
    dominatedPoints_efficient1 = np.array(dominatedPoints)
    # stop = timeit.default_timer()
    print("Number of Pareto Optimal Points", len(paretoPoints_efficient1))
    print("Pareto Optimal Points:",paretoPoints_efficient1)
    
    final_results['Pareto'] = result1
    final_results_pareto = final_results[final_results['Pareto']==True]
    
    final_results_pareto.to_csv('AP_Final_2_Pareto.csv')


if __name__ == "__main__":
    main()


