import numpy as np
import pandas as pd
from scipy.stats import linregress

from fffit.utils import values_scaled_to_real, values_real_to_scaled


def prepare_df(df_csv, AP):
    """Prepare a pandas dataframe for fitting GP models to AP data"""

    # Drop simulations that failed
    df_all = df_csv.dropna(subset=['uc_mean_distance'])
    # Add scaled_temperature
    df_all["scaled_temperature"] = values_real_to_scaled(
        df_all["temperature"].values, AP.temperature_bounds
    )

    params = df_all[list(AP.param_names)]
    scaled_params = values_real_to_scaled(params, AP.param_bounds)
    df_all[list(AP.param_names)] = scaled_params

    return df_all


