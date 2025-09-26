#!/usr/bin/env python
# coding: utf-8

# # Code to apply harmonization models to GOSAT data collocated to either TCCON or OCO-2
#
# # Laurel Hopkins Manella 9/18/25


import glob
import os
import datetime
import math
import numpy as np
import pandas as pd
import h5py as h5
import xarray as xr
import time
import utils
import xgboost as xgb
from xgboost import XGBRegressor
import glob
from multiprocessing import Pool


# define directories
#fn = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\match_gosat_v9_oco2_v11.1_20140906_20200630_time2_lat2_lon3_min3_dist300_qf=all_parallel_airmass_sza.nc'
fn = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\match_gosat_v9_tccon_ggg2020_time1_lat2.5_lon5.0_min15.nc'
model_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\models\\'

landH_fn = 'landH_qf=all_train2014-2017_eval2018.json'
landM_fn = 'landM_qf=all_train2014-2017_eval2018.json'
oceanH_fn = 'oceanH_qf=all_train2014-2017_eval2018.json'

out_fn = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\harmonized_gosat_v9_to_oco2_v11.1_gosat_oco_collocations_version1a.nc'


# Define if GOSAT is collocated to TCCON or OCO-2 - this matters when performing ak
# correction.
colloc_to_tccon = True  # Collocated to TCCON, otherwise OCO-2

# read collocated data
colloc_data_xr = xr.open_dataset(fn)
colloc_data = colloc_data_xr.to_dataframe()
colloc_data['gosat_sounding_id'] = colloc_data.index
colloc_data.shape

model_settings = {'landH': {'gain': 'H', 'surface_type': 1},
                'landM': {'gain': 'M', 'surface_type': 1},
                'oceanH':{'gain': 'H', 'surface_type': 0}}

# Read in the three models
model_landH = xgb.XGBRegressor()
model_landH.load_model(os.path.join(model_dir, landH_fn))

model_landM = xgb.XGBRegressor()
model_landM.load_model(os.path.join(model_dir, landM_fn))

model_oceanH = xgb.XGBRegressor()
model_oceanH.load_model(os.path.join(model_dir, oceanH_fn))

models = {'landH': model_landH,
          'landM': model_landM,
          'oceanH': model_oceanH,
         }

# define model features
# these features were determined to be the most informative (feature_importance_analysis.ipynb)

features = {'landM': ['gosat_retrieval_dws',
                      'gosat_retrieval_dust_height',
                      'gosat_retrieval_co2_grad_del',
                      'gosat_retrieval_aod_strataer',
                      'gosat_retrieval_albedo_sco2',
                      'gosat_retrieval_aod_sulfate',
                      'gosat_retrieval_dpfrac',
                      'gosat_retrieval_eof2_3',
                      'gosat_retrieval_eof1_3',
                      'gosat_retrieval_s32',
                      'gosat_retrieval_aod_ice',
                      'gosat_h2o_ratio',
                      'gosat_co2_ratio',
                      'gosat_retrieval_ice_height',
                      'gosat_retrieval_eof3_2',
                      'gosat_retrieval_psurf',
                      'gosat_retrieval_albedo_slope_wco2',
                      'gosat_retrieval_tcwv_uncertainty'],
            'landH': ['gosat_retrieval_dws',
                      'gosat_retrieval_co2_grad_del',
                      'gosat_retrieval_fs',
                      'gosat_retrieval_aod_fine',
                      'gosat_dp_abp',
                      'gosat_h2o_ratio',
                      'gosat_retrieval_albedo_wco2',
                      'gosat_retrieval_albedo_slope_sco2',
                      'gosat_retrieval_aod_strataer',
                      'gosat_retrieval_eof1_3',
                      'gosat_retrieval_albedo_slope_o2a',
                      'gosat_retrieval_water_height',
                      'gosat_retrieval_rms_rel_wco2',
                      'gosat_altitude_stddev',
                      'gosat_retrieval_h2o_scale',
                      'gosat_retrieval_dust_height',
                      'gosat_retrieval_aod_ice',
                      'gosat_retrieval_eof2_3',
                      'gosat_retrieval_deltaT',
                      'gosat_retrieval_offset_o2a_rel'],
            'oceanH': ['gosat_retrieval_albedo_slope_sco2',
                       'gosat_retrieval_co2_grad_del',
                       'gosat_retrieval_eof3_3',
                       'gosat_retrieval_albedo_slope_o2a',
                       'gosat_retrieval_aod_fine',
                       'gosat_dp_abp',
                       'gosat_snr_o2a',
                       'gosat_h2o_ratio',
                       'gosat_retrieval_aod_water',
                       'gosat_retrieval_albedo_sco2',
                       'gosat_retrieval_aod_strataer',
                       'gosat_retrieval_aod_dust',
                       'gosat_co2_ratio',
                       'gosat_retrieval_deltaT',
                       'gosat_retrieval_dws',
                       'gosat_retrieval_tcwv']}


corrected_list = []

colloc_data
# iterate over the three 'settings': landH, landM, oceanH
for setting in model_settings:
    _gain = model_settings[setting]['gain']
    _surface_type = model_settings[setting]['surface_type']
    # subset data to that in which model was trained
    subset_df = colloc_data[((colloc_data['gosat_retrieval_surface_type'] == _surface_type) & (colloc_data['gosat_gain'] == _gain))]
    if _surface_type == 0:
        print(f'ocean{_gain}: {subset_df.shape[0]}')
    else:
        print(f'land{_gain}: {subset_df.shape[0]}')
    x_eval = subset_df[features[setting]]
    # apply model - alpha = [0.16, 0.25, 0.5, 0.75, 0.84]
    predictions = models[setting].predict(x_eval)  # apply the model for the given setting
    # extract predictions
    subset_df = subset_df.assign(y_16 = predictions[:, 0])
    subset_df = subset_df.assign(y_25 = predictions[:, 1])
    subset_df = subset_df.assign(y_50 = predictions[:, 2])
    subset_df = subset_df.assign(y_75 = predictions[:, 3])
    subset_df = subset_df.assign(y_84 = predictions[:, 4])
    subset_df = subset_df.assign(gosat_xco2_harmonized = subset_df['gosat_retrieval_xco2_raw'] - subset_df['y_50'])
    # 84th - 16th percentile is 2 sigma (for a normal distribution)
    subset_df = subset_df.assign(gosat_ml_corrected_uncertainty = (subset_df['y_84'] - subset_df['y_16']) / 2)
    corrected_list += [subset_df]

full_corrected_df = pd.concat(corrected_list)
full_corrected_df['gosat_gain'] = full_corrected_df['gosat_gain'].astype(str)

print('\nDONE PERFORMING HARMONIZATION')
print(f'Harmonized {full_corrected_df.shape[0]} GOSAT soundings')


# Perform AK/prioir correction to enable comparison across different instruments
if colloc_to_tccon:
    full_corrected_df['tccon_xco2_ak_corrected_mean'], full_corrected_df['tccon_xco2_mean_ak_adj'] = utils.correct_tccon_ak(full_corrected_df, 'tccon_xco2_mean')
    full_corrected_df['tccon_xco2_ak_corrected_median'], full_corrected_df['tccon_xco2_median_ak_adj'] = utils.correct_tccon_ak(full_corrected_df, 'tccon_xco2_median')

else:
    # Correct for different CO2 priors - based on
    gosat_prior = full_corrected_df[full_corrected_df.columns[full_corrected_df.columns.str.startswith('gosat_co2_profile_apriori')]].values
    gosat_h = full_corrected_df[full_corrected_df.columns[full_corrected_df.columns.str.startswith('gosat_pressure_weight')]].values
    gosat_a = full_corrected_df[full_corrected_df.columns[full_corrected_df.columns.str.startswith('gosat_xco2_averaging_kernel')]].values
    oco2_prior = full_corrected_df[full_corrected_df.columns[full_corrected_df.columns.str.startswith('oco2_co2_profile_apriori')]].values
    oco2_h = full_corrected_df[full_corrected_df.columns[full_corrected_df.columns.str.startswith('oco2_pressure_weight')]].values
    oco2_a = full_corrected_df[full_corrected_df.columns[full_corrected_df.columns.str.startswith('oco2_xco2_averaging_kernel')]].values

    # CO2 prior and AK corrections
    # Former as in Eq A10 in [Wunch, ACP, 2010]. Also equivelant to Eq 3 in [Taylor, AMT, 2023].
    # The sum of the two terms yields the total correction to GOSAT XCO2 to harmonize it to OCO-2.
    true_profile = oco2_prior
    full_corrected_df.loc[:,['gosat_xco2_ak_adj_a']] = np.sum(gosat_h * (gosat_a - 1) * (gosat_prior - true_profile), axis=1)
    full_corrected_df.loc[:,['gosat_xco2_ak_adj_b']] = np.sum(gosat_h * (gosat_a - oco2_a) * (gosat_prior - true_profile), axis=1)

    # Apply correction
    full_corrected_df.loc[:,['gosat_xco2_harmonized_ak_corrected']] = full_corrected_df['gosat_xco2_harmonized'] + full_corrected_df['gosat_xco2_ak_adj_a'] + full_corrected_df['gosat_xco2_ak_adj_b']
    full_corrected_df.loc[:,['gosat_xco2_ak_corrected']] = full_corrected_df['gosat_xco2'] + full_corrected_df['gosat_xco2_ak_adj_a'] + full_corrected_df['gosat_xco2_ak_adj_b']


# Save harmonized collocations
full_corrected_xr = xr.Dataset.from_dataframe(full_corrected_df)

full_corrected_xr.to_netcdf(path=out_fn, mode='w', format='NETCDF4', engine='netcdf4') #, encoding=type_dict)
print(f'Saved harmonized collocations to: {out_fn}')
