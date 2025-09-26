#!/usr/bin/env python
# coding: utf-8

# Code to apply harmonization models to the GOSAT record

# Laurel Hopkins Manella 9/18/25


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


# define directories
gosat_dir = '/data6/GOSAT/product/Lite/B9/'
#gosat_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\gosat_test\\' 
#model_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\models\\'
model_dir = '/home/laurel/models/'

# define harmonization models
andH_fn = 'landH_qf=all_train2014-2017_eval2018.json'
landM_fn = 'landM_qf=all_train2014-2017_eval2018.json'
oceanH_fn = 'oceanH_qf=all_train2014-2017_eval2018.json'

#out_fn = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\harmonized_gosat_v9_to_oco2_v11.1.nc'
out_fn = '/home/laurel/harmonized_gosat_v9_to_oco2_v11.1_version1.nc


# model and parameter settings
model_settings = {'landH': {'gain': 'H', 'surface_type': 1},
                'landM': {'gain': 'M', 'surface_type': 1},
                'oceanH':{'gain': 'H', 'surface_type': 0}}


# Read in the three models
model_landH = xgb.XGBRegressor()
model_landH.load_model(os.path.join(model_dir, landH_fn))

model_landM = xgb.XGBRegressor()
model_landM.load_model(os.path.join(model_dir, landM_fn))

model_oceanH = xgb.XGBRegressor()
model_oceanH.load_model(os.path.join(model_dir, ocean_fn))

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


# Read in GOSAT files
gosat_full_path = glob.glob(gosat_dir + '*.nc4')  # acos Lite files
gosat_files = [os.path.basename(gosat_file) for gosat_file in gosat_full_path]
gosat_files.sort()
print(f'{len(gosat_files)} GOSAT files found')

# Begin harmonization 
corrected_list = []
for gosat_file in gosat_files:
    gosat = h5.File(os.path.join(gosat_dir, gosat_file),'r')
    gosat_date = '20' + gosat_file[11:17]
    gosat_df = utils.gosat_to_dataframe(gosat)
    print(f'Correcting GOSAT date: {gosat_date}')
    print(gosat_df.shape)
    
  # iterate over the three 'settings': landH, landM, oceanH
    for setting in model_settings:
        _gain = model_settings[setting]['gain']
        _surface_type = model_settings[setting]['surface_type']
        # subset data to that in which model was trained
        subset_df = gosat_df[((gosat_df['gosat_retrieval_surface_type'] == _surface_type) & (gosat_df['gosat_gain'] == _gain))]
        x_eval = subset_df[features[setting]]
        # apply model - alpha = [0.16, 0.25, 0.5, 0.75, 0.84]
        predictions = models[setting].predict(x_eval)  # apply the model for the given setting
        # extract predictions
        subset_df = subset_df.assign(y_16 = predictions[:, 0]) 
        subset_df = subset_df.assign(y_25 = predictions[:, 1]) 
        subset_df = subset_df.assign(y_50 = predictions[:, 2])
        subset_df = subset_df.assign(y_75 = predictions[:, 3]) 
        subset_df = subset_df.assign(y_84 = predictions[:, 4])
        subset_df = subset_df.assign(gosat_ml_corrected = subset_df['gosat_retrieval_xco2_raw'] - subset_df['y_50'])
        # 16th - 84th percentile is 2 sigma (for a normal distribution)
        subset_df = subset_df.assign(gosat_xco2_harmonized_uncertainty = (subset_df['y_84'] - subset_df['y_16']) / 2) 
        corrected_list += [subset_df]
full_corrected_df = pd.concat(corrected_list) 
full_corrected_df['gosat_gain'] = full_corrected_df['gosat_gain'].astype(str) 

print('\nDONE PERFORMING HARMONIZATION')
print(f'Harmonized {full_corrected_df.shape[0]} GOSAT soundings')

full_corrected_xr = xr.Dataset.from_dataframe(full_corrected_df)

# Save TCCON collocations
full_corrected_xr.to_netcdf(path=out_fn, mode='w', format='NETCDF4', engine='netcdf4') #, encoding=type_dict)
print(f'Saved corrected GOSAT to: {out_fn}')

