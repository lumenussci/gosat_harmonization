#!/usr/bin/env python
# coding: utf-8

# ###Python implementation of main_1_match_gosat_to_oco2.pro code on ocomaster (Tommy Taylor's IDL code): /home/ttaylor/analysis_utilities/tropical_iav/code/
# This code takes a set of GOSAT soundings and OCO-2 soundings and collocates the OCO-2 soundings with the GOSAT soundings based on conditions set below.

# Laurel Hopkins Manella 9/4/2025

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
from multiprocessing import Pool
from pandarallel import pandarallel


# define data paths -- PC paths
#oco2_dir  = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\oco2_B11.1\\'  # on ocomaster: /data11/OCO2/product/Lite/B11.1/LtCO2/
#gosat_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\gosat_B9\\'  # on ocomaster: /data6/GOSAT/product/Lite/B9/
#tccon_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\tccon\\'
#out_fn =  'C:\\Users\\hopki\\Projects\\gosat_oco2\\TEST.nc' # '/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_2hrs_300km_2degLat_3degLon_3-100soundings_closest.hdf'

# define ocomaster paths
oco2_dir = '/data11/OCO2/product/Lite/B11.1/LtCO2/'
gosat_dir = '/data6/GOSAT/product/Lite/B9/'
#out_fn = '/home/laurel/match_gosat_v9_oco2_v11.1_20140906_20200630_time2_lat2_lon3_min3_dist300_qf=all_updated_logic_ocoFiltered.nc'
out_fn = '/home/laurel/match_gosat_v9_oco2_v11.1_20140906_20200630_time2_lat2_lon3_min3_dist300_qf=all_parallel.nc'
#out_fn = '/home/laurel/match_gosat_v9_oco2_v11.1_20140906_20141001_time2_lat2_lon3_min3_dist300_qf=all_parallel.nc'

# Define parameters for what is considered a collocation
start_date = 20140906  # define date range for matching
end_date = 20200630

dtime_max = 2.0  # time difference in hours
dlat_max  = 2.0  # latitude difference in deg
dlon_max  = 3.0  # longitude difference in deg
dist_max = 300.0  # max distance in km
min_soundings_oco2 = 3  # min number of coc2 soundings per GOSAT sounding
#max_soundings_oco2 = 100  # max number of oco2 soundings per GOSAT sounding

# Initialize parallelization
# nb_workers can be increased for faster runtime; selected 1/3 of OCO ocomaster
# to not hog all of the resources
pandarallel.initialize(nb_workers=int(os.cpu_count()/3), progress_bar=False,
verbose=2)  # 0: no logs, 1: only warning logs, 2: all logs

# Collect GOSAT files
gosat_full_path = glob.glob(gosat_dir + '*.nc4')  # acos Lite files
gosat_files = [os.path.basename(gosat_file) for gosat_file in gosat_full_path]
gosat_files.sort()
print(f'{len(gosat_files)} GOSAT files found')
gosat_dates = [int('20' + f[11:17]) for f in gosat_files]  # extract date and prepend '20' for year
gosat_files = [gosat_file for gosat_file, gosat_date in list(zip(gosat_files, gosat_dates)) if gosat_date >= start_date and gosat_date <= end_date]
gosat_dates = [int('20' + f[11:17]) for f in gosat_files]  # extract dateis again from filtered files
print(f'{len(gosat_files)} GOSAT files match the observation period of OCO-2')


def load_oco_single_date(date, oco2_dir):
    date6 = ''.join(date.split())[2:8]
    oco_files = glob.glob(f'{oco2_dir}/*{date6}*')
    oco_files = [of for of in oco_files if of.endswith('.nc4')]  # only consider nc4 files

    if len(oco_files) == 0:  # if an oco file is NOT found
        print(f'No OCO Lite file found for the current date ({date6})')
        return pd.DataFrame()

    # If an oco file is found
    oco_current = h5.File(oco_files[-1],'r')
    oco_quality_mask = np.array(oco_current['xco2_quality_flag']) == 0  # True if sounding is good (quality_flag=0)
    n_oco_current = oco_quality_mask.sum()

    if n_oco_current > 0:  # if at least one good quality OCO sounding exists
        oco_current_df = utils.oco_to_dataframe(oco_current)
        oco_current_df = oco_current_df[o:wq
          co_quality_mask]  # only keep good soundings"""
        #oco_current_df = utils.oco_to_dataframe(oco_current)
        print(f'Loaded OCO date {date} with {oco_current_df.shape[0]} soundings')

        return oco_current_df

    else:
        print(f'No good soundings for OCO {date}')
        return pd.DataFrame()


def load_oco_files_in_parallel(dates_need, oco2_dir):
    # run load_oco_single_date in parallel for all three files
    with Pool(processes=3) as pool:
        _args = [(date, oco2_dir) for date in dates_need]
        results = pool.starmap(load_oco_single_date, _args)
    return pd.concat(results)


def perform_collocation(gosat_sounding, oco_df, oco_cols_to_ave, dtime_max,
                        dlat_max, dlon_max, dist_max, min_soundings_oco2):

    aid = gosat_sounding.gosat_sounding_id

    # Only keep OCO soundings that are +/- (dtime + 1hr) of GOSAT sounding
    gosat_min_time = gosat_sounding.timestamp + pd.Timedelta(hours=(-dtime_max - 1))
    gosat_max_time = gosat_sounding.timestamp + pd.Timedelta(hours=(dtime_max + 1))
    oco_soundings = oco_df[(oco_df['timestamp'] <= gosat_max_time) & (oco_df['timestamp'] >= gosat_min_time)]

    # Only collocate soundings with the same surface type
    surface_type_mask = oco_soundings.oco2_retrieval_surface_type == gosat_sounding.gosat_retrieval_surface_type
    if len(surface_type_mask) == 0:  # continue if there aren't any matches w/same surface type
        return pd.DataFrame()
    oco_soundings = oco_soundings[surface_type_mask]

    # Time, lat, lon, distance masks
    dtime = oco_soundings['timestamp'] - gosat_sounding.timestamp
    dlat  = oco_soundings.oco2_latitude.to_numpy() - gosat_sounding.gosat_latitude
    dlon  = utils.longitude_difference(oco_soundings.oco2_longitude.to_numpy(), gosat_sounding.gosat_longitude)
    ddist = utils.co_sphdist(oco_soundings.oco2_longitude.to_numpy(), oco_soundings.oco2_latitude.to_numpy(), gosat_sounding.gosat_longitude, gosat_sounding.gosat_latitude,
                       degrees=True, alternate=False, approximate=True) * 111.1  # distnace in km

    mask = np.logical_and.reduce((abs(dtime) <= datetime.timedelta(hours=dtime_max),
                              np.abs(dlat) <= dlat_max,
                              np.abs(dlon) <= dlon_max,
                              ddist <= dist_max))

    # Continue to next sounding if mask is empty
    if ~np.any(mask):
        return pd.DataFrame()
    # Otherwise, apply mask
    collocated_oco = oco_soundings[mask] # only keep valid collocations

    # Basic statistics over xco2
    xco2_mean = collocated_oco['oco2_xco2'].mean()
    xco2_std = collocated_oco['oco2_xco2'].std()

    # Filter out oco2 values <= 300 and >=900 and beyond 3 standard deviations from the mean
    collocated_oco = collocated_oco[(collocated_oco['oco2_xco2'] >= 300) & (collocated_oco['oco2_xco2'] <= 900)]
    collocated_oco[(collocated_oco['oco2_xco2'] < (xco2_mean + 3 * xco2_std)) & (collocated_oco['oco2_xco2'] > (xco2_mean - 3 * xco2_std))]

    # Recalculate statistics after filtering
    collocated_oco_mean = collocated_oco['oco2_xco2'].mean()
    collocated_oco_median = collocated_oco['oco2_xco2'].median()
    collocated_oco_std = collocated_oco['oco2_xco2'].std()

    num_matches = collocated_oco.shape[0]

    # Not enough oco2 soundings
    if num_matches < min_soundings_oco2:
        return pd.DataFrame()
    else:
        saved_matches = num_matches

    # Average over all collocated soundings and then append statistics
    ave_collocated_oco = utils.average_oco2_struct_df(collocated_oco, oco_cols_to_ave)
    ave_collocated_oco.loc[:,['oco2_xco2_mean']] = collocated_oco_mean
    ave_collocated_oco.loc[:,['oco2_xco2_median']] = collocated_oco_median
    ave_collocated_oco.loc[:,['oco2_xco2_std']] = collocated_oco_std

    # Merge collocations into a single dataframe and append to running list
    gosat_single_df = gosat_sounding.to_frame().T
    gosat_single_df = gosat_single_df.drop('timestamp', axis=1)
    merged = gosat_single_df.reset_index(drop=True).join(ave_collocated_oco.reset_index(drop=True))
    merged['num_total_collocations'] = num_matches
    merged['num_saved_collocations'] = saved_matches
    #merged['gosat_sounding_id'] = gosat_sounding.gosat_sounding_id
    merged.index = merged.gosat_sounding_id

    #list_of_collocations.append(merged)
    return merged

# OCO columns to keep for collocated data
oco_cols_to_ave = ['oco2_co2_profile_apriori0',
        'oco2_co2_profile_apriori1', 'oco2_co2_profile_apriori2',
        'oco2_co2_profile_apriori3', 'oco2_co2_profile_apriori4',
        'oco2_co2_profile_apriori5', 'oco2_co2_profile_apriori6',
        'oco2_co2_profile_apriori7', 'oco2_co2_profile_apriori8',
        'oco2_co2_profile_apriori9', 'oco2_co2_profile_apriori10',
        'oco2_co2_profile_apriori11', 'oco2_co2_profile_apriori12',
        'oco2_co2_profile_apriori13', 'oco2_co2_profile_apriori14',
        'oco2_co2_profile_apriori15', 'oco2_co2_profile_apriori16',
        'oco2_co2_profile_apriori17', 'oco2_co2_profile_apriori18',
        'oco2_co2_profile_apriori19', 'oco2_pressure_weight0',
        'oco2_pressure_weight1', 'oco2_pressure_weight2',
        'oco2_pressure_weight3', 'oco2_pressure_weight4',
        'oco2_pressure_weight5', 'oco2_pressure_weight6',
        'oco2_pressure_weight7', 'oco2_pressure_weight8',
        'oco2_pressure_weight9', 'oco2_pressure_weight10',
        'oco2_pressure_weight11', 'oco2_pressure_weight12',
        'oco2_pressure_weight13', 'oco2_pressure_weight14',
        'oco2_pressure_weight15', 'oco2_pressure_weight16',
        'oco2_pressure_weight17', 'oco2_pressure_weight18',
        'oco2_pressure_weight19', 'oco2_xco2_averaging_kernel0',
        'oco2_xco2_averaging_kernel1', 'oco2_xco2_averaging_kernel2',
        'oco2_xco2_averaging_kernel3', 'oco2_xco2_averaging_kernel4',
        'oco2_xco2_averaging_kernel5', 'oco2_xco2_averaging_kernel6',
        'oco2_xco2_averaging_kernel7', 'oco2_xco2_averaging_kernel8',
        'oco2_xco2_averaging_kernel9', 'oco2_xco2_averaging_kernel10',
        'oco2_xco2_averaging_kernel11', 'oco2_xco2_averaging_kernel12',
        'oco2_xco2_averaging_kernel13', 'oco2_xco2_averaging_kernel14',
        'oco2_xco2_averaging_kernel15', 'oco2_xco2_averaging_kernel16',
        'oco2_xco2_averaging_kernel17', 'oco2_xco2_averaging_kernel18',
        'oco2_xco2_averaging_kernel19', 'oco2_latitude',
        'oco2_longitude', 'oco2_xco2_quality_flag',
        'oco2_retrieval_surface_type', 'oco2_operation_mode', 'oco2_orbit']

# --------------------------------------------------- #
# Part 1. Read in the GOSAT data and the corresponding
# OCO data.
# --------------------------------------------------- #
list_of_collocations = []  # stores all collocated data
n = 0  # counter for collocations
c = np.pi/180  # for converting degrees to radians; needed for calculating distances

# Loop over GOSAT files
start_time = time.time()
for gosat_file in gosat_files:
    gosat = h5.File(os.path.join(gosat_dir, gosat_file),'r')
    gosat_date = '20' + gosat_file[11:17]
    gosat_df = utils.gosat_to_dataframe(gosat)

    # Only keep good soundings
    """gosat_quality_mask = np.array(gosat['xco2_quality_flag']) == 0  # find good soundings (quality_flag=0)
    if len(gosat_quality_mask) <= 0:  # None of the soundings are good
        print(f'No good GOSAT soundings for ({date})')
        continue
    gosat_df = gosat_df[gosat_quality_mask] """

    gosat_df = gosat_df.sort_values(by='gosat_sounding_id')
    gosat_dt_df = gosat_df[['gosat_year', 'gosat_month', 'gosat_day', 'gosat_hour', 'gosat_minute', 'gosat_second']].copy()
    gosat_dt_df.columns = ['year', 'month', 'day', 'hour', 'minute', 'second']
    gosat_df.loc[:,'timestamp'] = pd.to_datetime(gosat_dt_df[['year', 'month', 'day', 'hour', 'minute', 'second']])

    sid_min, sid_max = str(min(gosat_df['gosat_sounding_id'])), str(max(gosat_df['gosat_sounding_id']))

    # Convert sids to datetime objects
    date_min = datetime.datetime(int(sid_min[0:4]), int(sid_min[4:6]), int(sid_min[6:8]), 0, 0, 0)
    date_max = datetime.datetime(int(sid_max[0:4]), int(sid_max[4:6]), int(sid_max[6:8]), 0, 0, 0)
    # Always get one day before and after the GOSAT target date
    dates_need = np.unique(np.array([date_min - datetime.timedelta(days=1), date_min, date_max, date_max + datetime.timedelta(days=1)]))
    dates_need = [dt.strftime('%Y%m%d') for dt in dates_need]

    print(f'Current date of GOSAT file: {gosat_date}')
    print(f'Dates of OCO data that will be loaded: {dates_need}')

    # Load OCO data for matching dates
    oco_df = load_oco_files_in_parallel(dates_need, oco2_dir)

    if oco_df.empty:
        print(f'No matching OCO-2 soundings for ({dates_need})')
        continue  # no matching OCO soundings, move onto next GOSAT file

    # --------------------------------------------------- #
    # Part 2. Match the OCO data to the GOSAT sounding.
    # --------------------------------------------------- #
    print('Matching OCO2 to GOSAT')
    oco_dt_df = oco_df[['oco2_year', 'oco2_month', 'oco2_day', 'oco2_hour', 'oco2_minute', 'oco2_second']].copy()
    oco_dt_df.columns = ['year', 'month', 'day', 'hour', 'minute', 'second']
    _timestamp = pd.to_datetime(oco_dt_df[['year', 'month', 'day', 'hour', 'minute', 'second']])
    _timestamp.name = 'timestamp'
    oco_df = pd.concat([oco_df, _timestamp], axis=1)

    # Find gosat soundings that are between the min and max oco2 sounding ids
    min_oco_sid, max_oco_sid = min(oco_df['oco2_sounding_id'].to_numpy() // 100), max(oco_df['oco2_sounding_id'].to_numpy() // 100)
    gosat_df = gosat_df[(gosat_df['gosat_sounding_id'] >= min_oco_sid) & (gosat_df['gosat_sounding_id'] <= max_oco_sid)]  # cutting off times here?
    if gosat_df.shape[0] == 0:
        print(f'{gosat_date}: 0 soundings collocated')
        continue
    # Iterate over each sounding in a GOSAT file
    _collocated = gosat_df.parallel_apply(perform_collocation, axis=1,
                           args = (oco_df, oco_cols_to_ave, dtime_max, dlat_max,
                           dlon_max, dist_max, min_soundings_oco2))

    # Remove empty dataframes (GOSAT soundings with no OCO collocations)
    non_empty_mask = _collocated.apply(lambda df: not df.empty)
    _collocated = _collocated[non_empty_mask].tolist()

    list_of_collocations.extend(_collocated)
    n_oco_current = len(_collocated)
    n += len(_collocated)

    print(f'{gosat_date}: {n_oco_current} soundings collocated out of {n} total matches')

end_time = time.time()
# Save GOSAT-OCO2 collocations
# Merge all collocation dataframes into a single dataframe
full_collocation_df = pd.concat(list_of_collocations)

print('\nDONE PERFORMING COLLOCATIONS')
print(f'Time elapsed: {(end_time - start_time):.2f} seconds')

# convert to xarray in order to save file as a netcdf
full_collocation_xr = xr.Dataset.from_dataframe(full_collocation_df)

# get variable types
type_dict = {df_col: {'dtype': str(df_dtype)} for df_col, df_dtype in zip(full_collocation_df.columns.to_list(), full_collocation_df.dtypes.to_list())}
full_collocation_xr.to_netcdf(path=out_fn, mode='w', format='NETCDF4', engine='netcdf4', encoding=type_dict)
print(f'Saved collocations to: {out_fn}')

