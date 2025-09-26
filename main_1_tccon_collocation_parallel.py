#!/usr/bin/env python
# coding: utf-8

# Code to collocate GOSAT or OCO-2 soundings to an aggregation of TCCON soundings. Collocation requirements  
# based on Das et al. (2025) https://doi.org/10.1029/2024EA003935 and Chris O'Dell's TCCON collocations 
# as defined in tccon_acos_match.pro
#
# Laurel Hopkins Manella 9/1/2025


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
import datetime
from pandarallel import pandarallel


# define PC paths
#oco2_dir  = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\oco2_B11.1\\'  
#gosat_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\gosat_B9\\'  
#tccon_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\tccon_mini\\'
#out_fn =  'C:\\Users\\hopki\\Projects\\gosat_oco2\\tccon_collocation_test.nc'


# define ocomaster paths
oco2_dir = '/data11/OCO2/product/Lite/B11.1/LtCO2/'
gosat_dir = '/data6/GOSAT/product/Lite/B9/'
tccon_dir = '/home/laurel/data/tccon/'
out_fn = '/home/laurel/match_gosat_v9_tccon_ggg2020_time1_lat2.5_lon5.0_min15.nc'

# define if we are collocating to GOSAT, otherwise OCO-2
instrument = 'GOSAT'  # ['GOSAT' or 'OCO2']

# Initialize parallelization - nb_workers can be increased for faster runtime; selected 1/3 
# of OCO ocomaster's cpus to not hog all of the resources
pandarallel.initialize(nb_workers=int(os.cpu_count()/3), progress_bar=False, 
                       verbose=2)  # 0: no logs, 1: only warning logs, 2: all logs


# Define parameters for what is considered a collocation
qf_all = True  # If True, keep all data, otherwise only QF=0 data
start_date = 20090106  # earliset TCCON record is 1/6/2009
#end_date = 20110101

dtime_max = 1.0  # +/- time difference in hours -- CHANGE TO 2?? 
dlat_max  = 1.25  # latitude difference in deg -- CHANGE to 7.5??
dlon_max  = 2.5  # longitude difference in deg -- CHANGE TO 20???
#min_soundings_oco2 = 3  # min number of good quality oco2 soundings
min_soundings_tccon = 15  # min number of TCCON soundings to collocate 

verbose = False  # print extra statements; set to True to match idl outputs

# get GOSAT files
if instrument == 'GOSAT':
    instrmunet_full_path = glob.glob(gosat_dir + '*.nc4')
else:
    instrmunet_full_path = glob.glob(oco2_dir + '*.nc4')
instrument_files = [os.path.basename(instrument_file) for instrument_file in instrmunet_full_path]
instrument_files.sort()
print(f'{len(instrument_files)} {instrument} files found')
instrument_dates = [int('20' + f[11:17]) for f in instrument_files]  # extract date and prepend '20' for year
#instrument_files = [instrument_file for instrument_file, instrument_date in list(zip(instrument_files, instrument_dates)) if instrument_date >= start_date and instrument_date <= end_date]
instrument_files = [instrument_file for instrument_file, instrument_date in list(zip(instrument_files, instrument_dates)) if instrument_date >= start_date]

print(f'{len(instrument_files)} {instrument} files match TCCON observation period')


def match_tccon_to_instrument(instrument_sounding, tccon_single_day_df, dtime_max, 
                         dlat_max, dlon_max, min_soundings_tccon, instrument): 
    """Function responsible for perfomring collocation based on site-specific collocation
       criteria. Function takes in a single GOSAT/OCO-2 sounding and averages the TCCON 
       soundings that fall within the criteria."""
                           
    if instrument == 'GOSAT':
        instrument_lat = instrument_sounding.gosat_latitude
        instrument_lon = instrument_sounding.gosat_longitude
        instrument_sid = instrument_sounding.gosat_sounding_id
        instrument_surface_type = instrument_sounding.gosat_retrieval_surface_type
    else:
        instrument_lat = instrument_sounding.oco2_latitude
        instrument_lon = instrument_sounding.oco2_longitude  
        instrument_sid = instrument_sounding.oco2_sounding_id     
        instrument_surface_type = instrument_sounding.oco2_retrieval_surface_type

    
    # Only keep TCCON soundings that are +/- (dtime) of instrument sounding
    instrument_min_time = instrument_sounding.timestamp + pd.Timedelta(hours=(-dtime_max))
    instrument_max_time = instrument_sounding.timestamp + pd.Timedelta(hours=(dtime_max))
    tccon_single_day_df = tccon_single_day_df[(tccon_single_day_df['datetime'] <= instrument_max_time) & (tccon_single_day_df['datetime'] >= instrument_min_time)]
    if tccon_single_day_df.shape[0] == 0:
        return pd.DataFrame()
    # iterate over each unique TCCON site since difference sites have different requirements 
    list_of_collocations = []
    for site in tccon_single_day_df['tccon_site'].unique():
        tccon_site_df = tccon_single_day_df[tccon_single_day_df['tccon_site']==site]
        
        site_lat = tccon_site_df.iloc[0]['lat']
        site_lon = tccon_site_df.iloc[0]['lon']

        dlat = site_lat - instrument_lat
        dlon = utils.longitude_difference(site_lon, instrument_lon) 
        dist = utils.co_sphdist(site_lon, site_lat, instrument_lon, instrument_lat, degrees=True, 
                                alternate=False, approximate=True)

        # TCCON site-specific collocation requirements
        if site == 'Caltech':
            mask = np.logical_and.reduce((np.abs(dlat) <= 0.25,
                                          np.abs(dlon) <= 0.25,
                                          instrument_surface_type == 1)) # Must be land 
        elif site == 'Armstrong' or site == 'Dryden' or site == 'Edwards':
            mask = np.logical_and.reduce((np.abs(dlat) <= 0.5,
                                          np.abs(dlon) <=2.5))
        elif site =='Orléans':
            mask = np.logical_and.reduce((np.abs(dlat) <=0.75,
                                          np.abs(dlon) <= 2.5))
        elif site == 'Tsukuba':
            # can optionally exclude these two cities
            # Tokyo: [35.27,35.84,139.25,140.14] and Nagoya: [34.74,35.35,136.64,137.18]
            mask = np.logical_and.reduce((#~utils.in_box(instrument_lat, instrument_lon, [35.27,35.84,139.25,140.14]),
                                          #~utils.in_box(instrument_lat, instrument_lon, [34.74,35.35,136.64,137.18]),
                                          np.abs(dlat) <= 1.25,
                                          np.abs(dlon) <= 0.5)) 
        elif site_lat <= -25.0:  # Southern hemisphere sites south of 25S
            mask = np.logical_and.reduce((np.abs(dlat) <= 5.0,  
                                          np.abs(dlon) <= 30.0))
        else:
            mask = np.logical_and.reduce((np.abs(dlat) <= dlat_max,  
                                          np.abs(dlon) <= dlon_max))
        
        if ~mask:
            continue # if mask is False, skip to next site  
        if tccon_site_df.shape[0] > min_soundings_tccon:
            print(f'{tccon_site_df.shape[0]} TCCON soundings collocated to {instrument} \
                    {str(instrument_sid)}')
    
            # average TCCON 
            ave_tccon_df = utils.average_tccon_struct_df(tccon_site_df)
            ave_tccon_df['num_tccon_soundings'] = tccon_site_df.shape[0]
    
            # combine TCCON and instrument dataframes
            instrument_single_df = instrument_sounding.to_frame().T
            merged = instrument_single_df.reset_index(drop=True).join(ave_tccon_df.reset_index(drop=True))
            
            list_of_collocations.append(merged)
    
    if len(list_of_collocations) <= 0:
        return pd.DataFrame()
    else:
        collocated_tccon_df = pd.concat(list_of_collocations)     
        return collocated_tccon_df



tccon_df = utils.tccon_to_dataframe(tccon_dir)
print(f'Loaded TCCON data - {len(tccon_df['tccon_site'].unique())} sites: {list(tccon_df['tccon_site'].unique())}')

instrument_str = instrument.lower()
start_time = time.time()
list_of_collocations = []
for instrument_file in instrument_files:
    instrument_date = '20' + instrument_file[11:17]
    print(f'{instrument} date: {instrument_date}')

    # Check if there are any TCCON soundings for the given date
    if not tccon_df['match_date'].isin([instrument_date]).any():  
        print(f'No TCCON soundings for {instrument_date}')
        continue
      
    # Otherwise, subset TCCON data and open instrument file
    tccon_single_day_df = tccon_df[tccon_df['match_date']==str(instrument_date)]

    if instrument == 'GOSAT':
        instrument_h5 = h5.File(os.path.join(gosat_dir, instrument_file),'r')
        instrument_df = utils.gosat_to_dataframe(instrument_h5)
    else:
        instrument_h5 = h5.File(os.path.join(oco2_dir, instrument_file),'r')
        instrument_df = utils.oco_to_dataframe(instrument_h5)
        
    # Only keep QF=0 (good) soundings
    if not qf_all:
        instrument_df = instrument_df[instrument_df[f'{instrument_str}_xco2_quality_flag']==0]  
        #instrument_quality_mask = np.where(np.array(instrument_h5['xco2_quality_flag']) == 0)[0]  # find good soundings
        #if len(instrument_quality_mask) <= 0:  # None of the soundings are good
        if instrument_df.shape[0] == 0:  # None of the soundings are good
            print(f'No good {instrument} soundings for ({instrument_date})')
            continue
        # If there are any good soundings
        #instrument_df = instrument_df.iloc[instrument_quality_mask]  # only keep good soundings
   
    instrument_df = instrument_df.sort_values(by=f'{instrument_str}_sounding_id')
    instrument_dt_df = instrument_df[[f'{instrument_str}_year', f'{instrument_str}_month', 
                                      f'{instrument_str}_day', f'{instrument_str}_hour', 
                                      f'{instrument_str}_minute', f'{instrument_str}_second']].copy()
    instrument_dt_df.columns = ['year', 'month', 'day', 'hour', 'minute', 'second']
    instrument_df.loc[:,'timestamp'] = pd.to_datetime(instrument_dt_df[['year', 'month', 'day', 'hour', 'minute', 'second']])

    _collocated = instrument_df.parallel_apply(match_tccon_to_instrument, axis=1,
                           args = (tccon_single_day_df, dtime_max, dlat_max, 
                                   dlon_max, min_soundings_tccon, instrument))
     
    # Remove empty dataframes (GOSAT/OCO soundings with no TCCON collocations)
    non_empty_mask = _collocated.apply(lambda df: not df.empty)
    _collocated = _collocated[non_empty_mask].tolist()
    list_of_collocations.extend(_collocated)

end_time = time.time()
print('\nDONE PERFORMING TCCON COLLOCATIONS')
print(f'Time elapsed: {(end_time - start_time):.2f} seconds')


if len(list_of_collocations) > 0:
    full_collocation_df = pd.concat(list_of_collocations) 
    
    print(f'{full_collocation_df.shape[0]} total TCCON-{instrument} collocations')
    
    full_collocation_xr = xr.Dataset.from_dataframe(full_collocation_df)
    
    # Save TCCON collocations
    #type_dict = {df_col: {'dtype': str(df_dtype)} for df_col, df_dtype in zip(full_collocation_df.columns.to_list(), full_collocation_df.dtypes.to_list())}
    full_collocation_xr.to_netcdf(path=out_fn, mode='w', format='NETCDF4', engine='netcdf4') #, encoding=type_dict)
    print(f'Saved collocations to: {out_fn}')
else:
    print(f'0 total TCCON-{instrument} collocations')




