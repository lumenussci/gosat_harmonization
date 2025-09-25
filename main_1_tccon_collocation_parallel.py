#!/usr/bin/env python
# coding: utf-8

# Code to collocate GOSAT soundings to TCCON soundings. Collocation requirements based on 
# Das et al. (2025) https://doi.org/10.1029/2024EA003935 and Chris O'Dell's TCCON collocations 
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


# define data paths -- PC paths
#oco2_dir  = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\oco2_B11.1\\'  # on ocomaster: /data11/OCO2/product/Lite/B11.1/LtCO2/
#gosat_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\gosat_B9\\'  # on ocomaster: /data6/GOSAT/product/Lite/B9/
#tccon_dir = 'C:\\Users\\hopki\\Projects\\gosat_oco2\\data\\tccon_mini\\'
#out_fn =  'C:\\Users\\hopki\\Projects\\gosat_oco2\\TCCON_COD_9-12-25.nc' # '/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_2hrs_300km_2degLat_3degLon_3-100soundings_closest.hdf'


# define ocomaster paths
oco2_dir = '/data11/OCO2/product/Lite/B11.1/LtCO2/'
gosat_dir = '/data6/GOSAT/product/Lite/B9/'
tccon_dir = '/home/laurel/data/tccon/'
out_fn = '/home/laurel/match_gosat_v9_tccon_ggg2020_time1_lat2.5_lon5.0_min15.nc'


# Initialize parallelization - nb_workers can be increased for faster runtime; selected 1/3 
# of OCO ocomaster's cpus to not hog all of the resources
pandarallel.initialize(nb_workers=int(os.cpu_count()/3), progress_bar=False, 
                       verbose=2)  # 0: no logs, 1: only warning logs, 2: all logs


# Define parameters for what is considered a collocation
qf_all = True  # If True, keep all data, otherwise only QF=0 data
start_date = 20090106  # earliset TCCON record is 1/6/2009
#end_date = 20110101

dtime_max = 1.0  # +/- time difference in hours  
dlat_max  = 1.25  # latitude difference in deg
dlon_max  = 2.5  # longitude difference in deg
dist_max = 20.0  # max distance in degrees 
#min_soundings_oco2 = 3  # min number of good quality oco2 soundings
min_soundings_tccon = 15  # min number of TCCON soundings to collocate 

verbose = False  # print extra statements; set to True to match idl outputs

# get GOSAT files
gosat_full_path = glob.glob(gosat_dir + '*.nc4')  # acos Lite files
gosat_files = [os.path.basename(gosat_file) for gosat_file in gosat_full_path]
gosat_files.sort()
print(f'{len(gosat_files)} GOSAT files found')
gosat_dates = [int('20' + f[11:17]) for f in gosat_files]  # extract date and prepend '20' for year
#gosat_files = [gosat_file for gosat_file, gosat_date in list(zip(gosat_files, gosat_dates)) if gosat_date >= start_date and gosat_date <= end_date]
gosat_files = [gosat_file for gosat_file, gosat_date in list(zip(gosat_files, gosat_dates)) if gosat_date >= start_date]

print(f'{len(gosat_files)} GOSAT files match TCCON observation period')


def match_tccon_to_gosat(gosat_sounding, tccon_single_day_df, dtime_max, dist_max, 
                         dlat_max, dlon_max, min_soundings_tccon): 
    # Chris' collocation requirements in tccon_acos_match.pro (lines 137 - 165): 

    gosat_lat = gosat_sounding.gosat_latitude
    gosat_lon = gosat_sounding.gosat_longitude  
    
    # only keep TCCON soundings that are within the allowed time difference
    gosat_sid = gosat_sounding.gosat_sounding_id #.astype(int)
    
    # Only keep TCCON soundings that are +/- (dtime) of GOSAT sounding
    gosat_min_time = gosat_sounding.timestamp + pd.Timedelta(hours=(-dtime_max))
    gosat_max_time = gosat_sounding.timestamp + pd.Timedelta(hours=(dtime_max))
    tccon_single_day_df = tccon_single_day_df[(tccon_single_day_df['datetime'] <= gosat_max_time) & (tccon_single_day_df['datetime'] >= gosat_min_time)]
    if tccon_single_day_df.shape[0] == 0:
        #print(f'No TCCON soundings within the time window for {str(gosat_sounding.gosat_sounding_id)}')
        return pd.DataFrame()
    # iterate over each unique TCCON site since difference sites have different requirements 
    list_of_collocations = []
    for site in tccon_single_day_df['tccon_site'].unique():
        tccon_site_df = tccon_single_day_df[tccon_single_day_df['tccon_site']==site]
        
        site_lat = tccon_site_df.iloc[0]['lat']
        site_lon = tccon_site_df.iloc[0]['lon']

        dlat = site_lat - gosat_lat
        dlon = utils.longitude_difference(site_lon, gosat_lon) 
        dist = utils.co_sphdist(site_lon, site_lat, gosat_lon, gosat_lat, degrees=True, alternate=False, approximate=True)
        
        # distance requirements based on TCCON sites
        if site == 'Caltech':
            # Land only and must be within [33.38, 34.27] lat and [-118.49, -117.55] lon 
            mask = np.logical_and.reduce((np.abs(dlat) <= 0.25,
                                          np.abs(dlon) <= 0.25,
                                          gosat_sounding.gosat_retrieval_surface_type == 1)) # Must be land 
        elif site == 'Armstrong' or site == 'Dryden' or site == 'Edwards':
            # must be within [34.68, 37.47] lat and [-127.88, -112.88] lon 
            mask = np.logical_and.reduce((np.abs(dlat) <= 0.5,
                                          np.abs(dlon) <=2.5))
        elif site =='Orléans':
            mask = np.logical_and.reduce((np.abs(dlat) <=0.75,
                                          np.abs(dlon) <= 2.5))
        elif site == 'Tsukuba':
            # must be within +/- 1.25 deg. lat and 2.5 deg. lon of TCCON site and must exclude 
            # Tokyo: [35.27,35.84,139.25,140.14] and Nagoya: [34.74,35.35,136.64,137.18]
            mask = np.logical_and.reduce((#np.abs(dlat) <= dlat_max,
                                          #np.abs(dlon) <= dlon_max,
                                          #~utils.in_box(gosat_lat, gosat_lon, [35.27,35.84,139.25,140.14]),
                                          #~utils.in_box(gosat_lat, gosat_lon, [34.74,35.35,136.64,137.18]),
                                          #dist <= dist_max
                                          np.abs(dlat) <= 1.25,
                                          np.abs(dlon) <= 0.5)) 
        elif site_lat <= -25.0:  # Southern hemisphere sites south of 25S
            # must be within +/- 5.0 deg. lat and +/- 30.0 deg. lon of TCCON site
            mask = np.logical_and.reduce((np.abs(dlat) <= 5.0,  #### Double for COD
                                          np.abs(dlon) <= 30.0))
        else:
            # all opther sites must be within +/- 1.25 deg. lat and 2.5 deg. lon of TCCON site
            mask = np.logical_and.reduce((np.abs(dlat) <= dlat_max,  
                                          np.abs(dlon) <= dlon_max))
        
        # apply mask
        #if ~np.any(mask):
        if ~mask:
            continue # mask is False - skip to next site   #none of the masks are True, go to next site 
        #tccon_site_df = tccon_site_df[mask]
        if tccon_site_df.shape[0] > min_soundings_tccon:
            print(f'{tccon_site_df.shape[0]} TCCON soundings collocated to GOSAT {str(gosat_sounding.gosat_sounding_id)}')
    
            # average TCCON 
            ave_tccon_df = utils.average_tccon_struct_df(tccon_site_df)
            ave_tccon_df['num_tccon_soundings'] = tccon_site_df.shape[0]
    
            # combine TCCON and GOSAT dataframes
            gosat_single_df = gosat_sounding.to_frame().T
            merged = gosat_single_df.reset_index(drop=True).join(ave_tccon_df.reset_index(drop=True))
            
            list_of_collocations.append(merged)
    
    if len(list_of_collocations) <= 0:
        return pd.DataFrame()
    else:
        collocated_tccon_df = pd.concat(list_of_collocations)     
        return collocated_tccon_df



tccon_df = utils.tccon_to_dataframe(tccon_dir)
print(f'Loaded TCCON data - {len(tccon_df['tccon_site'].unique())} sites: {list(tccon_df['tccon_site'].unique())}')

start_time = time.time()
list_of_collocations = []
for gosat_file in gosat_files:
    gosat_date = '20' + gosat_file[11:17]
    print(f'GOSAT date: {gosat_date}')

    # Check if there are any TCCON soundings for the given date
    if not tccon_df['match_date'].isin([gosat_date]).any():  
        print(f'No TCCON soundings for {gosat_date}')
        continue
      
    # Otherwise, subset TCCON data and open GOSAT file
    tccon_single_day_df = tccon_df[tccon_df['match_date']==str(gosat_date)]
    
    gosat = h5.File(os.path.join(gosat_dir, gosat_file),'r')
    gosat_df = utils.gosat_to_dataframe(gosat)
    # Filter to only keep QF=0
    if not qf_all:
        gosat_quality_mask = np.where(np.array(gosat['xco2_quality_flag']) == 0)[0]  # find good soundings (quality_flag=0)
        if len(gosat_quality_mask) <= 0:  # None of the soundings are good
            print(f'No good GOSAT soundings for ({gosat_date})')
            continue
        # If there are any good soundings
        gosat_df = gosat_df.iloc[gosat_quality_mask]  # only keep good soundings
   
    gosat_df = gosat_df.sort_values(by='gosat_sounding_id')
    gosat_dt_df = gosat_df[['gosat_year', 'gosat_month', 'gosat_day', 'gosat_hour', 'gosat_minute', 'gosat_second']].copy()
    gosat_dt_df.columns = ['year', 'month', 'day', 'hour', 'minute', 'second']
    gosat_df.loc[:,'timestamp'] = pd.to_datetime(gosat_dt_df[['year', 'month', 'day', 'hour', 'minute', 'second']])

    _collocated = gosat_df.parallel_apply(match_tccon_to_gosat, axis=1,
                           args = (tccon_single_day_df, dtime_max, 
                                   dist_max, dlat_max, dlon_max, min_soundings_tccon))
     
    # Remove empty dataframes (GOSAT soundings with no OCO collocations)
    non_empty_mask = _collocated.apply(lambda df: not df.empty)
    _collocated = _collocated[non_empty_mask].tolist()
    list_of_collocations.extend(_collocated)

end_time = time.time()
print('\nDONE PERFORMING TCCON COLLOCATIONS')
print(f'Time elapsed: {(end_time - start_time):.2f} seconds')


if len(list_of_collocations) > 0:
    full_collocation_df = pd.concat(list_of_collocations) 
    
    print(f'{full_collocation_df.shape[0]} total TCCON-GOSAT collocations')
    
    full_collocation_xr = xr.Dataset.from_dataframe(full_collocation_df)
    
    # Save TCCON collocations
    #type_dict = {df_col: {'dtype': str(df_dtype)} for df_col, df_dtype in zip(full_collocation_df.columns.to_list(), full_collocation_df.dtypes.to_list())}
    full_collocation_xr.to_netcdf(path=out_fn, mode='w', format='NETCDF4', engine='netcdf4') #, encoding=type_dict)
    print(f'Saved collocations to: {out_fn}')
else:
    print(f'0 total TCCON-GOSAT collocations')




