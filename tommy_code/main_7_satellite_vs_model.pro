;;----- This procedure borrowed heavily from CWO's /home/codell/idlprogs/aos/models/OCO/l2/b9_acos/plot_acos_vs_models_b9.pro
;;----- Modifications by T.E.T Autumn, 2020 to adapt and generate analysis/plots for the GOSAT v9 data paper.
;;----- New version <main_satellite_vs_model.pro> developed for the OCO ROSES Tropical IAV work to harmonize GOSAT v9 and OCO-2 v10.
;;-----
;;----- Keep in mind that the matchup of models to OCO uses the ten second averages. So the "sounding densities" will be 
;;----- roughly order 10sec*3Hz*8sound/frame=240 times smaller than native. So instead of approx 100k good quality flagged soundings per day
;;----- there will be more like 400 matchups per day. Translates to 150k matchups per year.

;;----- Select data/plots to generate.
GridMap=0
GridHovmoller=1
PlotHistogram=0

;;----- Set the sensor and version
IF (1) THEN BEGIN
   ;VerName='gosat_v9a'
   VerName='gosat_v9c3'
   mode_list=['landH','landM','oceanH']
   ;mode_list=['landH']
   ;mode_list=['oceanH']
ENDIF
IF (0) THEN BEGIN
   ;VerName='oco2_v9'
   VerName='oco2_v10'
   mode_list=['land','oceanG']
ENDIF

;;----- Set latitude range for analysis.
;;----- Tropical IAV focuses on Land between -30/+30 latitude.
;;----- Sean will assimilate all latitudes, but for analysis it
;;----- is useful to highlight tropical land only.
IF (1) THEN BEGIN
   LatText='Global'
   LatRange=[-90.,90.]
   LonRange=[-180.,180.]
   DeltaLat_1=2.5 & DeltaLon_1=5.0
   BinMinCount=2
ENDIF
IF (0) THEN BEGIN
   LatText='Tropics'
   LatRange=[-30.,30.]
   LonRange=[-180.,180.]
   DeltaLat_1=2.5 & DeltaLon_1=5.0
   BinMinCount=2
ENDIF
IF (0) THEN BEGIN
   LatText='NTA'
   LatRange=[0.,20.]
   LonRange=[-20.,45.]
   DeltaLat_1=1. & DeltaLon_1=2.0
   BinMinCount=2
ENDIF
IF (0) THEN BEGIN
   LatText='Africa'
   LatRange=[-35.,35.] 
   LonRange=[-20.,50.]
   DeltaLat_1=1. & DeltaLon_1=2.0
   BinMinCount=2
ENDIF
IF (0) THEN BEGIN
   LatText='SouthAmerica'
   LatRange=[-65.,15.] 
   LonRange=[-85.,-35.]
   DeltaLat_1=1. & DeltaLon_1=2.0
   BinMinCount=2
ENDIF
IF (0) THEN BEGIN
   LatText='CONUS'
   LatRange=[25.,45.] 
   LonRange=[-130.,-65.]
   DeltaLat_1=1. & DeltaLon_1=2.0
   BinMinCount=2
ENDIF
BinSizeString_1=STRING(DeltaLat_1,FORMAT='(I2)')+' lat x '+STRING(DeltaLon_1,FORMAT='(I2)')+' lon'

;LatText='lat_'+STRTRIM(STRING(LatRange[0],FORMAT='(I3)'),2)+'_'+STRTRIM(STRING(LatRange[1],FORMAT='(I2)'),2)

;;----- Set the list of season names for analysis.
IF (GridMap OR PlotHistogram) THEN season_names = ['MAM','JJA','SON','DJF']
IF (GridHovmoller) THEN season_names = ['Annual']

;;----- Set start and end dates for analysis.
;;----- When working with the GOSAT/OCO-2 overlapping record,
;;----- it is useful to truncate the GOSAT data to start of OCO-2.
;time_range_oco2=1 & time_range_gosat=0
;IF (time_range_oco2)  THEN StartYear=2014 & StartMonth=9 & StartDay=6
;IF (time_range_gosat) THEN StartYear=2014 & StartMonth=9 & StartDay=6

;;----- Model agreement criteria.
max_model_std = 1.0
min_nmodels = 2
max_model_maxdiff = 1.5

;;----- Set suffix values.
;;----- Disabled if set to empty string.
bc_suffix='_bc'
ak_suffix='_ak'

;;-----
;;----- For GOSAT v9 processing.
;;-----
IF (STRMID(VerName,0,5) EQ 'gosat') THEN BEGIN

   IF (STRMID(VerName,6,3) EQ 'v9a') THEN $
      satellite_file = '/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9/acos_b9_lite_20090420_20200630_small.sav'
   IF (STRMID(VerName,6,3) EQ 'v9b') THEN $
      satellite_file = '/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9b/gosat_v9b_lite_20090420_20200630_filtered_superlite.sav'
   ;IF (STRMID(VerName,6,3) EQ 'v9c') THEN $
   ;   satellite_file = '/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c_r20240326/gosat_v9c_r20240326_lite_20090420_20200630_filtered_xmedium.sav'
   IF (STRMID(VerName,6,4) EQ 'v9c2') THEN $
      satellite_file = '/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c2_r20240329/gosat_v9c2_r20240329_lite_20090420_20200630_filtered_xmedium.sav'
   IF (STRMID(VerName,6,4) EQ 'v9c3') THEN $
      satellite_file = '/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c3_r20240422/gosat_v9c3_r20240422_lite_20090420_20200630_filtered_xmedium.sav'

   ;;-----
   ;;----- Model file used for GOSAT v9 ESSD publication.
   ;;-----
   ;;----- Model file=gosat_b9acos_4models20191119_ak_NO_OFFSET_REMOVED_stats.h5
   ;;----- StartDate=2009-04-20
   ;;----- EndDate=2019-03-28
   ;;----- Model_Names= CT_2017+NRTv2 Jena_s04oc-v4 MACC_2018-V2 UnivEd_v4.0
   ;;model_file='/home/codell/GOSAT/acos_results/b9_acos/models/b9test/results/gosat_b9acos_4models20191119'+ak_suffix+'_NO_OFFSET_REMOVED_stats.h5'

   ;;-----
   ;;----- Model file used Tropical IAV harmonization project.
   ;;-----
   ;;----- Model file=gosat_b9acos_4models20200820_ak_NO_OFFSET_REMOVED_stats.h5
   ;;----- StartDate=2009-04-20
   ;;----- EndDate=2020-05-25
   ;;----- Model_Names= Jena_s10oc-v2020_2009prefilled UnivEd_v4.0a MACC_v20r1 CT_2019
   ;model_file='/home/codell/GOSAT/acos_results/b9_acos/models/b9/results/gosat_b9acos_4models20200820_NO_OFFSET_REMOVED_stats.h5'

   ;;-----
   ;;----- Model file used Tropical IAV harmonization project.
   ;;-----
   ;;----- Model file=gosat_b9acos_4models_20210322.h5
   ;;----- StartDate=2009-04-20
   ;;----- EndDate=2020-05-25
   ;;----- Model_Names= Jena_s10oc-v2020_2009prefilled UnivEd_v4.0a MACC_v20r1 CT_2019
   ;model_file='/home/codell/GOSAT/acos_results/b9_acos/models/b9/results/gosat_b9acos_4models_20210322.h5'

   ;;-----
   ;;----- Model file used Tropical IAV harmonization project.
   ;;-----
   ;;----- Model file=gosat_b9acos_4models_20220223.h5
   ;;----- StartDate=2009-04-20
   ;;----- EndDate=2020-06-30
   ;;----- Model_Names= Jena_s10oc-v2021 UnivEd_v5 MACC_v20r2 CT_2019B+NRT2021-3
   model_file='/home/codell/GOSAT/acos_results/b9_acos/models/b9/results/gosat_b9acos_4models_20220223.h5'

   OutDir = '/home/ttaylor/analysis_utilities/tropical_iav/plots/satellite_vs_models/'+VerName+'/' & FILE_MKDIR,OutDir
   Pdf_File_Root=OutDir+VerName+'_vs_4models'

   ;;----- Set select data range using SID format.
   ;;----- This is because different model files have different date ranges.
   ;;----- So this hard cutoff allows for consistency in different versions.
   ;;DateRange=[20140906000000ULL, 20181231000000ULL]
   ;DateRange=[20140906000000ULL, 20191231000000ULL]
   ;DateRange=[20140906000000ULL, 20200630000000ULL]
   ;DateRange=[20180101000000ULL, 20191231000000ULL]
   IF (GridHovmoller) THEN DateRange=[20090420000000ULL, 20200630000000ULL]
   IF (GridMap) THEN BEGIN
      DateRange=[20091201000000ULL, 20191130000000ULL] ;; For DJF/MAM/JJA/SON analysis
      ;DateRange=[20091201000000ULL, 20101130000000ULL] ;; 2010 DJF/MAM/JJA/SON analysis
      ;DateRange=[20181201000000ULL, 20191130000000ULL] ;; 2019 DJF/MAM/JJA/SON analysis
      ;STOP
   ENDIF
   DateRangeString=STRMID(NUM2STR(DateRange[0],0),0,8)+'-'+STRMID(NUM2STR(DateRange[1],0),0,8)
   PRINT,'DateRangeString=',DateRangeString
   STOP
ENDIF 

;;-----
;;----- For OCO-2 v9 processing.
;;-----
IF (VerName EQ 'oco2_v9') THEN BEGIN
   ;satellite_file = '/home/codell/OCO2_results/b9/models/10sec/dave_files/OCO2_b91_10sec_GOOD_r23_newak.nc4'
   satellite_file = '/home/codell/OCO2_results/b9/models/10sec/dave_files/OCO2_b91_10sec_GOOD_r27.nc4'

   ;;-----
   ;;----- Model file used Tropical IAV harmonization project.
   ;;-----
   ;;----- Model file=oco2_b9_10sv27_4models_20230802.h5
   ;;----- StartDate=2014-09-06
   ;;----- EndDate=2020-01-21
   ;;----- Model_Names= CT_2019B+NRT2021-3 Jena_s10oc-v2021   MACC_v20r2         UnivEd_v5
   ;model_file='/home/codell/OCO2_results/b9/models/10sec/results/oco2_b9_10sv23_4modelRef_20200520_NO_OFFSET_REMOVED_stats.h5'
   ;model_file='/home/codell/OCO2_results/b9/models/10sec/results/oco2_b9_10sv27_4models_20220913.h5'
   model_file='/home/codell/OCO2_results/b9/models/10sec/results/oco2_b9_10sv27_4models_20230802.h5

   OutDir = '/home/ttaylor/analysis_utilities/tropical_iav/plots/'+SubDir+'/'+VerName+'/' & FILE_MKDIR,OutDir
   Pdf_File_Root=OutDir+VerName+'_vs_4models'

   ;;----- Set select data range using SID format.
   ;;----- This is because different model files have different date ranges.
   ;;----- So this hard cutoff allows for consistency in different versions.
   IF (GridHovmoller) THEN DateRange=[20140906000000ULL, 20200121000000ULL]
   IF (GridMap) THEN DateRange=[20141201000000ULL, 20191130000000ULL] ;; For DJF/MAM/JJA/SON analysis
ENDIF 

;;-----
;;----- For OCO-2 v10 processing.
;;-----
IF (VerName EQ 'oco2_v10') THEN BEGIN
   satellite_file = '/home/codell/OCO2_results/b10_tests/models/10sec/dave_files/OCO2_b10c_10sec_GOOD_r6.nc4'

   ;;-----
   ;;----- Model file used for Tropical IAV harmonization project.
   ;;-----
   ;;----- Model file=oco2_b10_10sv6_4models_20220219.h5
   ;;----- StartDate=2014-09-06
   ;;----- EndDate=2020-12-31
   ;;----- Model_Names= CT_2019B+NRT2021-3 Jena_s10oc-v2021 MACC_v20r2 UnivEd_v5
   ;;-----
   model_file='/home/codell/OCO2_results/b10_tests/models/10sec/results/oco2_b10_10sv6_4models_20220219.h5'

   OutDir = '/home/ttaylor/analysis_utilities/tropical_iav/plots/'+SubDir+'/'+VerName+'/' & FILE_MKDIR,OutDir
   Pdf_File_Root=OutDir+VerName+'_vs_4models'

   ;;----- Set select data range using SID format.
   ;;----- This is because different model files have different date ranges.
   ;;----- So this hard cutoff allows for consistency in different versions.
   ;DateRange=[20180101000000ULL, 20191231000000ULL]
   DateRange=[20140906000000ULL, 20190629000000ULL]
   IF (GridHovmoller) THEN DateRange=[20140906000000ULL, 201912311000000ULL]
   IF (GridMap) THEN DateRange=[20141201000000ULL, 20191130000000ULL] ;; For DJF/MAM/JJA/SON analysis

ENDIF 

;;----- Build output file root based on settings.
IF (bc_suffix) NE '' THEN pdf_file_root += bc_suffix
IF (ak_suffix) NE '' THEN pdf_file_root += ak_suffix
PRINT,'pdf_file_root=',pdf_file_root
;;STOP

;;-----
;;----- Read the model data H5 file.
;;-----
IF (N_ELEMENTS(loaded_model_file) EQ 0) THEN loaded_model_file=''
IF (loaded_model_file NE model_file) THEN BEGIN
   PRINT, 'Reading ' + model_file
   ;models = READ_H5_FILE(model_file, READ=['id','StatsAK/median','StatsAK/std','StatsAK/nmodels','StatsAK/maxdiff']) ;, 'lat','lon'], /cond)
   models = READ_H5_FILE(model_file)
   model_names = READ_H5_FIELD(model_file, 'model_names')
   loaded_model_file=model_file

   ID_models=SC(models.id)
   jd_models = ACOS_ID_TO_JD(ID_models)
   CALDAT, jd_models, month, day, year, hour, min, sec
   DOY_models = DATE2DOY(year, month, day)+0.
   Stats,DOY_models
   N=N_ELEMENTS(Year)
   ModelStartDateString=STRING(Year[0],FORMAT='(I4)')+'-'+STRING(month[0],FORMAT='(I2.2)')+'-'+STRING(Day[0],FORMAT='(I2.2)')
   ModelEndDateString=STRING(Year[N-1],FORMAT='(I4)')+'-'+STRING(month[N-1],FORMAT='(I2.2)')+'-'+STRING(Day[N-1],FORMAT='(I2.2)')
   PRINT,'VerName=',VerName
   PRINT,'Model file=',STRTRIM(FILE_BASENAME(model_file),2)
   PRINT,'ModelStartDate=',ModelStartDateString & PRINT,'ModelEndDate=',ModelEndDateString
   PRINT,'Model_Names=',Model_Names
   IF (1) THEN BEGIN
      STATS,models.statsak.nmodels
      STATS,models.statsak.maxdiff
      STATS,models.statsak.median
      STATS,models.statsak.std
   ENDIF

   ;;----- Do some plotting of model fields.
   ;;WIS_IMAGE_TET,models.nmodels,models.lat,models.lon,div=3
   ;STOP
ENDIF
;;STOP 

;;-----
;;----- Load the ACOS GOSAT data.
;;-----
IF (N_ELEMENTS(loaded_satellite_file) EQ 0) THEN loaded_satellite_file=''
IF (STRMID(VerName,0,5) EQ 'gosat') AND (loaded_satellite_file NE satellite_file) THEN BEGIN
   PRINT, 'Loading ' + satellite_file
   RESTORE, satellite_file
   so = sort(lite.sounding_id)
   lite=lite[so]
   rd = remove_duplicates(lite.sounding_id)
   lite=lite[rd] & rd=!NULL
   all_id = lite.sounding_id
   all_lat = lite.latitude
   all_lon = lite.longitude

   ;;----- Select the XCO2 from Lite file based on the version.
   IF (STRMID(VerName,6,3) EQ 'v9a') THEN all_xco2_bc = lite.xco2
   IF (STRMID(VerName,6,3) EQ 'v9b') THEN BEGIN
      all_xco2_bc = lite.xco2_bc_to_oco2_v10
      xco2_adjust=lite.xco2_bc_to_oco2_v10-lite.xco2
      ;STOP
   ENDIF
   IF (STRMID(VerName,6,3) EQ 'v9c') THEN BEGIN
      all_xco2_bc = lite.xco2_harmonized_to_oco2_v11_1
      all_xco2_adjust=lite.xco2_additive_adjustment_to_gosat_v9
      ;STOP
   ENDIF
   IF (VerName EQ 'gosat_v9a') THEN all_xco2_adjust=MAKE_ARRAY(N_ELEMENTS(all_id),/FLOAT,VALUE=0.0)

   all_xco2_uncert = lite.xco2_uncertainty
   all_xco2 = lite.retrieval.xco2_raw
   all_xco2_ap = lite.xco2_apriori
   all_gain = lite.sounding.gain eq 'H'
   all_land = lite.retrieval.surface_type eq 1 
   all_qf_mask = lite.xco2_quality_flag eq 0

   ;;----- Build "basic_mask"
   basic_mask = (lite.xco2_uncertainty GE 0.1) AND (lite.retrieval.aod_total LE 20.)
   lat_mask = INRANGE(all_lat,LatRange)
   lon_mask = INRANGE(all_lon,LonRange)
   date_mask=INRANGE(all_id,DateRange)

   ym = long(ulong64(all_id)/100000000)
   all_year = long(ym/100)
   all_month = fix((ym-all_year*100))
   all_jd=ACOS_ID_TO_JD(all_id)

   loaded_satellite_file = satellite_file

   MGAIN=0b
   HGAIN=1b
   ;STOP
ENDIF
;STOP

;;-----
;;----- Load the OCO-2 data.
;;-----
IF (N_ELEMENTS(loaded_satellite_file) EQ 0) THEN loaded_satellite_file=''
IF (STRMID(VerName,0,4) EQ 'oco2') AND (loaded_satellite_file NE satellite_file) THEN BEGIN
   PRINT, 'Loading ' + satellite_file

   lite = READ_H5_FILE(satellite_file, read=['sounding_id','latitude','longitude', $
                                             'xco2','xco2_quality_flag','xco2_uncertainty','xco2_raw',$     
                                             'N_total_shots','surface_type','operation_mode'], /cond)

   PRINT, 'Matching Soundings IDs'
   PRINT,'Number of satellite soundings prior to match =',N_ELEMENTS(Lite)
   PRINT,'Number of model soundings prior to match =',N_ELEMENTS(Models)
   wm = CMSET_OP(models.id, 'AND', ULONG64(lite.sounding_id), /index)
   wa = CMSET_OP(ULONG64(lite.sounding_id), 'AND', models.id, /index)
   IF (wm[0] eq  -1) then begin
      wm = CMSET_OP(models.id/10, 'AND', ULONG64(lite.sounding_id/10), /index)
      wa = CMSET_OP(ULONG64(lite.sounding_id/10), 'AND', models.id/10, /index)
   ENDIF   
   models=models[wm]
   lite=lite[wa]
   PRINT,'Number of satellite soundings after match =',N_ELEMENTS(Lite)
   PRINT,'Number of model soundings after match =',N_ELEMENTS(Models)
   ;STOP
   
   ;;----- Check the time span of the satellite data.
   all_id = lite.sounding_id
   ID_sat=SC(all_id)
   jd_sat = ACOS_ID_TO_JD(ID_sat)
   CALDAT, jd_sat, month, day, year, hour, min, sec
   N=N_ELEMENTS(Year)
   SatStartDateString=STRING(Year[0],FORMAT='(I4)')+'-'+STRING(month[0],FORMAT='(I2.2)')+'-'+STRING(Day[0],FORMAT='(I2.2)')
   SatEndDateString=STRING(Year[N-1],FORMAT='(I4)')+'-'+STRING(month[N-1],FORMAT='(I2.2)')+'-'+STRING(Day[N-1],FORMAT='(I2.2)')
   PRINT,'SatStartDate=',SatStartDateString & PRINT,'SatEndDate=',SatEndDateString
   ID_sat=!NULL & jd_sat=!NULL & month=!NULL & day=!NULL & year=!NULL & hour=!NULL & min=!NULL & sec=!NULL
   ;STOP

   ;;----- Extract satellite data fields.   
   all_lat = lite.latitude
   all_lon = lite.longitude
   all_xco2_bc = lite.xco2
   all_xco2_uncert = lite.xco2_uncertainty
   all_xco2 = lite.xco2_raw
   ;all_xco2_ap = lite.xco2_apriori
   all_n_total_shots = lite.n_total_shots
   last_digit = FIX(lite.sounding_id MOD 10)
   all_qf_mask = lite.xco2_quality_flag eq 0
   basic_mask = (lite.xco2_uncertainty GE 0.1) AND (lite.N_total_shots GE 10)
   lat_mask = INRANGE(all_lat,LatRange)
   lon_mask = INRANGE(all_lon,LonRange)
   date_mask=INRANGE(all_id,DateRange)
   all_land = lite.surface_type EQ 1

   ym = long(ulong64(all_id)/100000000)
   all_year = long(ym/100)
   all_month = fix((ym-all_year*100))
   all_jd=ACOS_ID_TO_JD(all_id)

   loaded_satellite_file = satellite_file

ENDIF
;STOP

;;--------------------------------------------------------------------------------
;;-----
;;--------------------------------------------------------------------------------

; NOTE on save files last digit in sounding_id
;------------------------------------------------
; DDDDHHMMSo, where DDDD = days starting 20140901, 
;   S = 10-sec span, o = measurement type: 
;   1=LnNa, 2=LnGl, 3=LnTarg, 4=LnTrans, 5=OcNa, 6=OcGl, 
;   7=OcTarg, 8=OcTrans, 9=other (mostly non-ocean water)

;;-----
;;----- Get the match-to-good-models criterion
;;-----
model_scrit = NUM2STR(max_model_std,2) + '_' + SC(min_nmodels) + '_' + NUM2STR(max_model_maxdiff,2)
this_match_criteria = STRJOIN(FILE_BASENAME([satellite_file,model_file]),';') + '_' + model_scrit
IF (N_ELEMENTS(match_criteria) EQ 0) THEN match_criteria =''
IF (match_criteria NE this_match_criteria) THEN BEGIN
   PRINT, 'Getting model_match_mask'

   ;;----- Determine the soundings for which we have truth AND they satisfy the model agreement condition.
   i_models = CMSET_OP(Models.id, 'AND', Lite.sounding_id, /index)
   N_models=N_ELEMENTS(i_models)
   PRINT,'Number of model matches='+NUM2STR(n_models,0)     
   P_models= FLOAT(n_models)/FLOAT(N_ELEMENTS(i_models))*100.
   PRINT,'Percent of model matches='+NUM2STR(P_models,2)  
   ;STOP   

   ;;----- Subset to the matched satellite/model elements
   IF (P_models NE 100.00) THEN BEGIN
      StatsAK=Models.StatsAK[i_models]
      Model_id=Models.id[i_models]

      ;;----- Determine satellite matches to subsetted model file.
      i_sat = CMSET_OP(all_id, 'AND', Model_id_good, /index)
      Lite=Lite[i_sat]
      PRINT,'Need to reset all the "all_xxx" data fields...'
      STOP
   ENDIF ELSE BEGIN
      StatsAK=Models.StatsAK
      Model_id=Models.id
   ENDELSE

   ;;----- Subset to the model points that meet selection criteria.
   where_models_criteria = $
      WHERE( $
             (StatsAK.std LE max_model_std) $
              AND (StatsAK.nmodels GE min_nmodels) $
              AND (ABS(StatsAK.maxdiff) LE max_model_maxdiff) $
              , N_Good_Models) 
   P_Good_Models=FLOAT(N_Good_Models)/FLOAT(N_ELEMENTS(StatsAK))*100.
   PRINT,'N_Good_Models=',N_Good_Models
   PRINT,'P_Good_Models=',P_Good_Models

   Model_id_good=Model_id[where_models_criteria]
   StatsAK_good=StatsAK[where_models_criteria]
   HELP,Model_id_good, StatsAK_good
   ;STOP

   ;;----- Subset to the matched satellite/model elements
   IF (P_Good_models NE 100.00) THEN BEGIN

      ;;----- Determine the satellite soundings that have a corresponding model value.
      i_sat = CMSET_OP(all_id, 'AND', Model_id_good, /index)
      Lite=Lite[i_sat]

      N_sat = N_ELEMENTS(i_sat)
      P_sat = FLOAT(n_sat)/FLOAT(N_ELEMENTS(all_id))*100.
      PRINT,'Number of satellite matches='+NUM2STR(n_sat,0)     
      PRINT,'Percent of satellite matches='+NUM2STR(P_sat,2) 
      ;STOP

   ENDIF

   model_match_mask = MAKE_ARRAY(N_ELEMENTS(all_id),/BYTE,VALUE=0b)
   model_match_mask[i_sat] = 1b   
   ;mmm=WHERE(model_match_mask EQ 1b, n_mmm)
   ;PRINT,'Number of model/sat matches='+NUM2STR(n_mmm,0)     
   ;PRINT,'Percent of model/sat matches='+NUM2STR(FLOAT(n_mmm)/FLOAT(N_ELEMENTS(model_match_mask))*100.,2)     
   ;STOP
   match_criteria = this_match_criteria

ENDIF
;;STOP

;;-----
;;----- Set some plotting parameters.
;;-----

;;-----
;;----- Loop over observation modes in list.
;;-----
for m_ = 0, n_elements(mode_list)-1 do begin
;;for m_ = 5, 5 do begin
   do_mode=mode_list[m_]

   ;;-----
   ;;----- Build mask for current observation mode.
   ;;-----
   print, 'Developing full mask for mode = ' + do_mode
   IF (STRMID(VerName,0,5) EQ 'gosat') THEN BEGIN
      CASE do_mode OF
         'landH' : $
            BEGIN
               type_mask = all_land AND all_gain eq HGAIN
               plot_color='green'
            END
         'landM' : $
            BEGIN
               type_mask = all_land AND all_gain eq MGAIN
               plot_color='red'
            END
         'land' : type_mask = all_land
         'oceanH' : $
            BEGIN
               type_mask = ~all_land AND all_gain eq HGAIN
               plot_color='blue'
            END
         'allH'   : type_mask = all_gain eq HGAIN
         'all'   : type_mask = (~all_land AND all_gain eq HGAIN) OR all_land
      ENDCASE
   ENDIF
   IF ( STRMID(VerName,0,4) EQ 'oco2') THEN BEGIN
      CASE do_mode OF
         'landN' : type_mask = last_digit eq 1
         'landG' : type_mask = last_digit eq 2
         'landNG' : type_mask = last_digit eq 1 OR last_digit eq 2
         'landT' : type_mask = last_digit eq 3
         'land' : type_mask = elt(last_digit, [1,2,3])
         'landNG_oceanG': type_mask = elt(last_digit, [1,2,6])
         'all' : type_mask=elt(last_digit,[1,2,3,6])
         'oceanG' : type_mask = last_digit eq 6
         'oceanN' : type_mask = last_digit eq 5
         'ocean' : type_mask = elt(last_digit, [5,6,7,8])
      ENDCASE
   ENDIF

   ;;-----
   ;;----- Gather full mask.
   ;;----- basic_mask; convergence of L2FP
   ;;----- type_mask; observation mode of GOSAT
   ;;----- model_match_mask; valid model data
   ;;----- lat_mask; Satellite latitude in desired range
   ;;----- lon_mask; Satellite longitude in desired range
   ;;----- date_mask; Satellite date in desired range
   ;;----- all_qf_mask; Good quality flag
   ;;-----
   HELP,basic_mask, type_mask, model_match_mask, lat_mask, lon_mask, date_mask, all_qf_mask
   main_mask = basic_mask AND type_mask AND model_match_mask AND lat_mask AND lon_mask AND date_mask AND all_qf_mask
   mm = WHERE(main_mask, N_Mask)
   P_mask = FLOAT(N_Mask) /FLOAT(N_ELEMENTS(main_mask))*100.
   PRINT, 'Using '+NUM2STR(N_Mask,0)+' of '+NUM2STR(N_sat,0)+' satellite/model matched ('+ NUM2STR(P_mask)+'%)'
   ;STOP

   ;;----- Subset model data to main mask.
   ;;----- First match model_id to all_id  (11-April-2024. This seems unnecessary?) 
   PRINT,'Subsetting model data to main_mask'
   i_models = CMSET_OP(Model_id_good, 'AND', all_id[mm], /index)
   model_xco2 = StatsAK_good[i_models].median 
   model_xco2_uncert=StatsAK_good[i_models].std
   model_id=Model_id_good[i_models]

   ;;----- Subset satellite data to main_mask.
   PRINT,'Subsetting satellite data to main_mask'
   ;i_sat = CMSET_OP(all_id[mm], 'AND', model_id, /index)
   ;mm_sat=mm[i_sat]
   mm_sat=mm
   acos_xco2_raw = all_xco2[mm_sat]
   acos_xco2 = all_xco2_bc[mm_sat]
   acos_xco2_adjust = all_xco2_adjust[mm_sat]
   acos_xco2_ap = all_xco2_ap[mm_sat]
   acos_xco2_uncert = all_xco2_uncert[mm_sat]

   HELP,model_xco2,acos_xco2
   
   form = '(a, i8, 2f10.2)'
   PRINT, 'N, Mean, Std dxco2 (Filtered  ) = ', N_Mask, MEAN((acos_xco2-model_xco2)), STDDEV( (acos_xco2-model_xco2)), form=form
   

   ;;----- Set end time. Usually through end of model data.
   ;;----- As of Jan 2022, have model data running through end of 2020. (there is always a long lag)
   ;;----- TO-DO. CHECK THAT THERE IS CONSISTENCY BETWEEN THESE START/STOP AND THE ACTUAL DATA.
   ;;----- LOOKS LIKE THESE START/STOP ARE USED TO MASK DATA FOR SEASONAL ANALYSIS (MAPS).
   StartYear=MIN(all_year[mm_sat])
   EndYear=MAX(all_year[mm_sat])
   PRINT,'Start/end year of satellite data after all filtering = '+NUM2STR(StartYear,0)+'/'+NUM2STR(EndYear,0)
   ;STOP

   SubTitle=VerName+' - MMM :: '
   ;STOP

   ;;----- Create initial scatter plot of the XCO2 adjustment vs time.
   IF (0) THEN BEGIN
      BINXY,all_jd,all_xco2_adjust,XB_adjust_all,YB_adjust_all,DX=30.
      p=PLOT(XB_adjust_all,YB_adjust_all,Symbol='o',Sym_Size=1,Sym_Filled=1,Sym_color='black',Sym_Fill_Color='black',Title=SubTitle+do_mode,XTitle='Julian Day',YTitle='XCO$_2$ adjustment (ppm)',YRange=[-1.5,1.0])

   ENDIF

   ;;-----
   ;;----- Loop over seasons.
   ;;-----
   ii=-1
   for season=0, N_ELEMENTS(Season_Names)-1 do begin
   ;;for season=0, 0 do begin
      ii+=1

      sname=season_names[season]

      CASE STRMID(sname,0,3) OF
         'Ann' : season_mask = ELT(all_month[mm_sat], [1,2,3,4,5,6,7,8,9,10,11,12]) AND ELT(all_year[mm_sat], RANGE(StartYear,EndYear))
         'MAM' : season_mask = ELT(all_month[mm_sat], [3,4,5]) AND ELT(all_year[mm_sat], RANGE(StartYear,EndYear))
         'JJA' : season_mask = ELT(all_month[mm_sat], [6,7,8]) AND ELT(all_year[mm_sat], RANGE(StartYear,EndYear)) 
         'SON' : season_mask = ELT(all_month[mm_sat], [9,10,11]) AND ELT(all_year[mm_sat], RANGE(StartYear,EndYear))  
         'DJF' : season_mask = (all_month[mm_sat] EQ 12 AND ELT(all_year[mm_sat], RANGE(StartYear-1,EndYear-1))) $
                               OR (ELT(all_month[mm_sat], [1,2]) AND ELT(all_year[mm_sat], RANGE(StartYear,EndYear)) )
      ENDCASE
    
      ;;-----
      ;;----- Combine season_mask with main_mask.
      ;;----- Use these on data that was already subsetted by main_mask.
      ;;----- sm=SeasonMask
      ;;-----
      sm = WHERE(season_mask, N_Season)
      P_season = FLOAT(N_Season)/FLOAT(N_Mask)*100. ; percent passing mask for this month/season  
      PRINT, SName + '   ; N='+SC(N_Season)+', fraction=',+P_season 
      ;STOP

      ;;-----
      ;;----- Another set of masks to be used on data arrays
      ;;----- that have not yet been subsetted by main_mask.
      ;;----- smf=SeasonMaskFull
      ;;-----
      smf = mm_sat[sm]
      ;STOP

      lat=all_lat[smf]
      lon=all_lon[smf]
      id=all_id[smf]
      jd=all_jd[smf]
      dxco2=acos_xco2[sm] - model_xco2[sm]
      aod=lite[smf].retrieval.aod_total
      adjust=acos_xco2_adjust[sm]
      HELP,lat,lon,id,jd,dxco2,aod,adjust
      ;;STOP

      IF (0) THEN BEGIN
         DX=30.
         BINXY,jd,adjust,XB_adjust_mode,YB_adjust_mode,DX=DX
         p=PLOT(/OVERPLOT,XB_adjust_mode,YB_adjust_mode,Symbol='d',Sym_Size=1,Sym_Filled=1,Sym_color=plot_color,Sym_Fill_Color=plot_color,Color=plot_color)
         xco2_adjust_bin={mode:do_mode,XB_all:XB_adjust_all,YB_all:YB_adjust_all,XB:XB_adjust_mode,YB:YB_adjust_mode,DX:DX}

         BINXY,jd,dxco2,XB_dxco2_mode,YB_dxco2_mode,DX=30.
         dxco2_bin={mode:do_mode,XB:XB_dxco2_mode,YB:YB_dxco2_mode,DX:DX}

         SaveFileName_binned=Pdf_File_Root+'_binned_'+LatText+'_'+do_mode+'_'+DateRangeString + '_'+sname+'.sav'
         PRINT,SaveFileName_binned ;& STOP
         SAVE,FILENAME=SaveFileName_binned,xco2_adjust_bin, dxco2_bin
      ENDIF

      ;;-----    
      ;;----- Select which data variable is to be plotted.
      ;;-----
      IF (STRMID(VerName,6,3) EQ 'v9') THEN BEGIN
         ZData=[[dxco2],[aod]] 
         VarList=['dxco2','aod_total']
      ENDIF ELSE BEGIN
         ZData=[[dxco2],[aod],[adjust]] 
         VarList=['dxco2','aod_total','xco2_adjust']
      ENDELSE

      ;;-----
      ;;----- Grid data by lat/lon.
      ;;-----
      IF (GridMap) THEN BEGIN

         PythonCT='bwr'
         Python_XLabel='Longitude'
         Python_YLabel='Latitude'

         FillFloat=-999.9 & FillInteger=-999
         Pdf_File =Pdf_File_Root+'_map_'+LatText+'_'+do_mode+'_'+DateRangeString + '_'+sname+'.pdf'
         OutNCFile=Pdf_File_Root+'_data_'+LatText+'_'+do_mode+'_'+DateRangeString + '_'+sname+'.nc'
         OutInfo=MAKE_ARRAY(18,/STRING,VALUE='')
         OutInfo[0]=Pdf_File
         OutInfo[1]='$\Delta$XCO$_2$ [ppm]'
         OutInfo[2]=SubTitle+SName+' :: '+Do_Mode
         OutInfo[3]='N='+STRING(FLOAT(N_Season)/1000.,FORMAT='(F6.1)')+'k (SS)'
         OutInfo[4]='G='+num2str(P_season,1)+'%'
         OutInfo[10]=PythonCT
         OutInfo[11]=Python_XLabel
         OutInfo[12]=Python_YLabel
         OutInfo[13]=STRING(FillFloat,FORMAT='(F6.1)')
         OutInfo[14]=STRING(FillInteger,FORMAT='(I4)')

         ;;----- Determine DX/DY for gridding.
         NX=LONG((LonRange[1]-LonRange[0])/DeltaLon_1)
         NY=LONG((LatRange[1]-LatRange[0])/DeltaLat_1)
         
         ;STOP
         PRINT,'GRID_DATA_WRITE_NETCDF for spatial map of dxco2...'
         GRID_DATA_WRITE_NETCDF_MULTI, $
            Lon, $
            Lat, $
            ZData, $   ;; inputs
            VarList, $
            ;X_Array,Y_Array,Z_Binned, $ ;; output gridded X/Y/Z data
            'XVar',DeltaLon_1, NX, LonRange,'lon','longitude (degrees)', $
            'YVar',DeltaLat_1, NY, LatRange,'lat','latitude (degrees)', $
            'NVar','N','Number of soundings per lat/lon grid box', $
            ;'ZVar','ppm','Delta XCO2 (measured-MMM)', $
            Pdf_File,OutNcFile,OutInfo,BinMinCount=BinMinCount,Verbose=0
         ;STOP
      ENDIF


      ;;-----
      ;;----- Plot histograms.
      ;;-----
      IF (PlotHistogram) THEN BEGIN
         pdf_file =pdf_file_root+'_hist_'+do_mode+'_'+sname+'.eps'
         PRINT,'Pdf_File=',Pdf_File
         SETUP_STANDARD_DEVICE, PDF_File, PCharThick=4., Thick=8, PCharSize=2., XCharSize=2., YCharSize=2.

         ;;----- Search for various surface types.
         ;;----- Use mask [smf]=SeasonMaskFull, which combines the main_mask+season_mask
         LandNH=WHERE( (all_land[smf] EQ 1) AND (all_lat[smf] GE 0.),NLandNH) & PRINT,'NLandNH=',NLandNH
         IF (NLandNH LT 1) THEN LandNH=0
         LandSH=WHERE( (all_land[smf] EQ 1) AND (all_lat[smf] LT 0.),NLandSH) & PRINT,'NLandSH=',NLandSH
         IF (NLandSH LT 1) THEN LandSH=0
         OceanNH=WHERE( (all_land[smf] EQ 0) AND (all_lat[smf] GE 0.), NOceanNH) & PRINT,'NOceanNH=',NOceanNH
         OceanSH=WHERE( (all_land[smf] EQ 0) AND (all_lat[smf] LT 0.), NOceanSH) & PRINT,'NOceanSH=',NOceanSH
         PRINT,'Total LandNH/LandSH/OceanNH/OceanSH=',NLandNH+NLandSH+NOceanNH+NOceanSH
         PRINT,'N Elements in ZData=',N_ELEMENTS(ZData)
         ;;STOP

         Cmap=50 ;; 50=CB Blue/Green
         XMin=-5. & XMax=5. & XTickLen=0.05
         loc = RANGE(XMin, XMax, 41)
         charsize=2.5 & charthick=12 & XCharSize=1.5 & YCharSize=1.5
         XTitle='dXCO2 (ppm)'
         Thick=12 & BChar=2.
         LOADCT,39

         ;;----- Generate histogram plot
         HIST, ZData, loc=loc, $
               Position=[0.55,0.20,0.95,0.92], $ 
               Title=SubTitle+SName+' :: '+Do_Mode, $
               XTitle=XTitle, xr=[-5,5], XTickLen=XTickLen, $
               YTitle='Number of Soundings', YMinor=-1, YTickLen=0.05, YRange=[10,9000], $
               XCharSize=XCharSize, YCharSize=YCharSize, $
               CharSize=CharSize,CharThick=CharThick, $
               Font=1
         LOADCT,Cmap
         Color_L1=100 & Color_L2=225
         Color_O1=25 & Color_O2=75
         IF (NLandNH GT 0) THEN HIST, ZData[LandNH], loc=loc, /oplot, col=Color_L1, th=Thick
         IF (NLandSH GT 0) THEN HIST, ZData[LandSH], loc=loc, /oplot, col=Color_L2, th=Thick
         IF (NOceanNH GT 0) THEN HIST, ZData[OceanNH], loc=loc, /oplot, col=Color_O1, th=Thick
         IF (NLandSH GT 0) THEN HIST, ZData[OceanSH], loc=loc, /oplot, col=Color_O2, th=Thick

         ;;----- Calculate Statistics
         IF (NLandNH GT 0) THEN STATS_TET,ZData[LandNH],ZData_landNH_stats
         IF (NLandSH GT 0) THEN STATS_TET,ZData[LandSH],ZData_landSH_stats
         IF (NOceanNH GT 0) THEN STATS_TET,ZData[OceanNH],ZData_OceanNH_stats
         IF (NOceanSH GT 0) THEN STATS_TET,ZData[OceanSH],ZData_OceanSH_stats
         LOADCT,39
         VLine,0.0,Color=0,LineStyle=2
         IF (0) THEN BEGIN
            VLine,ZData_landH_stats.median,Color=150,LineStyle=2
            VLine,ZData_landM_stats.median,Color=250,LineStyle=2
            VLine,ZData_Ocean_stats.median,Color=50,LineStyle=2
         ENDIF

         ;;----- Write statistics to plot legend.
         LOADCT,Cmap
         X0=0.025 & Y0=0.85 & DX=.15 & DY=-.05
         Y1=Y0-ABS(DY*6.)
         Y2=Y1-ABS(DY*6.)
         LOADCT,39
         co_xyouts_matrix, X0, Y0, DX, DY, ['Land (All)'], /norm, $
            ZData, color=[0], digits=2, $
            /N,/mean,/std,CharSize=3.,CharThick=10,Font=1
         LOADCT,Cmap
         IF (NLandNH GT 0) THEN $
            co_xyouts_matrix, X0, Y1, DX, DY, ['Land (NH)'], /norm, $
               ZData[LandNH], color=[Color_L1], digits=2, $
               /N,/mean,/std,CharSize=3.,CharThick=10,Font=1
         IF (NLandSH GT 0) THEN $
            co_xyouts_matrix, X0, Y2, DX, DY, ['Land (SH)'], /norm, $
               ZData[LandSH], color=[Color_L2], digits=2, $
               /N,/mean,/std,CharSize=3.,CharThick=10,Font=1

         X0=0.65 & Y0=0.85
         LOADCT,Cmap 
         IF (NOceanNH GT 0) AND (NOceanSH GT 0) THEN $
            co_xyouts_matrix, X0, Y0, DX, DY, ['OceanNH','OceanSH'], /norm, $
               ZData[OceanNH],ZData[OceanSH], color=[0,Color_ONH,Color_OSH], digits=2, $
               /N,/mean,/median,/std,CharSize=2

         ;;----- Close output hist file.
         ;pause
         DEVICE, /close_file
         PlotPng='ps2png -n -l '+Pdf_File
         SPAWN, PlotPng
         ;STOP

      ENDIF

      ;;-----
      ;;----- Hovmoller plots.
      ;;----- X=DOY, Y=Lat, Z=ZData, e.g., Delta XCO2
      ;;-----
      IF (GridHovmoller AND STRMID(sname,0,3) EQ 'Ann') THEN BEGIN

         SurfaceLabel=do_mode
         NumSound=N_ELEMENTS(Lat)
         DSS=jd-MIN(jd)+1L & MinDSS=MIN(DSS) & MaxDSS=MAX(DSS)
         CALDAT,jd,mon2,day2,year2,hour2,min2,sec2
         N=N_ELEMENTS(jd)
         StartDateString=STRING(Year2[0],FORMAT='(I4)')+'-'+STRING(Mon2[0],FORMAT='(I2.2)')+'-'+STRING(Day2[0],FORMAT='(I2.2)')
         EndDateString=STRING(Year2[N-1],FORMAT='(I4)')+'-'+STRING(Mon2[N-1],FORMAT='(I2.2)')+'-'+STRING(Day2[N-1],FORMAT='(I2.2)')
         PRINT,'StartDateString=',StartDateString
         PRINT,'EndDateString=',EndDateString
         PythonCT='bwr'
         Python_XLabel='Date (YYYY-MM)'
         Python_YLabel='Latitude'

         NDaysRange=[0,MAX(DSS)]
         DX=30 & NX=CEIL( (FLOAT(MaxDSS)-FLOAT(MinDSS) ) / FLOAT(DX) ) & PRINT,'NX=',NX 
         DY=10 & NY=CEIL(ABS(LatRange[1]-LatRange[0])/FLOAT(DY)) & PRINT,'NY=',NY

         ;;GOSAT StartDateString=2014-05-01 EndDateString=2018-12-31
         ;;OCO-2 StartDateString=2014-09-06 EndDateString=2020-12-30

         ;;----- Calculate mean signal for time range 2009 through June 2016 for direct v7 to v9 comparison.
         IF (0) THEN BEGIN
            xel=WHERE(INRANGE(year2,[2009,2016]) AND INRANGE(mon2,[1,6]),nfit) & PRINT,'nfit=',nfit ;;& STOP
            MeanComp=MEAN(ZData[xel]) & PRINT,'MeanComp=',MeanComp
            SDComp=STDDEV(ZData[xel]) & PRINT,'SDComp=',SDComp
            ;STOP
         ENDIF

         ;;-----
         ;;----- Grid data by time/lat.
         ;;-----
         IF (1) THEN BEGIN
            OutPlotFile=Pdf_File_Root+'_hovmoller_'+do_mode+'_'+sname+'.pdf'
            OutNCFile=Pdf_File_Root+'_hovmoller_'+do_mode+'_'+sname+'.nc'
            OutInfo=MAKE_ARRAY(15,/STRING,VALUE='')
            OutInfo[0]=OutPlotFile
            OutInfo[1]='$\Delta$XCO$_2$ [ppm]'
            OutInfo[2]=VerName+', '+SName+', '+Do_Mode
            OutInfo[3]='N='+STRING(FLOAT(N_Season)/1000.,FORMAT='(F6.1)')+'k'
            OutInfo[4]='G='+num2str(P_season,1)+'%'
            OutInfo[10]=PythonCT
            OutInfo[11]=Python_XLabel
            OutInfo[12]=Python_YLabel
         GRID_DATA_WRITE_NETCDF_MULTI, $
            DSS, $
            Lat, $
            ZData, $   ;; inputs
            VarList, $
            ;X_Array,Y_Array,Z_Binned, $ ;; output gridded X/Y/Z data
            'XVar',DX, NX, NDaysRange,'days','days since start', $
            'YVar',DY, NY, LatRange,'lat','latitude (degrees)', $
            'NVar','N','Number of soundings per day/lat grid box', $
            ;'ZVar','ppm','Delta XCO2 (measured-MMM)', $
            Pdf_File,OutNcFile,OutInfo,BinMinCount=BinMinCount,Verbose=0

            ;STOP
         ENDIF
         ;DSS=!NULL
         
      ENDIF


      ;;----- End loop over Season.
   ENDFOR

   ;;----- End loop over Mode_List.
   ;;STOP
ENDFOR

STOP          
END
