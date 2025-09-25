;;PRO MAIN_COLLOCATED_GOSAT_VS_OCO2
;;----- This procedure borrowed heavily from CWO's /home/codell/idlprogs/aos/models/OCO/l2/b9_acos/compare_acos_oco2_b9.pro
;;----- Modifications by T.E.T Autumn, 2020 to adapt and generate analysis/plots for the GOSAT v9 data paper.
;;----- Modifications by T.E.T late Winter, 2022 to adapt to Tropical IAV work. 
;;-----                        Following the logic of the OCO-3 vs OCO-2 v10 collocation work, why not apply the methodology
;;-----                        to GOSAT v9 versus OCO-2 v10? So use the collocations to calculate dxco2 (already done),
;;-----                        troll for some explanatory variable (geometry?).

;;----- Highlevel settings.
;VerNum_OCO='9'
;VerNum_OCO='10'
VerNum_OCO='11.1'
Min_oco=3 ;; minimum number of OCO soundings collocated to GOSAT.

;VerNum_GOSAT='v9'
;VerNum_GOSAT='v9b'
;VerNum_GOSAT='v9c_20240326' ;; land-H=XCO2_RAW, land-M=AOD_STRATAER, ocean-H=AOD_FINE
;VerNum_GOSAT='v9c3_r20240422' ;; skewed-sine fit to day for all modes
;VerNum_GOSAT='v9_ml_r20240503' ;; Will Keely ML correction
VerNum_GOSAT='v9_ml_r20240626' ;; Will Keely ML correction

Correct_AK=1 & Two_Term_AK=0

;;----- List of years to process
YearList=[2014,2015,2016,2017,2018,2019,2020]
;;YearList=[2016]
;;YearList=[2016,2017,2018,2019,2020]

;;----- Latitude range.
LatRange=[-90.,90.]
Region='global'
;LatRange=[-30.,30.]
;Region='tropics'

ID_Start=20140906 & ID_Stop =20200601

;;----- Select plots.
PlotAdjXco2TimeSeries=0
PlotAvgKer=0
PlotCollocationSummary=0

;;----- Grid for maps or for Hovmollers.
GridData=1
GridMap=0
GridHovmoller=1
PlotTimeSeries=0
WriteStats=0

;;----- Viewing mode.
;view_mode='landH' & view_type=1.
;view_mode='landM' & view_type=1.
view_mode='oceanH' & view_type=0.
;view_mode='land'

;;----- Use this for gridded spatial maps.
IF (GridMap) THEN BEGIN
   plot_type='seasonal' & independent_years=0
   ;plot_type='annual' & independent_years=0
ENDIF

;;----- Use this for gridded Hovmoller data and lat bin time series.
IF (GridHovmoller OR PlotTimeSeries) THEN BEGIN
   plot_type='annual' & independent_years=0
ENDIF

IF (plot_type EQ 'seasonal') THEN seaslist=['DJF','MAM','JJA','SON'] ELSE seaslist=['annual']
NS=N_ELEMENTS(seaslist)

;;----- Other settings.
;skip_calc=0

;;-----
;;----- Set some plotting parameters.
;;-----
OutDirBase='/home/ttaylor/analysis_utilities/tropical_iav/plots/'
OutDir=OutDirBase+'gosat_'+VerNum_gosat+'_vs_oco2_v'+VerNum_oco+'/' & FILE_MKDIR,OutDir
DeltaLat_1=2.5 & DeltaLon_1=5.0
;DeltaLat_1=15. & DeltaLon_1=20.
BinSizeString_1=STRING(DeltaLat_1,FORMAT='(I2)')+' lat x '+STRING(DeltaLon_1,FORMAT='(I2)')+' lon'
BinMinCount=1
FillFloat=-999.9 & FillInteger=-999L
HovInc_Day=30 & HovInc_Lat=10

;;----- Build structure of region information.
T={Name:'',LatBin:FLTARR(2),LonBin:FLTARR(2)}
NLocations=1
Location=REPLICATE(T,NLocations)
;Location[0].Name='Tropical-Americas' & Location[0].LatBin=[-35,35.] & Location[0].LonBin=[-120.,-30.]
;Location[1].Name='Tropical-Africa' & Location[1].LatBin=[-35,35.] & Location[1].LonBin=[-20.,60.]
;Location[2].Name='Tropical-Asia-Australia' & Location[2].LatBin=[-35,35.] & Location[2].LonBin=[60.,180.]
;Location[0].Name='Tropics' & Location[0].LatBin=[-35.,35.] & Location[0].LonBin=[-180.,180.]
Location[0].Name='global' & Location[0].LatBin=[-90.,90.] & Location[0].LonBin=[-180.,180.]

IF (0) THEN BEGIN
   ;;-----
   ;;----- 11-March-2021.
   ;;----- Use CWO collocated file.
   ;;----- max_dlat=2, max_dlon=3, max_dtime =2, max_dist = 300 km, and random when above 100.
   ;;----- Chris recommends that to change dlat/dlon/ddist/dtime, the collocation code should be rerun
   ;;----- instead of subsetting within this analysis.
   ;;-----
   file = '/home/codell/GOSAT/acos_results/b9_acos/match_oco2/match_acos_oco2_b10_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'
   max_dist=300. & max_dlat=2. & max_dlon=3. & dtime_allowed_range=[-2.,2.]
   oco_max_xco2_variability=1.5 ;; variablity of the 1 to 100 OCO-2 soundings collocated to GOSAT.
ENDIF

IF (0) THEN BEGIN
   ;;-----
   ;;----- 13-Nov-2023.
   ;;----- Use TET collocated file.
   ;;----- max_dlat=2, max_dlon=3, max_dtime =2, max_dist = 300 km, and random when above 100.
   ;;----- Generated from code <analysis_utilties/gosat_b9_analysis/code/acos_b9_paper/match_b9_acos_b10_oco2_v2.pro
   ;;----- The output save file now includes OCO2 co2_profile. Also, change "acos" to "gosat".
   ;;----- Current output directory is "/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/"
   ;;----- 
   ;;-----
   collocate_dir='/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/'
   file = collocate_dir+'match_gosat_v9_oco2_v10_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'
   max_dist=300. & max_dlat=2. & max_dlon=3. & dtime_allowed_range=[-2.,2.]
   oco_max_xco2_variability=1.5 ;; variablity of the 1 to 100 OCO-2 soundings collocated to GOSAT.
ENDIF

IF (VerNum_GOSAT EQ 'v9c_20240326') THEN BEGIN
   ;;-----
   ;;----- 15-Feb-2024.
   ;;----- Use TET collocated file for OCO-2 v11.1.
   ;;----- max_dlat=2, max_dlon=3, max_dtime =2, max_dist = 300 km, and random when above 100.
   ;;----- Generated from code <analysis_utilties/tropical_iav/code/main_1_match_gosat_to_oco2.pro
   ;;-----
   collocate_dir='/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/'
   CollocatedFile = collocate_dir+'match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'
   max_dist=300. & max_dlat=2. & max_dlon=3. & dtime_allowed_range=[-2.,2.]
   oco_max_xco2_variability=1.5 ;; variablity of the 1 to 100 OCO-2 soundings collocated to GOSAT.

   GosatAggSaveFile='/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c_r20240326/gosat_v9c_r20240326_lite_20090420_20200630_filtered_xmedium.sav'
ENDIF

IF (VerNum_GOSAT EQ 'v9c2_20240329') THEN BEGIN
   collocate_dir='/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/'
   CollocatedFile = collocate_dir+'match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'
   max_dist=300. & max_dlat=2. & max_dlon=3. & dtime_allowed_range=[-2.,2.]
   oco_max_xco2_variability=1.5 ;; variablity of the 1 to 100 OCO-2 soundings collocated to GOSAT.

   GosatAggSaveFile='/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c2_r20240329/gosat_v9c2_r20240329_lite_20090420_20200630_filtered_xmedium.sav'
ENDIF

IF (VerNum_GOSAT EQ 'v9c3_r20240422') THEN BEGIN
   collocate_dir='/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/'
   CollocatedFile = collocate_dir+'match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'
   max_dist=300. & max_dlat=2. & max_dlon=3. & dtime_allowed_range=[-2.,2.]
   oco_max_xco2_variability=1.5 ;; variablity of the 1 to 100 OCO-2 soundings collocated to GOSAT.

   GosatAggSaveFile='/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c3_r20240422/gosat_v9c3_r20240422_lite_20090420_20200630_filtered_xmedium.sav'
ENDIF

IF (VerNum_GOSAT EQ 'v9_ml_r20240503') THEN BEGIN
   collocate_dir='/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/'
   CollocatedFile = collocate_dir+'match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'
   max_dist=300. & max_dlat=2. & max_dlon=3. & dtime_allowed_range=[-2.,2.]
   oco_max_xco2_variability=1.5 ;; variablity of the 1 to 100 OCO-2 soundings collocated to GOSAT.

   GosatAggSaveFile='/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9_ml_r20240503/gosat_v9_ml_r20240503_lite_20090420_20200630_filtered_xmedium.sav'
   CorrectionVariable='Machine Learning'
   ;;STOP
ENDIF

IF (VerNum_GOSAT EQ 'v9_ml_r20240626') THEN BEGIN
   collocate_dir='/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/'
   CollocatedFile = collocate_dir+'match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'
   max_dist=300. & max_dlat=2. & max_dlon=3. & dtime_allowed_range=[-2.,2.]
   oco_max_xco2_variability=1.5 ;; variablity of the 1 to 100 OCO-2 soundings collocated to GOSAT.

   GosatAggSaveFile='/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9_ml_r20240626/gosat_v9_ml_r20240626_lite_20090420_20200630_filtered_xmedium.sav'
   CorrectionVariable='Machine Learning'

   ;STOP
ENDIF


IF (N_ELEMENTS(loaded_file) EQ 0) then loaded_file=''
;IF loaded_File NE CollocatedFile then skip_calc=0
;IF skip_calc then goto, skip_calc_

IF (loaded_file NE CollocatedFile) THEN BEGIN

   PRINT, 'Restoring satellite matched file... ' + FILE_BASENAME(CollocatedFile)

   RESTORE_COLLOCATED_DATA_FROM_MATCH_FILE, $
      CollocatedFile, $
      GosatAggSaveFile, $
      VerNum_GOSAT, $
      Match, $
      Lite, $
      Gosat, $
      Oco, $
      Metrics, $
      loaded_file 
ENDIF
;STOP

;;-----
;;----- 1. Match at the oco2 individual sounding level
;;----- Chris recommends that to change dlat/dlon/ddist/dtime, the collocation code should be rerun
;;----- instead of subsetting within this analysis. Not sure why.
;;----- Current settings pass 100% of OCO data.
;;-----
oco2_good =  $
   1 EQ 1 $
   AND INRANGE(match.oco2.xco2, [300,900.]) $
   AND (match.oco2.xco2_quality_flag EQ 0b) $
   AND (ABS(Metrics.dlat) LE max_dlat) $
   AND (ABS(Metrics.dlon) LE max_dlon) $ 
   AND (ABS(Metrics.dist) LE max_dist) $
   AND INRANGE(Metrics.dtime, dtime_allowed_range) ;$
   ;AND match.oco2.sounding.land_water_indicator EQ 0

N_oco = LONG(TOTAL(oco2_good, 1))
w1 = WHERE(N_oco GE Min_oco, N_w1)
PRINT,'Number of good OCO2 collocations=',N_w1
PRINT,'Percent of good OCO2 collocations=',FLOAT(N_w1)/FLOAT(N_ELEMENTS(N_oco))*100.
Match_metrics=Match[w1] 
Gosat_metrics=Gosat[w1]
Oco_metrics=Oco[w1]
;;STOP

;;-----
;;----- 2. Compute mean OCO-2 sounding statistics per match
;;-----
;oco_stddev_xco2=!NULL
IF (N_ELEMENTS(oco_stddev_xco2) NE N_w1) THEN BEGIN

   PRINT, 'Averaging the OCO soundings collocated to each GOSAT sounding.'

   ;;----- Generate data arrays.
   oco_mean_struct = REPLICATE(Match_metrics[0].oco2[0], N_w1)
   oco_mean_xco2 = FLTARR(N_w1) + 9999.9
   oco_median_xco2 = FLTARR(N_w1) + 9999.9
   oco_stddev_xco2 = FLTARR(N_w1) + 9999.9
   oco_mean_jd = DBLARR(N_w1)
   mean_dtime=oco_stddev_xco2*0.
   mean_ddist=oco_stddev_xco2*0.

   xco2_prior_adj_a = FLTARR(N_w1)
   xco2_prior_adj_b = FLTARR(N_w1)

   ;;----- Loop over number of GOSAT/OCO collocations.
   FOR i=0,N_w1-1 DO BEGIN

      PRINT,'i='+NUM2STR(i+1,0)+' of '+NUM2STR(N_w1,0)
      q = WHERE(oco2_good[*,w1[i]], N_q)
      oco_mean_struct[i] = AVERAGE_STRUCT(Match_metrics[i].oco2[q])

      ;----- fix the longitude
      this_lon =  Match_metrics[i].oco2[q].longitude
      mnlon = this_lon[0] + MEAN(longitude_difference(this_lon, this_lon[0]))
      IF (mnlon GT 180.) THEN mnlon=mnlon-360. 
      IF (mnlon LT -180.) THEN mnlon=mnlon+360.
      oco_mean_struct[i].longitude=mnlon

      ;;----- Calculate average OCO stats from the 1-100 collocated soundings.
      xco2 = Match_metrics[i].oco2[q].xco2
      oco_mean_xco2[i] = MEAN(xco2)       
      oco_median_xco2[i] = MEDIAN(xco2)       
      oco_stddev_xco2[i] = STDDEV(xco2)
      oco_mean_jd[i]= MEAN(oco_metrics[i].jd[q])

      ;;----- Calculate time and distance differences between OCO and GOSAT.
      mean_dtime[i] = (oco_mean_jd[i]-gosat_metrics[i].jd)*24.
      mean_ddist[i] = CO_SPHDIST(oco_mean_struct[i].longitude, oco_mean_struct[i].latitude, gosat_metrics[i].lon, gosat_metrics[i].lat, /deg, /approx)*111.1
      IF ( (mean_ddist[i] GT max_dist) OR (~INRANGE(mean_dtime[i], dtime_allowed_range)) ) THEN STOP

      ;;----- Correct for different CO2 priors between versions.
      ;;----- This provides the adjustment to XCO2 based on the difference in priors, i.e., order 0.5 ppm.
      ;;----- 9-Nov-2023. Reran analysis_utilities/gosat_b9_analysis/code/acos_b9_paper/match_b9_acos_b10_oco2_v2.pro
      ;;----- to regenerate the gosat-oco2 matched file to include oco2 pressure weight and averaging kernel profiles.
      ;;----- Neither GOSAT nor OCO2 Lite files contain actual retrieved CO2 profiles, but it seems like it should
      ;;----- be calculable from what we have in hand? Essentially reverse the XCO2 calculation.
      ;;-----
      ;IF (VerNum_GOSAT EQ 'v9c_20240214') $
      ;   THEN gosat_xco2=gosat_xco2_harm_w[i] $
      ;   ELSE gosat_xco2 = Match_metrics[i].gosat.xco2
      gosat_prior = Match_metrics[i].gosat.co2_profile_apriori
      gosat_h = Match_metrics[i].gosat.pressure_weight
      gosat_a = Match_metrics[i].gosat.xco2_averaging_kernel
      oco2_prior = oco_mean_struct[i].co2_profile_apriori
      oco2_h = oco_mean_struct[i].pressure_weight
      oco2_a = oco_mean_struct[i].xco2_averaging_kernel

      ;;----- 29-Nov-2023. 
      ;;----- Implement CO2 prior and AK corrections
      ;;----- Former as in Eq A10 in [Wunch, ACP, 2010].
      ;;----- Also equivelant to Eq 3 in [Taylor, AMT, 2023].
      ;;----- The sum of the two terms yields the total correction to GOSAT XCO2 to harmonize it to OCO-2. 
      ;;-----
      true_profile=oco2_prior
      xco2_prior_adj_a[i] = TOTAL( gosat_h * (gosat_a - 1.    ) * (gosat_prior - true_profile) )
      xco2_prior_adj_a_text='Term 1 = $\Sigma$ { gosat_pwf * (gosat_ak - 1.    ) * (gosat_prior - true_profile) }'  
      xco2_prior_adj_b[i] = TOTAL( gosat_h * (gosat_a - oco2_a) * (gosat_prior - true_profile) )
      xco2_prior_adj_b_text='Term 2 = $\Sigma$ { gosat_pwf * (gosat_ak - oco2_ak ) * (gosat_prior - true_profile) }'  
      
   ENDFOR  

   STATS,xco2_prior_adj_a
   PRINT,' '
   PRINT,' '
   STATS,xco2_prior_adj_b
   PRINT,' '
   PRINT,' '
   ;;STOP 
ENDIF

;;-----
;;----- 7-Nov-2023.
;;----- Write matched data to H5 output file.
;;-----
WriteH5=0
IF (WriteH5) THEN BEGIN

   OutH5Dir=collocate_dir
   OutH5File=OutH5Dir+FILE_BASENAME(CollocatedFile,'.sav')+'.h5'
   IF (FILE_TEST(OutH5File)) THEN FILE_DELETE,OutH5File
   PRINT,'Writing collocated data to file ',OutH5File
   STRUCT_TO_H5_SIMPLE,Match_metrics,OutH5File

   ;;----- Add XCO2 CO2 prior and AK adjustment arrays to H5 output file.
   ADD_H5_FIELD,OutH5File,xco2_prior_adj_a,'xco2_ak_term_1'
   ADD_H5_FIELD,OutH5File,xco2_prior_adj_b,'xco2_ak_term_2'

   ;;---- Add OCO-2 aggregated fields to H5 output file.
   ADD_H5_FIELD,OutH5File,oco_mean_struct,'oco_aggregated_mean_fields'
   ADD_H5_FIELD,OutH5File,oco_mean_jd,'oco_aggregated_julian_day_mean'
   ADD_H5_FIELD,OutH5File,oco_mean_xco2,'oco_aggregated_xco2_mean_ppm'
   ADD_H5_FIELD,OutH5File,oco_median_xco2,'oco_aggregated_xco2_median_ppm'
   ADD_H5_FIELD,OutH5File,oco_stddev_xco2,'oco_aggregated_xco2_stddev_ppm'
   ADD_H5_FIELD,OutH5File,oco_max_xco2_variability,'oco_max_allowed_xco2_variability'

   ;;---- Add mean dtime and dist fields.
   ADD_H5_FIELD,OutH5File,mean_dtime,'mean_delta_time_oco_minus_gosat_hours'
   ADD_H5_FIELD,OutH5File,mean_ddist,'mean_delta_distance_km'

   ;;----- Add GOSAT harmonized XCO2 and XCO2 adjustment to H5 output file.
   ADD_H5_FIELD,OutH5File,gosat_metrics.xco2_harm,'gosat_xco2_harmonized_to_oco2'
   ADD_H5_FIELD,OutH5File,gosat_metrics.xco2_adj,'gosat_xco2_adjustment'

   ;STOP
ENDIF


;;-----
;;----- Sanity check plot of the calculated xco2 before and after adjustment.
;;-----
VarName='JULIAN_DAY'
DataDir='bc_vars_time'+'/'+VarName+'/'
IF (PlotAdjXco2TimeSeries) THEN BEGIN

   IF (1) THEN BEGIN
      ViewModes=['All Observation Modes', 'Land-H','Land-M','Ocean-H'] & NumViewModes=N_ELEMENTS(ViewModes)
   ENDIF

   PlotSingles=0
   PlotSinglesFit=1
   PlotBins=1

   ;;----- Plot settings.
   Margin=[0.1,0.2,0.1,0.2]
   NCols=1 & NRows=4
   FontSize_1=10
   FitThick=2
   BinSymSize=0.75
   Symbol='Circle'
   SymIncrement=20
   DisplayText='Displaying 1 out of every '+NUM2STR(SymIncrement,0)+' points per series'

   ;;----- Output file for plot.
   FILE_MKDIR,OutDir+DataDir
   OutPlotFile=OutDir+DataDir+'dxco2_vs_'+VarName+'_harmonized.png'

IF (1) THEN BEGIN
   ;;----- XData is GOSAT Julian Day.
   XData=gosat_metrics.jd

   ;;----- YData is delta XCO2.
   ;;----- Calculated using original BC XCO2.
   IF (Correct_AK) THEN BEGIN
      IF (Two_Term_AK) THEN BEGIN
         YData=(Gosat_metrics.xco2_bc + xco2_prior_adj_a + xco2_prior_adj_b) - oco_mean_struct.xco2
         YData_2=(Gosat_metrics.xco2_harm + xco2_prior_adj_a + xco2_prior_adj_b) - oco_mean_struct.xco2
      ENDIF ELSE BEGIN
         YData=(Gosat_metrics.xco2_bc + xco2_prior_adj_b) - oco_mean_struct.xco2
         YData_2=(Gosat_metrics.xco2_harm + xco2_prior_adj_b) - oco_mean_struct.xco2
      ENDELSE

   ENDIF ELSE BEGIN
      YData=Gosat_metrics.xco2_bc - oco_mean_struct.xco2
      YData_2=Gosat_metrics.xco2_harm - oco_mean_struct.xco2
   ENDELSE

   ;;----- Index by observation modes.
   landh = WHERE( (Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '72') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)
   landm = WHERE( (Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '77') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)
   oceanh = WHERE( (~Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '72') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)
 


ENDIF ELSE BEGIN
   ;;----- collocated values saved in <main_2_fit_and_plot_collocated_gosat_to_oco2.pro>
   SaveFileName='./main_2_fit_x_y_data.sav'
   RESTORE,FILENAME=SaveFileName,/VERBOSE
   ;STOP
ENDELSE

   ;;----- Check for finite elements in the input data set.
   x_good=WHERE(FINITE(XData), nx_good)
   y_good=WHERE(FINITE(YData), ny_good)
   IF (nx_good NE ny_good) THEN BEGIN
      MATCH,x_good,y_good,x_final,y_final
      XData=XData[x_good[x_final]]
      YData_1=YData[y_good[y_final]]
      YData=!NULL
      YData_2=YData_2[y_good[y_final]]
      ;STOP
   ENDIF

   ;;----- X axis setup.
   Dummy = LABEL_DATE(DATE_FORMAT='%N/%Z')
   XTitle='Time [MM/YY]'
   ;XRange=[MIN(gosat_metrics.jd)-60L,MAX(gosat_metrics.jd)+60L]
   stretch_lower= ( MAX(XData) - MIN(XData) ) * 0.05
   stretch_upper= ( MAX(XData) - MIN(XData) ) * 0.55
   xrange_plot=[MIN(XData)-stretch_lower,MAX(XData)+stretch_upper]

   ;;----- Y axis setup
   YRange=[-1.0, 1.0]

   ;;-----
   ;;----- Loop over GOSAT viewing modes.
   ;;-----
   FOR Mode=0,NumViewModes-1 DO BEGIN

      CASE Mode OF
         0: BEGIN
               PlotPos=1
               CURRENT=0
               PlotLabel='(a)'
               XTitle=''
               YTitle=''
               ws=[landh,landm,oceanh]
               Color_1='black'
               Color_2='gray'
            END
         1: BEGIN
               PlotPos=2
               CURRENT=1
               PlotLabel='(b)'
               XTitle=''
               YTitle=''
               ws=landh
               Color_1='green'
               Color_2='lime'
               ;CorrectionVariable='XCO2_RAW'
               ;CorrectionVariable='AOD_TOTAL'
               ;CorrectionVariable='Skewed-sine vs time'
            END
         2: BEGIN
               PlotPos=3
               CURRENT=1
               PlotLabel='(c)'
               XTitle=''
               YTitle='$\Delta$XCO$_2$ [ppm] (GOSAT $-$ OCO)'
               ws=landm
               Color_1='orange'
               Color_2='red'
               ;CorrectionVariable='AOD_STRATAER'
               ;CorrectionVariable='AOD_TOTAL'
               ;CorrectionVariable='Skewed-sine vs time'
            END
         3: BEGIN
               PlotPos=4
               CURRENT=1
               PlotLabel='(d)'
               XTitle='Time [MM/YY]'
               YTitle=''
               ws=oceanh
               Color_1='blue'
               Color_2='aqua'
               ;CorrectionVariable='AOD_FINE'
               ;CorrectionVariable='AOD_TOTAL'
               ;CorrectionVariable='Skewed-sine vs time'
            END
      ENDCASE

      N_mode=N_ELEMENTS(ws)
      Title=PlotLabel+' '+ViewModes[Mode]+ ' :: Corrected via '+CorrectionVariable

      ;;----- Subset XData and YData for current observation mode.
      ;;----- This defines XD and YD.
      XD=XData[ws] & YD=YData[ws] & YD_2=YData_2[ws]

      ;;-----
      ;;----- Initial plot with no data.
      ;;-----
      PRINT,'Initializing plot....'
      p=SCATTERPLOT( $
                    /NODATA, $
                    CURRENT=Current, $
                    LAYOUT=[NCols,NRows,PlotPos], $
                    MARGIN=Margin, $
                    XD, $
                    YD, $
                    AXIS_STYLE=1, $
                    Title=Title, $
                    XRange=xrange_plot, XTickFormat='LABEL_DATE', XTitle=XTitle, $
                    YRange=yrange_plot, YTitle=YTitle, $
                    _EXTRA=PlotExtras, $
                    FONT_SIZE=FontSize_1, FONT_STYLE=1 $
                   )
      p.title.font_size = FontSize_1*1.4

      ;;----- Bin individual collocations into about 100 bins with at least 10 soundings each.
      ;;----- This generates XB1, YB1.
      BINXY, XD,YD, XB1, YB1, BINS=100, MINCOUNT=10
      N_bins=N_ELEMENTS(XB1)
      PRINT,'After binning [XD, YD], determined '+NUM2STR(N_bins,0)+' bins.'



      ;;----- Overplot original dXCO2
      text_1a='$\Delta$XCO$_2$ (GOSAT v9 BC - OCO-2 v'+VerNum_oco+', w/AKC)'
      IF (PlotSingles) THEN BEGIN   
         PRINT,'Scatter plot data...'
         p=SCATTERPLOT(OVERPLOT=p,XD,YD)
         p.SYM_INCREMENT = SymIncrement
         p.SYM_COLOR = Color_1
         p.SYM_SIZE=SymSize
         p.SYMBOL=Symbol
      ENDIF

         ;;----- Linear fit to data to get trend.
         PRINT,'Fit data...'
         Fit=LINFIT(XD,YD,YFit=YFit)
         M_a1=Fit[1]*365. ;; per year
         Mean=Mean(YD)
         Sigma=STDDEV(YD)
         R=CORRELATE(XD,YD)
         Text_1b='Linear trend: '+NUM2STR(M_a1,3)+' ppm/year, Bias='+NUM2STR(Mean,2)+' $\pm$ '+NUM2STR(Sigma,2)+' ppm'
      IF (PlotSinglesFit) THEN BEGIN
         p=PLOT(OVERPLOT=p,XD,YFit,COLOR=Color_1,Thick=5)
      ENDIF

      ;;------ Overplot binned data.
      IF (PlotBins) THEN BEGIN   
         p=SCATTERPLOT(OVERPLOT=p,XB1,YB1)
         p.SYM_INCREMENT = 1
         p.SYM_COLOR=Color_1
         p.SYMBOL='o'
         p.Sym_Size=BinSymSize
         p.SYM_FILLED=0
         p.SYM_FILL_COLOR=Color_1
         Fit=LINFIT(XB1,YB1,YFit=YFit)
         M_a2=Fit[1]*365. ;; per year
      ENDIF

      ;;----- Bin individual collocations into about 100 bins with at least 10 soundings each.
      ;;----- This generates XB1, YB1.
      BINXY, XD,YD_2, XB2, YB2, BINS=100, MINCOUNT=10
   
      ;;----- Overplot prior and AK adjusted dXCO2
      text_2a='$\Delta$XCO$_2$ (GOSAT '+VerNum_gosat+' harmonized - OCO-2 v'+VerNum_oco+', w/AKC)'
      IF (PlotSingles) THEN BEGIN   
        PRINT,'Scatter plot data...'
        p=SCATTERPLOT(OVERPLOT=p,XD,YD_2)
        p.SYM_INCREMENT = SymIncrement
        p.SYM_COLOR = Color_2
        p.SYM_SIZE=SymSize
        p.SYMBOL=Symbol
      ENDIF

      IF (PlotSinglesFit) THEN BEGIN
         ;;----- Linear fit to data to get trend.
         PRINT,'Fit data...'
         Fit=LINFIT(XD,YD_2,YFit=YFit)
         M_b1=Fit[1]*365. ;; per year
         Mean=MEAN(YD_2)
         Sigma=STDDEV(YD_2)
         R=CORRELATE(XD,YD_2)
         Text_2b='Linear trend: '+NUM2STR(M_b1,3)+' ppm/year, Bias='+NUM2STR(Mean,2)+' $\pm$ '+NUM2STR(Sigma,2)+' ppm'
         p=PLOT(OVERPLOT=p,XD,YFit,COLOR=Color_2,Thick=5)
      ENDIF

      IF (PlotBins) THEN BEGIN
         ;;------ Overplot binned data.
         p=SCATTERPLOT(OVERPLOT=p,XB2,YB2)
         p.SYM_INCREMENT = 1
         p.SYM_COLOR=Color_2
         p.SYMBOL='d'
         p.Sym_Size=BinSymSize
         p.SYM_FILLED=1
         p.SYM_FILL_COLOR=Color_2
         Fit=LINFIT(XB2,YB2,YFit=YFit)
         M_b2=Fit[1]*365. ;; per year

      ENDIF

      ;;----- Overplot horizontal zero lines
      IF (0) THEN BEGIN
         p=PLOT(OVERPLOT=p,XRange,[0.,0.],Thick=2,Color=Black)
      ENDIF
               ;;-----
               ;;----- Write stats to plot.
               ;;-----
               IF (1) THEN BEGIN
                  PrintFmt='(A)'
                  FontSize=6 & FontStyle=1 ; 1=bold
                  Xt1=0.075
                  Xt2=0.68
                  Yt1=0.85
                  Yt2=0.35
                  DY=0.1
            
                  ;;----- Print fit type to plot
                  t=TEXT(TARGET=p,/RELATIVE,Xt1,Yt1,'N collocations='+NUM2STR(N_mode,0), FONT_SIZE=FontSize,FONT_COLOR='black',FONT_STYLE=FontStyle)
                  t=TEXT(TARGET=p,/RELATIVE,Xt1,Yt1-DY,'N bins='+NUM2STR(N_bins,0), FONT_SIZE=FontSize,FONT_COLOR='black',FONT_STYLE=FontStyle)

                  ;;----- Header "Prior to adjustment (collocations)".
                  t=TEXT(TARGET=p,/RELATIVE,Xt2,Yt1,'Prior to adjustment (collocations)', FONT_SIZE=FontSize,FONT_COLOR=Color_1,FONT_STYLE=FontStyle)
                  
                  ;;----- Single collocations prior to adjustment (YD).
                  StatsString='$\mu$='+NUM2STR(MEAN(YD),2) +' $\pm$ '+NUM2STR(STDDEV(YD),2)+' ppm, m='+NUM2STR(m_a1,3)+' ppm/yr'
                     
                  t=TEXT(TARGET=p,/RELATIVE,Xt2,Yt1-DY, StatsString, FONT_SIZE=FontSize,FONT_COLOR=Color_1,FONT_STYLE=FontStyle)

                  ;;----- Header "Prior to adjustment (binned)".
                  t=TEXT(TARGET=p,/RELATIVE,Xt2,Yt1-DY*2,'Prior to adjustment (binned)', FONT_SIZE=FontSize,FONT_COLOR=Color_1,FONT_STYLE=FontStyle)

                  ;;----- Binned data prior to adjustment (YB)
                  StatsString='$\mu$='+NUM2STR(MEAN(YB1),2) +' $\pm$ '+NUM2STR(STDDEV(YB1),2)+' ppm, m='+NUM2STR(m_a2,3) ; $
                  t=TEXT(TARGET=p,/RELATIVE,Xt2,Yt1-DY*3, StatsString, FONT_SIZE=FontSize,FONT_COLOR=Color_1,FONT_STYLE=FontStyle)


                  ;;----- Header "After adjustment (collocations)".
                  t=TEXT(TARGET=p,/RELATIVE,Xt2,Yt2,'After adjustment (collocations)', FONT_SIZE=FontSize,FONT_COLOR=Color_2,FONT_STYLE=FontStyle)

                  ;;----- Single collocations after adjustment (dxco2_2).
                  StatsString='$\mu$='+NUM2STR(MEAN(YD_2),2) +' $\pm$ '+NUM2STR(STDDEV(YD_2),2)+' ppm, m='+NUM2STR(m_b1,3) 
                  t=TEXT(TARGET=p,/RELATIVE,Xt2,Yt2-DY, StatsString, FONT_SIZE=FontSize,FONT_COLOR=Color_2,FONT_STYLE=FontStyle)

                  ;;----- Header "After adjustment (binned)".
                  t=TEXT(TARGET=p,/RELATIVE,Xt2,Yt2-DY*2,'After adjustment (binned)', FONT_SIZE=FontSize,FONT_COLOR=Color_2,FONT_STYLE=FontStyle)

                  ;;----- Binned data after adjustment (YB3)
                  StatsString='$\mu$='+NUM2STR(MEAN(YB2),2) +' $\pm$ '+NUM2STR(STDDEV(YB2),2)+' ppm, m='+NUM2STR(m_b2,3)
                  t=TEXT(TARGET=p,/RELATIVE,Xt2,Yt2-DY*3, StatsString, FONT_SIZE=FontSize,FONT_COLOR=Color_2,FONT_STYLE=FontStyle)


                ENDIF

      ;;STOP
   ;;----- End loop over GOSAT viewing modes.
   ENDFOR


   ;;----- Add the legend.
   IF (0) THEN BEGIN
      FontSize=12
      IF (PlotSingles) THEN t=TEXT(TARGET=p,0.17,0.82,DisplayText,FONT_SIZE=FontSize+2,FONT_COLOR='black')

      t=TEXT(TARGET=p,0.17,0.77,text_1a,FONT_SIZE=FontSize,FONT_COLOR=Color_1a)
      t=TEXT(TARGET=p,0.17,0.72,text_2a,FONT_SIZE=FontSize,FONT_COLOR=Color_2a)

      t=TEXT(TARGET=p,0.17,0.22,Text_1b,FONT_SIZE=FontSize,FONT_COLOR=Color_1a)
      t=TEXT(TARGET=p,0.17,0.17,Text_2b,FONT_SIZE=FontSize,FONT_COLOR=Color_2a)
   ENDIF

   ;;----- Save plot to output file.
   PRINT,'Saving plot to file ',OutPlotFile
   p.save,OutPlotFile

   STOP

ENDIF

;;-----
;;----- Sanity check plot of the calculated xco2 adjustment due to prior.
;;-----
IF (PlotAvgKer) THEN BEGIN

   IF (1) THEN BEGIN
      ViewModes=['All Observation Modes', 'Land-H','Land-M','Ocean-H'] & NumViewModes=N_ELEMENTS(ViewModes)
   ENDIF

   ;;----- Plot settings.
   Symbol='Circle'
   Symbol_2='Diamond'
   Symbol_3='Asterisk'
   Col_1a='deep pink'
   Col_1b='medium violet red'
   Col_2a='blue'
   Col_2b='medium blue'
   Col_3a='olive'
   Col_3b='olive drab'
   Col_4a='orange'
   Col_4b='gold'
   Col_5a='dark red'
   Col_5b='maroon'
   Col_6a='powder blue'
   Col_6b='light blue'
   
   ;;----- Plot settings.
   Margin=[0.2,0.2,0.1,0.1]
   NCols=2 & NRows=2
   FontSize_1=10
   FitThick=2
   SymSize=0.75
   DisplayText='Points represent monthly binned median values.'

   ;;----- Index by observation modes.
   landh = WHERE( (Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '72') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)
   landm = WHERE( (Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '77') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)
   oceanh = WHERE( (~Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '72') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)

   ;;----- Define X and Y Data.
   XData=gosat.jd
   YData_1=xco2_prior_adj_a
   YData_2=xco2_prior_adj_b
   YData_3=xco2_prior_adj_a+xco2_prior_adj_b


   Dummy = LABEL_DATE(DATE_FORMAT='%N/%Z')
   XTitle='Time [MM/YY]'
   XRange=[MIN(XData)-60,MAX(XData)+60]
   YRange=[-0.75,1.25]

   ;;----- Output plot filename.
   DataDir='ak_plots/'
   FILE_MKDIR,OutDir+DataDir
   OutPlotFile=OutDir+DataDir+'xco2_ak_correction.png'

   ;;-----
   ;;----- Loop over GOSAT viewing modes.
   ;;-----
   FOR Mode=0,NumViewModes-1 DO BEGIN

      CASE Mode OF
         0: BEGIN
               PlotPos=1
               CURRENT=0
               PlotLabel='(a)'
               XTitle=XTitle
               YTitle='XCO$_2$ correction [ppm]'
               ws=[landh,landm,oceanh]
               Color_1='black'
               Color_2='gray'
            END
         1: BEGIN
               PlotPos=2
               CURRENT=1
               PlotLabel='(b)'
               XTitle=XTitle
               YTitle='XCO$_2$ correction [ppm]'
               ws=landh
               Color_1='green'
               Color_2='lime'
            END
         2: BEGIN
               PlotPos=3
               CURRENT=1
               PlotLabel='(c)'
               XTitle=XTitle
               YTitle='XCO$_2$ correction [ppm]'
               ws=landm
               Color_1='orange'
               Color_2='red'
            END
         3: BEGIN
               PlotPos=4
               CURRENT=1
               PlotLabel='(d)'
               XTitle=XTitle
               YTitle='XCO$_2$ correction [ppm]'
               ws=oceanh
               Color_1='blue'
               Color_2='aqua'
            END
      ENDCASE

      N_mode=N_ELEMENTS(ws)
      Title=PlotLabel+' '+ViewModes[Mode]

      ;;----- Subset XData and YData for current observation mode.
      ;;----- This defines XD and YD.
      XD=XData[ws] & YD_1=YData_1[ws] & YD_2=YData_2[ws] & YD_3=YData_3[ws]

      ;;----- Bin X/Y Data.
      BINXY,XD,YD_1,XB1,YB1,DX=30,MEDIAN=1,MINCOUNT=100
      BINXY,XD,YD_2,XB2,YB2,DX=30,MEDIAN=1,MINCOUNT=100
      BINXY,XD,YD_3,XB3,YB3,DX=30,MEDIAN=1,MINCOUNT=100
      ;STOP

      ;;-----
      ;;----- Initial plot with no data.
      ;;-----
      PRINT,'Initializing plot....'
      p=SCATTERPLOT( $
                    /NODATA, $
                    CURRENT=Current, $
                    LAYOUT=[NCols,NRows,PlotPos], $
                    MARGIN=Margin, $
                    XB1, $
                    YB1, $
                    Title=Title, $
                    XTitle=XTitle, XStyle=1, XRange=XRange, XTickFormat='LABEL_DATE', XMajor=4, XMinor=0, $
                    YTitle=YTitle, YRange=YRange, $
                    FONT_SIZE=FontSize_1, FONT_STYLE=1 $
                   )

      ;;----- Overplot xco2_prior_adj_a
      IF (1) THEN BEGIN   
         text_1a='Term 1: CO$_2$ prior adj. via Eq. A10 [Wunch, ACP, 2010]'
         text_1aa=' (equivalent to Eq. 3 [Taylor, ESSD, 2023])'
         PRINT,'Scatter plot data (a)...'
         p=SCATTERPLOT(OVERPLOT=p,XB1,YB1)
         p.SYM_INCREMENT = 1
         p.SYM_COLOR = Col_1a
         p.SYM_SIZE=SymSize
         p.SYMBOL=Symbol
         p.SYM_FILLED=1
         p.SYM_FILL_COLOR=Col_1a

         ;;----- Linear fit to data to get trend.
         PRINT,'Fit data (a)...'
         Fit=LINFIT(XB1,YB1,YFit=YFit1)
         Slope=Fit[1]*365. ;; per year
         Sigma=STDDEV(YB1)
         R=CORRELATE(XB1,YB1)
         Text_1b='Linear trend: '+NUM2STR(Slope,3)+' ppm/year, 1-$\sigma$='+NUM2STR(Sigma,2)+' ppm'
         p=PLOT(OVERPLOT=p,XB1,YFit1,COLOR=Col_1a,Thick=FitThick)

      ENDIF

      ;;----- Overplot xco2_prior_adj_b   
      IF (1) THEN BEGIN   
         text_2a='Term 2: AK "smoothing" term'
         PRINT,'Scatter plot data (b)...'
         p=SCATTERPLOT(OVERPLOT=p,XB2,YB2)
         p.SYM_INCREMENT = 1
         p.SYM_COLOR = Col_2a
         p.SYM_SIZE=SymSize
         p.SYMBOL=Symbol_2
         p.SYM_FILLED=1
         p.SYM_FILL_COLOR=Col_2a

         ;;----- Linear fit to data to get trend.
         PRINT,'Fit data (b)...'
         Fit=LINFIT(XB1,YB2,YFit=YFit2)
         Slope=Fit[1]*365. ;; per year
         R=CORRELATE(XB2,YB2)
         Sigma=STDDEV(YB2)
         Text_2b='Linear trend: '+NUM2STR(Slope,3)+' ppm/year, 1-$\sigma$='+NUM2STR(Sigma,2)+' ppm'
         p=PLOT(OVERPLOT=p,XB2,YFit2,COLOR=Col_2a,Thick=FitThick)

      ENDIF

      ;;----- Overplot sum of terms  
      IF (1) THEN BEGIN   
         text_3a='Sum of terms = T1 + T2'
         PRINT,'Scatter plot data (sum)...'
         p=SCATTERPLOT(OVERPLOT=p,XB3,YB3)
         p.SYM_INCREMENT = 1
         p.SYM_COLOR = Col_3a
         p.SYM_SIZE=SymSize
         p.SYMBOL=Symbol_3
         p.SYM_FILLED=1
         p.SYM_FILL_COLOR=Col_3a

         ;;----- Linear fit to data to get trend.
         PRINT,'Fit data (sum)...'
         Fit=LINFIT(XB3,YB3,YFit=YFit3)
         Slope=Fit[1]*365. ;; per year
         R=CORRELATE(XB3,YB3)
         Sigma=STDDEV(YB3)
         Text_3b='Linear trend: '+NUM2STR(Slope,3)+' ppm/year, 1-$\sigma$='+NUM2STR(Sigma,2)+' ppm'
         p=PLOT(OVERPLOT=p,XB3,YFit3,COLOR=Col_3a,Thick=FitThick)

      ENDIF

      ;;----- Overplot horizontal zero lines
      p=PLOT(OVERPLOT=p,XRange,[0.,0.],Thick=2,Color=Black)

      ;;----- Add the legend.
      IF (1) THEN BEGIN   
         FontSize=4
         FontSize_2=6
         tx0=0.1
         ty0=0.87 & ty1=0.2
         tdy=0.04 & tdy2=0.06
         t=TEXT(TARGET=p,/RELATIVE,tx0,ty0,DisplayText,FONT_SIZE=FontSize+2,FONT_COLOR='black')

         t=TEXT(TARGET=p,/RELATIVE,tx0,ty0-(tdy),text_1a,FONT_SIZE=FontSize,FONT_COLOR=Col_1a)
         t=TEXT(TARGET=p,/RELATIVE,tx0,ty0-(tdy*2),text_1aa,FONT_SIZE=FontSize,FONT_COLOR=Col_1a)
         t=TEXT(TARGET=p,/RELATIVE,tx0,ty0-(tdy*3),xco2_prior_adj_a_text,FONT_SIZE=FontSize,FONT_COLOR=Col_1a)
         t=TEXT(TARGET=p,/RELATIVE,tx0,ty0-(tdy*4),text_2a,FONT_SIZE=FontSize,FONT_COLOR=Col_2a)
         t=TEXT(TARGET=p,/RELATIVE,tx0,ty0-(tdy*5),xco2_prior_adj_b_text,FONT_SIZE=FontSize,FONT_COLOR=Col_2a)
         t=TEXT(TARGET=p,/RELATIVE,tx0,ty0-(tdy*6),text_3a,FONT_SIZE=FontSize,FONT_COLOR=Col_3a)

         t=TEXT(TARGET=p,/RELATIVE,tx0,ty1,Text_1b,FONT_SIZE=FontSize_2,FONT_COLOR=Col_1a)
         t=TEXT(TARGET=p,/RELATIVE,tx0,ty1-tdy2,Text_2b,FONT_SIZE=FontSize_2,FONT_COLOR=Col_2a)
         t=TEXT(TARGET=p,/RELATIVE,tx0,ty1-tdy2*2,Text_3b,FONT_SIZE=FontSize_2,FONT_COLOR=Col_3a)
      ENDIF

      ;;----- End loop over observation modes.
   ENDFOR

   ;;----- Save plot to output file.
   PRINT,'Saving plot to file ',OutPlotFile
   p.save,OutPlotFile
   STOP

ENDIF







;;---- Additional filter on OCO-2 data for variability over the scene.
;;---- Makes sense to move this to early in the code?...
PRINT, 'Determining overall (per match) mask.'
match_mask =  (oco_stddev_xco2 LE oco_max_xco2_variability) ;AND (mean_dist_g LE max_dist) AND (gosat_qf[w])
Mask_match=WHERE(match_mask, NMask_variability)
PRINT,'NMask_variability=',NMask_variability
PRINT,'Percent of total retained by OCO-2 XCO2 variability check=',FLOAT(NMask_variability)/FLOAT(N_ELEMENTS(match_mask))*100. 
;;STOP










;;-----
;;----- Build mask for current observation mode.
;;-----
PRINT, 'Developing mask for mode ' + view_mode
MGAIN=77 & HGAIN=72
;MGAIN='M' & HGAIN='H'
CASE view_mode OF
   'landH'  : type_mask = gosat.land AND (gosat.gain EQ HGAIN)
   'landM'  : type_mask = gosat.land AND (gosat.gain EQ MGAIN)
   'land'   : type_mask = gosat.land
   'oceanH' : type_mask = ~gosat.land AND (gosat.gain EQ HGAIN)
   'allH'   : type_mask = gosat.gain eq HGAIN
   'all'    : type_mask = (~gosat.land AND gosat.gain eq HGAIN) OR gosat.land
ENDCASE
Mask_type=WHERE(type_mask, NMask_surface)
PRINT,'Percent of total retained by surface/observation mask =',FLOAT(NMask_surface)/FLOAT(N_ELEMENTS(type_mask))*100. 

;;----- Mask for latitude range.
lat_mask=INRANGE(gosat.lat,LatRange)
Mask_last=WHERE(lat_mask, NMask_lat)
PRINT,'Percent of total retained by latitude mask =',FLOAT(NMask_lat)/FLOAT(N_ELEMENTS(lat_mask))*100. 

;;----- Combine submasks
main_mask = match_mask AND type_mask AND lat_mask
g = WHERE(main_mask, nsound)
wg=w1[g]
PRINT, 'Using ' + sc(nsound) + ' matched soundings.'
PRINT,'Percent of total retained by main mask =',FLOAT(NSound)/FLOAT(N_ELEMENTS(main_mask))*100. 
;STOP

;;-----
;;----- Extract the needed data variables.
;;-----
matchg = Match_metrics[g]
;mean_dist_g=mean_ddist[g]
;mean_dtime_g=mean_dtime[g]
oco_meang = oco_mean_struct[g]
dxco2_bc   = gosat_metrics[g].xco2_bc - oco_meang.xco2 
dxco2_harm = gosat_metrics[g].xco2_harm - oco_meang.xco2
IF (Correct_AK) THEN BEGIN
   IF (two_term_ak) THEN BEGIN
      dxco2_bc_ak   = dxco2_bc + xco2_prior_adj_a[g] + xco2_prior_adj_b[g]
      dxco2_harm_ak = dxco2_harm + xco2_prior_adj_a[g] + xco2_prior_adj_b[g] 
   ENDIF ELSE BEGIN
      dxco2_bc_ak   = dxco2_bc + xco2_prior_adj_b[g]
      dxco2_harm_ak = dxco2_harm + xco2_prior_adj_b[g] 
   ENDELSE

ENDIF
;STOP

tdiff = mean_dtime[g]
ddiff = mean_ddist[g]
lat=gosat_metrics[g].lat
lon=gosat_metrics[g].lon
year = fix(strmid(sc(gosat_metrics[g].id),0,4))
month = fix(strmid(sc(gosat_metrics[g].id),4,2))
season = (month mod 12)/3
;;;;l = where(gosat.land, comp=o)
julian_day=gosat_metrics[g].jd
land_flag=gosat_metrics[g].land

;;----- If skip_calc is set or code has not been reset.
;skip_calc_:

;;-----
;;----- Write Summary file of collocated GOSAT/OCO2 soundings.
;;-----
;STOP
IF (PlotCollocationSummary) THEN BEGIN

   OutDir=OutDirBase+'gosat_'+VerNum_gosat+'_vs_oco2_v'+VerNum_oco+'/' & FILE_MKDIR,OutDir
   SummaryFileName=OutDir+'collocated_summary_gosat_'+VerNum_gosat+'_vs_oco2_v'+VerNum_oco+'.txt'
   SummaryHeader=['SoundingID_gosat', $
                  +'Latitude_gosat', $
                  +'Longitude_gosat', $
                  +'Gain Setting', $
                  +'Land Fraction', $
                  +'Solar Zenith Angle', $
                  +'Sensor Zenith Angle', $
                  +'Mean Separation (km)', $
                  +'Mean Delta Time (hr)', $
                  +'SoundingID_oco2 (10 closest soundings)']
                  ;;01234567890123456789012345678901234567
                  ;;0         1         2         3   
   OPENW,SummaryLUN,SummaryFileName,/GET_LUN,WIDTH=250
   PRINTF,SummaryLUN,STRTRIM(SummaryHeader,2),FORMAT='( 9(A-20), 1(A-40) )'

   SummaryFormat='( 1(I14,6X), 2(F8.3,12X), 1(A1,19X), 1(F5.1,15X), 4(F8.3,12X), 10(I16,4X) )'
   FOR III=0L,NSound-1L DO BEGIN
      IF (0) THEN BEGIN
         PRINTF,SummaryLUN, $
         matchg[III].gosat.sounding_ID, $
         matchg[III].gosat.latitude, matchg[III].gosat.longitude, $
         matchg[III].gosat.sounding.gain, $
         ;;matchg[III].gosat.sounding.land_fraction, $
         matchg[III].gosat.solar_zenith_angle, matchg[III].gosat.sensor_zenith_angle, $
         ddiff[III], $
         tdiff[III], $
         matchg[III].oco2[0:9].sounding_id, $
         FORMAT=SummaryFormat
      ENDIF
   ENDFOR
   FREE_LUN,SummaryLUN

   ;;----- Index by observation modes.
   landh = WHERE( (Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '72') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)
   landm = WHERE( (Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '77') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)
   oceanh = WHERE( (~Gosat.land) $
                   AND (STRTRIM(Gosat.gain,2) EQ '72') $
                   ;AND (INRANGE(Gosat.lat,CurrentLatRange)) $
                   ;AND (INRANGE(lon,CurrentLonRange)) $
                   , N_mode)

   ;;----- Output plot filename.
   DataDir='collocation_plots/'
   FILE_MKDIR,OutDir+DataDir
   OutPlotFile=OutDir+DataDir+'gosat_oco2_collocations_dtime.png'

   ;;-----
   ;;----- Plot histograms of collocation time.
   ;;-----
   Title=' GOSAT vs OCO-2 collocations
   ;;Title=SubTitle+STRING(CurrentYear,FORMAT='(I4)')+', '+SeasList[s]
   XTitle='Mean collocation time difference (minutes)'
   Thick=3

   YD_hist_time=mean_dtime*60.
   binsize=10. ;; time in minutes

   pdf_yd=HISTOGRAM(YD_hist_time, LOCATIONS=xbin_hist_time, BINSIZE=binsize)
   p=PLOT(xbin_hist_time, pdf_yd, /STAIRSTEP, THICK=5, COLOR='black', Title=Title, xtit=XTitle, ytit='Number of Collocations',XRange=[-150,150],NAME='All obs modes')
   pdf_yd=HISTOGRAM(YD_hist_time[landH], LOCATIONS=xbin_hist_time, BINSIZE=binsize)

   p=PLOT(OVERPLOT=p, $
          /STAIRSTEP, $
           xbin_hist_time, pdf_yd, $
          THICK=Thick, COLOR='green', $
          NAME='Land-H')
   pdf_yd=HISTOGRAM(YD_hist_time[landM], LOCATIONS=xbin_hist_time, BINSIZE=binsize)
   p=PLOT(OVERPLOT=p, $
          /STAIRSTEP, $
           xbin_hist_time, pdf_yd, $
          THICK=Thick, COLOR='red', $
          NAME='Land-M')
   pdf_yd=HISTOGRAM(YD_hist_time[oceanH], LOCATIONS=xbin_hist_time, BINSIZE=binsize)
   p=PLOT(OVERPLOT=p, $
          /STAIRSTEP, $
           xbin_hist_time, pdf_yd, $
          THICK=Thick, COLOR='blue', $
          NAME='Ocean-H')
   ;VLINE,MEAN(YData),LineStyle=2,Thick=5

   ;;----- Save plot to output file.
   PRINT,'Saving plot to file ',OutPlotFile
   p.save,OutPlotFile



   ;;----- Output plot filename.
   DataDir='collocation_plots/'
   FILE_MKDIR,OutDir+DataDir
   OutPlotFile=OutDir+DataDir+'gosat_oco2_collocations_ddist.png'

   ;;-----
   ;;----- Plot histograms of collocation time.
   ;;-----
   Title=' GOSAT vs OCO-2 collocations
   ;;Title=SubTitle+STRING(CurrentYear,FORMAT='(I4)')+', '+SeasList[s]
   XTitle='Mean collocation distance difference (km)'
   Thick=3

   YD_hist_dist=mean_ddist
   binsize=10. ;; distance in km

   pdf_yd=HISTOGRAM(YD_hist_dist, LOCATIONS=xbin_hist_dist, BINSIZE=binsize)
   p=PLOT(xbin_hist_dist, pdf_yd, /STAIRSTEP, THICK=5, COLOR='black', Title=Title, xtit=XTitle, ytit='Number of Collocations',XRange=[0,300],NAME='All obs modes')
   pdf_yd=HISTOGRAM(YD_hist_dist[landH], LOCATIONS=xbin_hist_dist, BINSIZE=binsize)

   p=PLOT(OVERPLOT=p, $
          /STAIRSTEP, $
           xbin_hist_dist, pdf_yd, $
          THICK=Thick, COLOR='green', $
          NAME='Land-H')
   pdf_yd=HISTOGRAM(YD_hist_dist[landM], LOCATIONS=xbin_hist_dist, BINSIZE=binsize)
   p=PLOT(OVERPLOT=p, $
          /STAIRSTEP, $
           xbin_hist_dist, pdf_yd, $
          THICK=Thick, COLOR='red', $
          NAME='Land-M')
   pdf_yd=HISTOGRAM(YD_hist_dist[oceanH], LOCATIONS=xbin_hist_dist, BINSIZE=binsize)
   p=PLOT(OVERPLOT=p, $
          /STAIRSTEP, $
           xbin_hist_dist, pdf_yd, $
          THICK=Thick, COLOR='blue', $
          NAME='Ocean-H')
   ;VLINE,MEAN(YData),LineStyle=2,Thick=5

   ;;----- Save plot to output file.
   PRINT,'Saving plot to file ',OutPlotFile
   p.save,OutPlotFile

   ;;-----
   ;;----- TO-DO. Add histogram of difference in surface elevation.
   ;;----- Probably need to read in altitude variable either in Save file generation or here.
   ;;----- Plot histograms of collocation distance.
   ;;-----

ENDIF
;STOP

;;---------------------------
;;----- Grid data.
;;---------------------------
IF (GridData) THEN BEGIN


   BinMinCount=10
   IF plot_type EQ 'seasonal' THEN BEGIN
      nx=2 & ny=2 
      seaslist=['DJF','MAM','JJA','SON']
      Label_1='year_'
      Label_2='season_'
   ENDIF else begin
      nx=1 & ny=1
      seaslist=['annual']
      Label='annual_'
   endelse
   ns=nx*ny

   ;;----- Prepare the ZData array to be gridded.
   IF (VerNum_GOSAT EQ 'v9b') THEN BEGIN
         ZData=[[dxco2_bc],[gosat_xco2_adj_w[g]]]
         var_list=[var_list,'xco2_adjust']
   ENDIF
   IF (STRMID(VerNum_GOSAT,0,4) EQ 'v9c2') THEN BEGIN
         ZData=[ [dxco2_bc], [gosat_xco2_adj_w[g]], [dxco2_harm_ak] ]
         var_list=[var_list,'xco2_adjust','dxco2_harm_ak']
   ENDIF
   IF (STRMID(VerNum_GOSAT,0,4) EQ 'v9c3') THEN BEGIN
      ZData=[ [dxco2_bc], [dxco2_bc_ak], [gosat_metrics[g].xco2_adj], [dxco2_harm_ak] ]
      var_list=['dxco2_bc','dxco2_bc_ak','xco2_adjust','dxco2_harm_ak']
      ;STOP
   ENDIF
   IF (STRMID(VerNum_GOSAT,0,5) EQ 'v9_ml') THEN BEGIN
      IF (correct_ak) THEN BEGIN
         ZData=[ [dxco2_bc], [dxco2_bc_ak], [gosat_metrics[g].xco2_adj], [dxco2_harm_ak] ]
         var_list=['dxco2_bc','dxco2_bc_ak','xco2_adjust','dxco2_harm_ak']
      ENDIF ELSE BEGIN
         ZData=[ [dxco2_bc], [gosat_metrics[g].xco2_adj], [dxco2_harm] ]
         var_list=['dxco2_bc','xco2_adjust','dxco2_harm']
      ENDELSE
      STOP
   ENDIF
   ;;STOP

   ;;-----
   ;;----- Loop over number of years (or once for NOT Independent_Years).
   ;;-----
   IF (independent_years) THEN NYears=N_ELEMENTS(YearList) ELSE NYears=1

   ;;----- Create storage for summary statistics.
   ss={Year:0,Season:'', n_ss:0, mu_ss:0.0, sd_ss:0.0, n_grid:0, mu_grid:0.0, sd_grid:0.0}
   SummaryStats=REPLICATE(ss, [NYears,NS])

   FOR YearIndex=0,NYears-1 DO BEGIN      

      ;;-----
      ;;----- Loop over number of seasons (or once for annual).
      ;;-----
      FOR s=0,ns-1 DO BEGIN

         IF (independent_years) THEN CurrentYear=YearList[YearIndex] ELSE CurrentYear=YearList

         ;;----- Search for months in the season.
         IF (plot_type EQ 'seasonal') THEN BEGIN
            CASE SeasList[s] OF
                  'DJF' : ws=WHERE( ((Month EQ 12) AND INRANGE(Year,[MIN(CurrentYear-1),MAX(CurrentYear-1)] ) ) $
                                  OR ( INRANGE(Month, [1,2])  AND INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]) ), N_season )
                  'MAM' : ws=WHERE( INRANGE(Month, [3,5])   AND INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]), N_season) 
                  'JJA' : ws=WHERE( INRANGE(Month, [6,8])   AND INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]), N_season) 
                  'SON' : ws=WHERE( INRANGE(Month, [9,11]) AND INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]), N_season) 
            ENDCASE            

         ENDIF ELSE BEGIN
            ws=WHERE( INRANGE(Year, [MIN(CurrentYear),MAX(CurrentYear)]) AND (Season LE 2) , N_season)
         ENDELSE
         PRINT,'N Elements in Season/Year='+SC(N_season)
         ;STOP

         ;;-----
         ;;----- Grid data density for spatial maps.
         ;;-----
         IF (GridMap) THEN BEGIN

            FillFloat=-999.9 & FillInteger=-999
            IF (independent_years) $
               THEN OutFileBase = $
                     OutDir+'gosat_'+VerNum_GOSAT+'_vs_oco2_v'+VerNum_OCO+ $
                     '_lite-vars-vs-lat-lon_'+STRING(CurrentYear,FORMAT='(I4)')+'_'+SeasList[s]+'_'+view_mode+'-'+Region $
               ELSE OutFileBase = $
                     OutDir+'gosat_'+VerNum_GOSAT+'_vs_oco2_v'+VerNum_OCO+ $
                     '_lite-vars-vs-lat-lon_'+STRING(MIN(CurrentYear),FORMAT='(I4)') $
                     +'-'+STRING(MAX(CurrentYear),FORMAT='(I4)')+'_'+SeasList[s]+'_'+view_mode+'-'+Region

            OutPdfFile=OutFileBase+'.pdf'
            OutNCFile=OutFileBase+'.nc'
            OutInfo=MAKE_ARRAY(15,/STRING,VALUE='')
            OutInfo[0]=OutPdfFile
            OutInfo[1]='lite vars'
            OutInfo[2]=' GOSAT (v' + VerNum_GOSAT +') - OCO2 (v' + VerNum_OCO +') :: '+STRING(CurrentYear,FORMAT='(I4)')+' :: '+SeasList[s]
            OutInfo[3]='N='+STRING(FLOAT(N_season)/1000.,FORMAT='(F6.1)')+'k (SS)'
            OutInfo[4]=''
            OutInfo[10]='plasma'
            OutInfo[11]='Longitude'
            OutInfo[12]='Latitude'
            OutInfo[13]=STRING(FillFloat,FORMAT='(F6.1)')
            OutInfo[14]=STRING(FillInteger,FORMAT='(I4)')
            ;STOP
            
            ;;----- Subset the ZData to 'ws' elements.
            ;NZ=N_ELEMENTS(ZData[0,*])
            ;ZData_sub=MAKE_ARRAY([N_season,NZ], /FLOAT)
            ;FOR ZZZ=0,NZ-1 DO BEGIN
            ;   ZT=ZData[ws,ZZZ] 
            ;   ZData_sub[*,ZZZ]=ZT
            ;ENDFOR
            ;ZData_sub

            NX=LONG(360./DeltaLon_1)
            NY=LONG((LatRange[1]-LatRange[0])/DeltaLat_1)
            ;STOP
            PRINT,'GRID_DATA_WRITE_NETCDF_MULTI for spatial map of collocated sounding density...'
            GRID_DATA_WRITE_NETCDF_MULTI, $
                  Lon[ws], $
                  Lat[ws], $
                  ZData[ws,*], $  
                  var_list, $
                  LWI=0, DataDensity=0, $ 
                  'XVar',DeltaLon_1, NX, [-180,180],'lon','longitude (degrees)', $
                  'YVar',DeltaLat_1, NY, LatRange,'lat','latitude (degrees)', $
                  'NVar','N','Number of soundings per lat/lon grid box', $
                  OutPdfFile,OutNcFile,OutInfo, $
                  BinMinCount=BinMinCount, $
                  FillFloat=FillFloat, $
                  FillInteger=FillInteger, $
                  Verbose=0

            ;STOP
         ENDIF


         ;;-----
         ;;----- Grid data in time/lat/dxco2.
         ;;-----
         IF (GridHovmoller) THEN BEGIN
            ;OutDir=OutDirBase+'/gosat_l2fp_hovmoller/' & FILE_MKDIR,OutDir

            DSS=julian_day-MIN(julian_day)+1L & MinDSS=MIN(DSS) & MaxDSS=MAX(DSS)
            CALDAT,julian_day,mon2,day2,year2,hour2,min2,sec2
            StartDateString=STRING(Year2[0],FORMAT='(I4)')+'-'+STRING(Mon2[0],FORMAT='(I2.2)')+'-'+STRING(Day2[0],FORMAT='(I2.2)')
            PRINT,'StartDateString=',StartDateString

            NLoc=N_ELEMENTS(Location)
            FOR Loc=0,NLoc-1 DO BEGIN

               CurrentLatRange=Location[Loc].LatBin
               CurrentLonRange=Location[Loc].LonBin
               LocString=Location[Loc].Name

               Title=View_Mode+ ' :: '+LocString  

               w_hov=WHERE( INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]) $
                          AND (INRANGE(lat,CurrentLatRange)) $
                          AND (INRANGE(lon,CurrentLonRange)), NHov)
               NumSound=NHov
               Python_XLabel='Date (YYYY-MM)'
               Python_YLabel='Latitude'
               NDaysRange=[0,MAX(DSS[w_hov])]
               DX=HovInc_Day & NX=ROUND( (FLOAT(MAX(DSS))-FLOAT(MIN(DSS)) ) / FLOAT(DX) ) & PRINT,'NX=',NX 
               DY=HovInc_Lat & NY=ROUND(180./FLOAT(DY)) & PRINT,'NY=',NY
               HovText=STRING(HovInc_Day,FORMAT='(I2)') +'day_'+ STRING(HovInc_Lat,FORMAT='(I2)')+'lat'

               ;;----- Subset the ZData to w_hov
               ;z1=zdata[w_hov,0] & z2=zdata[w_hov,1] & zdata=[[z1],[z2]]

               OutFileBase = $
                     OutDir+'gosat_'+VerNum_GOSAT+'_vs_oco2_v'+VerNum_OCO+ $
                     '_lite-vars-vs-time-lat_' $
                     +HovText+'_' $
                     +STRMID(STRTRIM(STRING(ID_Start),2),0,8) $
                     +'-'+STRMID(STRTRIM(STRING(ID_Stop),2),0,8)+'_'+SeasList[s]+'_'+view_mode+'-'+LocString
               OutPdfFile=OutFileBase+'.pdf'
               OutNCFile=OutFileBase+'.nc'

               OutInfo=MAKE_ARRAY(15,/STRING,VALUE='')
               OutInfo[0]=OutPdfFile
               OutInfo[1]='$\Delta$XCO$_2$ [ppm]'
               OutInfo[2]=VerNum_OCO
               OutInfo[3]='N='+STRING(FLOAT(nhov)/1000L,FORMAT='(I4)')+'k (SS)'
               OutInfo[4]='G='+num2str(float(nhov)/float(n_elements(year))*100.,1)+'%'
               OutInfo[10]='coolwarm'
               OutInfo[11]=Python_XLabel
               OutInfo[12]=Python_YLabel
               PRINT,'GRID_DATA_WRITE_NETCDF_MULTI for dxco2 Hovmoller for region '+LocString
               ;STOP
               GRID_DATA_WRITE_NETCDF_MULTI, $
                  DSS[w_hov], $
                  Lat[w_hov], $
                  ZData[w_hov,*], $ ;; inputs
                     var_list, $
                     ;X_Array,Y_Array,Z_Binned, $ ;; output gridded X/Y/Z data
                    'XVar',DX, NX, NDaysRange,'days','days since start', $
                    'YVar',DY, NY, LatRange,'lat','latitude (degrees)', $
                    'NVar','N','Number of soundings per lat/lon grid box', $
                    ;'ZVar','ppm','Number of soundings per time/lat bin', $
                    OutPdfFile,OutNcFile,OutInfo,BinMinCount=BinMinCount, $
                    FillFloat=FillFloat,FillInteger=FillInteger, $
                    Verbose=0


                  ;;----- End loop over number of regions.
               ENDFOR
               DSS=!NULL

            ENDIF
            ;;STOP



         ;;-----
         ;;----- Plot timeseries.
         ;;-----
         IF (PlotTimeseries) THEN BEGIN
            PLOT_XCO2_TIMESERIES, $
               plot_type,Season,gosat_land,MatchG,lat,lon,wg,s, $
               Julian_Day, gosat_xco2_bc,oco_meang, $
               OutDir,Var,Label,SeasList,Vernum_oco,Location=Location
         ENDIF




         ;;----- End loop over seasons.
      ENDFOR

         ;;----- End loop over years.
      ENDFOR

      ;;----- Write stats to ascii file.
      IF (WriteStats) THEN BEGIN
         ViewModes=['Land-H','Land-M'] & NumViewModes=N_ELEMENTS(ViewModes)

         FOR Loc=0,NLocations-1 DO BEGIN
            CurrentLatRange=Location[Loc].LatBin
            CurrentLonRange=Location[Loc].LonBin
            LocString=Location[Loc].Name
            FOR Mode=0,NumViewModes-1 DO BEGIN
               Title=ViewModes[Mode]+ ' :: '+LocString  

               FOR s=0,ns-1 DO BEGIN
                  ;;----- Search for various surface types.
                  IF ( ViewModes[Mode] EQ 'Land-H') $
                  THEN ws = WHERE( (Season EQ s) AND (gosat_land) AND (MatchG.gosat.Sounding.gain EQ 'H') $
                                AND (INRANGE(lat,CurrentLatRange)) $
                                AND (INRANGE(lon,CurrentLonRange)) $
                                , N_season) 
                  IF ( ViewModes[Mode] EQ 'Land-M') $
                  THEN ws = WHERE( (Season EQ s) AND (gosat_land) AND (MatchG.gosat.Sounding.gain EQ 'M') $
                               AND (INRANGE(lat,CurrentLatRange)) $
                               AND (INRANGE(lon,CurrentLonRange)) $
                                , N_season) 
                  IF ( ViewModes[Mode] EQ 'Ocean' ) $
                  THEN ws = WHERE( (Season EQ s) AND (INRANGE(lat,CurrentLatRange)) AND (INRANGE(lon,CurrentLonRange)) AND (~gosat_land), N_season) 
                  
                  IF (N_season GT 2) THEN BEGIN
                     STATS_TET,gosat_xco2_bc[wg[ws]] - oco_meang[ws].xco2,dxco2_stats,VERBOSE=0
                     PrintFormat='(A20,A2,A8,A2,A4,A2,I5,A2,F5.2,A2,F4.2,A3)'
                     PRINT,LocString, ' &', ViewModes[Mode], ' &', SeasList[s], ' &', dxco2_stats.n,' &',dxco2_stats.mean,' &',dxco2_stats.stddev,' \\', FORMAT=PrintFormat
                  ENDIF ELSE BEGIN
                     PRINT,LocString, ' &', ViewModes[Mode], ' &', SeasList[s], ' &', 0,' &',-999.,' &',-999.,' \\', FORMAT=PrintFormat
                  ENDELSE
               ENDFOR
            ENDFOR
         ENDFOR
         STOP
      ENDIF


    
   IF 0 then l2_filter_plot, gosat_sza[o], xdiff[o], psy=csym(17), nbin=30, $
    xr=[0,80], yr=[-1.5, 1.5], xtit='GOSAT solar zenith angle [deg]', $
    ytit=Title,min=30,sd_p=0,gray=200

   ;;----- Complete plotting routines.
ENDIF

STOP
END
    
