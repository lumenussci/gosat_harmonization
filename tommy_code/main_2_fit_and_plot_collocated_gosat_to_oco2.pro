;;PRO MAIN_2_FIT_AND_PLOT_COLLOCATED_GOSAT_TO_OCO2

;;----- This procedure borrowed heavily from CWO's /home/codell/idlprogs/aos/models/OCO/l2/b9_acos/compare_acos_oco2_b9.pro
;;----- Modifications by T.E.T Autumn, 2020 to adapt and generate analysis/plots for the GOSAT v9 data paper.
;;----- Modifications by T.E.T late Winter, 2021 to adapt to Tropical IAV work.
;;-----                        Following the logic of the OCO-3 vs OCO-2 v10 collocation work, why not apply the methodology
;;-----                        to GOSAT v9 versus OCO-2 v10? So use the collocations to calculate dxco2 (already done),
;;-----                        troll for some explanatory variable (geometry?).
;;----- Mods by T.E.T. Summer 2022. Broke away from <main_gosat_vs_oco2.pro> to focus on plotting dxco2 vs. BC variables.
;;-----                        The previous code was becoming too unweildy.
;;----- Mods by TET Winter 2024. Revisit code to apply GOSAT v9 vs OCO-2 v11.1 and to attempt a multi-variate regression.

;;----- Highlevel settings.
;VerNum_OCO='9'
;VerNum_OCO='10'
VerNum_OCO='11.1'

;;----- Output version number
VerNum_out='v9c3_r20240422'

;;----- List of years to process
YearList=[2014,2015,2016,2017,2018,2019,2020]
;;YearList=[2016]
;;YearList=[2016,2017,2018,2019,2020]

;;----- Latitude range.
LatRange=[-90.,90.]
;LatRange=[-30.,30.]

;;----- Select plotting.
PlotBc_Time=0
PlotBc_LiteRetrievalVars=1
PlotBc_LiteSoundingVars=0
PlotBc_LiteAncillaryVars=0
PlotBc_LitePreprocessorsVars=0

max_dist=300. & max_dlat=2. & max_dlon=3. & dtime_allowed_range=[-2.,2.]

;;----- Other settings.
skip_calc=0
pdf=0
correct_ak = 1 & two_term_ak=0

WriteStats=0
FitBinned=0

plot_type='annual'
;view_mode='land' & view_type=1

IF (plot_type EQ 'seasonal') THEN seaslist=['DJF','MAM','JJA','SON'] ELSE seaslist=['annual']
NS=N_ELEMENTS(seaslist)

;;-----
;;----- Set some plotting parameters.
;;-----
OutDirBase='/home/ttaylor/analysis_utilities/tropical_iav/plots/'
OutDir_data='/home/ttaylor/analysis_utilities/tropical_iav/data/gosat_vs_oco2_v'+VerNum_OCO+'/bc_fit_coefficients/'+VerNum_out+'/' & FILE_MKDIR,OutDir_data
DeltaLat_1=2.5 & DeltaLon_1=5.0
;DeltaLat_1=15. & DeltaLon_1=20.
BinSizeString_1=STRING(DeltaLat_1,FORMAT='(I2)')+' lat x '+STRING(DeltaLon_1,FORMAT='(I2)')+' lon'
BinMinCount=1
FillFloat=-999.9 & FillInteger=-999L

;;----- Build structure of region information.
T={Name:'',LatBin:FLTARR(2),LonBin:FLTARR(2)}

IF (0) THEN BEGIN
   NLocations=4
   Location=REPLICATE(T,NLocations)
   Location[0].Name='Tropical-Americas'
   Location[0].LatBin=[-35,35.]
   Location[0].LonBin=[-120.,-30.]
   Location[1].Name='Tropical-Africa'
   Location[1].LatBin=[-35,35.]
   Location[1].LonBin=[-20.,60.]
   Location[2].Name='Tropical-Asia-Australia'
   Location[2].LatBin=[-35,35.]
   Location[2].LonBin=[60.,180.]
   Location[3].Name='Global'
   Location[3].LatBin=[-90.,90.]
   Location[3].LonBin=[-180.,180.]
ENDIF
IF (1) THEN BEGIN
   NLocations=1
   Location=REPLICATE(T,NLocations)
   Location[0].Name='Global'
   Location[0].LatBin=[-90.,90.]
   Location[0].LonBin=[-180.,180.]
ENDIF


;;-----
;;----- Read the full GOSAT L2Lite record.
;;-----
IF (N_ELEMENTS(Gosat) LT 1) THEN BEGIN

   ID_Start=20140906000000ULL
   ID_Stop =20200601000000ULL
   GosatFile='/home/codell/GOSAT/acos_results/b9_acos/lite/save/acos_b9_lite_20090420_20200630_small.sav'
   PRINT,'Restoring GOSAT Lite file ',GosatFile
   RESTORE,FILENAME=GosatFile,VERBOSE=1
   Gosat=TEMPORARY(Lite)
   PRINT,'Initial N soundings=',N_ELEMENTS(Gosat)
   Mask=WHERE(INRANGE(Gosat.sounding_id,[ID_Start,ID_Stop]), NMask)
   PRINT,'NMask=',NMask
   Gosat=Gosat[Mask]
   STATS,Gosat.sounding_id
   Mask=!NULL

   ;;----- Build masks.
   PRINT,'Size of Gosat prior to lat/type masking =',N_ELEMENTS(Gosat)
   qf_mask=Gosat.xco2_quality_flag EQ 0
   lat_mask=INRANGE(Gosat.latitude,LatRange)
   ;type_mask=Gosat.retrieval.surface_type EQ view_type
   ;w=WHERE(qf_mask AND lat_mask AND type_mask)
   w=WHERE(qf_mask AND lat_mask)
   Gosat=Gosat[w]
   PRINT,'Size of Gosat after qf/lat/type masking =',N_ELEMENTS(Gosat)

   ;;----- Relinquish storage.
   qf_mask=!NULL & lat_mask=!NULL & type_mask=!NULL & w=!NULL
   ;STOP

ENDIF
;STOP

;;-----
;;----- 11-March-2021.
;;----- Use new CWO collocated file.
;;----- max_dlat=2, max_dlon=3, max_dtime =2, max_dist = 300 km, and random when above 100.
;;----- Chris recommends that to change dlat/dlon/ddist/dtime, the collocation code should be rerun
;;----- instead of subsetting within this analysis.
;;-----
;file = '/home/codell/GOSAT/acos_results/b9_acos/match_oco2/match_acos_oco2_b10_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'

;;-----
;;----- 13-Feb-2024.
;;----- Use <match_gosat_to_oco2.pro> to generate versions of collocation files in directory:
;;----- /data8/ttaylor/data_ttaylor/gosat_oco2_collocations/
;;-----
IF (VerNum_OCO EQ '10') THEN $
   file='/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/match_gosat_v9_oco2_v10_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'
IF (VerNum_OCO EQ '11.1') THEN $
   file='/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'

IF (N_ELEMENTS(loaded_file) EQ 0) then loaded_file=''
IF loaded_File NE file then skip_calc=0
IF skip_calc then goto, skip_calc_

   ;;---------------------------------------------------------------
   ;;----- IF skip_calc=0, execute this block.
   ;;---------------------------------------------------------------
   PRINT,'Entering primary calculation block for the Match structure...'

   IF (loaded_file NE file) THEN BEGIN
      PRINT, 'Reading satellite matched file... ' + file_basename(file)
      RESTORE, file
      PRINT, 'Determining OCO2 Nsound'
      nall = N_ELEMENTS(match)

      gosat_qf = 1b-match.gosat.xco2_quality_flag
      gosat_xco2_bc= match.gosat.xco2
      gosat_lat=match.gosat.latitude & gosat_lon=match.gosat.longitude
      gosat_jd = acos_id_to_jd(match.gosat.sounding_id)
      gosat_id = match.gosat.sounding_id
      gosat_sza=match.gosat.solar_zenith_angle
      gosat_gain=match.gosat.Sounding.gain
      gosat_land = match.gosat.retrieval.surface_type EQ 1
      oco_jd = acos_id_to_jd(match.oco2.sounding_id)

      ;;----- Extract variables used in the collocation matching.
      nper=N_ELEMENTS(match[0].oco2)
      dlat = (fltarr(nper)+1.)#gosat_lat - match.oco2.latitude
      dlon = longitude_difference((fltarr(nper)+1.)#gosat_lon, match.oco2.longitude)
      dist = sqrt(dlat^2 + cosd((fltarr(nper)+1.)#gosat_lat)*cosd(match.oco2.latitude) * dlon^2) * 111.1 ; distance in km
      dtime = (oco_jd - (dblarr(nper)+1.d0)#gosat_jd)*24.
      del_var, oco_stddev
      loaded_file=file
   ENDIF
   ;;STOP

   ;;-----
   ;;----- 1. Match at the oco2 individual sounding level
   ;;-----
   oco2_good =  $
      1 EQ 1 $
      AND inrange(match.oco2.xco2, [300,900.]) $
      AND match.oco2.xco2_quality_flag EQ 0b $
      AND abs(dlat) LE max_dlat $
      AND abs(dlon) LE max_dlon $
      AND abs(dist) LE max_dist $
      AND inrange(dtime, dtime_allowed_range)

   Noco2 = long(total(oco2_good, 1))
   w = where(Noco2 GE 3, nw)
   PRINT,'Number of good OCO2 collocations=',nw
   PRINT,'Percent of good OCO2 collocations=',FLOAT(nw)/FLOAT(N_ELEMENTS(Noco2))*100.
   matchw=match[w]
   ;;STOP

   ;;-----
   ;;----- 2. Compute mean OCO-2 sounding statistics per match
   ;;-----
   ;oco_stddev_xco2=!NULL
   IF N_ELEMENTS(oco_stddev_xco2) NE nw THEN BEGIN
      PRINT, 'Averaging OCO-2 soundings within each GOSAT sounding.'

      ;;----- Generate data arrays.
      oco_mean_struct = REPLICATE(matchw[0].oco2[0], nall)
      oco_mean_xco2 = FLTARR(nw) + 9999.9
      oco_median_xco2 = FLTARR(nw) + 9999.9
      oco_stddev_xco2 = FLTARR(nw) + 9999.9
      oco_mean_jd = DBLARR(nw)
      mean_dtime=oco_stddev_xco2*0.
      mean_dist=oco_stddev_xco2*0.
      xco2_prior_adj_a = FLTARR(nw)
      xco2_prior_adj_b = FLTARR(nw)
      ;delta_xco2_due_to_prior = fltarr(nw)

      ;;----- Loop over number of GOSAT/OCO collocations.
      PRINT,'Looping over GOSAT/OCO collocations...'
      FOR i=0,nw-1 DO BEGIN

         ;;PRINT,'i='+NUM2STR(i+1,0)+' of '+NUM2STR(nw,0)

         q = WHERE(oco2_good[*,w[i]], nq)
         oco_mean_struct[i] = AVERAGE_STRUCT(matchw[i].oco2[q])

         ;----- fix the longitude
         this_lon =  matchw[i].oco2[q].longitude
         mnlon = this_lon[0] + MEAN(longitude_difference(this_lon, this_lon[0]))
         IF (mnlon GT 180.) THEN mnlon=mnlon-360.
         IF (mnlon LT -180.) THEN mnlon=mnlon+360.
         oco_mean_struct[i].longitude=mnlon

         ;;----- Calculate average OCO stats from the 1-100 collocated soundings.
         xco2 = matchw[i].oco2[q].xco2
         oco_mean_xco2[i] = MEAN(xco2)
         oco_median_xco2[i] = MEDIAN(xco2)
         oco_stddev_xco2[i] = STDDEV(xco2)
         oco_mean_jd[i]= MEAN(oco_jd[q,w[i]])

         ;;----- Calculate time and distance differences between OCO and GOSAT.
         mean_dtime[i] = (oco_mean_jd[i]-gosat_jd[w[i]])*24.
         mean_dist[i] = CO_SPHDIST(oco_mean_struct[i].longitude, oco_mean_struct[i].latitude, gosat_lon[w[i]], gosat_lat[w[i]], /deg, /approx)*111.1
         IF ( (mean_dist[i] GT max_dist) OR (~INRANGE(mean_dtime[i], dtime_allowed_range)) ) THEN STOP

         ;;----- Correct for different CO2 priors between versions.
         ;;----- TET. 16-Feb-2024. Added from <main_4_collocated_gosat_vs_oco2.pro>
         gosat_prior = matchw[i].gosat.co2_profile_apriori
         gosat_h = matchw[i].gosat.pressure_weight
         gosat_a = matchw[i].gosat.xco2_averaging_kernel
         oco2_prior = oco_mean_struct[i].co2_profile_apriori
         oco2_h = oco_mean_struct[i].pressure_weight
         oco2_a = oco_mean_struct[i].xco2_averaging_kernel

         ;;----- Implement CO2 prior and AK corrections
         ;;----- Former as in Eq A10 in [Wunch, ACP, 2010].
         ;;----- Also equivelant to Eq 3 in [Taylor, AMT, 2023].
         ;;----- The sum of the two terms yields the total correction to GOSAT XCO2 to harmonize it to OCO-2.
         ;;-----
         true_profile=oco2_prior
         xco2_prior_adj_a[i] = TOTAL( gosat_h * (gosat_a - 1.    ) * (gosat_prior - true_profile) )
         xco2_prior_adj_a_text='t1 = SUM( gosat_pwf * (gosat_ak - 1.    ) * (gosat_prior - true_profile) )'
         xco2_prior_adj_b[i] = TOTAL( gosat_h * (gosat_a - oco2_a) * (gosat_prior - true_profile) )
         xco2_prior_adj_b_text='t2 = SUM( gosat_pwf * (gosat_ak - oco2_ak ) * (gosat_prior - true_profile) )'

      ENDFOR
      STATS,xco2_prior_adj_a
      PRINT,' '
      PRINT,' '
      STATS,xco2_prior_adj_b
      PRINT,' '
      PRINT,' '

   ENDIF
   ;;STOP
   PRINT, 'Determining overall (per match) mask.'
   match_mask =  (oco_stddev_xco2 LE 1.5) AND (mean_dist LE max_dist) AND (gosat_qf[w])

   ;;----- Mask for latitude range.
   lat_mask=INRANGE(gosat_lat,LatRange)

   ;;-----
   ;;----- Build mask for current observation mode.
   ;;-----
   IF (N_ELEMENTS(view_mode) NE 0) THEN BEGIN
      PRINT, 'Developing mask for mode ' + view_mode
      ;MGAIN=0b & HGAIN=1b
      MGAIN='M' & HGAIN='H'
      CASE view_mode OF
         'landH'  : type_mask = gosat_land AND (gosat_gain EQ HGAIN)
         'landM'  : type_mask = gosat_land AND (gosat_gain EQ MGAIN)
         'land'   : type_mask = gosat_land
         'oceanH' : type_mask = ~gosat_land AND (gosat_gain EQ HGAIN)
         'allH'   : type_mask = gosat_gain eq HGAIN
         'all'    : type_mask = (~gosat_land AND gosat_gain eq HGAIN) OR gosat_land
      ENDCASE
      main_mask = match_mask AND type_mask AND lat_mask
   ENDIF ELSE BEGIN
      main_mask = match_mask AND lat_mask
   ENDELSE

   ;;----- Apply main_mask filter.
   g = WHERE(main_mask, nsound)
   wg=w[g]
   PRINT, 'Using ' + sc(nsound) + ' matched soundings.'
   ;STOP
   matchg = matchw[g]
   mean_dist_g=mean_dist[g]
   mean_dtime_g=mean_dtime[g]
   ;;STOP

   ;;----- Yet another layer of filtering...
   matchZ=matchG

   ;;----- SID match the full GOSAT record to matchg structure.
   PRINT,'N elements of matchg prior to GOSAT SID matching = ',N_ELEMENTS(Matchg)
   MATCH_DATA_SETS_BY_SOUNDING_ID, $
      matchZ.gosat.sounding_id, $
      matchZ, $
      Gosat.sounding_id, $
      Gosat, $
      NumMatchedSoundings, $
      Index_1=z, $
      Index_2=Index_2
   PRINT,'N elements of matchZ afer GOSAT SID matching = ',N_ELEMENTS(MatchZ)
   wgz=w[g[z]]
   ;;STOP

   ;;-----
   ;;----- Extract the needed data variables.
   ;;----- Filter variable arrays that were previously extracted from structures.
   ;;-----
   ;oco_meang = oco_mean[wgz]

   xco2_bc_diff = gosat_xco2_bc[wgz] - oco_mean_struct[wgz].xco2

   tdiff = mean_dtime[wgz]

   lat=gosat_lat[wgz] & lon=gosat_lon[wgz]

   julian_day=gosat_jd[wgz]

   ;; NOT USING IN THIS CODE
   ;year = fix(strmid(sc(gosat_id[wg]),0,4))
   ;month = fix(strmid(sc(gosat_id[wg]),4,2))
   ;season = (month mod 12)/3
   ;l = where(gosat_land, comp=o)
   land_flag=gosat_land[wgz]

   ;;----- If skip_calc is set or code has not been reset.
skip_calc_:
IF (skip_calc) $
   THEN PRINT,'Skipped the calculation block for the Match structure based on skip_calc setting.' $
   ELSE PRINT,'Completed primary calculation block for the Match structure.'
;STOP

;;----- Set the dxco2 array (y variable).
IF (correct_ak) THEN BEGIN

   IF (0) THEN HELP,xco2_bc_diff,xco2_prior_adj_a[wgz],xco2_prior_adj_b[wgz]

   IF (two_term_ak) THEN BEGIN
      PRINT,'Implement the AK correction to XCO2 to account for different CO2 priors and smoothing.'

      ;;----- Two-term AK correction.
      dxco2=xco2_bc_diff + xco2_prior_adj_a[wgz] + xco2_prior_adj_b[wgz]

   ENDIF ELSE BEGIN

      ;;----- single-terrm AK correction.
      dxco2=xco2_bc_diff + xco2_prior_adj_b[wgz]

   ENDELSE

   ;;STOP

ENDIF ELSE BEGIN

   dxco2=xco2_bc_diff

ENDELSE


STATS_TET,dxco2,dxco2_Stats
YTitle='$\Delta$XCO$_2$ [ppm]'
YRange=[-1.5,1.5] & YTicks=5 & YMinor=-1 & YTickLen=0.1 & YTickInterval=1
SubTitle=' GOSAT (v9) - OCO2 (v' + VerNum_OCO +')'
;STOP


;;--------------------------------------------------
;;----- Peform plotting on Lite.Sounding variables.
;;--------------------------------------------------
IF (1 AND PlotBc_Time) THEN BEGIN

   LegendText_x1=0.075 & LegendText_x2=0.68
   LegendText_y1=0.85 & LegendText_y2=0.35 & LegendText_dy=0.1

   JulianDay=1
   IF (JulianDay) THEN BEGIN
      XData=julian_day & HELP,XData
      full_list=['JULIAN_DAY']
      var_list=full_list
      var_index=WHERE(full_list EQ var_list)
      YRange=[[-0.5,0.5],[-1.,1.],[-1.5,1.5]]

      ;;----- dummy settings
      ;XRange=-1
      HistBinWidth=-1

   ENDIF ELSE BEGIN
      STOP
   ENDELSE

   AllFit=!NULL
   FitTypes=['skewed-sine','skewed-sine','skewed-sine']

   ;;-----
   ;;----- Loop over the list of select variables.
   ;;-----
   FOR xxx=0, N_ELEMENTS(var_list)-1 DO BEGIN
   ;;FOR xxx=16,18 DO BEGIN

      var=var_list[xxx]
      index=var_index[xxx]
      PRINT,'var='+var
      PRINT,'index='+NUM2STR(index,0)
      ;STOP

      IF (1) THEN BEGIN

         PRINT,'Current variable is ',STRTRIM(var,2)
         OutDir_plot=OutDirBase+'/gosat_vs_oco2_v'+VerNum_OCO+'/bc_vars_time/'+var+'/' & FILE_MKDIR,OutDir_plot

         ;;----- Plot settings.
         XTickLen=0.075 & XMinor=0
         YTickLen=0.025 & YMinor=0
         Dummy = LABEL_DATE(DATE_FORMAT='%N/%Z')
         XTitle='Time [MM/YY]'
         ;XTitle='GOSAT '+var
         XTickLen=0.1 & XMinor=0
         HistBinWidth= (MAX(XData)-MIN(XData)) / 20. & PRINT,'HistBinWidth=',HistBinWidth

         ;;----- Plot dxco2 vs critical variables such as SZA, or retrieval variables
         FuncExtras={OutDir:OutDir_plot,OutDir_data:OutDir_data, OutPlotFileType:'.png', $
                  VarLabel:var,XTitle:XTitle, $
                  LegendText_x1:LegendText_x1, LegendText_x2:LegendText_x2, $
                  LegendText_y1:LegendText_y1, LegendText_y2:LegendText_y2, LegendText_dy:LegendText_dy, $
                  BackHist:0,HistBinWidth:HistBinWidth,FullScatter:0,FitBinned:FitBinned,FitType:FitTypes, ZeroLine:1, $
                  YRange:YRange,ExtendXRange:1}
         PlotExtras={XTitle:XTitle, XTickFormat:'LABEL_DATE', YTitle:YTitle, XTickLen:XTickLen, YTickLen:YTickLen, XMinor:XMinor, YMinor:YMinor}
         ;PlotExtras={XTitle:XTitle, XTickFormat:'', YTitle:YTitle, XTickLen:XTickLen, YTickLen:YTickLen, XMinor:XMinor, YMinor:YMinor}

         ;;STOP
         PLOT_DXCO2_VS_BC_VARS, $
               FuncExtras, $
               PlotExtras, $
               XData, dxco2, $
               land_flag,MatchZ,lat,lon, $
               Location=Location, $
               Fit
         IF (N_ELEMENTS(AllFit) EQ 0) THEN BEGIN
            AllFit=Fit
         ENDIF ELSE BEGIN
            AllFit=[AllFit,Fit]
         ENDELSE

      ENDIF

      ;;----- End loop over select variables.
   ENDFOR
   STOP

   ;;----- Complete plotting routines.
ENDIF



;;--------------------------------------------------
;;----- Peform plotting on Lite.Retrieval variables.
;;--------------------------------------------------
IF (PlotBc_LiteRetrievalVars) THEN BEGIN

   LegendText_x2=0.075 & LegendText_x1=0.68
   LegendText_y1=0.85 & LegendText_y2=0.35 & LegendText_dy=0.1
   ;YRange=[[-0.5,0.5],[-1.,1.],[-1.5,1.5]]

   WhichVars='kitchen_sink'
   ;WhichVars='aod_total'
   ;WhichVars='xco2_raw'
   ;WhichVars='fs'
   ;WhichVars='iterations' ;; var_index=33
   ;WhichVars='windspeed_apriori
   ;WhichVars='s31'
   SkipList=['ALBEDO_SLOPE_O2A','ALBEDO_SLOPE_SCO2','ALBEDO_SLOPE_WCO2','DIVERGING_STEPS','ITERATIONS', $
             'PSURF','PSURF_APRIORI','RMS_REL_SCO2','RMS_REL_WCO2','SURFACE_TYPE']

   CASE WhichVars OF

      'kitchen_sink': $
          BEGIN
             full_list=TAG_NAMES(Gosat.Retrieval)
             var_list=TAG_NAMES(Gosat.Retrieval)
             var_index=WHERE(full_list EQ var_list)

             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])

             ;;----- dummy settings
             ;XRange=-1
             HistBinWidth=-1

             FitTypes=['quadratic','quadratic','quadratic']

          END

      'aod_total': $
          BEGIN
             full_list=TAG_NAMES(Gosat.Retrieval)
             var_list='AOD_TOTAL'
             var_index=WHERE(full_list EQ var_list)
             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])
             XRange=[0.0,0.52]
             HistBinWidth=0.02
             FitTypes=['quadratic','quadratic','quadratic']
          END

      'xco2_raw': $
          BEGIN
             full_list=TAG_NAMES(Gosat.Retrieval)
             var_list='XCO2_RAW'
             var_index=WHERE(full_list EQ var_list)
             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])
             XRange=[380.,420.]
             HistBinWidth=5.
             FitTypes=['quadratic','quadratic','quadratic']
          END
      'fs': $
          BEGIN
             full_list=TAG_NAMES(Gosat.Retrieval)
             var_list='FS'
             var_index=WHERE(full_list EQ var_list)
             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])
             ;XRange=[-]
             HistBinWidth=-1.
             FitTypes=['linear','linear','linear']
             ;;FitTypes=['quadratic','quadratic','quadratic']
          END
      'iterations': $
          BEGIN
             full_list=TAG_NAMES(Gosat.Retrieval)
             var_list='ITERATIONS'
             var_index=WHERE(full_list EQ var_list)
             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])
             ;XRange=[-]
             HistBinWidth=-1.
             ;;FitTypes=['linear','linear','linear']
             FitTypes=['quadratic','quadratic','quadratic']
          END
      'windspeed_apriori': $
          BEGIN
             full_list=TAG_NAMES(Gosat.Retrieval)
             var_list='WINDSPEED_APRIORI'
             var_index=WHERE(full_list EQ var_list)
             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])
             ;XRange=[-]
             HistBinWidth=-1.
             ;;FitTypes=['linear','linear','linear']
             FitTypes=['quadratic','quadratic','quadratic']
          END
      's31': $
          BEGIN
             full_list=TAG_NAMES(Gosat.Retrieval)
             var_list='S31'
             var_index=WHERE(full_list EQ var_list)
             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])
             ;XRange=[-]
             HistBinWidth=-1.
             ;;FitTypes=['linear','linear','linear']
             FitTypes=['quadratic','quadratic','quadratic']
          END
      ' ': $
         BEGIN
            var_list=['DWS']
         END

   ENDCASE
   ;;STOP

   AllFit=!NULL

   ;;-----
   ;;----- Loop over the list of select variables.
   ;;-----
   PRINT,'About to loop over '+NUM2STR(N_ELEMENTS(var_list),0)+' variables...'
   ;STOP
   FOR xxx=0, N_ELEMENTS(var_list)-1 DO BEGIN
   ;;FOR xxx=34,49 DO BEGIN
   ;;FOR xxx=45,49 DO BEGIN

      var=var_list[xxx]
      index=var_index[xxx]
      PRINT,'var='+var
      PRINT,'index='+NUM2STR(index,0)
      ;;STOP

      IF (var EQ 'FS') THEN FitTypes=['linear','linear','linear'] ELSE FitTypes=['quadratic','quadratic','quadratic']
      w=WHERE(STRUPCASE(var) EQ SkipList,Skip)

      IF (~Skip) THEN BEGIN

         PRINT,'Current L2Lite variable is ',STRTRIM(var,2)
         OutDir_plot=OutDirBase+'/gosat_vs_oco2_v'+VerNum_OCO+'/bc_vars_retrieval/' & FILE_MKDIR,OutDir_plot

         ;;----- Build x-data array.
         ;;var_index=
         XData=Gosat.retrieval.(index) & HELP,XData ;& STOP
         XTitle='GOSAT '+var
         XTickLen=0.075 & XMinor=0
         YTickLen=0.025 & YMinor=0
         ;HistBinWidth= (MAX(XData)-MIN(XData)) / 20. & PRINT,'HistBinWidth=',HistBinWidth
         OutPlotFileType='.png'

         ;;----- Plot dxco2 vs critical variables such as SZA, or retrieval variables
         FuncExtras={OutDir:OutDir_plot,OutDir_data:OutDir_data, OutPlotFileType:OutPlotFileType, $
                  VarLabel:var,XTitle:XTitle, $
                  LegendText_x1:LegendText_x1, LegendText_x2:LegendText_x2, $
                  LegendText_y1:LegendText_y1, LegendText_y2:LegendText_y2, LegendText_dy:LegendText_dy, $
                  BackHist:1,HistBinWidth:HistBinWidth,FullScatter:0, $
                  FitType:FitTypes,FitBinned:FitBinned,ZeroLine:1, $
                  YRange:YRange,ExtendXRange:0}

         PlotExtras={XTitle:XTitle, YTitle:YTitle, XTickLen:XTickLen, YTickLen:YTickLen, XMinor:XMinor, YMinor:YMinor}

         ;;STOP
         PLOT_DXCO2_VS_BC_VARS, $
               FuncExtras, $
               PlotExtras, $
               XData, dxco2, $
               land_flag,MatchZ,lat,lon, $
               Location=Location, $
               Fit
         IF (N_ELEMENTS(AllFit) EQ 0) THEN BEGIN
            AllFit=Fit & Fit=!NULL
         ENDIF ELSE BEGIN
            AllFit=[AllFit,Fit] & Fit=!NULL
         ENDELSE
         XData=!NULL

         ;STOP
      ENDIF ELSE BEGIN
         PRINT,'Current variable is being skipped because it matches an entry in SkipList.'
         ;STOP
      ENDELSE

      ;;----- End loop over select variables.
   ENDFOR
   ;;STOP

   ;;----- Complete plotting routines.
ENDIF

;;--------------------------------------------------
;;----- Peform plotting on Lite.Sounding variables.
;;--------------------------------------------------
IF (PlotBc_LiteSoundingVars) THEN BEGIN

   AllFit=!NULL

   WhichVars='kitchen_sink'
   CASE WhichVars OF
      'kitchen_sink': $
          BEGIN
             full_list=TAG_NAMES(Gosat.Sounding)
             var_list=full_list
             var_index=WHERE(STRLOWCASE(var_list) NE 'gain' AND STRLOWCASE(var_list) NE 'l1b_type')
             var_list=var_list[var_index]

             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])

             ;;----- dummy settings
             ;XRange=-1
             HistBinWidth=-1

             FitTypes=['quadratic','quadratic','quadratic']
          END

   ENDCASE
   ;STOP

   PRINT,'About to loop over '+NUM2STR(N_ELEMENTS(var_list),0)+' variables...'
   ;STOP
   FOR xxx=0, N_ELEMENTS(var_list)-1 DO BEGIN
   ;;FOR xxx=16,18 DO BEGIN

      var=var_list[xxx]
      index=var_index[xxx]
      PRINT,'var='+var
      PRINT,'index='+NUM2STR(index,0)
      ;;STOP

      IF (1) THEN BEGIN

         PRINT,'Current L2Lite variable is ',STRTRIM(var,2)
         OutDir_plot=OutDirBase+'/gosat_vs_oco2_v'+VerNum_OCO+'/bc_vars_sounding/' & FILE_MKDIR,OutDir_plot

         ;;----- Build x-data array.
         ;;var_index=
         XData=Gosat.Sounding.(index) & HELP,XData
         XTitle='GOSAT '+var
         XTickLen=0.075 & XMinor=0
         YTickLen=0.025 & YMinor=0
         ;HistBinWidth= (MAX(XData)-MIN(XData)) / 20. & PRINT,'HistBinWidth=',HistBinWidth

         ;;----- Plot dxco2 vs critical variables such as SZA, or retrieval variables
         FuncExtras={OutDir:OutDir_plot,OutDir_data:OutDir_data,VarLabel:var,XTitle:XTitle, $
                  LegendText_x1:LegendText_x1, LegendText_x2:LegendText_x2, $
                  LegendText_y1:LegendText_y1, LegendText_y2:LegendText_y2, LegendText_dy:LegendText_dy, $
                  BackHist:1,HistBinWidth:HistBinWidth,FullScatter:0, $
                  FitType:FitTypes,FitBinned:FitBinned,ZeroLine:1, $
                  YRange:YRange,ExtendXRange:0}

         PlotExtras={XTitle:XTitle, YTitle:YTitle, XTickLen:XTickLen, YTickLen:YTickLen, XMinor:XMinor, YMinor:YMinor}

         ;;STOP
         PLOT_DXCO2_VS_BC_VARS, $
               FuncExtras, $
               PlotExtras, $
               XData, dxco2, $
               land_flag,MatchZ,lat,lon, $
               Location=Location, $
               Fit
         IF (N_ELEMENTS(AllFit) EQ 0) THEN BEGIN
            AllFit=Fit
         ENDIF ELSE BEGIN
            AllFit=[AllFit,Fit]
         ENDELSE

      ENDIF
      ;;----- End loop over select variables.
   ENDFOR
   ;STOP

   ;;----- Complete plotting routines.
ENDIF

;;--------------------------------------------------
;;----- Peform plotting on Lite.Retrieval variables.
;;--------------------------------------------------
IF (PlotBc_LiteAncillaryVars) THEN BEGIN

   AllFit=!NULL

   WhichVars='kitchen_sink'
   CASE WhichVars OF
      'kitchen_sink': $
          BEGIN
             full_list=TAG_NAMES(Gosat)
             var_list=full_list
             var_index=WHERE( $
                             STRLOWCASE(var_list) NE 'meteorology' $
                             AND STRLOWCASE(var_list) NE 'preprocessors' $
                             AND STRLOWCASE(var_list) NE 'retrieval' $
                             AND STRLOWCASE(var_list) NE 'sounding' $
                             AND STRLOWCASE(var_list) NE 'sounding_id' $
                             AND STRLOWCASE(var_list) NE 'xco2_qf_bitflag' $
                             AND STRLOWCASE(var_list) NE 'xco2_quality_flag')

             var_list=var_list[var_index]

             ;i_1=FIX(var_index[0]) & i_2=FIX(var_index[0])

             ;;----- dummy settings
             ;XRange=-1
             HistBinWidth=-1

             FitTypes=['quadratic','quadratic','quadratic']
          END

   ENDCASE
   ;STOP

   PRINT,'About to loop over '+NUM2STR(N_ELEMENTS(var_list),0)+' variables...'
   ;STOP
   FOR xxx=0, N_ELEMENTS(var_list)-1 DO BEGIN
   ;;FOR xxx=16,18 DO BEGIN

      var=var_list[xxx]
      index=var_index[xxx]
      PRINT,'var='+var
      PRINT,'index='+NUM2STR(index,0)
      ;;STOP

      IF (1) THEN BEGIN

         PRINT,'Current L2Lite variable is ',STRTRIM(var,2)
         OutDir_plot=OutDirBase+'/gosat_vs_oco2_v'+VerNum_OCO+'/bc_vars_ancillary/' & FILE_MKDIR,OutDir_plot

         ;;----- Build x-data array.
         ;;var_index=
         XData=Gosat.(index) & HELP,XData
         XTitle='GOSAT '+var
         XTickLen=0.075 & XMinor=0
         YTickLen=0.025 & YMinor=0
         ;HistBinWidth= (MAX(XData)-MIN(XData)) / 20. & PRINT,'HistBinWidth=',HistBinWidth

         ;;----- Plot dxco2 vs critical variables such as SZA, or retrieval variables
         FuncExtras={OutDir:OutDir_plot,OutDir_data:OutDir_data,VarLabel:var,XTitle:XTitle, $
                  LegendText_x1:LegendText_x1, LegendText_x2:LegendText_x2, $
                  LegendText_y1:LegendText_y1, LegendText_y2:LegendText_y2, LegendText_dy:LegendText_dy, $
                  BackHist:1,HistBinWidth:HistBinWidth,FullScatter:0, $
                  FitType:FitTypes,FitBinned:FitBinned,ZeroLine:1, $
                  YRange:YRange,ExtendXRange:0}

         PlotExtras={XTitle:XTitle, YTitle:YTitle, XTickLen:XTickLen, YTickLen:YTickLen, XMinor:XMinor, YMinor:YMinor}

         ;;STOP
         PLOT_DXCO2_VS_BC_VARS, $
               FuncExtras, $
               PlotExtras, $
               XData, dxco2, $
               land_flag,MatchZ,lat,lon, $
               Location=Location, $
               Fit
         IF (N_ELEMENTS(AllFit) EQ 0) THEN BEGIN
            AllFit=Fit
         ENDIF ELSE BEGIN
            AllFit=[AllFit,Fit]
         ENDELSE

      ENDIF
      ;;----- End loop over select variables.
   ENDFOR
   ;STOP

   ;;----- Complete plotting routines.
ENDIF


;;--------------------------------------------------
;;----- Peform plotting on Lite.Preprocessor variables.
;;--------------------------------------------------
IF (0 AND PlotBc_LitePreprocessorsVars) THEN BEGIN

   KitchenSink=1
   IF (KitchenSink) THEN BEGIN
      var_list=TAG_NAMES(Gosat.preprocessors)
   ENDIF ELSE BEGIN
      ;STOP
      var_list=['airmass']
   ENDELSE

   ;;-----
   ;;----- Loop over the list of select variables.
   ;;-----
   for i_v=0, N_ELEMENTS(var_list)-1 do begin
   ;for i_v=0,2 do begin

      var=var_list[i_v]
      IF (1) THEN BEGIN

         PRINT,'Current L2Lite variable is ',STRTRIM(var,2)
         OutDir=OutDirBase+'/gosat_vs_oco2_v'+VerNum_OCO+'/bc_vars_preprocessors/'+var+'/' & FILE_MKDIR,OutDir_plot

         ;;----- Build x-data array.
         XData=Gosat.preprocessors.(i_v) & HELP,XData
         XTitle='GOSAT '+var
         XTickLen=0.1 & XMinor=0
         HistBinWidth= (MAX(XData)-MIN(XData)) / 20. & PRINT,'HistBinWidth=',HistBinWidth

         ;;----- Plot dxco2 vs critical variables such as SZA, or retrieval variables
         PlotExtras={OutDir:OutDir_plot,OutDir_data:OutDir_data,VarLabel:var,XTitle:XTitle, $
                  XMinor:XMinor,$
                  XTickLen:XTickLen,  $
                  YRange:YRange,YTickInterval:YTickInterval,YTicks:YTicks,YMinor:YMinor,YTickLen:YTickLen, $
                  BackHist:1,HistBinWidth:HistBinWidth,FullScatter:1,FitBinned:FitBinned}

         PLOT_DXCO2_VS_BC_VARS, $
               PlotExtras, $
               XData, dxco2, $
               land_flag,MatchZ,lat,lon, $
               Location=Location
            ;STOP
      ENDIF

      ;;----- End loop over select variables.
   ENDFOR
   ;STOP

   ;;----- Complete plotting routines.
ENDIF
;STOP

;;-----
;;----- Plot ranked order metrics based on saved AllFit structure.
;;-----
IF (1) THEN BEGIN
   ViewModes=['Land-H','Land-M','Ocean-H'] & NumViewModes=N_ELEMENTS(ViewModes)
ENDIF
FOR mode=0,NumViewModes-1 DO BEGIN
   CurrentMode=ViewModes[mode]
   w=WHERE(AllFit.viewmode EQ CurrentMode, NW)
   IF (NW GT 0) THEN BEGIN
      PRINT, ' '
      PRINT,'***** RESULTS FROM FITTING *****'
      PRINT,'Viewing mode ',CurrentMode
      metric=ABS(AllFit[w].reduction_upper_abs)
      variable=AllFit[w].variable
      so=SORT(metric)
      N_metrics=MIN([100,N_ELEMENTS(variable)])
      variable_ranked=REVERSE(variable[so])
      metric_ranked=REVERSE(metric[so])
      PRINT,'Ranked list of variables in descending order using the absolute correction to XCO2 as the sorting metric:'
      FOR var=0,N_metrics-1 DO BEGIN
         PRINT, variable_ranked[var]+' ('+NUM2STR(metric_ranked[var],1)+' ppm)'
      ENDFOR




      ;STOP
   ENDIF

         ;;----- Plot dxco2 vs critical variables such as SZA, or retrieval variables
   XTitle=''
   YTitle='$\delta$$\Delta$XCO$_2$'
         XTickLen=0.075 & XMinor=0 & XTickInterval=1
         YTickLen=0.025 & YMinor=0
         FuncExtras={OutDir:OutDir_plot,OutDir_data:OutDir_data,VarLabel:var, $
                  LegendText_x1:LegendText_x1, LegendText_x2:LegendText_x2, $
                  LegendText_y1:LegendText_y1, LegendText_y2:LegendText_y2, LegendText_dy:LegendText_dy}

         PlotExtras={XTitle:XTitle, YTitle:YTitle, XTickLen:XTickLen, YTickLen:YTickLen, XMinor:XMinor, YMinor:YMinor, $
              XTickInterval:XTickInterval}

   PLOT_DXCO2_METRICS, $
      AllFit, $
      FuncExtras=FuncExtras, $
      PlotExtras=PlotExtras, $
      Location=Location, $
      CurrentMode, $
      variable_ranked[0:N_metrics-1], $
      metric_ranked[0:N_metrics-1]

   ;STOP
ENDFOR


STOP
END
