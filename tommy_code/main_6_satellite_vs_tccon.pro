
;;----- 
;;----- High level settings.
;;-----
IF (1) THEN BEGIN
   Sensor='gosat'
   ;Version='v9'
   Version='v9c'
   ;Group = 'acos'
   Retrieval = strupcase(Sensor)

   Ver_tccon='ggg2020'

   ;mode = 'landH' & ModeLong='Land H-gain (dark)' 
   ;mode = 'landM' & ModeLong='Land M-gain (bright)'  
   mode = 'oceanH' & ModeLong='OceanH'  

   ;mode = 'land' & ModeLong='Land (H+M gains)' 
   ;mode = 'all' & ModeLong='Combined Ocean/LandH/LandM' 

   OneToOneTitle='(a)'

   ;;----- Set Julian Date range.
   jdrange=[julday(4,1,2009, 0.), julday(6,31,2020,0.)] ; GOSAT v9
   OcoTimeRange=0

ENDIF

IF (0) THEN BEGIN
   Sensor='oco2'
   Version='v10'
   ;Group='oco'
   Retrieval = strupcase(Sensor)

   mode = 'oceanG' & ModeLong='Ocean Glint' 
   ;mode = 'land' & ModeLong='Land (ND+GL+TG)'
   ;mode = 'landT' & ModeLong='Land Target'

   OneToOneTitle='(c)'

   ;;----- Set Julian Date range.
   jdrange=[julday(9,1,2014, 0.), julday(2,28,2022,0.)] ; OCO2 v10
   OcoTimeRange=1
   
ENDIF

;;----- Set the data aggregation level.
overpass_mean=1
daily = 0
weekly = 0
monthly = 0
annual = 0 ; this is for the last panel of the per-site plots.  Typically it is monthly.

;;----- Define output data directory
OutDir='/home/ttaylor/analysis_utilities/tropical_iav/plots/satellite_vs_tccon/'+Sensor+'_'+Version+'/'
FILE_MKDIR,OutDir

;;----- Plot settings.
Plot11Combined=1
PlotTimeDiff=1 ;; strange problem with the SYMBOL function in this procedure
AnalyzePerSite=0
IF (AnalyzePerSite) THEN BEGIN
   Plot11PerSite=0
   PlotTimePerSite=1
ENDIF

;;-----
;;----- Original v9 file from CWO copied to data11 drive.
;;----- file = '/home/codell/GOSAT/acos_results/b9_acos/tccon_match/acos_b9_lite_20090420_20191231_small_tccon_match20200605_Debra.sav' 
;;-----
IF (Sensor EQ 'gosat') THEN BEGIN

   IF (Ver_tccon EQ 'ggg2014') THEN BEGIN
      MatchFile='/home/codell/GOSAT/acos_results/b9_acos/tccon_match/acos_b9_lite_20090420_20200630_small_tccon_match20210419_Debra.sav'
   ENDIF
   IF (Ver_tccon EQ 'ggg2020') THEN BEGIN
      MatchFile='/home/codell/GOSAT/acos_results/b9_acos/tccon_match/ggg2020/acos_b9_lite_20090420_20200630_small_tccon_match20240410_Debra.sav'
   ENDIF

   IF (Version EQ 'v9c') THEN BEGIN
      ;LiteSaveFile_v9x='/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c2_r20240329/gosat_v9c2_r20240329_lite_20090420_20200630_filtered_xmedium.sav'
      LiteSaveFile_v9x='/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c3_r20240422/gosat_v9c3_r20240422_lite_20090420_20200630_filtered_xmedium.sav'
   ENDIF
   ;STOP
ENDIF

;;-----
;;----- OCO-2 v10 files from Mattheus Kiel (used in the OCO v10 data paper).
;;-----
IF (Sensor EQ 'oco2') THEN BEGIN

   IF (Version EQ 'v10') THEN BEGIN
      MatchFile='/home/codell/OCO2_results/b10_tests/tccon_match/ggg2014/oco2_b10v1_lite_20140906_20220228_filtered_tiny_tccon_match20221010_Debra.sav'

   ENDIF
ENDIF

;IF (retrieval EQ 'ACOS') THEN retrieval += '/' + Build

PRINT,'Processing observation mode ',mode
SubTitle=mode
;STOP

;;-----
;;----- Set collocation criteria.
;;----- Really the collocation is performed in a separate code that generates the input save file.
;;----- Keep very large lat/lon box so that everything in collocation save file is used in analysis.
;;-----
max_lat =25. ;; degrees
max_lon =50. ;; degrees
max_time=2.0 ;; hours

;;----- Set land/ocean multiplier. Not really sure how these are derived?
land_mult=1.0 & ocean_mult=1.0
;land_mult=0.9958/0.9959 & ocean_mult=0.9956/0.9950

;;----- Set mon_av value.
;;----- 1=
;;----- 2=
;;----- 3=
IF (daily OR weekly) THEN mon_av=1 ELSE mon_av=2
IF (annual) THEN mon_av=3

;;----- Set the TCCON and GOSAT averaging kernel.
correct_ak_tccon = 1 ; Correct TCCON for the satellite AK

;;----- 
;;----- Optional settings.
;;-----
bc=1  ;; not clear to me difference in "bc" and "bc_only"
bc_only=1 ;& IF (bc_only OR bc) THEN SubTitle+='_bc'
fit_type=2 ;; 0=linfit, 1=ladfit, 2=linfit_york
seasonal_fit=1
median=1
min_per_overpass=2 & min_per_site=10 
new_no_top=1
Range_xco2=[375.,430.]
acos_b9_tighter=0
fix_yrange_site_plots=1
yranges_sites=[[380,420.],[-6,6.],[-3,3.]]
IF (annual) THEN yranges_sites[*,2] = [-1,1.]

;;----- Load fossil fuel map. Used to prescreen TCCON data.
FFmap=READ_H5_FILE('/data1/ODIAC/odiac2013a_5km_excl_intl_2012_smoothed1.0.h5') 
;FFlimit = 300.0 ; seems a decent limit for japan
FFlimit = 100000.0 ; effectively no limit

;;----- Set mode information.
ND_mode = 0 & GL_mode = 1 & TG_mode=2 & XS_mode=3
gainH = 'H' & gainM = 'M'
land_mode = 1 & ocean_mode=0

;;-----
;;----- Restore the satellite/TCCON collocated data from IDL save.
;;-----
IF ( N_ELEMENTS(restored_file) EQ 0) THEN restored_file=''
IF (restored_file NE MatchFile) THEN BEGIN

   PRINT, 'Reading ' + MatchFile
   RESTORE, MatchFile, VERBOSE=1
   ;STOP

   ;;----- If GOSAT version other than v9 (harmonized), then read appropriate Lite Save file
   ;;----- from which to obtain modified XCO2 that can be SID matched to collocations.
   IF (Sensor EQ 'gosat') THEN BEGIN 

      IF (Version EQ 'v9c') THEN BEGIN
 
         ;;----- This section is under construction 8-April-2024
         IF (N_ELEMENTS(Lite) EQ 0) THEN RESTORE,FILENAME=LiteSaveFile_v9x,VERBOSE=1
         ;STOP

         ;;----- Obtain the GOSAT v9X XCO2 values from Lite Save.
         ID_Start=20140906000000ULL
         ID_Stop =20200601000000ULL

         PRINT,'Initial number of GOSAT v9b soundings=',N_ELEMENTS(Lite)
         PRINT,'Sounding match GOSAT Lite record to GOSAT v9 (collocated to TCCON).
         MATCH_DATA_SETS_BY_SOUNDING_ID, $
            Match.ACOS.sounding_id, $
            Match, $
            Lite.sounding_id, $
            Lite, $
            NumMatchedSoundings, $
            Index_1=Index_1, $
            Index_2=Index_2
         STATS,Lite.sounding_id
         STATS,Lite.xco2_quality_flag
         ;STOP

      ENDIF

      ;;----- Subset Lite data that is specific to GOSAT sensor.
      acos_operation_mode = match.acos.sounding.gain ; H or M
      acos_quality_flag = match.acos.xco2_quality_flag


      ;STOP
      ;;----- End block for restoring ACOS GOSAT data.
   ENDIF

   IF ( (Sensor EQ 'oco2') AND (Version EQ 'v10')) THEN BEGIN 

      ;;----- Subset Lite data that is specific to OCO sensor.
      acos_operation_mode = match.acos.sounding.operation_mode ; 0=nd, 1=gl, 2=tg, 3=xs

      ;;----- For OCO, only good quality flagged data was included in collocation.
      acos_quality_flag=MAKE_ARRAY(N_ELEMENTS(acos_operation_mode),/BYTE,VALUE=0b)

      ;;----- End block for restoring ACOS GOSAT data.
   ENDIF

   ;;----- Gather Lite data that is invariant with sensor.
   acos_lon = match.acos.longitude
   acos_lat = match.acos.latitude
   acos_surface_type = match.acos.retrieval.surface_type  ; 0=ocean, 1=land
   acos_jd = acos_id_to_jd(match.acos.sounding_id)
   tccon_jd = match.tccon.jd
   dtime = (tccon_jd - acos_jd)*24.

   IF (Version EQ 'v9c') THEN BEGIN
      acos_xco2_bc = Lite.xco2_harmonized_to_oco2_v11_1
      acos_xco2_raw = acos_xco2_bc
   ENDIF ELSE BEGIN
      acos_xco2_bc = match.acos.xco2
      acos_xco2_raw = match.acos.retrieval.xco2_raw
   ENDELSE

   restored_file =MatchFile
ENDIF
PRINT,'Total input data size for sensor/version '+Sensor+' '+Version+'='+NUM2STR(N_ELEMENTS(acos_xco2_bc),0 ) 
;STOP      

;;----- Reassign the original ACOS GOSAT BC XCO2, while preserving the original.
;;----- Modify the values as necessary using the land/ocean multipliers.
acos_xco2_bc = acos_xco2_bc
IF (ocean_mult NE 1.0) THEN acos_xco2_bc *= (acos_surface_type eq 0)*ocean_mult + (acos_surface_type eq 1)*1.0
IF (land_mult  NE 1.0) THEN acos_xco2_bc *= (acos_surface_type eq 1)*land_mult + (acos_surface_type eq 0)*1.0

;;----- Extract fossil fuel variables. Used to help filter the TCCON data.
IF (N_ELEMENTS(FFmatch) EQ 0) THEN BEGIN
   ; which gridbox is each observation in?
   nlon_ff = N_ELEMENTS(ffmap.lon)
   nlat_ff = N_ELEMENTS(ffmap.lat)
   ilon = FIX( ((acos_lon+360.) mod 360.) / (360.0/nlon_ff) )
   ilat = FIX(  (acos_lat+90.) / (180./nlat_ff) )
   FFmatch = ffmap.ff[ilon + nlon_ff*ilat]
ENDIF

;;-----
;;----- Set spatial requirements for collocations. 
;;-----
spatial_requirement = num2str(2*max_lon,0)+'x'+num2str(2*max_lat,0)+' deg box, '
time_requirement=tex('+- 2 hrs') ;& SubTitle+='_temporal-2hrs'

tg = n_tags(match)-1

;;-----
;;----- Not sure this code will do anything if Group ~= 'acos'.
;;-----
;IF (Group EQ 'acos') THEN BEGIN
IF (1) THEN BEGIN

    IF (N_ELEMENTS(distance) NE N_ELEMENTS(match) ) THEN distance= CO_SPHDIST(acos_lon, acos_lat, match.tccon.lon, match.tccon.lat, /deg)    
    IF (N_ELEMENTS(dlat) NE N_ELEMENTS(match) )     THEN dlat = ABS(match.tccon.lat-acos_lat)
    IF (N_ELEMENTS(dlon) NE N_ELEMENTS(match) )     THEN dlon = ABS(longitude_difference(match.tccon.lon, acos_lon))
    ;;STOP

    ;;----- Set prefilter.
    prefilter = acos_quality_flag EQ 0
    IF ( (acos_b9_tighter) AND (Sensor EQ 'gosat') AND (Version EQ 'v9') ) $
       THEN prefilter= (acos_quality_flag EQ 0) AND (acos_b9_tighter_quality_flag_lite(match.acos) )

    suffix = '_' + mode

    CASE mode OF
       'landN': group_filter =  acos_operation_mode eq ND_mode and $
                                acos_surface_type eq land_mode
       'landG': group_filter =  acos_operation_mode eq GL_mode and $
                                acos_surface_type eq land_mode
  
       'landNG' : group_filter =  acos_surface_type eq land_mode $
                              AND acos_operation_mode NE TG_mode

       'land' : group_filter =  acos_surface_type eq land_mode
   
       'landT':group_filter =  acos_surface_type eq land_mode AND $
                               acos_operation_mode eq TG_mode
       'landH': group_filter = acos_surface_type eq land_mode AND $
                               acos_operation_mode eq gainH                     
       'landM': group_filter = acos_surface_type eq land_mode AND $
                               acos_operation_mode eq gainM                      
       'oceanG': group_filter =  acos_operation_mode eq GL_mode and $
                                 acos_surface_type ne land_mode
       'oceanN': group_filter =  acos_operation_mode eq ND_mode and $
                                 acos_surface_type ne land_mode      
       'oceanH': group_filter = acos_surface_type NE land_mode AND $
                                acos_operation_mode eq gainH  
       'all'   : group_filter = bytarr(n_elements(match))+1b
       'all_tropics' : group_filter = abs(acos_lat) LE 30.
   ENDCASE
ENDIF
;;STOP

;;-------------------------------------------------
;;----- 1) Filter the data.
;;-------------------------------------------------
good1 =  prefilter AND group_filter 
w1 = WHERE(good1, N_good_1)
PRINT,'Fraction of data selected for observation mode and prefilter =',FLOAT(N_good_1)/FLOAT(N_ELEMENTS(prefilter))*100.0 
;STOP

match1 = match[w1]
IF (Sensor EQ 'gosat') AND (Version EQ 'v9c') THEN Lite_1=Lite[w1]

distance1 = distance[w1]
dlat1 = dlat[w1]
dlon1 = dlon[w1]
dtime1=dtime[w1]
ffmatch1 = ffmatch[w1]

good =  1 eq 1 $
        AND INRANGE(match1.tccon.jd, jdrange) $          
        AND (dlat1 LE max_lat) $
        AND (dlon1 LE max_lon) $
        AND (ABS(dtime1) LE max_time)
        ;;AND ~(elt(match1.tccon.station, exclude_ocean_sites) AND acos_surface_type[w1] eq 0) $   
        ;;AND ~elt(match1.tccon.station, exclude_sites) $   
        ;AND match1.tccon.n GT 1 $
        ;AND ffmatch1 LE FFlimit $
w2 = WHERE(good, N_good_2)
PRINT,'Fraction of data within selected time/lat/lon boundaries =',FLOAT(N_good_2)/FLOAT(N_ELEMENTS(Match1))*100.0 
;STOP

;;-----
;;----- Impose requirement that # matched soundings per overpass exceed some minimum value
;;-----
flagminN = BYTARR(N_good_2)
IF (STRMID(sensor,0,3) EQ 'oco') THEN BEGIN 
   overpass=match1.tccon.station + '_' + string(match1.acos.sounding.orbit, form='(i5.5)')
ENDIF ELSE BEGIN
   overpass=match1.tccon.station + '_'+ string(match1.acos.sounding_id/1000000,form='(i8.8)') + '_'+string(match1.acos.sounding.path,form='(i2.2)')
ENDELSE

overpass_names = DIFFERENT(overpass[w2])
so=SORT(overpass[w2])
nover = N_ELEMENTS(overpass_names)
v=VALUE_LOCATE(overpass[w2[so]],overpass_names)
bb=-1L
FOR s=0,nover-1L DO BEGIN
   aa=bb+1 & bb=v[s] & g=so[aa:bb] & ng=bb-aa+1 
   IF (ng GE min_per_overpass) THEN flagminN[g] = 1b
ENDFOR
N_bad_1=LONG(TOTAL(flagminN EQ 0))
PRINT, 'Removed ' + sc(N_bad_1) + ' soundings that failed min_per_overpass='+sc(min_per_overpass)+' test.'
good[w2] = good[w2] AND flagminN
w2 = WHERE(good, N_good_3)
PRINT,'Fraction of good data meeting minimum number of overpasses =',FLOAT(N_good_3)/FLOAT(N_ELEMENTS(FlagMinN))*100.0 
;STOP

;;-----
;;----- Impose requirement that # matched soundings per TCCON station exceed some minimum value
;;-----
flagminN = BYTARR(N_good_3)

;;----- Rename "cyprus" to "nicosia"
nicosia=WHERE(match1.tccon.station EQ 'cyprus',NNic)
IF (NNic GT 0) THEN match1[nicosia].tccon.station='nicosia'
 
site=match1[w2].tccon.station
site_names = different(site)
so=sort(site)
nsites = n_elements(site_names)
v=value_locate(site[so],site_names)
bb=-1L
FOR s=0,nsites-1 DO BEGIN
   aa=bb+1 & bb=v[s] & g=so[aa:bb] & ng=bb-aa+1    
   IF (ng GE min_per_site) THEN flagminN[g] = 1b
ENDFOR
N_bad_2=long(total(flagminN eq 0))
PRINT, 'Removed ' + sc(N_bad_2) + ' soundings that failed min_per_site='+sc(min_per_site)+' test.' ;& STOP
good[w2] = good[w2] AND flagminN
w2 = WHERE(good, N_good_4)
PRINT,'Fraction of good data meeting minimum number of soundings per site =',FLOAT(N_good_4)/FLOAT(N_ELEMENTS(FlagMinN))*100.0 
;;STOP

;;-----
;;----- Find all matched soundings
;;-----
w = w2
ww = w1[w]
matchw = match1[w] ; this still includes all sites
IF (Sensor EQ 'gosat') AND (Version EQ 'v9c') THEN Lite_w=Lite_1[w]

overpass=overpass[w]
IF (N_good_4 LE 2) THEN STOP

;;-----
;;----- Select either BC or RAW satellite XCO2.
;;-----
IF (bc) THEN BEGIN
    xco2_satellite = acos_xco2_bc[ww]
ENDIF ELSE BEGIN
   xco2_satellite=acos_xco2_raw[ww]
ENDELSE
;;STOP

;;-----
;;----- use AK-corrected XCO2 from TCCON
;;-----
IF (correct_ak_tccon) THEN BEGIN
   xco2_tccon = CORRECT_TCCON_ACOS_AK(matchw)
   PRINT,'Statistics of the TCCON averaging kernel correction: '
   STATS,(xco2_tccon - matchw.tccon.xco2)
ENDIF ELSE BEGIN
   xco2_tccon=matchw.tccon.xco2
ENDELSE
;STOP

;;-----
;;----- Average down the data. 
;;-----
agg_name='Single Sounding' ; default

IF (weekly) THEN BEGIN
    caldat, matchw.tccon.jd, mon, day, year
    doy = date2doy(year, mon, day)
    week = (year-2009)*52 + ((doy-1)<363) / 7
    sort_thing = sc(week) + '_' + matchw.tccon.station
    agg_name= 'Weekly'
ENDIF 
IF (daily AND ~weekly) THEN BEGIN   
    caldat, matchw.tccon.jd, mon, day, year
    sort_thing = matchw.tccon.station + '_' + sc(year) + string(mon, form='(i2.2)') + string(day, form='(i2.2)') 
    sort_thing = sort_thing + '_'+(['ocean','land'])[acos_surface_type[ww]]
    agg_name = 'Daily'
    IF (STRMID(sensor,0,3) EQ 'oco') THEN mincount=10 
    IF (sensor EQ 'gosat') THEN mincount=3
ENDIF
IF (overpass_mean) THEN BEGIN
    caldat, matchw.tccon.jd, mon, day, year
    sort_thing = matchw.tccon.station + '_' + sc(year) + string(mon, form='(i2.2)') + string(day, form='(i2.2)') 
    sort_thing = overpass + '_'+(['ocean','land'])[acos_surface_type[ww]]
    agg_name = 'Overpass-Mean'
    IF (STRMID(sensor,0,3) EQ 'oco') THEN mincount=10 
    IF (sensor EQ 'gosat') THEN mincount=3
ENDIF
IF (monthly) THEN BEGIN   
    caldat, matchw.tccon.jd, mon, day, year
    sort_thing = matchw.tccon.station + '_' + sc(year) + string(mon, form='(i2.2)') 
    sort_thing = sort_thing + '_'+(['ocean','land'])[acos_surface_type[ww]]
    agg_name = 'Monthly'
    IF (STRMID(sensor,0,3) EQ 'oco') THEN mincount=100 
    IF (sensor EQ 'gosat') THEN mincount=10
ENDIF
;;STOP

;;------------------------------------------------------------------
;;----- 2) Create plots of all sites combined; Satellite vs. TCCON.
;;------------------------------------------------------------------
IF (Sensor EQ 'gosat') THEN gosat_err = matchw.acos.xco2_uncertainty ELSE gosat_err=0.0
;; gosat_err = gosat_err*0.+1. ; uniform weighting
;STOP

;;-----
;;----- Generate 1:1 scatter plot for all stations combined for current observation mode.
;;-----
IF (Plot11Combined) THEN BEGIN
   IF ~(daily or weekly) THEN psym=9 ELSE psym=17
   ;Suffix='.png'
   Suffix='.pdf'
   SymSize=0.75

   ;;----- Ensure not too small tccon error.    
   matchw.tccon.sd = matchw.tccon.sd > 0.2

   ;;----- Generate one:to:one scatter plot.
   IF (1) THEN BEGIN
      XTitle='TCCON '+STRUPCASE(Ver_tccon)+' XCO$_2$ [ppm]'
      YTitle='ACOS '+STRUPCASE(Sensor)+' '+Version+' XCO$_2$ [ppm]
      ;;IF (bc_only OR bc) THEN YTitle=YTitle+' (BC applied)' ELSE YTitle=YTitle+' (No BC applied)'
      STANDARD_PLOT_SATELLITE_VS_TCCON_ONE_TO_ONE, $
         matchw.tccon, $
         xco2_tccon, $
         xco2_satellite, $
         error=gosat_err, $
         sort_thing=sort_thing, $
         subtit=file_basename(MatchFile), $
         psy=psym, sym_size=SymSize, $
         agg_name=agg_name, $
         xr=Range_xco2,yr=Range_xco2, mincount=mincount, fit_type=fit_type, no_top=new_no_top, $
         font_size=13, XTitle=XTitle, YTitle=YTitle, PlotTitle=OneToOneTitle, Mode=Mode,biascorr=bc,bc_only=bc_only, $
         OutDir=OutDir,LongMode=ModeLong,STitle=SubTitle, $
         RPlot, Info=Info
      
      ;;----- Save one:one plot to output file.
      OD=OutDir & FILE_MKDIR,OD
      OutFile=OD+Sensor+'_'+Version+'_vs_tccon_'+Ver_tccon+'_one-to-one_'+agg_name+'_'+SubTitle+Suffix 
      PRINT,'OutFile=',OutFile ; & STOP
      RPlot.save,Outfile
      ;STOP
   ENDIF

   ;;-----
   ;;----- Write out summary statistics to compare with Hannakaisa analysis.
   ;;-----
   PRINT,'Total N stations (v9) [CSU analysis] = ',N_ELEMENTS(Info)
   PRINT,'Total N days for all stations (v9) [CSU analysis] = ',TOTAL(Info.N)
   PRINT,'StationMeanBiasStdDev (v9) [CSU analysis] = ',STDDEV(Info.mean)
   
   ;;----- Write info to ascii file.
   ;;----- Used by Abhishek for latitude band plots.
   IF (0) THEN BEGIN
      OD=OutDir+mode+'/site_info/' & FILE_MKDIR,OD
      SiteInfoFile=OD+'gosat_vs_tccon_site_info_'+agg_name+'_'+SubTitle+'.txt'
      PRINT,'SiteInfoFile=',SiteInfoFile ;;& STOP
      OPENW,Lun1,SiteInfoFile,/GET_LUN
      NSites=N_ELEMENTS(Info) & PRINT,'NSites=',NSites
      FOR S=0,NSites-1 DO BEGIN
         FOR N=0,Info[S].N-1 DO BEGIN
            PRINTF,Lun1, Info[S].Station, ', ', Info[S].lat, ', ', Info[S].N, ', ' ,Info[S].Julian[N], ', ', Info[S].xco2_tccon[N], ', ', Info[S].xco2_gosat[N], $
            FORMAT='(A,A,F7.3,A,I3,A,E15.8,A, 2(F7.3 ,A))'
         ENDFOR
      ENDFOR
      FREE_LUN,Lun1
      STOP
   ENDIF
ENDIF



;;---------------------------------------------------------------------------------------------
;;----- Generate Delta XCO2 time trend for all stations combined (in the current viewing mode).
;;----------------------------------------------------------------------------------------------
IF (PlotTimeDiff) THEN BEGIN

   ;;----- Generate delta xco2 vs time plot.
   YTitle='$\Delta$ XCO$_2$' + ' ('+STRUPCASE(Sensor)+' '+Version+ ' $-$ TCCON '+STRUPCASE(Ver_tccon)+') [ppm]'
   ;BeginJD=JULDAY(6,1,2004,12,0,0) ;; pad start date to leave room on left for site legend
   BeginJD=JULDAY(6,1,2007,12,0,0) ;; pad start date to leave room on left for site legend
   EndJD=JULDAY(12,1,2022,12,0,0)
   xr=[BeginJD,EndJD] & XMajor=4 & XMinor=0
   yr=[-6,6] & YMajor=7 & YMinor=3
   SymSize=0.5
   FontSize=16
   LegendFontSize=10
   StatsFontSize=10
   PlotDim=[1000,700]

   PlotExtras={XRange:xr,XMajor:XMajor, XMinor:XMinor, $
               YRange:yr,YMajor:YMajor,YMinor:YMinor, YStyle:1, $
               Sym_Size:SymSize, $
               Font_Size:FontSize, $
               Title:OneToOneTitle,YTitle:YTitle, $
               Dimensions:PlotDim}

   ;STOP
   STANDARD_PLOT_SATELLITE_VS_TCCON_TIMESERIES, $
      matchw.tccon, $
      xco2_tccon, $
      xco2_satellite, $
      error=gosat_err, $
      sort_thing=sort_thing, $
      mincount=mincount, fit_type=fit_type, $
      agg_name=agg_name, $
      no_top=new_no_top, $
      Mode=Mode, $
      Legend_Font_Size=LegendFontSize, $
      Stats_Font_Size=StatsFontSize, $
      biascorr=bc, bc_only=bc_only, $
      OutDir=OutDir,LongMode=ModeLong,STitle=SubTitle, $
      DeltaPlot, $
      _extra=PlotExtras

   ;;----- Save delta xco2 vs time to output file.
   OD=OutDir+'time_diff/' & FILE_MKDIR,OD
   DeltaOutFile=OD+Sensor+'_'+Version+'_vs_tccon_'+Ver_tccon+'_time-diff_'+agg_name+'_'+SubTitle+'.pdf'
   PRINT,'DeltaOutFile=',DeltaOutFile ;& STOP
   DeltaPlot.save,DeltaOutfile   
   xco2_sat=!NULL

   STOP
ENDIF


;;-------------------------------------------------------------
;;----- 3) Create plots of individual sites; ACOS vs. TCCON.
;;-------------------------------------------------------------
jdr=jdrange
;site_names = ['all', different(matchw.tccon.station)]
site_names = [different(matchw.tccon.station)]
nsites = n_elements(site_names)
 
;;-----
;;----- Analyze individual TCCON stations.
;;-----
IF (AnalyzePerSite) THEN BEGIN

   IF (fix_yrange_site_plots) THEN yran_sites=yranges_sites ELSE del_var, yran_sites
   
   ;;----- Loop over individual sites.
   FOR s=0,nsites-1 DO BEGIN
   
      site_name = site_names[s]
      PRINT,'Analyzing data for TCCON station ',site_name
      IF (STRLOWCASE(STRMID(site_name,0,3)) eq 'all') THEN all_sites=1 ELSE all_sites=0
      IF (~all_sites) THEN g = WHERE(matchw.tccon.station EQ site_name, ng) ELSE g = WHERE(matchw.tccon.station NE '',ng)
         
      
      ;;----- Check minimum number of sounding/site requirement.
      IF (ng LT Min_Per_Site) THEN BEGIN
         PRINT,'Number of soundings at site ',site_name,' =',ng
         PRINT,'Below minimum requirement set to ',Min_Per_Site
         CONTINUE  ; too few matches for this site     
      ENDIF ELSE BEGIN
         PRINT,'Found ',ng,' good collocated matches for site ',site_name
         ;;PAUSE
      ENDELSE

      del_var, station_day_info

      ;;----- Subset data sets.     
      sort_thing_g = sort_thing[g]
      match_g = matchw[g]
      xco2_tccon_g=xco2_tccon[g]
      xco2_satellite_g=xco2_satellite[g]
      match_g.tccon.sd = match_g.tccon.sd > 0.2
      ;yerror=SQRT(match_g.acos.xco2_uncertainty^2. + match_g.tccon.sd^2.)
      gosat_err_g=gosat_err[g]
      HELP,sort_thing_g,match_g,xco2_tccon_g,xco2_satellite_g,gosat_err_g
      ;STOP

      IF (all_sites) THEN match_g.tccon.station = 'all'
     
      ;;-----
      ;;----- Generate 1:1 scatter plot per site.
      ;;-----
      IF ( (ng GT 2) AND ~(STRLOWCASE(STRMID(site_name,0,3)) EQ 'all') ) THEN BEGIN
         IF ~(daily or weekly) THEN psym=9 ELSE psym=17


         ;;----- Set time range.
         ;;----- For analysis against OCO time range, truncate data.
         IF (OcoTimeRange) THEN BEGIN
            BeginMonth=8 & BeginYear=2014
            EndMonth=6 & EndYear=2020
            BeginJD=JULDAY(BeginMonth,1,BeginYear,12,0,0) & BeginJDPlot=JULDAY(6,1,2014,12,0,0) 
            EndJD=JULDAY(EndMonth,1,EndYear,12,0,0)       & EndJDPlot=JULDAY(1,1,2021,12,0,0)
            Good_Time = WHERE(INRANGE(matchg.tccon.jd,[BeginJD,EndJD]),NGoodTime)
            PRINT,'NGoodTime=',NGoodTime
            IF (NGoodTime GT 1) THEN BEGIN
               match_g=matchg[Good_Time]
               xco2_satellite_g=xco2_satellite[Good_Time]
               sort_thing_g=sort_thingg[Good_Time]
               YError=YError[Good_Time]
               ;;STOP
            ENDIF ELSE BEGIN
               
            ENDELSE
         ENDIF ELSE BEGIN
            BeginMonth=1 & BeginYear=2009
            EndMonth=10 & EndYear=2020
            BeginJD=JULDAY(BeginMonth,1,BeginYear,12,0,0) & BeginJDPlot=JULDAY(1,1,2009,12,0,0)
            EndJD=JULDAY(EndMonth,1,EndYear,12,0,0)       & EndJDPlot=JULDAY(1,1,2021,12,0,0)
            NGoodTime=N_ELEMENTS(xco2_satellite_g)
         ENDELSE
         xr=[BeginJDPlot,EndJDPlot]
         yr=[-5,5]

         Range_xco2=[ FLOOR(MIN([[xco2_tccon_g],[xco2_satellite_g]])), CEIL(MAX([[xco2_tccon_g],[xco2_satellite_g]]))*1.01 ]

         ;;----- Generate one:to:one scatter plot (and delta xco2 vs time).
         ;;nsites = n_elements(site_names)
         IF ( Plot11PerSite AND NGoodTime GT 1) THEN BEGIN
            STANDARD_PLOT_SATELLITE_VS_TCCON_ONE_TO_ONE, $
               match_g.tccon, $
               xco2_tccon_g, $
               xco2_satellite_g, $
               err=gosat_err_g, $
               sort_thing=sort_thing_g, $
               agg_name=agg_name, $
               xr=Range_xco2,yr=Range_xco2, mincount=mincount, fit_type=fit_type, no_top=new_no_top, $
                  Site_Name=Site_Name,nallsites=nsites,siteindex=s,font_size=14, $
                  Title=site_name, YTitle='ACOS GOSAT ('+Version+') XCO$_2$ [ppm]', Mode=Mode,biascorr=bc,bc_only=bc_only, $
                  OutDir=OutDir,LongMode=ModeLong,STitle=SubTitle, $
                  RPlot,DeltaPlot

            ;;----- Save one:one plot to output file.
            OD=OutDir+mode+'/one_to_one/time_range_'+NUM2STR(BeginYear,0)+STRING(BeginMonth,FORMAT='(I2.2)')+'-'+NUM2STR(EndYear,0)+STRING(EndMonth,FORMAT='(I2.2)')+'/'
            FILE_MKDIR,OD
            OutFile=OD+'gosat_vs_tccon_'+Ver_tccon+'_one-to-one_site_'+Site_Name+'_'+SubTitle+'.png' 
            PRINT,'OutFile=',OutFile ; & STOP
            RPlot.save,Outfile
            ;STOP
         ENDIF

         IF (PlotTimePerSite) THEN BEGIN

            ;;----- Generate delta xco2 vs time plot.
            YTitle='$\Delta$ XCO$_2$ (GOSAT '+Version+' - TCCON) [ppm]'
            ;;IF (bc_only OR bc) THEN YTitle=YTitle+' (BC applied)' ELSE YTitle=YTitle+' (No BC applied)'

            xr=[BeginJD,EndJD] & XMajor=4 & XMinor=0
            yr=[-6,6] & YMajor=7 & YMinor=3
            SymSize=1.5
            FontSize=20
            LegendFontSize=10
            StatsFontSize=15
            PlotDim=[1000,700]

            PlotExtras={XRange:xr,XMajor:XMajor, XMinor:XMinor, $
               YRange:yr,YMajor:YMajor,YMinor:YMinor, YStyle:1, $
               Sym_Size:SymSize, $
               Font_Size:FontSize, $
               Title:site_name,YTitle:YTitle, $
               Dimensions:PlotDim}
            LegendFontSize=10
            StatsFontSize=10

            fit_type_time=0 & PlotPosTime=[0.10,0.50,0.95,0.80]
            STANDARD_PLOT_SATELLITE_VS_TCCON_TIMESERIES, $
               match_g.tccon, $
               xco2_tccon_g, $
               xco2_satellite_g, $
               error=gosat_err_g, $
               sort_thing=sort_thing_g, $
               fit_type=fit_type_time, $
               agg_name=agg_name, $
               _extra=PlotExtras, $
               no_top=new_no_top, $
               Site_Name=Site_Name, $
               nallsites=nsites,siteindex=s, $
               Mode=Mode, $
               Legend_Font_Size=LegendFontSize, $
               Stats_Font_Size=StatsFontSize, $
               biascorr=bc,bc_only=bc_only, $
               OutDir=OutDir,LongMode=ModeLong,STitle=SubTitle, $
               DeltaPlot, Position=PlotPosTime

            ;;----- Save delta xco2 vs time to output file.
            OD=OutDir+mode+'/time_diff/time_range_'+NUM2STR(BeginYear,0)+STRING(BeginMonth,FORMAT='(I2.2)')+'-'+NUM2STR(EndYear,0)+STRING(EndMonth,FORMAT='(I2.2)')+'/'
            FILE_MKDIR,OD
            DeltaOutFile=OD+'gosat_vs_tccon_'+Ver_tccon+'_time-diff_site_'+Site_Name+'_'+SubTitle+'.png'
            PRINT,'DeltaOutFile=',DeltaOutFile ;& STOP
            DeltaPlot.save,DeltaOutfile   

         ENDIF
      ENDIF
   ;;----- End loop over individual sites.
   ;STOP 
   ENDFOR
 
ENDIF
 
END
