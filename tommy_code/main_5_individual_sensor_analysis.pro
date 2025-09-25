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

;VerNum_GOSAT='v9'
;VerNum_GOSAT='v9b'
VerNum_GOSAT='v9c_20240214'

;;----- List of years to process
YearList=[2014,2015,2016,2017,2018,2019,2020]
;;YearList=[2016]
;;YearList=[2016,2017,2018,2019,2020]

;;----- Latitude range.
LatRange=[-90.,90.]
Region='global'
;LatRange=[-30.,30.]
;Region='tropics'

;;----- Individual sensor analysis.
;;----- Use this to get sounding density of Good QF 
;;----- for GOSAT/OCO indepedently across the overlapping time range.
;;-----
;IndividualSensors=1 & WhichSensor='oco_v10' & ReadTruncateOcoLite=0 & GridIndivSensor=1
;IndividualSensors=1 & WhichSensor='gosat_v9' & GridIndivSensor=1
;IndividualSensors=1 & WhichSensor='gosat_v9b' & GridIndivSensor=1
IndividualSensors=0 & ID_Start=20140906 & ID_Stop =20200601

;;----- Grid for maps or for Hovmollers.
GridMap=1
GridHovmoller=0
PlotTimeSeries=0
WriteStats=0

view_mode='land' & view_type=1.
;view_mode='oceanH' & view_type=0.
;view_mode='all'

;;----- Use this for gridded spatial maps.
IF (GridMap) THEN BEGIN
   ;plot_type='seasonal' & independent_years=0
   plot_type='annual' & independent_years=0
ENDIF

;;----- Use this for gridded Hovmoller data and lat bin time series.
IF (GridHovmoller OR PlotTimeSeries) THEN BEGIN
   plot_type='annual' & independent_years=0
ENDIF

IF (plot_type EQ 'seasonal') THEN seaslist=['DJF','MAM','JJA','SON'] ELSE seaslist=['annual']
NS=N_ELEMENTS(seaslist)

;;----- Other settings.
skip_calc=0
correct_different_priors = 1

;;-----
;;----- Set some plotting parameters.
;;-----
OutDirBase='/home/ttaylor/analysis_utilities/tropical_iav/plots/'
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
;Location[0].Name='Tropical-Americas'
;Location[0].LatBin=[-35,35.]
;Location[0].LonBin=[-120.,-30.]
;Location[1].Name='Tropical-Africa'
;Location[1].LatBin=[-35,35.]
;Location[1].LonBin=[-20.,60.]
;Location[2].Name='Tropical-Asia-Australia'
;Location[2].LatBin=[-35,35.]
;Location[2].LonBin=[60.,180.]
Location[0].Name='global'
Location[0].LatBin=[-90.,90.]
Location[0].LonBin=[-180.,180.]
;Location[0].Name='Tropics'
;Location[0].LatBin=[-35.,35.]
;Location[0].LonBin=[-180.,180.]

;;-----
;;----- Perform analysis on individual sensors.
IF (IndividualSensors) THEN BEGIN


   CASE WhichSensor OF

      'gosat_v9': BEGIN
                  ID_Start=20140906000000ULL
                  ID_Stop =20200601000000ULL
                  ;STOP
                  ;; Instead of reading the below save file, probably better to read in raw Lite. 
                  ;; Never know what filtering was done. 
                  GosatFile='/home/codell/GOSAT/acos_results/b9_acos/lite/save/acos_b9_lite_20090420_20200630_small.sav'
                  IF (N_ELEMENTS(SatData) EQ 0) THEN BEGIN
                     PRINT,'Restoring file ',GosatFile
                     RESTORE,FILENAME=GosatFile,VERBOSE=1
                     Gosat=TEMPORARY(Lite)
                     PRINT,'Initial N soundings=',N_ELEMENTS(Gosat)
                     Mask=WHERE(INRANGE(Gosat.sounding_id,[ID_Start,ID_Stop]), NMask)
                     PRINT,'NMask=',NMask
                     Gosat=Gosat[Mask]
                     STATS,Gosat.sounding_id
                     SatData=TEMPORARY(Gosat)

                    ;;----- Build masks.
                    PRINT,'Size of SatData prior to qf/lat/type masking =',N_ELEMENTS(SatData)
                    qf_mask=satdata.xco2_quality_flag EQ 0
                    lat_mask=INRANGE(satdata.latitude,LatRange)
                    type_mask=satdata.retrieval.surface_type EQ view_type
                    w=WHERE(qf_mask AND lat_mask AND type_mask)
                    SatData=SatData[w]
                    PRINT,'Size of SatData after qf/lat/type masking =',N_ELEMENTS(SatData)

                    PRINT,'CALL ACOS_ID_TO_JD...' 
                    JD = ACOS_ID_TO_JD(SC(SatData.sounding_id))
                    ;STOP

                  ENDIF
                  END


      'gosat_v9b': BEGIN
                      ID_Start=20140906000000ULL
                      ID_Stop =20200601000000ULL
                  
                      NcFile='../data/gosat_v9_to_oco2_v10/gosat_mod_xco2.nc'
                      SatData=READ_H5_FILE(NcFile)

                      PRINT,'Initial N soundings=',N_ELEMENTS(SatData)
                      Mask=WHERE(INRANGE(SatData.sid,[ID_Start,ID_Stop]), NMask)
                      PRINT,'NMask=',NMask
                      SatData=SatData[Mask]
                      STATS,SatData.sid
                    
                      ;;----- Build masks.
                      PRINT,'Size of SatData prior to lat/type masking =',N_ELEMENTS(SatData)
                      qf_mask=satdata.qf EQ 0
                      lat_mask=INRANGE(satdata.latitude,LatRange)
                      type_mask=satdata.surface_type EQ view_type
                      w=WHERE(qf_mask AND lat_mask AND type_mask)
                      SatData=SatData[w]
                      PRINT,'Size of SatData after qf/lat/type masking =',N_ELEMENTS(SatData)

                      PRINT,'CALL ACOS_ID_TO_JD...' 
                      JD = ACOS_ID_TO_JD(SC(SatData.sid))
                      ;STOP
                   END


      'oco_v10': BEGIN
                 ;ID_Start=2014090600000000ULL
                 ;ID_Stop =2020060100000000ULL
                         ;2014090602072531
                 ID_Start=20140906L
                 ID_Stop=20200601L
                 OcoSaveFile=WhichSensor+'_agg_super_lite_20140906-20200601.sav'

                 IF (ReadTruncateOcoLite) THEN BEGIN
                    OcoFileBase='/data6/OCO2/product/Lite/B10/LtCO2/*.nc4'
                    OcoFiles=FILE_SEARCH(OcoFileBase,COUNT=NFiles)
                    PRINT,'Found ',NFiles,' OCO L2Lite files in archive.'
                    Dates=LONG('20'+STRTRIM(STRMID(FILE_BASENAME(OcoFiles,'.nc4'),11,6 ),2))
                    Mask=WHERE(INRANGE(Dates,[ID_Start,ID_Stop]), NMask)
                    PRINT,'Truncating L2Lite list to Start/Stop data range, NFiles=',NMask
                    OcoFiles=OcoFiles[Mask]

                    ReadList=['sounding_id','latitude','longitude','xco2','xco2_quality_flag','Sounding/operation_mode','Retrieval/surface_type']
                    PRINT,'Reading OCO files'
                    oco=READ_H5_FILES(OcoFiles,READ=ReadList,VERBOSE=1)
                    PRINT,'Number of OCO soundings=',N_ELEMENTS(oco)

                    ;;----- Build masks.
                    STATS,Oco.sounding_id
                    SatData=TEMPORARY(oco)

                    SAVE,FILENAME=OcoSaveFile, SatData
                    ;STOP
                 ENDIF ELSE BEGIN
                    PRINT,'Restoring file ',OcoSaveFile
                    IF (N_ELEMENTS(SatData) EQ 0) THEN RESTORE,FILENAME=OcoSaveFile
                 ENDELSE

                    PRINT,'Size of SatData prior to lat/type masking =',N_ELEMENTS(SatData)
                    qf_mask=satdata.xco2_quality_flag EQ 0
                    lat_mask=INRANGE(satdata.latitude,LatRange)
                    type_mask=satdata.retrieval.surface_type EQ view_type
                    ;opmode_mask=INRANGE(SatData.Sounding.operation_mode,[0,2])
                    w=WHERE(qf_mask AND lat_mask AND type_mask)
                    SatData=SatData[w]
                    PRINT,'Size of SatData after qf/lat/type/opmode masking =',N_ELEMENTS(SatData)

                    PRINT,'CALL ACOS_ID_TO_JD...' 
                    JD = ACOS_ID_TO_JD(SC(SatData.sounding_id/10L))
                    ;;STOP
             END
   ENDCASE

   
   ;;----- Extract year/month/season.   
   IF (WhichSensor EQ 'gosat_v9b') THEN BEGIN  
      year = FIX(STRMID(SC(SatData.sid),0,4))
      month = FIX(STRMID(SC(SatData.sid),4,2))
      season = (month MOD 12)/3
   ENDIF ELSE BEGIN
      year = FIX(STRMID(SC(SatData.sounding_id),0,4))
      month = FIX(STRMID(SC(SatData.sounding_id),4,2))
      season = (month MOD 12)/3
   ENDELSE

   ;;-----
   ;;----- Grid data density for spatial maps.
   ;;-----
   IF (GridIndivSensor) THEN BEGIN

      ;;-----
      ;;----- Loop over number of seasons (or once for annual).
      ;;-----
      FOR s=0,ns-1 DO BEGIN
      
         ;;----- Search for various surface types.
         CASE SeasList[s] OF
            'DJF' : ws=WHERE( (Month EQ 12) OR ( INRANGE(Month, [1,2])), NThis )
            'MAM' : ws=WHERE( INRANGE(Month, [3,5]), NThis) 
            'JJA' : ws=WHERE( INRANGE(Month, [6,8]), NThis) 
            'SON' : ws=WHERE( INRANGE(Month, [9,11]), NThis) 
            'annual' : ws=WHERE( INRANGE(Month, [1,12]), NThis)
         ENDCASE            
         PRINT,'N Elements in Season/Year=',NThis


         ;;----- Determine the ZData for gridding.
         IF (WhichSensor EQ 'gosat_v9b') THEN BEGIN
            ZData=[[SatData[ws].xco2_mod], [SatData[ws].xco2_mod - SatData[ws].xco2_orig]] 
            VarNames=['xco2','xco2_adjust']
         ENDIF ELSE BEGIN
            ZData=SatData[ws].xco2
            VarNames=['xco2']
         ENDELSE

         IF (1) THEN BEGIN
            ;WIS_IMAGE_TET,SatData.xco2,SatData.latitude,SatData.longitude,ISOTROPIC=1,Mini=405,Maxi=415,Magnify=1,OcoCtIndex=1
            WIS_IMAGE_BIN_TET,SatData[ws].latitude,SatData[ws].longitude,ISOTROPIC=1,Magnify=4,OcoCtIndex=1, NLon=LONG(360./DeltaLon_1),NLat=LONG(180./DeltaLat_1);,Mini=0,Maxi=75E3
            ;STOP
         ENDIF

         IF (1) THEN BEGIN

            OutDir=OutDirBase+WhichSensor+'/' & FILE_MKDIR,OutDir
            ;STOP
            OutFileBase=OutDir+WhichSensor+'_lite-vars-vs-lon-lat_' $
               +STRMID(STRTRIM(STRING(ID_Start),2),0,8) $
               +'-'+STRMID(STRTRIM(STRING(ID_Stop),2),0,8)+'_'+SeasList[s]+'_'+view_mode+'-'+Region
            OutPdfFile=OutFileBase+'.pdf'
            OutNcFile=OutFileBase+'.nc'

            OutInfo=MAKE_ARRAY(15,/STRING,VALUE='')
            OutInfo[0]=OutPdfFile
            OutInfo[1]='N soundings per bin'
            OutInfo[2]=SeasList[s]
            OutInfo[3]='N='+STRTRIM(STRING(FLOAT(NThis)/1E6,FORMAT='(F6.2)'),2)+'M (SS)'
            OutInfo[4]=''
            OutInfo[10]='plasma'
            OutInfo[11]='longitude'
            OutInfo[12]='latitude'
            OutInfo[13]=STRING(FillFloat,FORMAT='(F6.1)')
            OutInfo[14]=STRING(FillInteger,FORMAT='(I4)')
            PRINT,'GRID_DATA_WRITE_NETCDF_MULTI for spatial map of ',WhichSensor,' sounding density for season '+STRTRIM(SeasList[s],2)
            STOP
            GRID_DATA_WRITE_NETCDF_MULTI, $
               SatData[ws].Longitude, $
               SatData[ws].Latitude, $
               ZData, $   ;; inputs
               VarNames, $
               ;X_Array,Y_Array,Z_Binned, $ ;; output gridded X/Y/Z data
               'XVar',DeltaLon_1, LONG(360./DeltaLon_1), [-180,180],'lon','longitude (degrees)', $
               'YVar',DeltaLat_1, LONG(180./DeltaLat_1), LatRange,'lat','latitude (degrees)', $
               'NVar','N','Number of soundings per lat/lon grid box', $
               ;'ZVar','ppm','Delta XCO2 (GOSAT-OCO)', $
               OutPdfFile,OutNcFile,OutInfo,BinMinCount=BinMinCount, $
               FillFloat=FillFloat,FillInteger=FillInteger, Verbose=1
         ENDIF

         IF (1) THEN BEGIN
            OutDir=OutDirBase+WhichSensor+'/' & FILE_MKDIR,OutDir

            PRINT,'Calculting DSS from JD...'
            DSS=JD[ws]-MIN(JD[ws])+1L 
            ;;STOP

            NDaysRange=[0,MAX(DSS)]
            DX=HovInc_Day & NX=ROUND( (FLOAT(MAX(DSS))-FLOAT(MIN(DSS)) ) / FLOAT(DX) ) & PRINT,'NX=',NX 
            DY=HovInc_Lat & NY=ROUND(180./FLOAT(DY)) & PRINT,'NY=',NY
            HovText=STRING(HovInc_Day,FORMAT='(I2)') +'day_'+ STRING(HovInc_Lat,FORMAT='(I2)')+'lat'

            OutFileBase=OutDir+WhichSensor+'_lite-vars-vs-time-lat_' $
               +HovText+'_' $
               +STRMID(STRTRIM(STRING(ID_Start),2),0,8) $
               +'-'+STRMID(STRTRIM(STRING(ID_Stop),2),0,8)+'_'+SeasList[s]+'_'+view_mode+'-'+Region
            OutPdfFile=OutFileBase+'.pdf'
            OutNcFile=OutFileBase+'.nc'

            OutInfo=MAKE_ARRAY(15,/STRING,VALUE='')
            OutInfo[0]=OutPdfFile
            OutInfo[1]='lite vars'
            OutInfo[2]=SeasList[s]
            OutInfo[3]='N='+STRTRIM(STRING(FLOAT(NThis)/1E6,FORMAT='(F6.2)'),2)+'M (SS)'
            OutInfo[4]=''
            OutInfo[10]='plasma'
            OutInfo[11]='Date (YYYY-MM)'
            OutInfo[12]='latitude'
            OutInfo[13]=STRING(FillFloat,FORMAT='(F6.1)')
            OutInfo[14]=STRING(FillInteger,FORMAT='(I4)')
            ;STOP
            PRINT,'GRID_DATA_WRITE_NETCDF_MULTI for spatial map of ',WhichSensor,' sounding density for season '+STRTRIM(SeasList[s],2)
            GRID_DATA_WRITE_NETCDF_MULTI, $
               DSS, $
               SatData[ws].Latitude, $
               ZData, $   ;; inputs
               VarNames, $
               ;X_Array,Y_Array,Z_Binned, $ ;; output gridded X/Y/Z data
               'XVar',DX, NX, NDaysRange,'days','days since start', $
               'YVar',DY, NY, [-90,90],'lat','latitude (degrees)', $
               'NVar','N','Number of soundings per lat/lon grid box', $
               ;'ZVar','ppm','Delta XCO2 (GOSAT-OCO)', $
               OutPdfFile,OutNcFile,OutInfo,BinMinCount=BinMinCount, $
               FillFloat=FillFloat,FillInteger=FillInteger, Verbose=1


         ENDIF

      ENDFOR
   ENDIF


   STOP
ENDIF

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

IF (VerNum_GOSAT EQ 'v9c_20240214') THEN BEGIN
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

   GosatAggSaveFile='/data10/data8/ttaylor/data_ttaylor/lite_save_files/gosat/v9c_r20240214/gosat_v9c_r20240214_lite_20090420_20200630_filtered_xmedium.sav'
ENDIF


IF (N_ELEMENTS(loaded_file) EQ 0) then loaded_file=''
IF loaded_File NE CollocatedFile then skip_calc=0
IF skip_calc then goto, skip_calc_

IF (loaded_file NE CollocatedFile) THEN BEGIN

   PRINT, 'Reading satellite matched file... ' + FILE_BASENAME(CollocatedFile)
   IF (N_ELEMENTS(Match) EQ 0) THEN RESTORE, CollocatedFile, VERBOSE=1

   ;;----- Either read GOSAT XCO2 from Match structure,
   ;;IF (VerNum_GOSAT EQ 'v9') THEN gosat_xco2_bc= match.gosat.xco2
   gosat_xco2_bc= match.gosat.xco2

   ;;----- or replace it with values from modified GOSAT XCO2 files...
   ;;----- 15-Feb-2024. No longer using the python code to create the "gosat_mod_xco2.nc" file
   ;;----- Replace this block with a simple substitution of the xco2 or harmonized xco2.
   ;;-----
   IF (VerNum_GOSAT EQ 'v9b') THEN BEGIN
      NcFile='../data/gosat_v9_to_oco2_v10/gosat_mod_xco2.nc'
      ModData=READ_H5_FILE(NcFile)
      MATCH_DATA_SETS_BY_SOUNDING_ID, $
         ModData.sid, $
         ModData, $
         match.gosat.sounding_id, $
         Match, $
         NMatched, $
         Index_1=Index_1, $
         Index_2=Index_2

      ;;----- Replace the original GOSAT v9 bias corrected XCO2
      ;;----- with the modified values.
      gosat_xco2_bc=ModData.xco2_mod
      gosat_xco2_adj=ModData.xco2_mod - ModData.xco2_orig
   ENDIF

   ;;----- 15-Feb-2024. New block to replace above one.
   IF (VerNum_GOSAT EQ 'v9c_20240214') THEN BEGIN
      IF (N_ELEMENTS(Lite) EQ 0) THEN RESTORE,FILENAME=GosatAggSaveFile,VERBOSE=1
      ;;STOP
      MATCH_DATA_SETS_BY_SOUNDING_ID, $
         Lite.sounding_id, $
         Lite, $
         match.gosat.sounding_id, $
         Match, $
         NMatched, $
         Index_1=Index_1, $
         Index_2=Index_2

      ;;----- Replace the original GOSAT v9 bias corrected XCO2
      ;;----- with the modified values.
      gosat_xco2_harm=Lite.xco2_harmonized_to_oco2_v11_1
      gosat_xco2_adj=Lite.xco2_adjustment_to_gosat_v9

   ENDIF   

   ;;----- Subset needed variables for GOSAT.
   gosat_qf = 1b-match.gosat.xco2_quality_flag
   gosat_lat=match.gosat.latitude & gosat_lon=match.gosat.longitude
   gosat_jd = acos_id_to_jd(match.gosat.sounding_id)
   gosat_id = match.gosat.sounding_id
   gosat_sza=match.gosat.solar_zenith_angle
   gosat_gain=match.gosat.Sounding.gain
   gosat_land = match.gosat.retrieval.surface_type EQ 1

   ;;----- Subset/calculate needed variables for OCO-2.
   oco_jd = ACOS_ID_TO_JD(match.oco2.sounding_id)
   PRINT, 'Determining OCO2 Nsound'               
   nall = N_ELEMENTS(match)
   nper=N_ELEMENTS(match[0].oco2)
   dlat = (FLTARR(nper)+1.)#gosat_lat - match.oco2.latitude
   dlon = LONGITUDE_DIFFERENCE( (FLTARR(nper)+1.)#gosat_lon, match.oco2.longitude )
   dist = SQRT(dlat^2 + COSD( (FLTARR(nper)+1.)#gosat_lat)*COSD(match.oco2.latitude) * dlon^2) * 111.1 ; distance in km
   dtime = (oco_jd - (DBLARR(nper)+1.d0)#gosat_jd)*24.
   del_var, oco_stddev_xco2
   loaded_file=CollocatedFile

ENDIF
;;STOP

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
   AND (ABS(dlat) LE max_dlat) $
   AND (ABS(dlon) LE max_dlon) $ 
   AND (ABS(dist) LE max_dist) $
   AND INRANGE(dtime, dtime_allowed_range) ;$
   ;AND match.oco2.sounding.land_water_indicator EQ 0

Noco2 = LONG(TOTAL(oco2_good, 1)) 
w = WHERE(Noco2 GE 3, nw)
PRINT,'Number of good OCO2 collocations=',nw
PRINT,'Percent of good OCO2 collocations=',FLOAT(nw)/FLOAT(N_ELEMENTS(Noco2))*100.
matchw=match[w] 
IF (VerNum_GOSAT EQ 'v9c_20240214') THEN BEGIN
   gosat_xco2_harm_w=gosat_xco2_harm[w]
   gosat_xco2_adj_w=gosat_xco2_adj[w]
ENDIF
;STOP

;;-----
;;----- 2. Compute mean OCO-2 sounding statistics per match
;;-----
;oco_stddev_xco2=!NULL
IF (N_ELEMENTS(oco_stddev_xco2) NE nw) THEN BEGIN

   PRINT, 'Averaging the OCO soundings collocated to each GOSAT sounding.'

   ;;----- Generate data arrays.
   oco_mean_struct = REPLICATE(matchw[0].oco2[0], nall)
   oco_mean_xco2 = FLTARR(nw) + 9999.9
   oco_median_xco2 = FLTARR(nw) + 9999.9
   oco_stddev_xco2 = FLTARR(nw) + 9999.9
   oco_mean_jd = DBLARR(nw)
   mean_dtime=oco_stddev_xco2*0.
   mean_ddist=oco_stddev_xco2*0.

   xco2_prior_adj_a = FLTARR(nw)
   xco2_prior_adj_b = FLTARR(nw)

   ;;----- Loop over number of GOSAT/OCO collocations.
   FOR i=0,nw-1 DO BEGIN

      PRINT,'i='+NUM2STR(i+1,0)+' of '+NUM2STR(nw,0)
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
      mean_ddist[i] = CO_SPHDIST(oco_mean_struct[i].longitude, oco_mean_struct[i].latitude, gosat_lon[w[i]], gosat_lat[w[i]], /deg, /approx)*111.1
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
      ;   ELSE gosat_xco2 = matchw[i].gosat.xco2
      gosat_prior = matchw[i].gosat.co2_profile_apriori
      gosat_h = matchw[i].gosat.pressure_weight
      gosat_a = matchw[i].gosat.xco2_averaging_kernel
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
   ;STOP 
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
   PRINT,'Writing gridded data to file ',OutH5File
   STRUCT_TO_H5_SIMPLE,MatchW,OutH5File

   ;;----- Add XCO2 CO2 prior and AK adjustment arrays to H5 output file.
   ADD_H5_FIELD,OutH5File,xco2_prior_adj_a,'xco2_prior_adj_term_1'
   ADD_H5_FIELD,OutH5File,xco2_prior_adj_b,'xco2_ak_adj_term_2'

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
   ADD_H5_FIELD,OutH5File,gosat_xco2_harm_w,'gosat_xco2_harmonized_to_oco2'
   ADD_H5_FIELD,OutH5File,gosat_xco2_adj_w,'gosat_xco2_adjustment'

   ;STOP
ENDIF

;;-----
;;----- Sanity check plot of the calculated xco2 adjustment due to prior.
;;-----
AdjXco2Plot=0
IF (AdjXco2Plot) THEN BEGIN

   ;;----- Plot settings.
   Symbol='Circle'
   Col_0a='black'
   Col_0b='black'
   Col_1a='deep pink'
   Col_1b='medium violet red'
   Col_2a='blue'
   Col_2b='medium blue'
   Col_3a='olive'
   Col_3b='olive drab'
   SymSize=0.02
   SymIncrement=20
   DisplayText='Displaying 1 out of every '+NUM2STR(SymIncrement,0)+' points per series'
   Dummy = LABEL_DATE(DATE_FORMAT='%N/%Z')
   XTitle='Time [MM/YY]'
   XRange=[MIN(gosat_jd)-60,MAX(gosat_jd)+60]
   OutPlotFile=collocate_dir+'dxco2_with_ak-prior_correction.pdf'
   YRange=[-2.0, 2.0]

   ;;----- XData is GOSAT Julian Day.
   XD=gosat_jd

   ;;----- YData is delta XCO2.
   ;;----- Can be using either unharmonized or harmonized XCO2.
   IF (VerNum_GOSAT EQ 'v9c_20240214') $
      THEN YD_1=gosat_xco2_harm_w - oco_mean_struct.xco2 $
      ELSE YD_1=Matchw.gosat.xco2 - oco_mean_struct.xco2
   
   ;;----- Create initial plot with no data.
   p=SCATTERPLOT(/NODATA,XD,YD_1, $
           Title='dXCO$_2$ with AK/prior correction', $
           XTitle=XTitle, XStyle=1, XRange=XRange, XTickFormat='LABEL_DATE', $
           YTitle='dXCO$_2$ (GOSAT $-$ OCO-2) [ppm]', YRange=YRange, $
           FONT_SIZE=15, FONT_STYLE=1 $
          )

   ;;----- Overplot original dXCO2
   IF (1) THEN BEGIN   
      text_1a='Unmodified dXCO$_2$'
      PRINT,'Scatter plot data...'
      p=SCATTERPLOT(OVERPLOT=p,XD,YD_1)
      p.SYM_INCREMENT = SymIncrement
      p.SYM_COLOR = Col_1a
      p.SYM_SIZE=SymSize
      p.SYMBOL=Symbol

      ;;----- Linear fit to data to get trend.
      PRINT,'Fit data...'
      Fit=LINFIT(XD,YD_1,YFit=YFit)
      Slope=Fit[1]*365. ;; per year
      Mean=Mean(YD_1)
      Sigma=STDDEV(YD_1)
      R=CORRELATE(XD,YD_1)
      Text_1b='Linear trend: '+NUM2STR(Slope,3)+' ppm/year, Bias='+NUM2STR(Mean,2)+' $\pm$ '+NUM2STR(Sigma,2)+' ppm'
      p=PLOT(OVERPLOT=p,XD,YFit,COLOR=Col_1b,Thick=5)

   ENDIF

   ;;----- Overplot prior and AK adjusted dXCO2
   IF (1) THEN BEGIN   
      text_2a='dXCO$_2$ with 2-term CO2 prior + AK adjustments'
      ;;YD_2=(Matchw.gosat.xco2 + xco2_prior_adj_a + xco2_prior_adj_b) - oco_mean_struct.xco2
      YD_2=YD_1 + xco2_prior_adj_a + xco2_prior_adj_b
      PRINT,'Scatter plot data...'
      p=SCATTERPLOT(OVERPLOT=p,XD,YD_2)
      p.SYM_INCREMENT = SymIncrement
      p.SYM_COLOR = Col_2a
      p.SYM_SIZE=SymSize
      p.SYMBOL=Symbol

      ;;----- Linear fit to data to get trend.
      PRINT,'Fit data...'
      Fit=LINFIT(XD,YD_2,YFit=YFit)
      Slope=Fit[1]*365. ;; per year
      Mean=MEAN(YD_2)
      Sigma=STDDEV(YD_2)
      R=CORRELATE(XD,YD_2)
      Text_2b='Linear trend: '+NUM2STR(Slope,3)+' ppm/year, Bias='+NUM2STR(Mean,2)+' $\pm$ '+NUM2STR(Sigma,2)+' ppm'
      p=PLOT(OVERPLOT=p,XD,YFit,COLOR=Col_2b,Thick=5)

   ENDIF

   ;;----- Add the legend.
   FontSize=14
   t=TEXT(TARGET=p,0.17,0.82,DisplayText,FONT_SIZE=FontSize+2,FONT_COLOR='black')

   t=TEXT(TARGET=p,0.17,0.77,text_1a,FONT_SIZE=FontSize,FONT_COLOR=Col_1a)
   t=TEXT(TARGET=p,0.17,0.72,text_2a,FONT_SIZE=FontSize,FONT_COLOR=Col_2a)

   t=TEXT(TARGET=p,0.17,0.22,Text_1b,FONT_SIZE=FontSize,FONT_COLOR=Col_1b)
   t=TEXT(TARGET=p,0.17,0.17,Text_2b,FONT_SIZE=FontSize,FONT_COLOR=Col_2b)

   ;;----- Overplot horizontal zero lines
   p=PLOT(OVERPLOT=p,XRange,[0.,0.],Thick=2,Color=Black)

   ;;----- Save plot to output file.
   PRINT,'Saving plot to file ',OutPlotFile
   p.save,OutPlotFile

   ;STOP

ENDIF

;;-----
;;----- Sanity check plot of the calculated xco2 adjustment due to prior.
;;-----
AvgKerPlot=0
IF (AvgKerPlot) THEN BEGIN

   ;;----- Plot settings.
   Symbol='Circle'
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
   

   SymSize=0.02
   SymIncrement=20
   DisplayText='Displaying 1 out of every '+NUM2STR(SymIncrement,0)+' points per series'
   Dummy = LABEL_DATE(DATE_FORMAT='%N/%Z')
   XTitle='Time [MM/YY]'
   XRange=[MIN(gosat_jd)-60,MAX(gosat_jd)+60]
   OutPlotFile=collocate_dir+'xco2_prior_adj.pdf'
   YRange=[-0.5,1.5]

   ;;----- Create initial plot with no data.
   p=SCATTERPLOT(/NODATA,gosat_jd,xco2_prior_adj_a, $
           Title='Calculated xco2 adjustment due to prior', $
           XTitle=XTitle, XStyle=1, XRange=XRange, XTickFormat='LABEL_DATE', $
           YTitle='XCO$_2$ correction [ppm]', YRange=YRange, $
           FONT_SIZE=15, FONT_STYLE=1 $
          )

   ;;----- Overplot xco2_prior_adj_a
   IF (1) THEN BEGIN   
      text_1a='Term 1: CO$_2$ prior adj. via Eq. A10 [Wunch, ACP, 2010]'
      text_1aa=' (equivalent to Eq. 3 [Taylor, ESSD, 2023])'
      XD=gosat_jd
      YD=xco2_prior_adj_a
      PRINT,'Scatter plot data (a)...'
      p=SCATTERPLOT(OVERPLOT=p,XD,YD)
      p.SYM_INCREMENT = SymIncrement
      p.SYM_COLOR = Col_1a
      p.SYM_SIZE=SymSize
      p.SYMBOL=Symbol

      ;;----- Linear fit to data to get trend.
      PRINT,'Fit data (a)...'
      Fit=LINFIT(XD,YD,YFit=YFit)
      Slope=Fit[1]*365. ;; per year
      Sigma=STDDEV(YD)
      R=CORRELATE(XD,YD)
      Text_1b='Linear trend: '+NUM2STR(Slope,3)+' ppm/year, 1-$\sigma$='+NUM2STR(Sigma,2)+' ppm'
      p=PLOT(OVERPLOT=p,XD,YFit,COLOR=Col_1a,Thick=5)

   ENDIF

   ;;----- Overplot xco2_prior_adj_b   
   IF (1) THEN BEGIN   
      text_2a='Term 2: AK "smoothing" term'
      XD=gosat_jd
      YD=xco2_prior_adj_b
      PRINT,'Scatter plot data (b)...'
      p=SCATTERPLOT(OVERPLOT=p,XD,YD)
      p.SYM_INCREMENT = SymIncrement
      p.SYM_COLOR = Col_2a
      p.SYM_SIZE=SymSize
      p.SYMBOL=Symbol

      ;;----- Linear fit to data to get trend.
      PRINT,'Fit data (b)...'
      Fit=LINFIT(XD,YD,YFit=YFit)
      Slope=Fit[1]*365. ;; per year
      R=CORRELATE(XD,YD)
      Sigma=STDDEV(YD)
      Text_2b='Linear trend: '+NUM2STR(Slope,3)+' ppm/year, 1-$\sigma$='+NUM2STR(Sigma,2)+' ppm'
      p=PLOT(OVERPLOT=p,XD,YFit,COLOR=Col_2a,Thick=5)

   ENDIF

   ;;----- Overplot sum of terms  
   IF (1) THEN BEGIN   
      text_3a='Sum of terms = T1 + T2'
      XD=gosat_jd
      YD=xco2_prior_adj_a+xco2_prior_adj_b
      PRINT,'Scatter plot data (sum)...'
      p=SCATTERPLOT(OVERPLOT=p,XD,YD)
      p.SYM_INCREMENT = SymIncrement
      p.SYM_COLOR = Col_3a
      p.SYM_SIZE=SymSize
      p.SYMBOL=Symbol

      ;;----- Linear fit to data to get trend.
      PRINT,'Fit data (sum)...'
      Fit=LINFIT(XD,YD,YFit=YFit)
      Slope=Fit[1]*365. ;; per year
      R=CORRELATE(XD,YD)
      Sigma=STDDEV(YD)
      Text_3b='Linear trend: '+NUM2STR(Slope,3)+' ppm/year, 1-$\sigma$='+NUM2STR(Sigma,2)+' ppm'
      p=PLOT(OVERPLOT=p,XD,YFit,COLOR=Col_3a,Thick=5)

   ENDIF

   ;;----- Add the legend.
   FontSize=12
   t=TEXT(TARGET=p,0.17,0.82,DisplayText,FONT_SIZE=FontSize+2,FONT_COLOR='black')

   t=TEXT(TARGET=p,0.17,0.77,text_1a,FONT_SIZE=FontSize,FONT_COLOR=Col_1a)
   t=TEXT(TARGET=p,0.17,0.72,text_1aa,FONT_SIZE=FontSize,FONT_COLOR=Col_1a)
   t=TEXT(TARGET=p,0.17,0.67,xco2_prior_adj_a_text,FONT_SIZE=FontSize*0.9,FONT_COLOR=Col_1a)
   t=TEXT(TARGET=p,0.17,0.62,text_2a,FONT_SIZE=FontSize,FONT_COLOR=Col_2a)
   t=TEXT(TARGET=p,0.17,0.57,xco2_prior_adj_b_text,FONT_SIZE=FontSize*0.9,FONT_COLOR=Col_2a)
   t=TEXT(TARGET=p,0.17,0.52,text_3a,FONT_SIZE=FontSize,FONT_COLOR=Col_3a)

   t=TEXT(TARGET=p,0.17,0.25,Text_1b,FONT_SIZE=FontSize,FONT_COLOR=Col_1a)
   t=TEXT(TARGET=p,0.17,0.20,Text_2b,FONT_SIZE=FontSize,FONT_COLOR=Col_2a)
   t=TEXT(TARGET=p,0.17,0.15,Text_3b,FONT_SIZE=FontSize,FONT_COLOR=Col_3a)

   ;;----- Overplot horizontal zero lines
   p=PLOT(OVERPLOT=p,XRange,[0.,0.],Thick=2,Color=Black)

   ;;----- Save plot to output file.
   PRINT,'Saving plot to file ',OutPlotFile
   p.save,OutPlotFile
   STOP

ENDIF

;;---- Additional filter on OCO-2 data for variability over the scene.
PRINT, 'Determining overall (per match) mask.'
match_mask =  (oco_stddev_xco2 LE oco_max_xco2_variability) ;AND (mean_dist_g LE max_dist) AND (gosat_qf[w])
Mask_match=WHERE(match_mask, NMask_variability)
PRINT,'NMask_variability=',NMask_variability
PRINT,'Percent of total retained by OCO-2 XCO2 variability check=',FLOAT(NMask_variability)/FLOAT(N_ELEMENTS(match_mask))*100. 
;STOP

;;-----
;;----- Build mask for current observation mode.
;;-----
PRINT, 'Developing mask for mode ' + view_mode
MGAIN=77 & HGAIN=72
;MGAIN='M' & HGAIN='H'
CASE view_mode OF
   'landH'  : type_mask = gosat_land AND (gosat_gain EQ HGAIN)
   'landM'  : type_mask = gosat_land AND (gosat_gain EQ MGAIN)
   'land'   : type_mask = gosat_land
   'oceanH' : type_mask = ~gosat_land AND (gosat_gain EQ HGAIN)
   'allH'   : type_mask = gosat_gain eq HGAIN
   'all'    : type_mask = (~gosat_land AND gosat_gain eq HGAIN) OR gosat_land
ENDCASE
Mask_type=WHERE(type_mask, NMask_surface)
PRINT,'Percent of total retained by surface/observation mask =',FLOAT(NMask_surface)/FLOAT(N_ELEMENTS(type_mask))*100. 

;;----- Mask for latitude range.
lat_mask=INRANGE(gosat_lat,LatRange)
Mask_last=WHERE(lat_mask, NMask_lat)
PRINT,'Percent of total retained by latitude mask =',FLOAT(NMask_lat)/FLOAT(N_ELEMENTS(lat_mask))*100. 

;;----- Combine submasks
main_mask = match_mask AND type_mask AND lat_mask
g = WHERE(main_mask, nsound)
wg=w[g]
PRINT, 'Using ' + sc(nsound) + ' matched soundings.'
PRINT,'Percent of total retained by main mask =',FLOAT(NSound)/FLOAT(N_ELEMENTS(main_mask))*100. 
;;STOP

;;-----
;;----- Extract the needed data variables.
;;-----
matchg = matchw[g]
mean_dist_g=mean_ddist[g]
mean_dtime_g=mean_dtime[g]
oco_meang = oco_mean_struct[g]
xco2_bc_diff = gosat_xco2_bc[wg] - oco_meang.xco2 
xco2_harm_diff = gosat_xco2_harm_w[g] - oco_meang.xco2
;IF (correct_different_priors) THEN BEGIN
;   xco2_bc_diff = xco2_bc_diff + xco2_prior_adj_a[g] + xco2_prior_adj_b[g]
;   xco2_harm_diff = xco2_harm_diff + xco2_prior_adj_a[g] + xco2_prior_adj_b[g] 
;ENDIF

tdiff = mean_dtime[g]
lat=gosat_lat[wg] & lon=gosat_lon[wg]
year = fix(strmid(sc(gosat_id[wg]),0,4))
month = fix(strmid(sc(gosat_id[wg]),4,2))
season = (month mod 12)/3
l = where(gosat_land, comp=o)
julian_day=gosat_jd[wg]
land_flag=gosat_land[wg]

;;----- If skip_calc is set or code has not been reset.
skip_calc_:

;;-----
;;----- Write Summary file of collocated GOSAT/OCO2 soundings.
;;-----
;;STOP
WriteSummary=0
IF (WriteSummary) THEN BEGIN

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
         mean_dist_g[III], $
         mean_dtime_g[III], $
         matchg[III].oco2[0:9].sounding_id, $
         FORMAT=SummaryFormat
      ENDIF
   ENDFOR
   FREE_LUN,SummaryLUN

   ;;-----
   ;;----- Plot histograms of collocation time.
   ;;-----
   XMin=-1.5*60. & XMax=2.*60. & XTickLen=0.05
   loc = RANGE(XMin, XMax, 41)
   charsize=2. & charthick=2 & XCharSize=1.75 & YCharSize=1.75
   Title=' GOSAT (v9) - OCO2 (v' + VerNum_OCO +') :: '+View_Mode+'-'+Region
   ;;Title=SubTitle+STRING(CurrentYear,FORMAT='(I4)')+', '+SeasList[s]
   XTitle='Mean collocation time difference (minutes)'
   Thick=12 & BChar=2.
   YData=mean_dtime_g*60.
   HIST, YData, loc=loc, Title=Title, xtit=XTitle, ytit='Number of Collocations', xr=[XMin,XMax], XTickLen=XTickLen, $
      Position=[0.25,0.20,0.95,0.92], CharSi=BChar,XCharSize=XCharSize,YCharSize=YCharSize,CharSize=CharSize,CharThick=CharThick,YMinor=0
   VLINE,MEAN(YData),LineStyle=2,Thick=5
   STOP

   ;;-----
   ;;----- Plot histograms of collocation distance.
   ;;-----
   XMin=-50. & XMax=350. & XTickLen=0.05
   loc = RANGE(XMin, XMax, 41)
   charsize=2. & charthick=2 & XCharSize=1.75 & YCharSize=1.75
   Title=' GOSAT (v9) - OCO2 (v' + VerNum_OCO +') :: '+View_Mode+'-'+Region
   ;;Title=SubTitle+STRING(CurrentYear,FORMAT='(I4)')+', '+SeasList[s]
   XTitle='Mean collocation distance (km)'
   Thick=12 & BChar=2.
   YData=mean_dist_g
   HIST, YData, loc=loc, Title=Title, xtit=XTitle, ytit='Number of Collocations', xr=[XMin,XMax], XTickLen=XTickLen, $
      Position=[0.25,0.20,0.95,0.92], CharSi=BChar,XCharSize=XCharSize,YCharSize=YCharSize,CharSize=CharSize,CharThick=CharThick,YMinor=0
   VLINE,MEAN(YData),LineStyle=2,Thick=5
   STOP

ENDIF



;;----- Peform plotting.
IF (1) THEN BEGIN

   OutDir=OutDirBase+'gosat_'+VerNum_gosat+'_vs_oco2_v'+VerNum_oco+'/' & FILE_MKDIR,OutDir

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

   var_list=['dxco2']
   IF (correct_different_priors) THEN BEGIN
      PRINT,'Implement the correction to XCO2 to account for different CO2 priors between L2FP builds and the AK.'
      ZData=xco2_bc_diff ;+ xco2_prior_adj_a[g] + xco2_prior_adj_b[g]
      IF (VerNum_GOSAT EQ 'v9b') THEN BEGIN
         ZData=[[ZData],[gosat_xco2_adj_w[g]]]
         var_list=[var_list,'xco2_adjust']
      ENDIF
      IF (VerNum_GOSAT EQ 'v9c_20240214') THEN BEGIN
         ZData=[ [ZData], [gosat_xco2_adj_w[g]], [xco2_harm_diff] ]
         var_list=[var_list,'xco2_adjust','dxco2_harm']
      ENDIF
      ;;STOP
   ENDIF ELSE BEGIN
      ZData=xco2_bc_diff
      IF (VerNum_GOSAT EQ 'v9b') THEN BEGIN
         ZData=[[ZData],[gosat_xco2_adj[g]]]
         var_list=[var_list,'xco2_adjust']
      ENDIF
   ENDELSE      

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
                                  OR ( INRANGE(Month, [1,2])  AND INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]) ), NThis )
                  'MAM' : ws=WHERE( INRANGE(Month, [3,5])   AND INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]), NThis) 
                  'JJA' : ws=WHERE( INRANGE(Month, [6,8])   AND INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]), NThis) 
                  'SON' : ws=WHERE( INRANGE(Month, [9,11]) AND INRANGE(Year,[MIN(CurrentYear),MAX(CurrentYear)]), NThis) 
            ENDCASE            

         ENDIF ELSE BEGIN
            ws=WHERE( INRANGE(Year, [MIN(CurrentYear),MAX(CurrentYear)]) AND (Season LE 3) , nthis)
         ENDELSE
         PRINT,'N Elements in Season/Year='+SC(NThis)
         ;;STOP

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
            OutInfo[3]='N='+STRING(FLOAT(NThis)/1000.,FORMAT='(F6.1)')+'k (SS)'
            OutInfo[4]=''
            OutInfo[10]='plasma'
            OutInfo[11]='Longitude'
            OutInfo[12]='Latitude'
            OutInfo[13]=STRING(FillFloat,FORMAT='(F6.1)')
            OutInfo[14]=STRING(FillInteger,FORMAT='(I4)')
            ;;STOP
            
            ;;----- Subset the ZData to 'ws' elements.
            NZ=N_ELEMENTS(ZData[0,*])
            ZData_sub=MAKE_ARRAY([NThis,NZ], /FLOAT)
            FOR ZZZ=0,NZ-1 DO BEGIN
               ZT=ZData[ws,ZZZ] 
               ZData_sub[*,ZZZ]=ZT
            ENDFOR

            NX=LONG(360./DeltaLon_1)
            NY=LONG((LatRange[1]-LatRange[0])/DeltaLat_1)
            ;;STOP
            PRINT,'GRID_DATA_WRITE_NETCDF_MULTI for spatial map of collocated sounding density...'
            GRID_DATA_WRITE_NETCDF_MULTI, $
                  Lon[ws], $
                  Lat[ws], $
                  ZData_sub, $  
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
         ENDIF
         STOP




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
                  z1=zdata[w_hov,0] & z2=zdata[w_hov,1] & zdata=[[z1],[z2]]

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
                     ZData, $ ;; inputs
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
                                , nthis) 
                  IF ( ViewModes[Mode] EQ 'Land-M') $
                  THEN ws = WHERE( (Season EQ s) AND (gosat_land) AND (MatchG.gosat.Sounding.gain EQ 'M') $
                               AND (INRANGE(lat,CurrentLatRange)) $
                               AND (INRANGE(lon,CurrentLonRange)) $
                                , nthis) 
                  IF ( ViewModes[Mode] EQ 'Ocean' ) $
                  THEN ws = WHERE( (Season EQ s) AND (INRANGE(lat,CurrentLatRange)) AND (INRANGE(lon,CurrentLonRange)) AND (~gosat_land), nthis) 
                  
                  IF (nthis GT 2) THEN BEGIN
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
    
