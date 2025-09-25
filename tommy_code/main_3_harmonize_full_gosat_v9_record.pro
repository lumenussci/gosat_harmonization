;; PRO MAIN_3_HARMONIZE_FULL_GOSAT_V9_RECORD

;;----- 
;;----- T.E. Taylor 4-Jan-2023.
;;----- Use the previously determined fit coefficients to bias correct the GOSAT v9 Lite file record.
;;----- The fit coefficients were determined by <main_2_fit_and_plot_collocated_gosat_to_oco2.pro>...<plot_dxco2_vs_bc_vars.pro>.
;;----- Fit coef saved to IDL SAVE files, e.g., 
;;----- /home/ttaylor/analysis_utilities/tropical_iav/data/gosat_vs_oco2_v10/bc_fit_coefficients/bc_fit_coefficients_AOD_TOTAL_Global.SAV
;;-----
;;----- The fit coef were determined from best fits to the dxco2 (GOSAT v9 minus OCO-2 v10) verses the GOSAT v9 L2FP TOTAL AOD.
;;----- For Land-H and Ocean-H, quadratic fits were used. For Land-M a linear fit was used.
;;-----
;;----- An example of the additional XCO2 bias correction to the GOSAT Lite file bias corrected XCO2 is:
;;-----    XD=XData[ws] & YD=YData[ws]
;;-----    Fit2 = POLY_FIT( XD, YD, 2)
;;-----    xco2_adjustment2= Fit2[0] + Fit2[1]*XD + Fit2[2]*XD^2.
;;-----    dxco2_2=(gosat_xco2-xco2_adjustment2) - oco_xco2
;;-----
;;----- T.E.T. 25-Mar-2024.
;;----- Many updates to <main_2_fit_and_plot_collocated_gosat_to_oco2.pro> over past weeks.
;;----- Generated a new set of fit coefficients and ranked the variables to determine most powerful correctors.
;;----- 
;;----- Land-H: xco2_raw, aod_fine, aod_strataer, aod_total, aod_ice
;;----- Land-M: aod_strataer, dust_height, aod_fine, aod_sulfate, eof2_3
;;----- Ocean-H: aod_fine, albedo_slope_sco2, albedo_slope_wco2, aod_sulfate, albedo_o2a
;;-----
;;----- The primary steps in this code are:
;;----- (0) make a copy of the original lite files to be modified (only needs to be done once),
;;----- (1) read in the full GOSAT v9 Lite file record granule by granule, 
;;----- (2) extract the bias corrected XCO2, operation mode (H/M), and surface type (land/water),
;;----- (3) apply the appropriate set of mode-dependent fit coefficients as: 
;;-----        xco2_adjustment2= Fit2[0] + Fit2[1]*XD + Fit[2]*XD^2 + Bias (for quadratic fits)
;;-----        xco2_adjustment2= Fit2[0] + Fit2[1]*XD + Bias (for linear fit)
;;-----     where XD is the XCO2
;;----- (4) bias correct the GOSAT v9 data as:
;;-----       new_gosat_xco2 = old_gosat_xco2 + xco2_adjustment (new 25-Mar-2024, the xco2_adjustment is an additive quantity)
;;----- (5) copy the old GOSAT v9 Lite file
;;----- (6) overwrite the xco2 data field with the newly bias corrected xco2 values
;;-----
;;----- After this procedure, run <analysis_utilities/generate_lite_save_files/generate_lite_save.pro> to aggregate Lite files into IDL.SAV.

;;----- Input version number
VerNum_in='v9c3_r20240422'

;;-------------------------------------------------------------------------
;;----- (0) make a copy of the original lite files to be modified (only needs to be done once),
;;-------------------------------------------------------------------------
CopyLiteFiles=0
IF (CopyLiteFiles) THEN BEGIN

   InLiteDir='/data6/GOSAT/product/Lite/B9/*.nc4'
   InLiteFiles=FILE_SEARCH(InLiteDir,COUNT=NLiteFiles)
   PRINT,'Found ',STRING(NLiteFiles,FORMAT='(I4)'),' Lite files in data directory ',InLiteDir
   PRINT,'Making a copy of the original Lite files...'

   ;;----- 10-May-2023. For unknown reason, I can no longer make directories, or copy files, or modify files on this disk.
   ;OutLiteDir='/data6/GOSAT/product/Lite/B9_bc_to_oco2_v10/'

   ;;----- 10-May-2023. So...copy them to /data8/ttaylor/data_ttaylor/...
   ;OutLiteDir='/data8/ttaylor/data_ttaylor/GOSAT/product/Lite/v9b_r20230725/'

   ;;----- 14-Feb-2024. Creating a second iteration of the GOSAT v9b data product.
   ;;----- Call it v9c to preserve the earlier version.
   ;OutLiteDir='/data8/ttaylor/data_ttaylor/GOSAT/product/Lite/v9c_r20240214/'

   ;;----- 26-Mar-2024. Creating a third iteration of the GOSAT v9 data product.
   ;;----- Lots of changes to <main_2_fit_and_plot_collocated_gosat_to_oco2.pro> to derive fit coefficients and rank variables.
   ;;----- Call it v9c to preserve the earlier version.
   ;OutLiteDir='/data8/ttaylor/data_ttaylor/GOSAT/product/Lite/v9c2_r20240329/'

   ;;----- 22-April-2024. Creating a fourth iteration of the GOSAT v9 data product.
   ;;----- Applied only the single-term AK correction (no CO2 prior correction) during fitting.
   ;;----- v9c3_r20240422 to preserve the earlier version.
   OutLiteDir='/data8/ttaylor/data_ttaylor/GOSAT/product/Lite/'+VerNum_in+'/'

   FILE_MKDIR,OutLiteDir
   FILE_COPY,InLiteFiles,OutLiteDir,/OVERWRITE,/VERBOSE
   ;STOP   

ENDIF

;;----- GOSAT observation modes.
ViewModes=['Land-H','Land-M','Ocean-H'] & NumViewModes=N_ELEMENTS(ViewModes)

;;----- Read in the fit coefficients from IDL SAVE files.
;FitCoefFile='/home/ttaylor/analysis_utilities/tropical_iav/data/gosat_vs_oco2_v10/bc_fit_coefficients/bc_fit_coefficients_AOD_TOTAL_Global.SAV'
;FitCoefFileDir='/home/ttaylor/analysis_utilities/tropical_iav/data/gosat_vs_oco2_v11.1/bc_fit_coefficients/'
FitCoefFileDir='/home/ttaylor/analysis_utilities/tropical_iav/data/gosat_vs_oco2_v11.1/bc_fit_coefficients/'+VerNum_in+'/'
FitCoefFiles=FitCoefFileDir + [ $
                'bc_fit_coefficients_JULIAN_DAY_Global.SAV', $     ;; Land-M
                'bc_fit_coefficients_JULIAN_DAY_Global.SAV', $ ;; Land-H
                'bc_fit_coefficients_JULIAN_DAY_Global.SAV' $      ;; Ocean-H
                ]

AllFits=!NULL
FOR mode=0,NumViewModes-1 DO BEGIN
   RESTORE,FILENAME=FitCoefFiles[Mode],VERBOSE=1
   CASE mode OF
      0: AllFits=FitStructure[mode]
      1: AllFits=[AllFits,FitStructure[mode]]
      2: AllFits=[AllFits,FitStructure[mode]]
   ENDCASE
   FitStructure=!NULL
ENDFOR
PRINT,AllFits[0].FitCoef
;STOP

;;-------------------------------------------------------------------------
;;----- (1) read in the full GOSAT v9 Lite file record granule by granule, 
;;-------------------------------------------------------------------------
BiasCorrectRecord=1
IF (BiasCorrectRecord) THEN BEGIN

   InLiteDir='/data8/ttaylor/data_ttaylor/GOSAT/product/Lite/'+VerNum_in+'/*.nc4'
   InLiteFiles=FILE_SEARCH(InLiteDir,COUNT=NLiteFiles)
   PRINT,'Found ',STRING(NLiteFiles,FORMAT='(I4)'),' Lite files in data directory ',InLiteDir
   ;;STOP

   FOR FF=0L,NLiteFiles-1L DO BEGIN
   ;;FOR FF=0L,0L DO BEGIN

      PRINT,'Processing file ',STRTRIM(STRING(FF+1,FORMAT='(I4)'),2),' of ',STRING(NLiteFiles,FORMAT='(I4)')

      ;;----- Read gain and surface data fields used to determine observation mode.
      surface=READ_H5_FIELD(InLiteFiles[FF],'Retrieval/surface_type') ;; integer 0=ocean, 1=land
      gain=READ_H5_FIELD(InLiteFiles[FF],'Sounding/gain') ;; string H/M

      ;;----- Read xco2 data field. Generate output harmonized and adjusted arrays.
      xco2=READ_H5_FIELD(InLiteFiles[FF],'xco2') ;; float [ppm]
      xco2_harmonized=MAKE_ARRAY(N_ELEMENTS(xco2),/FLOAT,VALUE=-999.999)
      xco2_adjustment=MAKE_ARRAY(N_ELEMENTS(xco2),/FLOAT,VALUE=-999.999)

      ;;----- Read the appropriate variable data to use in adjustment/harmonization procedure. 
      IF (AllFits[0].FitType EQ 'skewed-sine') THEN BEGIN

         ;;----- Assumption here that if the skewed-sine fit against julian day was
         ;;----- used for one observation mode, it was used for all 3 modes.
         sid=READ_H5_FIELD(InLiteFiles[FF],'sounding_id')
         jd=ACOS_ID_TO_JD(sid)
         FOR mode=0,NumViewModes-1 DO BEGIN
            ;xdata=jd 
            CASE mode OF
               0: $
                  BEGIN
                     j0=AllFits[0].JD0
                     jd_0=jd-j0
                     AllX=jd_0
                  END
               1: $
                  BEGIN
                     j0=AllFits[1].JD0
                     jd_0=jd-j0
                     AllX=[[AllX],[jd_0]]
                  END
               2: $
                  BEGIN
                     j0=AllFits[2].JD0
                     jd_0=jd-j0
                     AllX=[[AllX],[jd_0]]
                  END
            ENDCASE
         ENDFOR

      ENDIF ELSE BEGIN
         variable_names='Retrieval/'+STRLOWCASE(AllFits.variable)
         FOR mode=0,NumViewModes-1 DO BEGIN
            xdata=READ_H5_FIELD(InLiteFiles[FF],variable_names[mode]) 
            CASE mode OF
               0: AllX=xdata
               1: AllX=[[AllX],[xdata]]
               2: AllX=[[AllX],[xdata]]
            ENDCASE
         ENDFOR
      ENDELSE
      ;STOP
      ;aod_total=READ_H5_FIELD(InLiteFiles[FF],'Retrieval/aod_total')  ;; float

      ;;----- bias correct land-H.
      FitIndex=0
      Fit=AllFits[FitIndex]
      land_h=WHERE( (gain EQ 'H') AND (surface EQ 1), n_land_h)
      XD=AllX[land_h,FitIndex]

      IF (Fit.FitType EQ 'skewed-sine') THEN BEGIN
         SKEWED_SINE, XD, Fit.SkewSine_A, result
         xco2_adjustment[land_h]= (-1) * result
         xco2_harmonized[land_h]=xco2[land_h] + xco2_adjustment[land_h]
         ;STOP
      ENDIF ELSE BEGIN
         xco2_adjustment[land_h]= (-1) * (Fit.FitCoef[0] + Fit.FitCoef[1]*XD + Fit.FitCoef[2]*XD^2.)
         xco2_harmonized[land_h]=xco2[land_h] + xco2_adjustment[land_h]
      ENDELSE
      ;;STOP

      ;;----- bias correct land-M.
      FitIndex=1
      Fit=AllFits[FitIndex]
      land_m=WHERE( (gain EQ 'M') AND (surface EQ 1), n_land_m)
      XD=AllX[land_m,FitIndex]
      IF (Fit.FitType EQ 'skewed-sine') THEN BEGIN
         SKEWED_SINE, XD, Fit.SkewSine_A, result
         xco2_adjustment[land_m]= (-1) * result
         xco2_harmonized[land_m]=xco2[land_m] + xco2_adjustment[land_m]
         ;STOP
      ENDIF ELSE BEGIN
         xco2_adjustment[land_m]= (-1) * (Fit.FitCoef[0] + Fit.FitCoef[1]*XD + Fit.FitCoef[2]*XD^2.)
         xco2_harmonized[land_m]=xco2[land_m] + xco2_adjustment[land_m]
      ENDELSE

      ;;----- bias correct ocean-H.
      FitIndex=2
      Fit=AllFits[FitIndex]
      ocean_h=WHERE( (gain EQ 'H') AND (surface EQ 0), n_ocean_h)
      XD=AllX[ocean_h,FitIndex]
      IF (Fit.FitType EQ 'skewed-sine') THEN BEGIN
         SKEWED_SINE, XD, Fit.SkewSine_A, result
         xco2_adjustment[ocean_h]= (-1) * result
         xco2_harmonized[ocean_h]=xco2[ocean_h] + xco2_adjustment[ocean_h]
         ;STOP
      ENDIF ELSE BEGIN
         xco2_adjustment[ocean_h]= (-1) * (Fit.FitCoef[0] + Fit.FitCoef[1]*XD + Fit.FitCoef[2]*XD^2.)
         xco2_harmonized[ocean_h]=xco2[ocean_h] + xco2_adjustment[ocean_h]
      ENDELSE

      ;;----- .compile /usr/local/harris/idl/lib/graphics/legend.pro
      ;;----- Sanity plot of the XCO2 adjustment vs the orginal BC XCO2.
      IF (FF MOD 1000 EQ 0) THEN BEGIN
         PlotDim=[500,400] & Margin=[0.15,0.1,0.1,0.1]
         p = PLOT(xco2,xco2_adjustment, /NODATA, $
                  YTITLE='XCO$_2$ adjustment', $ 
                  Title=FILE_BASENAME(InLiteFiles[FF]), $
                  XTITLE='GOSAT v9 bias corrected XCO$_2$', $
                  DIM=PlotDim, MARGIN=Margin)
         p=PLOT(OVERPLOT=p,p.xrange,[0.,0.],Thick=2,Color='black')
         p1=PLOT(OVERPLOT=p,xco2[land_h],xco2_adjustment[land_h],'go',NAME='Land-H')
         p2=PLOT(OVERPLOT=p1,xco2[land_m],xco2_adjustment[land_m],'rd',NAME='Land-M')
         p3=PLOT(OVERPLOT=p2,xco2[ocean_h],xco2_adjustment[ocean_h],'b+',NAME='Ocean-H')
         leg=LEGEND(TARGET=[p1,p2,p3],/AUTO_TEXT_COLOR,POSITION=[0.45,0.85],/NORM)
         ;STOP
      ENDIF

      ;;----- Sanity plot of the harmonized XCO2 vs the orginal BC XCO2.
      IF (FF MOD 1000 EQ 0) THEN BEGIN
         p = PLOT(xco2,xco2_harmonized, /NODATA, $
                  YTITLE='Harmonized XCO$_2$', $ 
                  Title=FILE_BASENAME(InLiteFiles[FF]), $
                  XTITLE='GOSAT v9 bias corrected XCO$_2$', $
                  DIM=PlotDim, MARGIN=Margin)
         p=PLOT(OVERPLOT=p,p.xrange,p.yrange,Thick=2,Color='black')
         p1=PLOT(OVERPLOT=p,xco2[land_h],xco2_harmonized[land_h],'go',NAME='Land-H')
         p2=PLOT(OVERPLOT=p1,xco2[land_m],xco2_harmonized[land_m],'rd',NAME='Land-M')
         p3=PLOT(OVERPLOT=p2,xco2[ocean_h],xco2_harmonized[ocean_h],'b+',NAME='Ocean-H')
         leg=LEGEND(TARGET=[p1,p2,p3],/AUTO_TEXT_COLOR,POSITION=[0.45,0.85],/NORM)
         ;;STOP
      ENDIF

      ;;----- Create save file of the adjusted xco2 values
      ;FileSave='./adjusted_xco2_v10.sav'
      ;xco2_v10=xco2 & xco2_harmonized_v10=xco2_harmonized
      ;SAVE,FILENAME=FileSave,xco2_v10, xco2_harmonized_v10

      ;FileSave='./adjusted_xco2_v11.sav'
      ;xco2_v11=xco2 & xco2_harmonized_v11=xco2_harmonized
      ;SAVE,FILENAME=FileSave,xco2_v11, xco2_harmonized_v11
      ;STOP

      ;;----- Create a new data field in the Lite files called "xco2_harmonized_to_oco2_vxx".
      PRINT,'Writing new data fields (xco2_adjustment and xco2_harmonized) to L2Lite file.'
      var={xco2_harmonized_to_oco2_v11_1:xco2_harmonized, $
           xco2_additive_adjustment_to_gosat_v9:xco2_adjustment, $
           harmonization_fit_modes:AllFits.ViewMode, $
           harmonization_fit_type:AllFits.FitType, $
           harmonization_fit_variable:AllFits.Variable $
           }
      NCDF_PUT,InLiteFiles[FF],VARIABLES=var

      ;STOP

   ENDFOR
ENDIF

STOP
END