;; .compile main_1_match_gosat_to_oco2.pro

; This version reads daily lite files directly
;;----- Copied from /home/codell/idlprogs/aos/models/OCO/l2/b9_acos/match_b9_acos_b10_oco2_v2.pro
;;-----
;;----- 9-Nov-2023. Update to include OCO2 co2_profile in output. Also, change "acos" to "gosat".
;;----- Current output directory is "/data8/ttaylor/data_ttaylor/gosat_oco2_collocations"
;;----- Remove defunct data directory "/data8/ttaylor/data_ttaylor/acos/acos_results/b9_acos/match_oco2/"
;;-----
;;----- 7-Jan-2025. TET. Some beautification and commenting.
;;-----

oco2dir  = '/data11/OCO2/product/Lite/B11.1/LtCO2/'
gosatdir= '/data6/GOSAT/product/Lite/B9/'

outfile = '/data8/ttaylor/data_ttaylor/gosat_oco2_collocations/match_gosat_v9_oco2_v11.1_20140906_20200630_LiteDirect_time2_lat2_lon3_min3_dist300_random.sav'

;;---- Search for the input GOSAT Lite file listing and determine the dates (6 digit YYMMDD string values).
files_gosat = FILE_SEARCH(gosatdir + '*.nc4', COUNT=n_files_gosat_total)
dates_gosat = STRMID(FILE_BASENAME(files_gosat), 11,6)

;;----- Truncate the GOSAT file listing to the period of overlap with OCO-2.
match_start_date=20140906 ;;-- start date of OCO-2 science observations.
match_end_date=20200630   ;;-- termination date of ACOS GOSAT v9 XCO2 product
wa = WHERE(INRANGE(LONG('20'+dates_gosat), [match_start_date,match_end_date]), n_files_gosat)
files_gosat=files_gosat[wa]
dates_gosat=dates_gosat[wa]
PRINT,'n_files_gosat_total=',n_files_gosat_total
PRINT,'n_files_gosat (matched to OCO-2 observation period)=',n_files_gosat
;STOP

;;----- Set the matchup collocation criteria.
dtime_max = 2.0        ;;-- time difference in hours
dlat_max  = 2.0        ;;-- latitude difference in deg
dlon_max  = 3.0        ;;-- longitude difference in deg
dist_max = 300.        ;;-- max distance in km
min_soundings_oco2 = 3 ;;-- min number of matching OCO-2 soundings
nmax = 100L            ;;-- max number of oco2 soundings per GOSAT sounding
take_closest=0         ;;-- 0= take nmax random; 1= take nmax closest

;;---- Set the list of variables to read from the GOSAT Lite files (otherwise the data structures will be MUCH larger).
var_list_gosat = ['sounding_id','latitude','longitude','xco2', 'xco2_quality_flag', 'solar_zenith_angle', $
         'co2_profile_apriori','pressure_weight','xco2_averaging_kernel',$
         'Sounding/gain',  'Retrieval/surface_type']
;var_list_gosat = ['sounding_id','latitude','longitude','xco2', 'xco2_quality_flag', 'solar_zenith_angle']

;;---- Set the list of variables to read from the OCO Lite files (otherwise the data structures will be MUCH larger).
var_list_oco = ['sounding_id','latitude','longitude','xco2', 'xco2_quality_flag','solar_zenith_angle',$
         'co2_profile_apriori','pressure_weight','xco2_averaging_kernel', $
         'Sounding/operation_mode', 'Sounding/orbit', 'Retrieval/surface_type']
;var_list_oco = ['sounding_id','latitude','longitude','xco2', 'xco2_quality_flag','solar_zenith_angle']

;;-----------------------------------------------------------------
;;----- Part 1.
;;----- Read in the GOSAT data and the corresponding OCO data.
;;----- For each day of GOSAT data (a single Lite file), three
;;----- days of OCO data are read into memory (+/-1 day around the
;;----- GOSAT date. This is because the two satellites are in
;;----- different orbits and the data files are chunked by UTC day
;;----- and there is the international date line!
;;-----------------------------------------------------------------

;;----- Initialize the loaded_oco_dates string array.
loaded_oco_dates = MAKE_ARRAY(3,/STRING,VALUE='n/a')

n=0L
;;----- Loop over the number of input GOSAT Lite files.
n_files_gosat=2 ;; falsify the value for a quick test. Comment this line out to loop the full GOSAT file set.

FOR f=0, n_files_gosat-1 DO BEGIN

   ;;----- Read the current GOSAT Lite file into structure "gosat".
   gosat = H5_FASTREAD(files_gosat[f], READ=var_list_gosat)

   ;;----- Index GOSAT data for quality flag good (=0) soundings only.
   wa = WHERE(gosat.xco2_quality_flag EQ 0, nwa)
   IF (nwa EQ 0) THEN CONTINUE

   ;;----- If at least one good quality GOSAT sounding exists,
   ;;----- then subset the original gosat structure to only those with good QF.
   gosat=gosat[wa]

   ;;----- Set current date.
   current_date = '20'+dates_gosat[f]

   ;;----- Set the Julian Dates needed to poll the matching OCO data.
   ;;----- Always get one day before and after the GOSAT target date (even if we don't need them).
   jd_need = dtime_max*[-1,1.]/24. + ACOS_ID_TO_JD(MINMAX(gosat.sounding_id))
   jd_need = DIFFERENT([jd_need[0], jd_need[0]+1., jd_need[1]])

   ;;----- Convert the Julian Day values (long integers) to string values.
   dates_need = DIFFERENT(STRMID(ACOS_JD_TO_ID(jd_need),0,8))
   ndates_need = N_ELEMENTS(dates_need)
   PRINT, 'Current date of GOSAT file = ',current_date
   PRINT, 'Dates of OCO data that will be loaded = ' + STRJOIN(dates_need,', ')
   ;STOP

   ;;----- Initialize the counter for the number of OCO files loaded for the current GOSAT file.
   n_oco_files_loaded=0L

   ;;----- Loop over the needed dates (+/-1 of the current date).
   FOR j=0, ndates_need-1 DO BEGIN

      fo = (WHERE(loaded_oco_dates EQ dates_need[j]))[0]
      PRINT,'fo=',fo
      PRINT,'(fo LT 0): fo should = -1 on the first iteration of this loop...new OCO data will be loaded.'
      PRINT,'(fo GE 0): then OCO data for current date has already been loaded...'
      ;STOP

      IF (fo GE 0) THEN BEGIN

         ;;----- OCO data has already been loaded for this date.
         ;;----- Assign 'loaded_ocoX' to 'oco_current' structure.
         ;;----- 'loaded_ocoX' were previously saved later in this code block on a previous iteration.
         PRINT,'Assigning the values from previously loaded OCO data to the "oco_current" structure.'
         dummy = EXECUTE('oco_current = loaded_oco'+SC(fo))
         n_oco_current = N_ELEMENTS(oco_current)
         ;STOP

      ENDIF ELSE BEGIN
         ;;----- Load OCO data for this date.

         ;;----- Find the currently loaded dates that we do NOT need
         fo = (WHERE(~ELT(loaded_oco_dates, dates_need), nnot))[0]
         date6=STRMID(dates_need[j],2,6)
         PRINT, 'Need to load OCO date ' + date6
         oco_file=(FILE_SEARCH(oco2dir + '*_'+date6+'_*.nc4', COUNT=ofound))[ofound-1]  ; LAUREL: why do we only read in the last file if multiple files are found? Are we assuming the files are identical?

         ;;----- If an OCO files does not exists for this date, then skip...
         IF (ofound EQ 0) THEN BEGIN
            PRINT,'No OCO Lite file found for current date.'
            CONTINUE
         ENDIF

         ;;----- If an OCO file does exist, then read it into structure 'oco_current'.
         PRINT,'oco_file=',oco_file
         ;STOP
         oco_current = H5_FASTREAD(oco_file, READ=var_list_oco)

         ;;----- Index for quality flag good (=0).
         wg = WHERE(oco_current.xco2_quality_flag EQ 0, n_oco_current)
         PRINT, 'Loaded OCO date ' + date6 + ' with ' + SC(n_oco_current) + ' good soundings.'
         IF (n_oco_current EQ 0) THEN CONTINUE

         ;;----- If at least one good quality OCO sounding exists,
         ;;----- then subset the original 'oco_current' structure.
         oco_current=oco_current[wg]

         ;;----- Update the loaded_oco_dates array.
         loaded_oco_dates[fo] = dates_need[j]

         ;;----- Assign 'oco_current' to 'loaded_ocoX' structure.
         ;;----- This saves the already loaded data into memory so that on the
         ;;----- next iteration of reading OCO files, the data does not
         ;;----- have to be reloaded.
         dum=EXECUTE('loaded_oco'+SC(fo)+' = oco_current')

         ;;----- Terminate block that loads OCO data.
      ENDELSE

      ;;----- Append the currently read in 'oco_current' data to the 'oco' structure.
      IF (n_oco_files_loaded EQ 0) THEN oco=oco_current ELSE oco=[oco,oco_current]

      ;;----- Update the counter
      n_oco_files_loaded = n_oco_files_loaded + 1
      PRINT,'n_oco_files_loaded='+SC(n_oco_files_loaded)

   ;;----- Terminate the loop over the number of dates needed.
   ENDFOR

   PRINT, 'Loaded dates =' + STRJOIN(loaded_oco_dates,', ')
   ;STOP


   ;;-------------------------------------------------------
   ;;----- Part 2.
   ;;----- Match the OCO data to the GOSAT soundings.
   ;;-------------------------------------------------------
   oco_hour = oco.sounding_id / 1000000ULL ; YYYYMMDDHH
   oco_hours = DIFFERENT(oco_hour)
   oco_jd = ACOS_ID_TO_JD_SLOW(oco.sounding_id)

   w = WHERE(gosat.sounding_id GE MIN(oco.sounding_id/100ULL) AND gosat.sounding_id LE MAX(oco.sounding_id/100ULL))
   gosat = gosat[w]
   gosat_id=gosat.sounding_id
   gosat_hour = gosat_id / 10000ULL
   gosat_jd = acos_id_to_jd(gosat_id)

   nv = N_ELEMENTS(oco_hours)
   PRINT, 'Value Locating Hours in OCO2 & GOSAT'
   v = VALUE_LOCATE(oco_hour, oco_hours)
   vj = VALUE_LOCATE(oco_hours, gosat_hour)

   ; cycle through all GOSAT soundings?
   if f eq 0 then begin
        aelt = zero_structure(gosat[0])  # LAUREL: What is zero_structure doing? Setting the types within the structure?
        oelt = zero_structure(oco[0])
        outelt = {gosat:aelt, oco2:replicate(oelt,nmax), ntot:0, n:0}
        match = replicate(outelt, 250000L)
   endif

   ngosat = n_elements(gosat)
   n_oco_current=0
   for i=0L,ngosat-1 do begin
        ; match THIS sounding.

        ; which OCO2 soundings might pass the time filter?
      ;  print, 'gosat Sounding ' + sc(i+1) + ' / ' + sc(ngosat)
        thisa = gosat[i]
        ajd = gosat_jd[i]
        aid = gosat_id[i]
        j =  vj[i] ; centered hour
        va = v[(j-2)>0] & vb = v[(j+3)<(nv-1)] ; which OCO soundings to try  LAUREL: why j-2 and j+3
        thiso = oco[va:vb]

        w = where(thiso.retrieval.surface_type eq thisa.retrieval.surface_type, nw)
        if nw eq 0 then continue
        thiso = thiso[w]
        ojd = (oco_jd[va:vb])[w]

        dtime = (ojd-ajd)*24.
        dlat  = thiso.latitude-thisa.latitude
        dlon  = longitude_difference(thiso.longitude, thisa.longitude)
        ddist = co_sphdist(thiso.longitude, thiso.latitude, thisa.longitude, thisa.latitude, /approx,/deg)*111.1

        ###START HERE

        mask = abs(dtime) LE dtime_max $
           AND abs(dlat) LE dlat_max $
           AND abs(dlon) LE dlon_max $
           AND ddist LE dist_max


        w = where(mask, nw)
        if nw LT min_soundings_oco2 then continue

       ; print, sc(aid) + ': MATCH with ' + sc(nw) + ' OCO-2 soundings; nmatch='+sc(n)

        match[n].ntot = nw
        if nw GT nmax then begin
          if take_closest then begin
            ; take the 100 closest
            dist = abs(dlat[w]) + cosd(thisa.latitude) * abs(dlon[w])
            so = sort(dist)
            w = w[so[0:nmax-1]]
          endif else begin
            ; take 100 random soundings
            w = (shotgun(w))[0:nmax-1]
          endelse
          nw = nmax
        endif
        match[n].n = nw
        match[n].gosat = thisa
        match[n].oco2[0:nw-1] = thiso[w]
        n += 1
        n_oco_current +=1
   endfor
   print, current_date + ': Nmatch = ' + sc(n_oco_current)+'; Nmatch_tot = ' + sc(n)

   ;;----- End loop over the number of input GOSAT Lite files.
ENDFOR

;;----- Truncate output 'match' structure to total number of matched soundings.
match = match[0:n-1]

;;----- Save the 'match' structure to IDL save file.
checkdir, outfile
print, 'Saving ' + outfile
SAVE, FILENAME=outfile, match


END
