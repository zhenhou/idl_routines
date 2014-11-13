;;;
; NAME: make_lt_list.pro
; PURPOSE:
;   Make a list of lead and trail fields from a runlist
; 
; OUTPUTS:
;   Saves lead_list, trail_list in
;   '/data/kstory/projects/lps12/runlists/lt/list_'+field_name+'_150.sav'
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/29/2012: (KTS) Created
;  03/10/2012: (KTS) Add skip_missing, add explicit lead or trail
;                    check, by summing columns in the weight map
;  03/26/2012: (KTS) use get_lps12_fieldinfo instead of get_lps12_map_npix
;;;


;...................................................................
; Make nmiss arrays by field
pro make_lt_list, field_idx, band=band, start_index=start_index, stopit=stopit, $
                  nosave=nosave, skip_missing=skip_missing
compile_opt IDL2, HIDDEN

if ~keyword_set(band) then band = 150

;------------------------
; Get field information
;------------------------
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

; Check that this is a lt field
if ~fst.lead_trail then begin
    print, "field "+ field_name + " is not LT.  Returning."
    return
endif

;------------------------
; get the idf fits files
;------------------------
date_list = get_lps12_runlist(field_idx, /obs_dates)
map_list  = get_lps12_runlist(field_idx, /obs_maps)
ndates    = n_elements(map_list)

; If desired, skip observations that are in the runlist but do not have maps.
if keyword_set(skip_missing) then begin
    ndates = n_elements(map_list)
    date_list = strarr(ndates)
    for ii=0, ndates-1 do begin
        date_list[ii] = extract_date_from_filename(map_list[ii])
    endfor
endif

;sort by date
mjd      = date_string_to_mjd(date_list)
mysort   = sort(mjd)
mjd      = mjd[mysort]

map_list = map_list[mysort]

; Output file name
;pre = '_tmp'
pre = ''
listname = '/data/kstory/projects/lps12/runlists/lt/list_'+field_name+'_150'+pre+'.sav'
if keyword_set(skip_missing) then listname = '/data/kstory/projects/lps12/runlists/lt/list_skip_missing_'+field_name+'_150'+pre+'.sav'
;prelistname = 'pre'+listname+'.sav'

if n_elements(start_index) eq 0 then start_index=0
    
;------------------------
npix = (get_lps12_fieldinfo(field_idx)).npix
l_col = floor(npix[0]/4)
t_col = floor(3*npix[0]/4)


;------------------------
; lead-trail pairs will be separated by a fixed amount of time.  Find
; the pairs based on this separation.
case field_name of 
    'ra23h30dec-55_2008': delta_start_lt = 0.023854
    'ra21hdec-60': delta_start_lt = 0.043634
    'ra3h30dec-60': delta_start_lt = 0.0644
    'ra21hdec-50': delta_start_lt = 0.0431597
endcase

;------------------------
; find lead/trail pairs
;------------------------
lead_list = ['']
trail_list = ['']

for j=0,ndates-2 do begin

    print, map_list[j]
    ; make some noise
    if (j mod 50 eq 2) then print, 'processing file '+strtrim(string(j),2)+'/'+strtrim(string(ndates),2)

    d = expand_fits_struct(read_spt_fits(map_list[j]))
    is_lead = total(d.weight.map[l_col, *]) gt total(d.weight.map[t_col, *])
    
    ; check if it is lead
    if (is_lead) then begin

        ; check if it has another observation which lags the appropriate amount of time
        if abs(1.-(mjd[j+1]-mjd[j])/delta_start_lt) lt 0.05 then begin

            ; check that it is a trail field
            dt = expand_fits_struct(read_spt_fits(map_list[j+1]))
            is_trail = total(dt.weight.map[l_col, *]) lt total(dt.weight.map[t_col, *])
            if is_trail then begin
                lead_list = [lead_list,map_list[j]]
                trail_list = [trail_list,map_list[j+1]]
            endif
        endif
    endif
endfor
lead_list = lead_list[1:*]
trail_list = trail_list[1:*]
npair = n_elements(lead_list)


if ~keyword_set(nosave) then save,lead_list,trail_list,filename=listname

if keyword_set(stopit) then stop
end







;...................................................................
; Make nmiss arrays by field
pro make_lt_list_old, field_idx, band, start_index=start_index, stopit=stopit, $
                  nosave=nosave, skip_missing=skip_missing
compile_opt IDL2, HIDDEN

;------------------------
; Get field information
;------------------------
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

; Check that this is a lt field
if ~fst.lead_trail then begin
    print, "field "+ field_name + " is not LT.  Returning."
    return
endif

;------------------------
; get the idf fits files
;------------------------
date_list = get_lps12_runlist(field_idx, /obs_dates)
map_list  = get_lps12_runlist(field_idx, /obs_maps)
ndates    = n_elements(map_list)

; If desired, skip observations that are in the runlist but do not have maps.
if keyword_set(skip_missing) then begin
    ndates = n_elements(map_list)
    date_list = strarr(ndates)
    for ii=0, ndates-1 do begin
        date_list[ii] = extract_date_from_filename(map_list[ii])
    endfor
endif

;sort by date
mjd      = date_string_to_mjd(date_list)
mysort   = sort(mjd)
mjd      = mjd[mysort]
map_list = map_list[mysort]

; Output file name
listname = '/data/kstory/projects/lps12/runlists/lt/list_'+field_name+'_150.sav'
if keyword_set(skip_missing) then listname = '/data/kstory/projects/lps12/runlists/lt/list_skip_missing_'+field_name+'_150.sav'
;prelistname = 'pre'+listname+'.sav'

if n_elements(start_index) eq 0 then start_index=0
    
;------------------------
; lead-trail pairs will be separated by a fixed amount of time.  Find
; the pairs based on this separation.
case field_name of 
    'ra23h30dec-55_2008': delta_start_lt = 0.023854
    'ra21hdec-60': delta_start_lt = 0.043634
    'ra3h30dec-60': delta_start_lt = 0.0644
    'ra21hdec-50': delta_start_lt = 0.0431597
endcase

;------------------------
; find lead/trail pairs
;------------------------
lead_list = ['']
trail_list = ['']
for j=0,ndates-2 do begin
    if abs(1.-(mjd[j+1]-mjd[j])/delta_start_lt) lt 0.05 then begin
        lead_list = [lead_list,map_list[j]]
        trail_list = [trail_list,map_list[j+1]]
    endif
endfor
lead_list = lead_list[1:*]
trail_list = trail_list[1:*]
npair = n_elements(lead_list)


if ~keyword_set(nosave) then save,lead_list,trail_list,filename=listname

if keyword_set(stopit) then stop
end


;;;;;;;;;;;;;;;;;;;;
;; obsolete
;;;;;;;;;;;;;;;;;;;;
; Grab list of bad trail files from ra21hdec-60
; There are observations of this one field that get labled as 'trail,'
; but are actually lead observations.  These were found by hand.
; mjd_tmp = mjd
; if (field_name eq 'ra21hdec-60') then begin
;     exclude_list = '/home/kstory/lps12/runlists/lt/exclude_ra21hdec-60.txt'
;     readcol,/silent,exclude_list,c1,format='a'

;     c1_mjd = date_string_to_mjd(c1)

;     nexclude = n_elements(c1_mjd)
;     for jj=0, nexclude-1 do begin
;         wh = where(mjd eq c1_mjd[jj], nwh, complement=whgood)
;         mjd = mjd[whgood]
;     endfor

;     ; re-adjust the mjd arrays
;     mysort   = sort(mjd)
;     mjd      = mjd[mysort]
;     ndates = n_elements(mjd)

; endif    
    
