;;;
; NAME: get_1hz_bolo_info
; PURPOSE:
;   Get info for cutting bad 1hz bolos, set up to run study_1hz_bolo.pro
;
; CALLING SEQUENCE: get_1hz_bolo_info, 17
;
; INPUTS: 
;   field index,        index in lps12_fieldstruct()
;
; OUTPUTS:
;   sav file called /data/kstory/projects/lps12/1hzBolo/lines_[field_name].sav
;
; NOTES: (from RK, who got this code from TC)
;   1) Analyzed by study_1hz_bolos.pro
;
; MODIFICATION HISTORY:
;   02/07/2012: (KTS) Created from /data/rkeisler/lines/gather3.pro
;;;

pro get_1hz_bolo_info, field_idx
compile_opt IDL2, HIDDEN

; Get the IDF directory
field_arr = lps12_fieldstruct()
fitsdir = field_arr[field_idx].idf_dirs[0]
field_name = field_arr[field_idx].name

; get the low ell run list
get_runlist_by_field, field_idx, list=lowlist
list = fitsdir + 'field_scan_150_'+lowlist+'.fits'
;spawn,'ls '+fitsdir+'field_scan_150_*',list;save
nlist = n_elements(list);save
date = extract_date_from_filename(list);save
mjd = date_string_to_mjd(date);save

; sort by mjd
mysort = sort(mjd)
list = list[mysort]
date = date[mysort]
mjd = mjd[mysort]

; initialize arrays
k=0
keep_going=1
while keep_going do begin
    print, "read file ", list[k]
    d=read_spt_fits(list[k],include='OBSERVATION')
    CATCH, Error_status
    if Error_status NE 0 then begin
        print, " ---- caught error reading in file, ", list[k]
        print, "stopping..."
        stop
    endif
; make sure it's a structure
    if datatype(d) ne 'STC' then begin
        k++
        continue
    endif
; check for OBSERVATION tag
    tags = get_tags(d)
    wh=where(tags eq '.OBSERVATION',nwh)
    if nwh eq 0 then begin
        k++
        continue
    endif
    keep_going=0
endwhile
boloid = d.observation.bolo_readout_idx;save
nbolos = n_elements(boloid)

nhz = 40
center = fltarr(nlist,nhz,nbolos)
width = fltarr(nlist,nhz,nbolos)
height = fltarr(nlist,nhz,nbolos)
power = fltarr(nlist,nhz,nbolos)
flag = intarr(nlist,nbolos)
bolocal = fltarr(nlist,nbolos)
bolowts = fltarr(nlist,nbolos)
ptc = fltarr(nlist)

;;;
; Loop over files
;;;

for i=0,nlist-1 do begin
    timea = systime(/sec)
    print,strtrim(i,2),'/',strtrim(nlist-1,2)
    print,'restoring ',list[i]

    ; Read file
    d = read_spt_fits(list[i],include='OBSERVATION')
    print,'...restoring took ',systime(/sec)-timea,' seconds.'

    ; make sure it's a structure
    if datatype(d) ne 'STC' then continue
    ; check for OBSERVATION tag
    tags = get_tags(d)
    wh=where(tags eq '.OBSERVATION',nwh)
    if nwh eq 0 then continue


    flag[i,*]=d.observation.bolo_flags
    bolocal[i,*] = d.observation.bolo_relcal
    bolowts[i,*] = (d.observation.bolo_relcal*d.observation.BOLO_PSD_1TO3)^(-2.)
    ptc[i] = ptc_frequency(d.observation.start_mjd)
    center[i,*,*] = (d.observation.BOLO_LINE_CENTERS)
    width[i,*,*] = (d.observation.BOLO_LINE_WIDTHS)
    this_height = (d.observation.BOLO_LINE_HEIGHTS)
    for j=0,nbolos-1 do this_height[*,j] *= ((d.observation.bolo_relcal)[j])
    height[i,*,*] = this_height
    power[i,*,*] = this_height*(d.observation.BOLO_LINE_WIDTHS)
    
    d=0
    print,'...that obs took ',systime(/sec)-timea,' seconds.'
    if (i mod 50) eq 0 or (i eq nlist-1) then $
      save,boloid,nbolos,date,mjd, $
      bolocal,bolowts,ptc,flag,list,nlist,center,width,height,power, $
      filename='/data/kstory/projects/lps12/1hzBolo/lines_'+field_name+'.sav'
endfor

end


