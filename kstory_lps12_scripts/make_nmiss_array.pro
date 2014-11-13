;;;
; NAME: make_nmiss_array.pro
; PURPOSE:
;   Make nmiss arrays which are used to make the apodization mask
;
; INPUTS:
;   ptsrc_config files:     /data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_50mJy_
;
; OUTPUTS:
;   saves arrays in /data/kstory/projects/lps12/masks/nmiss_sav/
;
; NOTES:
; 1) Also have procedure make_nmiss_all, which is useful for checking
;   individual nmiss arrays
;
; MODIFICATION HISTORY:
;  03/01/2012: (KTS) Created
;  03/02/2012: (KTS) add make_nmiss_all
;  03/05/2012: (KTS) add special case for ra5h30dec-55
;;;


;...................................................................
; Calculate nmiss arrays for one field
pro make_nmiss_array, field_idx, stopit=stopit
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
outdir = '/data/kstory/projects/lps12/masks/nmiss_sav/'

;------------------------
; Get field information
;------------------------
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

;------------------------
; get the map files
;------------------------
maps = get_lps12_runlist(field_idx, /xspec_maps)

; Special case for ra5h30dec-55: Take the conservative map coverage
; from before 4/29/2008
if (field_idx eq 0) then begin
    print, "Special case for field", field_name

    dates = get_lps12_runlist(field_idx, /xspec_dates)
    ndates = n_elements(dates)
    wh_good_date = intarr(ndates)

    mjd_cutoff = date_toolkit('20080429_000000', 'mjd')
    for ii=0, ndates -1 do begin
        mjd_tmp = date_toolkit(dates[ii], 'mjd')
        if mjd_tmp le mjd_cutoff then wh_good_date[ii] = 1
    endfor
    
    maps = maps[ where(wh_good_date eq 1) ]
endif

nmaps = n_elements(maps)

;------------------------
; Get nmiss for this observation
;------------------------
info = expand_fits_struct(read_spt_fits(maps[0])) ; get map shape
nmiss = 0*info.map.map

print, "MAKE_NMISS_ARRAY, nmaps = ", nmaps
for ii=0, nmaps-1 do begin
;for ii=0, 9 do begin

    ; make some noise
    if (ii mod 50 eq 2) then print, "processing " + strtrim(string(ii),2) + "/" + strtrim(string(nmaps),2)

    d = expand_fits_struct(read_spt_fits(maps[ii]))
    weight = d.weight.map
    weight /=max(weight)

    ; find pixels that have zero obs
    wh_0 = where(weight eq 0, nwh, complement=wh_1)

    nmiss[wh_0] += 1
endfor

; Save nmiss
outname = outdir+'nmiss_'+field_name+'.sav'
print, "Save file here: ", outname
save, nmiss, filename=outname

if keyword_set(stopit) then stop
end


;...................................................................
; Make a test nmiss array for each obs of ra23h30dec-55_2008
pro make_nmiss_all, field_idx, plot=plot, saveit=saveit
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
outdir = '/data/kstory/projects/lps12/masks/nmiss_sav/'

;------------------------
; Get field information
;------------------------
;field_idx = 1
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

;------------------------
; get the map files
;------------------------
maps = get_lps12_runlist(field_idx, /xspec_maps)
nmaps = n_elements(maps)

;------------------------
; Get nmiss for this observation
;------------------------
info = expand_fits_struct(read_spt_fits(maps[0])) ; get map shape
sss = size(info.map.map)
nmiss_all = fltarr(nmaps, sss[1], sss[2])

print, "MAKE_NMISS_ARRAY, nmaps = ", nmaps
for ii=0, nmaps-1 do begin
;for ii=150, nmaps-1 do begin

    ; make some noise
    if (ii mod 50 eq 2) then print, "processing " + strtrim(string(ii),2) + "/" + strtrim(string(nmaps),2)

    d = expand_fits_struct(read_spt_fits(maps[ii]))
    weight = d.weight.map
    weight /=max(weight)
    nmiss     = 0*weight

    ; find pixels that have zero obs
    wh_0 = where(weight eq 0, nwh, complement=wh_1)

    nmiss[wh_0] = 1
    nmiss_all[ii,*,*] = nmiss

    if keyword_set(plot) then begin
        ;tv_spt_map, nmiss[*,200:*], /norms, scale=1, /forcesize, winnum=1
        tv_spt_map, nmiss[*,200:*], /norms, winnum=1
        com = ''
        read, com, prompt='mapnum = ' + strtrim(string(ii))
    endif

endfor

; Save nmiss
if keyword_set(saveit) then begin 
    outname = outdir+'nmiss_all_'+field_name+'.sav'
    print, "Save file here: ", outname
    save, nmiss_all, filename=outname
endif
end


