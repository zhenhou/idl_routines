;;;
; NAME: script_0301
; PURPOSE:
;   Investigate map overlap
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  03/01/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Investigating nmiss
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;...................................................................
; Make nmiss arrays by field
pro sss;, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
;map_dir = '/data/kstory/projects/lps12/maps/20120224/'
outdir = '/data/kstory/projects/lps12/masks/apod_0228/'

field_idx = 1 ; ra23h30dec-55_2008


;------------------------
; Get field information
;------------------------
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name
map_dir = fst.xspec_map_dir


;------------------------
; get the xspec maps
;------------------------
files = get_lps12_runlist(field_idx, /xspec_maps)
nfiles = n_elements(files)

;------------------------
; Get nmiss for this observation
;------------------------
info = expand_fits_struct(read_spt_fits(files[0])) ; get map shape
nmiss = 0*info.map.map

window, 2
;for ii=0, nfiles-1 do begin
for ii=0, 9 do begin
    tmp = 0*info.map.map
    print, files[ii]
    d = expand_fits_struct(read_spt_fits(files[ii]))
    weight = d.weight.map
    weight /=max(weight)

    ; find pixels that have zero obs
    wh_0 = where(weight eq 0, nwh, complement=wh_1)

    nmiss[wh_0] += 1

    ; plot it
    tmp[wh_0] = 1
    tv_spt_map, tmp, /norms, scale=0.75, /forcesize, winnum=2
    com = ''
    read, com, prompt = 'mapnum = ' + strtrim(string(ii),2)
endfor

stop    
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Making ra5h30dec-55_2011 maps
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;...................................................................
; Make a runlist of the ra5h30dec-55_2011 with all idf's in the
; directory
pro make_ra5h30_2011_runlist_allidfs
compile_opt IDL2, HIDDEN
idf_dir = '/home/kstory/lps12/lowellfits/ra5h30dec-55_2011/'
spawn, 'ls ' + idf_dir + 'field_scan_150_*.fits', list
nfiles = n_elements(list)

out_fname = '/data/kstory/projects/lps12/runlists/runlist_allidf_ra5h30dec-55_2011.txt'

; Write to a runlist file

get_lun, lun1
openw, lun1, out_fname
for ii=0, nfiles-1 do begin
    printf, lun1, list[ii]
endfor
close, lun1
free_lun,lun1

end

;...................................................................
; Make a runlist of the ra5h30dec-55_2011 with all idf's in the
; directory
pro run_mk_make_maps_script_0301_ra5h30_2011
MK_MAKE_MAPS_SCRIPT_0301, 20, runlist='/data/kstory/projects/lps12/runlists/runlist_allidf_ra5h30dec-55_2011.txt'
end

;...................................................................
; Investigate nmiss in observations
pro script_0228
compile_opt IDL2, HIDDEN

; ; Get field information
; field_arr = lps12_fieldstruct()
; fst = field_arr[field_idx]
; field_name = fst.name
; field_dir_name = fst.dir_name

; r5h30dec-55 [0]
;myfile = '/data/kstory/projects/lps12/maps/20120224/coadds/coadd_ra5h30dec-55_2008_50mJy.fits'a
;myfile = '/data/kstory/projects/lps12/maps/20120224/ra5h30dec-55_2008/map_ra5h30dec-55_150_20080529_091936.fits'
;myfile = '/data/kstory/projects/lps12/maps/20120224/ra5h30dec-55_2008/map_ra5h30dec-55_150_20080528_042504.fits'

; ra23h30dec-55_2008 [1]
;myfile = '/data/kstory/projects/lps12/maps/20120224/coadds/coadd_ra23h30dec-55_2008_50mJy.fits'
myfile = '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008/map_ra23h30dec-55_150_20080527_070222.fits'
myfile = '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008/map_ra23h30dec-55_150_20080527_073643.fits'
;myfile = '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008/map_ra23h30dec-55_150_20080527_153901.fits'
;myfile = '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008/map_ra23h30dec-55_150_20080610_014537.fits'

d0full = expand_fits_struct(read_spt_fits(myfile))
d0full_ss = d0full.weight.map
d0full_ss /=max(d0full_ss)

; find pixels that have zero obs
wh_0 = where(d0full_ss eq 0, nwh, complement=wh_1)

; set to either 0 or 1
tmp = d0full_ss
tmp[wh_1] = 1

tmp = -1*tmp + 1 ; invert
tv_spt_map, tmp, /norms, scale=1, /forcesize


stop
end




