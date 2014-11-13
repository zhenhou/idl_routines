;;;
; NAME: mk_make_maps_script_20130428
; PURPOSE:
;   Make the .tcsh scripts which will make maps
;
; INPUTS: 
;   runlist:                If this is specified, it must be a txt file with the names of the idfs that need maps made.  
;                           Otherwise, the runlist from lps12_fieldstruct is used
;
; NOTES:
; (called in the program, not arguments to the pro)
;   ptsrc_config files:     /data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_50mJy_
;   output map dir:         /data/kstory/projects/lps12/maps/20120224/
;   IDF dir:                defined in lps12_fieldstruct
;
; 1) Makes maps in the lps12/maps/20120224/*field/ directories
; 2) Uses preaz runlists, so we have the option of changing the azrms cut later.
;
; MODIFICATION HISTORY:
;  02/22/2012: (KTS) Created
;  02/24/2012: (KTS) Modified to use fixed poly order to 7,
;                    Add n1 and n2 for all fields
;  02/28/2012: (KTS) Use get_lps12_map_npix.pro 
;  03/01/2012: (KTS) Change expected ptsrc_config name to _20120305.txt
;  03/05/2012: (KTS) Change expected ptsrc_config name to
;                    _20120305.txt, output directory maps/20120305/
;  03/26/2012: (KTS) use get_lps12_fieldinfo instead of get_lps12_map_npix
;  04/16/2012: (KTS) Copy from mk_make_maps_script_0305.pro
;  06/23/2012: (KTS) Copy from mk_make_maps_script_0420.pro
;  04/30/2013: (KTS) Copy from mk_make_maps_script_0620.pro
;;;


;...................................................................
; test noise psds
pro run_0428
list = [1,2,7,9,12,15,18]
n = n_elements(list)


for i=0, n-1 do begin
    idx = list[i]
    mk_make_maps_script_0428, idx
endfor
end


;...................................................................
; Make the runscript for a given field
pro mk_make_maps_script_0428, idx, runlist=runlist
compile_opt IDL2, HIDDEN

;------------------------
; Get the field information
;------------------------
f = lps12_fieldstruct()
field_name     = f[idx].name
field_dir_name = f[idx].dir_name
idf_dir        = f[idx].idf_lpfds_dirs
ra0            = strtrim( string(f[idx].ra0), 2)
dec0           = strtrim( string(f[idx].dec0), 2)


;------------------------
; Get the inputs to the scanmap_sim.x call
;------------------------
map_dir = '/data/kstory/projects/lps12/maps/20130428/'

; script name
script_fname = map_dir + 'make_maps_'+field_name+'_150.txt'

; get the idf fits files

; get the runlist from lps12_fieldstruct()
if ~keyword_set(runlist) then begin
    files     = get_lps12_runlist(idx, /idfs, /preaz)
    date_list = get_lps12_runlist(idx, /obs_dates, /preaz)
    nfiles    = n_elements(files)

; If a specific runlist of idf names is passed, read that
endif else begin
    readcol,/silent,runlist,files,format='a'
    nfiles = n_elements(files)
    date_list = strarr(nfiles)
    for ii=0, nfiles-1 do begin
        date_list[ii] = extract_date_from_filename(files[ii])
    endfor
endelse

; make output file names
outfiles = strarr(nfiles) ; output maps
for ii=0, nfiles-1 do begin
    outfiles[ii] = map_dir+field_name+'/map_'+field_dir_name+'_150_'+date_list[ii]+'.fits'
endfor

; ptsrc_config: 50mJy cut
;ptsrc_config = '/data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_50mJy_'+field_name+'_20120305.txt'
ptsrc_config = '/data/kstory/projects/lps12/ptsrc_lists/20130428_corrected/ptsrc_config_50mJy_'+field_name+'_20130428.txt'

; get field properties
npix = (get_lps12_fieldinfo(idx)).npix
snpix = [ strtrim(string(npix[0]), 2), strtrim(string(npix[1]), 2) ]

; set poly order to 7
poly_order = '7'

; Sample rate (field-dependent)
; Calculate all values based on down-sampling rate and kll values:
; K11: ds=6, map_lpf=5.0 - > lpf = 5.0*(6/ds)
ds = double(4)
if ( intersect(idx, [0,1,3,4,5]) ne -1 ) then begin
    ds = double(6)
endif 
if idx eq 2 then begin
    ds = double(5)
endif

samplerate = double(100.)/ds
ssrate = strtrim(string(samplerate, format='(f16.13)'), 2)

;------------------------
; write the file
;------------------------
get_lun, lun1
print, idx, ": ", script_fname
stop
openw,lun1,script_fname

; printf,lun1, '#!/bin/tcsh'
; printf,lun1, 'setenv SPT_ANAL /home/sptdat/spt_analysis'
; printf,lun1, 'setenv LD_LIBRARY_PATH /usr/lib'
; printf,lun1, 'setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:$SPT_ANAL/lib'
; printf,lun1, 'setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:$SPT_ANAL/util'
; printf,lun1, 'setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:/home/rkeisler/libs/acml/acml4.4.0/gfortran64_int64/lib:/home/tcrawfor/lib:/home/cr/hpic-0.52.2/lib'

for ii=0, nfiles-1 do begin

    printf,lun1, '/home/kstory/sptcdevel/mapping/scanmap_sim.x ' + $
      ' -i ' + files[ii] + $
      ' -doreal ' + outfiles[ii] + $
      ' -ptsrc ' + ptsrc_config + $
      ' -poly ' + poly_order + $
      ' -n1 ' + snpix[0] + ' -n2 ' + snpix[1] + ' -ra0 ' + ra0 + ' -dec0 ' + dec0 + $
      ' -maskedellhpf 270 -o 1 -jacklr a.fits -reallr a.fits' + $
      ' -nmap 1 -wedgespatial 0' + $
      ' -reso_arcmin 1.0' + $
      ' -samplerate '+ ssrate +' -hpf 0.0 -elllpf 6600 -proj 5 -shrinkoutput'
endfor
close, lun1
free_lun,lun1

end
