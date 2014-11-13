;;;
; NAME: mk_make_maps_script_6p4mJy.pro
; PURPOSE:
;   Make the .tcsh scripts which will make maps with ptsrc cut 6.4mJy
;
; INPUTS: 
;   runlist:                If this is specified, it must be a txt file with the names of the idfs that need maps made.  
;                           Otherwise, the runlist from lps12_fieldstruct is used
;
; NOTES:
; (called in the program, not arguments to the pro)
;   ptsrc_config files:     /data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_6p4mJy_
;   output map dir:         /data/kstory/projects/lps12/maps/20120402_6p4mJy/
;   IDF dir:                defined in lps12_fieldstruct
;
; 1) Makes maps in the lps12/maps/20120224/*field/ directories
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
;  04/02/2012: (KTS) Re-name to mk_make_maps_script_6p4mJy.pro
;;;


;...................................................................
; test noise psds
pro run_6p4mJy
for ii=0, 19 do begin
    mk_make_maps_script_6p4mJy, ii
endfor

; treat ra5h30dec-55_2011 differently
MK_MAKE_MAPS_SCRIPT_6p4mJy, 20, runlist='/data/kstory/projects/lps12/runlists/runlist_allidf_ra5h30dec-55_2011.txt'
end


;...................................................................
; Make the runscript for a given field
pro mk_make_maps_script_6p4mJy, field_idx, runlist=runlist
compile_opt IDL2, HIDDEN

;------------------------
; Get the field information
;------------------------
field_arr = lps12_fieldstruct()
field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name
idf_dir = field_arr[field_idx].idf_lpfds_dirs
ra0 = strtrim( string(field_arr[field_idx].ra0), 2)
dec0 = strtrim( string(field_arr[field_idx].dec0), 2)


;------------------------
; Get the inputs to the scanmap_sim.x call
;------------------------
map_dir = '/data/kstory/projects/lps12/maps/20120402_6p4mJy/'

; script name
script_fname = map_dir + 'make_maps_'+field_name+'_150.tcsh'

; get the idf fits files

; get the runlist from lps12_fieldstruct()
if ~keyword_set(runlist) then begin
    files     = get_lps12_runlist(field_idx, /idfs)
    date_list = get_lps12_runlist(field_idx, /obs_dates)
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

; ptsrc_config: 6p4mJy cut
ptsrc_config = '/data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_6p4mJy_'+field_name+'_20120305.txt'

; get field properties
npix = (get_lps12_fieldinfo(field_idx)).npix
snpix = [ strtrim(string(npix[0]), 2), strtrim(string(npix[1]), 2) ]

; set poly order to 7
poly_order = '7'

;------------------------
; write the file
;------------------------
get_lun, lun1
print, field_idx, ": ", script_fname
openw,lun1,script_fname

printf,lun1, '#!/bin/tcsh'
printf,lun1, 'setenv SPT_ANAL /home/sptdat/spt_analysis'
printf,lun1, 'setenv LD_LIBRARY_PATH /usr/lib'
printf,lun1, 'setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:$SPT_ANAL/lib'
printf,lun1, 'setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:$SPT_ANAL/util'
printf,lun1, 'setenv LD_LIBRARY_PATH {$LD_LIBRARY_PATH}:/home/rkeisler/libs/acml/acml4.4.0/gfortran64_int64/lib:/home/tcrawfor/lib:/home/cr/hpic-0.52.2/lib'

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
      ' -samplerate 16.6666666666667 -hpf 0.0 -lpf 5.0 -proj 5 -shrinkoutput'
endfor
close, lun1
free_lun,lun1

end
