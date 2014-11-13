;;;
; NAME: mk_make_maps_run2_ra3h30dec-42.5.pro
; PURPOSE:
;   Make run2 maps for ra3h30dec42.5, for making the runlist
;
; INPUTS: 
;
; NOTES:
; (called in the program, not arguments to the pro)
;   ptsrc_config files:     /data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_50mJy_
;   output map dir:         /data/kstory/projects/lps12/maps/20120224/
;   IDF dir:                defined in lps12_fieldstruct
;
; 1) Makes maps in the lps12/maps/20120224/*field/ directories
;
; MODIFICATION HISTORY:
;  03/20/2012: (KTS) Created from mk_make_maps_script_0305.pro
;;;


;...................................................................
; Make the runscript for a given field
pro mk_make_maps_run2
compile_opt IDL2, HIDDEN

;------------------------
; Get the field information
;------------------------
field_idx = 17
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name
full_idf_dir = fst.idf_dirs

;------------------------
; Get the inputs to the scanmap_sim.x call
;------------------------

; the base directory where maps will be written out
map_dir = '/data/kstory/projects/lps12/maps/20120320_run2/'

; script name
script_fname = map_dir + 'make_maps_'+field_name+'_150.tcsh'

; get the idf fits files

; get the runlist from full idf dir
spawn, 'ls ' + fst.idf_dirs + 'field_scan_150_*.fits', files
date_list =  extract_date_from_filename(files)
nfiles    = n_elements(files)


; make output file names
outfiles = strarr(nfiles) ; output maps
for ii=0, nfiles-1 do begin
    outfiles[ii] = map_dir+field_name+'/map_'+field_dir_name+'_150_'+date_list[ii]+'.fits'
endfor

; ptsrc_config: 50mJy cut
ptsrc_config = '/data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_50mJy_'+field_name+'_20120305.txt'

; get field properties
npix = get_lps12_map_npix(field_idx)
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

;/home/cr/sptcdevel/mapping/scanmap_sim.x -i /data/sptdat/run1_2011/ra3h30dec-42.5/fits/field_scan_150_20100912_132651.fits -o 1 -doreal /data17/sptdat/run2_2011/ra3h30dec-42.5/maps/map_ra3h30dec-42.5_150_20100912_132651.fits -nmap 1 -wedgespatial 0 -poly 9 -reso_arcmin 0.25 -n1 9000 -n2 1800 -proj 0 -maskedellhpf 400 -lpf 30.0 -hpf 0.0 -ra0 52.5000 -dec0 -42.5000 -ptsrc /home/sptdat/spt_analysis/config/ptsrc_config_ra3h30dec-42.520110915_234944.txt -shrinkoutput

    printf,lun1, '/home/kstory/sptcdevel/mapping/scanmap_sim.x ' + $
      ' -i ' + files[ii] + $
      ' -o 1 -doreal ' + outfiles[ii] + $
      ' -nmap 1 -wedgespatial 0 -poly 9 -reso_arcmin 0.25 -n1 9000 -n2 1800 -proj 0 -maskedellhpf 400 -lpf 30.0 -hpf 0.0 -ra0 52.5000 -dec0 -42.5000 -ptsrc /home/sptdat/spt_analysis/config/ptsrc_config_ra3h30dec-42.520110915_234944.txt -shrinkoutput'
    
endfor
close, lun1
free_lun,lun1

end
