;;;
; NAME: mk_make_maps_script_0224
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
;
; MODIFICATION HISTORY:
;  02/22/2012: (KTS) Created
;  02/24/2012: (KTS) Modified to use fixed poly order to 7,
;                    Add n1 and n2 for all fields
;  02/28/2012: (KTS) Use get_lps12_map_npix.pro 
;  03/01/2012: (KTS) Change expected ptsrc_config name to _20120301.txt
;;;


;...................................................................
; test noise psds
pro run_0301
for ii=0, 19 do begin
    print, ii
endfor
end


;...................................................................
; Make the runscript for a given field
pro mk_make_maps_script_0301, field_idx, runlist=runlist
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
map_dir = '/data/kstory/projects/lps12/maps/20120224/'

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

; ptsrc_config: 50mJy cut
ptsrc_config = '/data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_50mJy_'+field_name+'_20120301.txt'

; get field properties
npix = get_lps12_map_npix(field_idx, /obs_dates)
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




;;;;;;;;;;;;;;;;;;;;;;
; Obsolete

; ;...................................................................
; ; Quickly get the poly order and number of pixels for each field
; ;   Note, currently the poly order is fixed to 7
; function get_field_properties, field_idx
; case field_idx of

;     0: begin ; ra5h30dec-55_2008
;         ret = {n1:'960', n2:'960'}
;     end
;     1: begin ; ra23h30dec-55_2008
;         ret= {n1:'960', n2:'960'}
;     end
;     2: begin ; ra23h30dec-55_2010
;         ret= {n1:'960', n2:'960'}
;     end
;     3: begin ; ra21hdec-60
;         ret= {n1:'1280', n2:'960'}
;      end
;     4: begin ; ra3h30dec-60
;         ret= {n1:'2048', n2:'1024'}
;     end
;     5: begin ; ra21hdec-50
;         ret= {n1:'1536', n2:'960'}
;     end
;     6: begin ; ra4h10dec-50
;         ret= {n1:'1536', n2:'960'}
;     end
;     7: begin ; ra0h50dec-50
;         ret= {n1:'1536', n2:'960'}
;     end
;     8: begin ; ra2h30dec-50
;         ret= {n1:'1536', n2:'960'}
;     end
;     9: begin ; ra1hdec-60
;         ret= {n1:'1280', n2:'960'}
;     end
;     10: begin ; ra5h30dec-45
;         ret= {n1:'1024', n2:'1024'}
;     end
;     11: begin ; ra6h30dec-55
;         ret= {n1:'960', n2:'960'}
;     end
;     12: begin ; ra23hdec-62.5
;         ret= {n1:'1152', n2:'640'}
;     end
;     13: begin ; ra21hdec-42.5
;         ret= {n1:'1728', n2:'640'}
;     end
;     14: begin ; ra22h30dec-55
;         ret= {n1:'960', n2:'960'}
;     end
;     15: begin ; ra23hdec-45
;         ret= {n1:'1728', n2:'960'}
;     end
;     16: begin ; ra6hdec-62.5
;         ret= {n1:'1152', n2:'640'}
;     end
;     17: begin ; ra3h30dec-42.5
;         ret= {n1:'2160', n2:'960'}
;     end
;     18: begin ; ra1hdec-42.5
;         ret= {n1:'1728', n2:'640'}
;     end
;     19: begin ; ra6h30dec-45
;         ret= {n1:'1024', n2:'1024'}
;     end
; endcase

; return, ret
; end

