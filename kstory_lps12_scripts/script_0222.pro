;;;
; NAME: script_0222
; PURPOSE:
;   Make some maps!
;
; INPUTS:
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/22/2012: (KTS) Created
;;;


;...................................................................
; Quickly get the poly order and number of pixels for each field
function get_field_properties, field_idx
case field_idx of

    ; ra5h30dec-55_2008
    0: begin
        ret = {poly:5, n1:960, n2:960}
    end
    1: begin
        ret= {poly:5, n1:0, n2:0}
    end
    2: begin
        ret= {poly:5, n1:0, n2:0}
    end
    3: begin
        ret= {poly:5, n1:0, n2:0}
    end
    4: begin
        ret= {poly:6, n1:0, n2:0}
    end
    5: begin
        ret= {poly:5, n1:0, n2:0}
    end
    6: begin
        ret= {poly:8, n1:0, n2:0}
    end
    7: begin
        ret= {poly:8, n1:0, n2:0}
    end
    8: begin
        ret= {poly:8, n1:0, n2:0}
    end
    9: begin
        ret= {poly:8, n1:0, n2:0}
    end
    10: begin
        ret= {poly:6, n1:0, n2:0}
    end
    11: begin
        ret= {poly:5, n1:0, n2:0}
    end
    12: begin
        ret= {poly:7, n1:0, n2:0}
    end
    13: begin
        ret= {poly:12, n1:0, n2:0}
    end
    14: begin
        ret= {poly:5, n1:0, n2:0}
    end
    15: begin
        ret= {poly:11, n1:0, n2:0}
    end
    16: begin
        ret= {poly:7, n1:0, n2:0}
    end
    17: begin
        ret= {poly:18, n1:0, n2:0}
    end
    18: begin
        ret= {poly:12, n1:0, n2:0}
    end
    19: begin
        ret= {poly:6, n1:0, n2:0}
    end
endcase

return, ret
end

;...................................................................
; test noise psds
pro run_0222
for ii=0, 19 do begin
    print, ii
endfor
end


;...................................................................
; Make the runscript for a given field
pro mkrun_script, field_idx

;------------------------
; Get the field information
;------------------------
field_arr = lps12_fieldstruct()
field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name
idf_dir = field_arr[field_idx].idf_lpfds_dirs
ra0 = field_arr[field_idx].ra0
dec0 = field_arr[field_idx].dec0


;------------------------
; Get the inputs to the scanmap_sim.x call
;------------------------
map_dir = '/data/kstory/projects/lps12/maps/20120222/'

; script name
script_fname = map_dir + 'make_maps_'+field_name+'_150.tcsh'

; get the idf fits files
spawn, 'ls ' + idf_dir + 'field_scan_150_*.fits', files
nfiles = n_elements(files)

; Make a list of the final map filenames
outfiles = strarr(nfiles)
for ii=0, n_elements(files)-1 do begin
    ss  = strsplit(files[0], '/', /extract)
    ss1 = ss[n_elements(ss)-1]
    ss2 = (strsplit(ss1, '.', /extract))[0]
    ss3 = strsplit(ss2, '_', /extract)
    ss_date = ss3[3]+'_'+ss3[4]

    outfiles[ii] = map_dir+field_name+'/map_'+field_dir_name+'_150_'+ss_date+'.fits'
endfor

; ptsrc_config: 50mJy cut
ptsrc_config = '/data/kstory/projects/lps12/ptsrc_lists/ptsrc_config_50mJy_'+field_name+'_fixed_20120220.txt'

; get field properties
field_prop = get_field_properties(field_idx)


stop
;------------------------
; write the file
;------------------------
if 0 then begin; tmp
get_lun, lun1
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
      ' -poly ' + field_prop.poly + $
      ' -n1 ' + field_prop.n1 + ' -n2 ' + field_prop.n2 + ' -ra0 ' + ra0 + ' -dec0 ' + dec0 + $
      '-maskedellhpf 270 + -o 1 -jacklr a.fits -reallr a.fits ' + $
      '-nmap 1 -wedgespatial 0 ' + $
      '-reso_arcmin 1.0 ' + $
      '-samplerate 16.6666666666667 -hpf 0.0 -lpf 5.0 -proj 5 -shrinkoutput'
endfor
close, lun1

endif ; tmp

end




