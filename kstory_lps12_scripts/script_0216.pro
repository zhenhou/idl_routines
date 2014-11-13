;;;
; NAME: script_0216
; PURPOSE:
;   mess around with point-source lists
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/16/2012: (KTS) Created
;;;

;;;
pro check_ps_lists, type=type

if ~keyword_set(type) then type = 1
case type of
    1: type_str = '0p5_loosepm.txt'
    2: type_str = 'brady_astrom_corrected.txt'
    3: type_str = 'brady.txt'
endcase
    

field_arr = lps12_fieldstruct()
dir_tc = '/data/sptdat/point_source_lists/current/'


for ii=0, 20-1 do begin
    field_dir_name = field_arr[ii].dir_name
    print, "field = ", field_dir_name, ":   type_str = ", type_str
    spawn, 'ls ' + dir_tc + field_dir_name + '_list_' + type_str , list
    help, list

endfor
end

