;;;
; NAME: script_0829
; PURPOSE:
;   Make cuts on total(data.weight.map)
;
; NOTES:
;   Uses output from script_0829
;
; MODIFICATION HISTORY:
;  08/29/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro script_0829, field_name
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'

field_idx = where(all_fields eq field_name, nwh)
if(nwh le 0) then begin
    print, "field_name ["+field_name+"] does not exist!  exit..."
    return
endif

restore, '/data/kstory/projects/lps12/runlists/cut_map_total/runlist_cut_map_total_'+field_name+'.sav'

; Make histogram
mymin = min(total_map)
mymax = max(total_map)
nbins = 100
hh = histogram(total_map, nbins=nbins)
xvec = mymin + findgen(100)*(mymax-mymin)/double(nbins)
plot, xvec, hh, title='total(map) for field '+field_name

stop

end

;...................................................................
; Tmp script to run script_0829
;
pro my_main
field_idx_list = [1,2,5,11]
restore, 'lps12_fieldnames.sav'

for ii=0, 3 do begin
    field_idx = field_idx_list[ii]
    field_name = all_fields[field_idx]
    script_0829, field_name
endfor

end
