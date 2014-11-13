;;;
; NAME: script_0830
; PURPOSE:
;   Make different runlists from total_weight cuts
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
pro script_0830, field_name
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'

field_idx = where(all_fields eq field_name, nwh)
if(nwh le 0) then begin
    print, "field_name ["+field_name+"] does not exist!  exit..."
    return
endif

restore, '/data/kstory/projects/lps12/runlists/cut_map_total/runlist_cut_map_total_'+field_name+'.sav'

; Apply cuts
cutSet = cut_set()
cuts = cutSet.default2 ; Grab default cut set

;; Default
tweight_high = 2. & tweight_low = 2.
wh = where( medwts lt median(medwts)*cuts.high_medwt and $
            medwts gt median(medwts)/cuts.low_medwt or $
            rms   lt median(rms)*cuts.high_rms and $
            rms   gt median(rms)/cuts.low_rms and $
            wtrms lt median(wtrms)*cuts.high_wtrms and $
            wtrms gt median(wtrms)/cuts.low_wtrms and $
            total_map lt median(total_map)*tweight_high and $
            total_map gt median(total_map)/tweight_low,  nwh)




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
