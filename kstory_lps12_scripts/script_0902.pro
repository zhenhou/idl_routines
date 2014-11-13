;;;
; NAME: script_0902
; PURPOSE:
;   make runlists for high_tweight = 1.25
;
; NOTES:
;
; MODIFICATION HISTORY:
;  09/01/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro script_0902
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'

exclude_list = [-1]
field_idx_list = setdifference(lindgen(n_elements(all_fields)), exclude_list)
;field_idx_list = [10]

for ii=0, n_elements(field_idx_list)-1 do begin
    idx = field_idx_list[ii]
    field_name = all_fields[idx]

    nobs = n_elements(read_runlist('/data/sptdat/autotools/runlists/scans_' + field_name + '.txt'))
    print, nobs
    
endfor
end




; restore, 'lps12_fieldnames.sav'
; ;exclude_list = [1,2,5,11]
; ;exclude_list = [0,1,2,5,11]
; exclude_list = [-1]
; ;field_idx_list = setdifference(lindgen(n_elements(all_fields)), exclude_list)
; field_idx_list = [10]

; for ii=0, n_elements(field_idx_list)-1 do begin
;     idx = field_idx_list[ii]
;     field_name = all_fields[idx]

;     base_dir = '/data/kstory/projects/lps12/raw_maps/'
;     pdf_name = base_dir + 'all_raw_maps_' + field_name + '.pdf'
    
;     command = 'convert ' + base_dir + field_name + '/*.png ' + pdf_name

;     print, '************* Execute: $ ' + command
;     spawn, command
    
; endfor
; stop
; end


; pro script_0902
; compile_opt IDL2, HIDDEN

; restore, 'lps12_fieldnames.sav'
; ;exclude_list = [1,2,5,11]
; ;exclude_list = [0,1,2,5,11]
; exclude_list = [-1]
; field_idx_list = setdifference(lindgen(n_elements(all_fields)), exclude_list)
; ;field_idx_list = [15]

; for ii=0, n_elements(field_idx_list)-1 do begin
;     idx = field_idx_list[ii]
;     field_name = all_fields[idx]
;     print, "************* make_final_runlists, "+field_name
;     make_final_runlists, field_name
; endfor
; stop
; end

