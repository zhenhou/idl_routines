;;;
; NAME: script_0314
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Much of the work from today ended up in script_0313
; 2) 
;
; MODIFICATION HISTORY:
;  03/14/2012: (KTS) Created
;;;


;...................................................................
; Run 12 jack
pro run_jack_12_0to10
compile_opt IDL2, HIDDEN
for ii=0, 10 do begin
    lps12_jack, ii, '12'
end
end

; Run fields 11 to 19
pro run_jack_12_11to19
compile_opt IDL2, HIDDEN
for ii=11, 19 do begin
    lps12_jack, ii, '12'
end
end

