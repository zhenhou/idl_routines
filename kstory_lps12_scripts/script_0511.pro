;;;
; NAME: script_0510
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  05/09/2012: (KTS) Created
;;;

PRO coadd_jack_sims_12
for idx=0, 19 do begin
    print, 'coadd jack sims 12, field ', idx
    coadd_lps12_jack_sims, idx, '12'
endfor
END

