;;;
; NAME: script_0424
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) run try_end2end.pro, first attempt
;
; MODIFICATION HISTORY:
;  04/24/2012: (KTS) Created
;;;

;...................................................................
; run_e2e
PRO run_e2e

f = lps12_fieldstruct()
idx = 3
fst = f[idx]

try_end2end, idx, run='0424', /resume

END


