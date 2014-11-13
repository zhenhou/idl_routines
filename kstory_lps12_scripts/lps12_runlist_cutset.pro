;;;
; NAME: lps12_runlist_cutset.pro
; PURPOSE:
;   Return a structure of the final runlist cut values
;
; CALLING SEQUENCE: cutSet = lps12_runlist_cutset()
;
;
; MODIFICATION HISTORY:
;   09/13/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
function lps12_runlist_cutset
compile_opt IDL2, HIDDEN

return,  {high_medwt:1.25,  $
          low_medwt:2.0,    $
          high_rms:2.0,     $
          low_rms:3.0,      $
          high_wtrms:2.0,   $
          low_wtrms:2.0,    $
          high_tweight:1.2, $
          low_tweight:1.5}
end
