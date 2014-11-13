;;;
; NAME: cut_set_names
; PURPOSE:
;   Give the set of cut names for wtrms cuts
;
; MODIFICATION HISTORY:
;  08/20/2011: (KTS) Created
;  08/30/2011: (KTS) Add tweight cut
;;;

;...................................................................
; Set of cut names
;
function cut_set_names
cut_set_names = ['default2', $
'high_medwt_test1p5', $
'low_medwt_test3p0', $
'low_medwt_test1p5', $
'low_medwt_test1p3', $
$
'high_rms_test1p5', $
'low_rms_test3p0', $
'low_rms_test1p5', $
$
'high_wtrms_test1p5', $
'high_wtrms_test1p3', $
'low_wtrms_test1p5', $
'low_wtrms_test1p3', $
$
'high_tweight_test1p5', $ 
'high_tweight_test1p3', $
'high_tweight_test1p2', $
'low_tweight_test1p5', $ 
'low_tweight_test1p3']
return, cut_set_names
end

