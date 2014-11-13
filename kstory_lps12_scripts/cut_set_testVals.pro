;;;
; NAME: cut_set_testVals
; PURPOSE:
;   Test different cut values for runlists
;
; MODIFICATION HISTORY:
;  08/20/2011: (KTS) Created
;  08/30/2011: (KTS) Add tweight cut
;;;

;...................................................................
; Set of cuts
;
function cut_set_testVals
compile_opt IDL2, HIDDEN
cuts = {default2:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        $
        high_medwt_test1p5:{high_medwt:1.5,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        low_medwt_test3p0:{high_medwt:2.0,low_medwt:3.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        low_medwt_test1p5:{high_medwt:2.0,low_medwt:1.5, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        low_medwt_test1p3:{high_medwt:2.0,low_medwt:1.3, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        $
        high_rms_test1p5:{high_medwt:2.0,low_medwt:2.0, high_rms:1.5,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        low_rms_test3p0:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:3.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        low_rms_test1p5:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:1.5,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        $
        high_wtrms_test1p5:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:1.5,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        high_wtrms_test1p3:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:1.3,low_wtrms:2.0, high_tweight:2.0, low_tweight:2.0}, $
        low_wtrms_test1p5:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:1.5, high_tweight:2.0, low_tweight:2.0}, $
        low_wtrms_test1p3:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:1.3, high_tweight:2.0, low_tweight:2.0}, $
        $
        high_tweight_test1p5:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:1.5, low_tweight:2.0}, $
        high_tweight_test1p3:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:1.3, low_tweight:2.0}, $
        high_tweight_test1p2:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:1.2, low_tweight:2.0}, $
        low_tweight_test1p5:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:1.5}, $ 
        low_tweight_test1p3:{high_medwt:2.0,low_medwt:2.0, high_rms:2.0,low_rms:2.0,high_wtrms:2.0,low_wtrms:2.0, high_tweight:2.0, low_tweight:1.3} }

return, cuts
end

