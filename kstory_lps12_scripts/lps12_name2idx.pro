;;;
; NAME: lps12_name2idx
; PURPOSE:
;   Convert a string name into an idx from lps12_fieldstruct
;
; CALLING SEQUENCE: field = 'ra4h10dec-50' & idx = lps12_name2idx(field)
;
; INPUTS: field_name [string]
;
; OUTPUTS:
;   idx [int]
;
; NOTES:
;
; MODIFICATION HISTORY:
;   
;   04/24/2012: (KTS) Created
;;;


;...................................................................
; Main function
FUNCTION lps12_name2idx, field_name

case field_name of
; 2008
    'ra5h30dec-55_2008'  : idx=0
    'ra5h30dec-55'       : idx=0
    'ra23h30dec-55_2008' : idx=1
    'ra23h30dec-55'      : idx=1
    'ra23h30dec-55_2010' : idx=2
    'ra23h30dec-55'      : idx=2
; 2009    
    'ra21hdec-60'        : idx=3
    'ra3h30dec-60'       : idx=4
    'ra21hdec-50'        : idx=5
; 2010
    'ra4h10dec-50'       : idx=6
    'ra0h50dec-50'       : idx=7
    'ra2h30dec-50'       : idx=8
    'ra1hdec-60'         : idx=9
    'ra5h30dec-45'       : idx=10
    'ra6h30dec-55'       : idx=11
    'ra23hdec-62.5'      : idx=12
    'ra21hdec-42.5'      : idx=13
    'ra22h30dec-55'      : idx=14
    'ra23hdec-45'        : idx=15
    'ra6hdec-62.5'       : idx=16
    'ra3h30dec-42.5'     : idx=17
    'ra1hdec-42.5'       : idx=18
    'ra6h30dec-45'       : idx=19
    'ra5h30dec-55_2011'  : idx=20
    'ra5h30dec-55'       : idx=20
    else : begin
        print, 'LPS12_NAME2IDX, field_name ('+field_name+') not recognized.  Returning -1'
        idx=-1
    endelse
endcase

RETURN, idx
END
