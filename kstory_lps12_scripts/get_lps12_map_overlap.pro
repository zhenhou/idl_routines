;;;
; NAME: get_lps12_map_overlap
; PURPOSE:
;   Return 4-element array with overlap in deg
;     [+x, -x, +y, -y]
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  03/06/2012: (KTS) created
;;;


function get_lps12_map_overlap, idx
compile_opt IDL2, HIDDEN

ret = -1
case idx of
    0:  ret = [0.65, 1.10, 0.48, 0.48] ; ra5h30dec-55_2008
    1:  ret = [0.80, 0.75, 0.48, 0.48] ; ra23h30dec-55_2008 ;
    2:  ret = [0.80, 0.75, 0.48, 0.48] ; ra23h30dec-55_2010
    3:  ret = [0.80,   -1,   -1, 0.48] ; ra21hdec-60
    4:  ret = [0.80, 0.85,   -1, 0.48] ; ra3h30dec-60
    5:  ret = [0.75,   -1, 0.48, 0.48] ; ra21hdec-50
    6:  ret = [0.65, 0.7,  0.48, 0.48] ; ra4h10dec-50
    7:  ret = [0.7,  0.75, 0.48, 0.48] ; ra0h50dec-50
    8:  ret = [ 0.7, 0.7,  0.48, 0.48] ; ra2h30dec-50
    9:  ret = [0.85, 0.80,   -1, 0.48] ; ra1hdec-60
    10: ret = [0.65, 0.65, 0.48,   -1] ; ra5h30dec-45
    11: ret = [  -1, 0.65, 0.48, 0.48] ; ra6h30dec-55
    12: ret = [1.10, 1.10,   -1, 0.48] ; ra23hdec-62.5
    13: ret = [0.65,   -1, 0.48,   -1] ; ra21hdec-42.5
    14: ret = [0.75, 0.75, 0.48, 0.48] ; ra22h30dec-55 ;
    15: ret = [0.65, 0.65, 0.48,   -1] ; ra23hdec-45
    16: ret = [  -1,  1.3,   -1, 0.48] ; ra6hdec-62.5
    17: ret = [0.65, 0.65, 0.48,   -1] ; ra3h30dec-42.5
    18: ret = [0.65, 0.65, 0.48,   -1] ; ra1hdec-42.5e
    19: ret = [  -1, 0.65, 0.48,   -1] ; ra6h30dec-45
    20: ret = [0.65, 0.65, 0.48, 0.48] ; ra5h30dec-55_2011
    else: print, "bad idx argument.  Return -1"
endcase

print, "overlap = ", ret

return, ret
end
