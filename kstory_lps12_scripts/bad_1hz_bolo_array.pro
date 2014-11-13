;;;
; NAME: bad_1hz_bolo_array
; PURPOSE:
;   Returs an array of the bad 1hz bolo indices
;
; CALLING SEQUENCE: bad_bolo_arr = bad_1hz_bolo_array(2009)
;
; INPUTS: year,       integer in the set [2008, 2009, 2010, 2011]
;
; OUTPUTS:
;   an array of bolo indices
;
; NOTES:
;   1) access data in the following way:
;      bad_bolo_arr = bad_1hz_bolo_array(2009)
;
; MODIFICATION HISTORY:
;   
;   02/16/2011: (KTS) Created
;;;


;...................................................................
; Main function
;
function bad_1hz_bolo_array, year

idx = 0
case year of
    2008: idx = 1
    2009: idx = 2
    2010: idx = 2
    2011: idx = 2

    ; Throw an error
    else: print, "bad_1hz_bolo_array.pro: Bad value for year argument.  Returning -1."
endcase

case idx of 
    ; 2008
    1: boloset =  [22,  58,  66,  99, 116, 127, 129, 137, 326, 341, 355, 356, 362, 370, 426, 472, 502, 529, 543, 639]

    ; 2009, 2010, 2011
    2: boloset =  [58, 151, 158, 166, 170, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 220, 235, 268, 362, 369, 389, 407, 465, 644, 768,  769]

    ; Error
    0: return, -1
endcase

return, boloset
    
end
