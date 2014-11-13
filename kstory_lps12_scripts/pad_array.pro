;;;
; function to pad an array
; 04/25/2012: (KTS) created
;;;

FUNCTION pad_array, array, nbig

; get the data type
type = datatype(array)

case type of
    'DOU' : array_padded = dblarr(nbig, nbig)
    'FLO' : array_padded = fltarr(nbig, nbig)
    else : print, 'datatype of array, '+type+', not recognized.  Returning -1'
endcase

nx = (size(array))[1]
ny = (size(array))[2]
array_padded[0:nx-1, 0:ny-1] = array

RETURN, array_padded
END
