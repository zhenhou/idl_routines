;;;
; function to read fits files
; 03/14/2012: (KTS) Created
;;;

function krf, map
return, expand_fits_struct(read_spt_fits(map))
end
