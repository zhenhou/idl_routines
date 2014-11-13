;;;
; NAME: get_lps12_fieldinfo.pro
; PURPOSE:
;   Return a struct with field information
;
; INPUTS: field_idx,      index in lps12_fieldstruct()
;
; OUTPUTS:
;   info,           struct
;
; NOTES:
;   The idea is to return a struct, and add to the struct as new
;   things I want to keep track of come up.
;
; MODIFICATION HISTORY:
;  03/07/2012: (KTS) Created
;  03/26/2012: (KTS) Change nbig from 2160 to 4320
;;;


;...................................................................
; Quickly get the number of pixels for each field
;   Note, currently the poly order is fixed to 7
function get_lps12_fieldinfo, field_idx

nbig=4320
case field_idx of

    0: begin ; ra5h30dec-55_2008
        info={npix:[960, 960], nbig:nbig, ds:6, lpf:7.5, dps_on_sky:0.275317, el:55}
    end
    1: begin ; ra23h30dec-55_2008
        info={npix:[960, 960], nbig:nbig, ds:6, lpf:7.5, dps_on_sky:0.275317, el:55}
    end
    2: begin ; ra23h30dec-55_2010
        info={npix:[960, 960], nbig:nbig, ds:5, lpf:9.0, dps_on_sky:0.338410, el:55}
    end
    3: begin ; ra21hdec-60
        info={npix:[1280, 960], nbig:nbig, ds:6, lpf:7.5, dps_on_sky:0.240000, el:60}
     end
    4: begin ; ra3h30dec-60
        info={npix:[2048, 1024], nbig:nbig, ds:6, lpf:7.5, dps_on_sky:0.240000, el:60}
    end
    5: begin ; ra21hdec-50
        info={npix:[1536, 960], nbig:nbig, ds:6, lpf:7.5, dps_on_sky:0.308538, el:50}
    end
    6: begin ; ra4h10dec-50
        info={npix:[1536, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:50}
    end
    7: begin ; ra0h50dec-50
        info={npix:[1536, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:50}
    end
    8: begin ; ra2h30dec-50
        info={npix:[1536, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:50}
    end
    9: begin ; ra1hdec-60
        info={npix:[1280, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:60}
    end
    10: begin ; ra5h30dec-45
        info={npix:[1024, 1024], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:45}
    end
    11: begin ; ra6h30dec-55
        info={npix:[960, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:55}
    end
    12: begin ; ra23hdec-62.5
        info={npix:[1152, 640], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:62.5}
    end
    13: begin ; ra21hdec-42.5
        info={npix:[1728, 640], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el: 42.5}
    end
    14: begin ; ra22h30dec-55
        info={npix:[960, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:55}
    end
    15: begin ; ra23hdec-45
        info={npix:[1728, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:45}
    end
    16: begin ; ra6hdec-62.5
        info={npix:[1152, 640], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:62.5}
    end
    17: begin ; ra3h30dec-42.5
        info={npix:[2160, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:42.5}
    end
    18: begin ; ra1hdec-42.5
        info={npix:[1728, 640], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:42.5}
    end
    19: begin ; ra6h30dec-45
        info={npix:[1024, 1024], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:45}
    end
    20: begin ; ra5h30dec-55_2011
        info={npix:[960, 960], nbig:nbig, ds:4, lpf:11.25, dps_on_sky:0.417812, el:55}
    end
endcase

return, info
end
