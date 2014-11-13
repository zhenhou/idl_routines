;;;
; NAME: script_0523
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) run tf_prep to make tf_prep04.sav
; 2) make_tf_prep04, make tf_prep04.sav
; 3) make_tfs
; 4) kern_0to6, etc, make coupling kernels
; 5) noise_0toblah, make noise psd's
;
; MODIFICATION HISTORY:
;  05/23/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Processing
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO init_tf_prep ; done
s0=-1 & s1=-1 & s2=-1 & s3=-1 & s4=-1 & s5=-1 & s6=-1 & s7=-1 & s8=-1 & s9=-1
s10=-1 & s11=-1 & s12=-1 & s13=-1 & s14=-1 & s15=-1 & s16=-1 & s17=-1 & s18=-1 & s19=-1

savfile = 'tf_prep04.sav'

save,s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19,$ ;l_mode_coupling,mode_couplings, $
  filename=savfile
END

PRO make_tf_prep04 ; done
for idx=0, 19 do begin
    print, 'make tf_prep, field ', idx
    t0 = systime(0, /seconds)

    tf_prep, idx, '04'

    print, 'That tf_prep took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO make_tfs ; done
for idx=0, 19 do begin
    print, 'make TF, field ', idx
    t0 = systime(0, /seconds)

    make_twod_tfs, idx, '04', /dosave

    print, 'That TF took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


;;;;;;;;;;;;;;;;;;
; Kernels
;;;;;;;;;;;;;;;;;;

PRO kern_0to5 ; done
for idx=0, 5 do begin
    print, 'make kern, field ', idx
    t0 = systime(0, /seconds)

    make_coupling_kernel, idx

    print, 'That kern took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO kern_6to11 ; done
for idx=6, 11 do begin
    print, 'make kern, field ', idx
    t0 = systime(0, /seconds)

    make_coupling_kernel, idx

    print, 'That kern took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO kern_12to19 ; done
for idx=12, 19 do begin
    print, 'make kern, field ', idx
    t0 = systime(0, /seconds)

    make_coupling_kernel, idx

    print, 'That kern took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

;;;;;;;;;;;;;;;;;;
; noise
;;;;;;;;;;;;;;;;;;

PRO noise_0to9 ; done
for idx=0, 9 do begin
    print, 'make noise, field ', idx
    t0 = systime(0, /seconds)

    make_noise_psd_lps12, idx

    print, 'That noise took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO noise_10to19 ; done
for idx=10, 19 do begin
    print, 'make noise, field ', idx
    t0 = systime(0, /seconds)

    make_noise_psd_lps12, idx

    print, 'That noise took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

;;;;;;;;;;;;;;;;;;
; kweights
;;;;;;;;;;;;;;;;;;

PRO kweight_0to19 ; done
for idx=0, 19 do begin
    print, 'make kweight, field ', idx
    t0 = systime(0, /seconds)

    make_twod_kweights_lps12, idx, /save_files

    print, 'That kweight took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END
