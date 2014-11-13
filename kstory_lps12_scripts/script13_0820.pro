;;;
; NOTES:
;  1) Look at lrange differences
;  2) Make ini files lrange, to make ML spectra
;;;

PRO comp_lr
;--------------------
; 
;--------------------
restore, '/home/kstory/lps12/scripts/sav_files/lrange_0819.sav'

np=6
for i=0, np-1 do begin
    print, params[i]
    err = sqrt(lmx15.err[i]^2. + lmx30.err[i]^2.)
    ss = lmx15.value[i] - lmx30.value[i]
    print, ss/err
endfor

print, (lmx30.value - lmx15.value) / sqrt(lmx15.err^2. + lmx30.err^2.)

stop
END



;;;;;;;;;;;;;;;;
; 2) Make ini files lrange, to make ML spectra
;;;;;;;;;;;;;;;;
PRO make_ML_lrange
chain = 'c2_lcdm_pico_w7s12_lmax1500'
chain_dir = '/data23/hou/lps12/paramfits/chains_0828/c2_lcdm_pico_w7s12_lmax1500/chains/'
make_ML_ini, chain, chain_dir=chain_dir

chain = 'c2_lcdm_pico_w7s12_1500ell3000'
chain_dir = '/data23/hou/lps12/paramfits/chains_0828/c2_lcdm_pico_w7s12_1500ell3000/chains/'
make_ML_ini, chain, chain_dir=chain_dir

END
