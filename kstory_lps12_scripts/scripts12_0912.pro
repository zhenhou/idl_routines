;;;
; NAME: script12_0912
; PURPOSE:
;   General script
;
; NOTES:
;  1) Calculate the most Gaussian power for A_L
;  2) make ML spectrum for WMAP7+SPT
;;;

;;;;;;;;;;;;;;;;;;
; 
; Get best-fit alens power
;
;;;;;;;;;;;;;;;;;;

PRO alens_pow, make_sav=make_sav


; Re-make sav file?
if keyword_set(make_sav) then begin
    cdir = '/data23/hou/lps12/paramfits/chains_0828/'
    dir = cdir+'c54_lcdm_alens_camb_w7s12/chains/'
    files = file_search(dir+'c54_lcdm_alens_camb_w7s12*.txt')
    pname = dir+'c54_lcdm_alens_camb_w7s12.paramnames'

    subsamp=10
    nskip=1000
    power=1.0

    npow=20
    pset= findgen(npow)*0.05 + 0.05
    ss = fltarr(npow)
    
    for ii=0, npow-1 do begin
        print, '*** power = ',pset[ii]
        plot_like1dname,files,pname,'alens',subsamp=subsamp,nskip=nskip,power=pset[ii],/get_chain,results=results
        histogauss, results.pp, a
        ss[ii] = a[4]
    endfor

    plot, pset, ss
    save, pset, ss, filename='tmp_alens_pow.sav'
    stop
endif else restore, 'tmp_alens_pow.sav'

; Get the min valke
plot, pset, ss, yr=[0.95, 1.]
xx = poly_fit(pset, ss, 2)
xvec = findgen(1000)/1000.

fit = xx[0] + xx[1]*xvec + xx[2]*(xvec^2.)
oplot, xvec, fit, color=!red
wh = where(fit eq min(fit))
print, '*** most-gaussian power = ', xvec[wh]

END



PRO make_w7s12_ML
readcol, '/home/cr/paramfits/cosmomc.lps12/suxp_cls', cl_tt, te, ee
nn = n_elements(cl_tt)
l = indgen(nn) + 2

dl = cl_tt * l * (l+1)/(2*!pi)
;readcol, '/home/kstory/lps12/best_fit/obsolete_0717/lcdm_w7s12_ML_dl.txt', l_old, dl_old

; write out new ML file
openw,lun,/get_lun,'/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt'
for i=0, nn-1 do begin
  printf, lun, l[i], dl[i]
endfor
close, lun

stop
END
