;;;
; NAME: script_0614
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Plots
; 2) Check iskip in end2end
;
; MODIFICATION HISTORY:
;  06/14/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;
; Plot error bars
;;;;;;;;;;;;;;;;;;;;;;

PRO plot_ps

; plotting setup
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10

;run_05
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_022240_kweight.sav' ; to be sure.
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_205636_kweight.sav' fix bug in cov_data*calib
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_213028_kweight.sav'
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120614_192137_kweight.sav' ; Should be final
dl_05  = dl_all
diag05 = diag_nobeam

; diag05_e2e = dblarr(58)
; cov_e2e = cov_all_nobeam_data + cov_all_nobeam_mc
; for ii=0, 57 do diag05_e2e[ii] = sqrt(cov_e2e[ii, ii])


;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_202229_kweight_calib.sav' ; e2e_calib
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_213316_kweight_calib.sav'
;dl05_calib = dl_all
;diag05_calib = diag_nobeam

; run_05; 0809 only
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_024104_kweight_0809.sav'; best guess
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_025500_kweight_0809.sav' ; no combined_calib applied to cov_data
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_040328_kweight_0809.sav' ; corrected cov_sv
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120613_222432_kweight_0809.sav' ; use cov_sv
;restore, '/home/kstory/lps12/end2end/run_05/obsolete/combined_spectrum_20120613_222844_kweight_0809.sav' ; use cov from e2e
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120614_193416_kweight_0809.sav' ; Should be right
dl_0809     = dl_all
diag05_0809 = diag_nobeam

; run_04
;restore, '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav'
;diag04 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam


; cl_theory
readcol,'/home/kstory/lps12/theory_spectrum/Cls_100sims_ave_from_combined_map.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

;----------------
; Plots
;----------------
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

; ps plot
window, 0
plot,l[whplot],dl_05[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_05'
oplot,l_vec,dl_th,color=!red
errplot,l[whplot],(dl_05-diag05)[whplot]*1d12,(dl_05+diag05)[whplot]*1d12
; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0614')

; window, 6
; plot,l[whplot],dl05_calib[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
;   xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_05'
; oplot,l_vec,dl_th,color=!red
; errplot,l[whplot],(dl05_calib-diag05_calib)[whplot]*1d12,(dl05_calib+diag05_calib)[whplot]*1d12
; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0614')

; ; compare against K11
window, 2
plot, l[whplot], ((err_rk )/diag05)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05', title='Error bars, K11 vs 05'
oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v05_0614')

; window, 3
; plot, l[whplot], ((err_rk )/diag05_e2e)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05', title='Error bars, K11 vs 05_e2e'
; oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
; ;err=tvread(/png,/nodialog,filename=fdir+'err_k11v05e2e_0614')

; window, 4
; plot, l[whplot], (diag05/diag05_e2e)[whplot], xtitle=xtitle, ytitle='diag05/diag05_e2e', title='Error bars, 05_sim vs 05_e2e'
; oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
; ;err=tvread(/png,/nodialog,filename=fdir+'err_05simv05e2e_0614')

window, 5
plot, l[whplot], (err_rk/diag05_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05_0809', title='Error bars, K11 vs 05_0809'


stop
END




;;;;;;;;;;;;;;;;;;;;;;
; Check iskip in end2end
;;;;;;;;;;;;;;;;;;;;;;

PRO check_iskip

;setup
f = lps12_fieldstruct()
idx = 4
field = f[idx].name

; output arrays
dls  = dblarr(9,58)
covs = dblarr(9,58,58)
invkerns = dblarr(9,58,58)


; loop over test values of iskip
for jskip=0, 8 do begin
print, 'CHECK_ISKIP: ', jskip
iskip_test = jskip+1
file = '/home/kstory/lps12/end2end/end_ra3h30dec-60_06_kweight_calib.sav'
restore, file

;;;;;;;;;;;;;;;;;;;;;;;;;
; RE-do beams
;;;;;;;;;;;;;;;;;;;;;;;;;

nsets=(size(setdef_data))[2]
beams=dblarr(n_elements(ellkern), nspectra)
simbeams=dblarr(n_elements(ellkern), nspectra)

if n_elements(beamfiles) ne nsets then begin
    message, 'Please specify one beamfile per set',n_elements(beamfiles),nsets
endif


beam_interp=dblarr(n_elements(ellkern), nsets)
simbeam_interp=dblarr(n_elements(ellkern), nsets)

for i=0, nsets-1 do begin
    beam=read_ascii(beamfiles[i])
; fold the effect of the lpf used prior to down-sampling, 
; which is a really small effect anyway, like <0.5% in power at all ells
    lpf_ds = az_avg_lpf_lps12(beam.field1[0,*], field, f0=f0)
    beam_interp[*, i]=interpol(beam.field1[1, *]*lpf_ds, $
                               beam.field1[0, *], $
                               ellkern)    
    simbeam=read_ascii(simbeamfiles[i])
    simbeam_interp[*, i]=interpol(simbeam.field1[1, *], $
                               simbeam.field1[0, *], $
                               ellkern)    

endfor

k=0
for i=0, nsets-1 do begin
    for j=i, nsets-1 do begin
        idx=where((beam_interp[*, i] ge 0) * $
                  (beam_interp[*, j] ge 0))
        beams[idx, k]=$
          sqrt(beam_interp[idx, i]*beam_interp[idx, j])
        idx=where((beam_interp[*, i] ge 0) * $
                  (beam_interp[*, j] ge 0))
        simbeams[idx, k]=$
          sqrt(simbeam_interp[idx, i]*simbeam_interp[idx, j])
        k++
    endfor
endfor


;;;;;;;;;;;;;;;;;;;;;;;;;
; Do iskip test
;;;;;;;;;;;;;;;;;;;;;;;;;

superkern=dblarr(nbands, nspectra, nbands, nspectra)
invkern=dblarr(nbands, nspectra, nbands, nspectra)
sim_superkern=dblarr(nbands, nspectra, nbands, nspectra)
sim_invkern=dblarr(nbands, nspectra, nbands, nspectra)
defaultskip=1
iskips = intarr(nspectra)
for i=0, nspectra-1 do begin 
    superkern[*, i, *, i]=rebin_coupling_matrix(kernel, ellkern, banddef, $
                                                transfer=transfer[*, i], $
                                                beam=beams[*, i])
    sim_superkern[*, i, *, i]=rebin_coupling_matrix(kernel, ellkern, banddef, $
                                                transfer=transfer[*, i], $
                                                beam=simbeams[*, i])

    wow = lindgen(nbands) + (1l*nbands)*nspectra*(lindgen(nbands))
    ;iskip = max([0,where(superkern(wow) eq 0)])+1
    iskip = iskip_test
    iskips[i]=iskip
    ;; leave the first (usually bogus) bin out of the inversion
    invkern[iskip:*, i, iskip:*, i]=invert(/double,reform(superkern[iskip:*, i, iskip:*, i]))
    sim_invkern[iskip:*, i, iskip:*, i]=invert(/double,reform(sim_superkern[iskip:*, i, iskip:*, i]))
endfor

invkernmat=double(reform(invkern, nbands*nspectra, nbands*nspectra))
siminvkernmat=double(reform(sim_invkern, nbands*nspectra, nbands*nspectra))

invkernmattr=double(transpose(invkernmat))
siminvkernmattr=double(transpose(siminvkernmat))

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 6: Apply the binned coupling kernel to the spectra
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

spectrum=reform(invkernmat##double(spectrum_data_raw), $
                nbands, nspectra)

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 7: Apply the binned coupling kernel to the covariances
;; 
;;;;;;;;;;;;;;;;;;;;;;;;;

sample_cov=reform(siminvkernmat##(double(cov_mc_raw)##siminvkernmattr),$
                  nbands, nspectra, nbands, nspectra)
meas_cov=reform(invkernmat##(double(cov_data_raw)##invkernmattr),$
                nbands, nspectra, nbands, nspectra)

;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; Step 8: Sum the covariance to get the total covariance
;;
;;;;;;;;;;;;;;;;;;;;;;;;;

cov=meas_cov+sample_cov

;store output:
dls[jskip,*] = spectrum
covs[jskip,*,*] = reform(cov)
invkerns[jskip,*,*] = reform(invkern)

endfor

save, dls, covs, invkerns, filename='iskip_test.sav'

stop
END


PRO anal_iskip
restore, 'iskip_test.sav'

diags = dblarr(9,58)
for j=0, 8 do begin
    for i=0, 57 do begin
        diags[j,i] = sqrt(covs[j,i,i])
    endfor
endfor

;;;;;;;;;;;
; plots
ca = kscolor_array()
l = indgen(58)*50 + 250

window, 1
plot, l, diags[1,*] / diags[0,*], ytitle='diags[ii] / diags[0]', title='ra3h30dec-60'
for ii=1, 7 do begin
    oplot, l, diags[ii,*] / diags[0,*], color=ca[ii]
endfor

window, 2
plot, l, dls[1,*] / dls[0,*], ytitle='dls[ii] / dls[0]', title='ra3h30dec-60'
for ii=1, 7 do begin
    oplot, l, dls[ii,*] / dls[0,*], color=ca[ii]
endfor

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
stop
END


