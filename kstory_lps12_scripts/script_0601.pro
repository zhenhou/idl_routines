;;;
; NAME: script_0601
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) first pipeline check.
;
; MODIFICATION HISTORY:
;  06/01/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; calculate sample variance from sims
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Compare the calculated sim-spec to the input.
PRO pipe_1

; setup
nbin=58
end_dir = '/home/kstory/lps12/end2end/'

; Get calculated sim-spec
restore, end_dir+'sv_err_04_kweight.sav' ; get dl_all
l = banddef-25
mndl = fltarr(nbin)
for ii=0, nbin-1 do mndl[ii] = mean(dl_all[ii,*])
mndl *= 1d12

; get theory
readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_uK2
dl_th = interpol(dl_uK2, l_vec, l)


; plot
v = indgen(48)+8

window, 0
plot, l[v], dl_th[v], /ylog, title='sim and theory spectrum', xtitle='ell', ytitle='Dl [uK^2]'
oplot, l[v], dl_th[v], color=!red
oplot, l[v], mndl[v]
legend, ['Dl_sim', 'Dl_th'], colors=[!white, !red], pos=[5000], linestyle=[1,1]

window, 2
plot, l[v], (mndl/dl_th)[v]-1., ytitle='(sim-dl_th) / dl_th', xtitle='ell'


fdir = '/home/kstory/public_html/notebook/spt_lps12/'
stop
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot sim error bars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_sv_error
nbin=58
end_dir = '/home/kstory/lps12/end2end/'

; Get calculated sim-spec
restore, end_dir+'sv_err_04_kweight.sav' ; get dl_all
l = banddef-25

restore, end_dir+'run_04/combined_spectrum_20120531_194921_kweight.sav'

; plot
v = indgen(48)+8

window, 0
plot, l[v], sv_err[v], /ylog, title='Error bars', xtitle='ell', ytitle='err [K^2]'
oplot, l[v], diag_nobeam[v], color=!red
legend, ['err_sim', 'err_total_e2e'], colors=[!white, !red], pos=[1500,7e-11], linestyle=[1,1]


window, 2
plot, l[v], ((sv_err - diag_nobeam)/diag_nobeam)[v], ytitle='(err_sim - err_e2e) / err_e2e', xtitle='ell'

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
stop
END

