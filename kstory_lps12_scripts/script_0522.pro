;;;
; NAME: script_0521
; PURPOSE:
;   General script for today
;
; NOTES:
;
;
; MODIFICATION HISTORY:
;  05/21/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;----------------------------------
;;; combined error bars
PRO comp_err_comb
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120520_001439_kweight.sav' & err_ks = diag_nobeam
;restore, '/home/kstory/lps12/end2end/run_03/obsolete/combined_spectrum_20120521_203118_kweight.sav'
;restore,
;'/home/kstory/lps12/end2end/run_03/obsolete/combined_spectrum_20120521_212127_kweight.sav' ; something wrong
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120521_221658_kweight.sav'
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120521_223459_kweight.sav'
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120521_224750_kweight.sav'
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120522_150814_kweight.sav' ; set cov10, 11 = 09
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120522_152534_kweight.sav' ; use run_03_badBeam, set cov10, 11 = 09
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120522_153713_kweight.sav' ; wrong beam, wrong cov
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120522_155024_kweight.sav' ; right beam, right cov
;restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120522_184200_kweight.sav' ; 0809


restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120522_184600_kweight.sav' ; all fields
err_ks = diag_nobeam
restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120522_184955.sav' ; all fields, no kweight
err_nokw_ks = diag_nobeam
restore, '/home/kstory/lps12/end2end/run_03/combined_spectrum_20120522_184732_kweight_0809.sav' ; 0809
err_0809_ks = diag_nobeam

;restore, '/home/kstory/lps12/end2end/run_03_badBeam/run_03/combined_spectrum_20120511_200438_kweight.sav'
;err_nobeam_ks_old = diag_nobeam & err_ks_old = diag

restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam

l_bdf = l
;x = (err_nobeam_rk - err_nobeam_ks)/err_nobeam_rk
x = err_rk / err_ks
v = 6+indgen(58-6)
res = linfit( l_bdf[v], x[v] )

; PLOTS
window, 1
plot, l_bdf[v], x[v], title='change in error lps12 to k11', xtitle='ell', ytitle='err_k11 / err_lps12'
oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=3

fdir='/home/kstory/public_html/notebook/spt_lps12/'
;err=tvread(/png,/nodialog,filename=fdir+'err_0522')
stop
END


PRO sss
for ii=0, 5 do begin
    if ii eq 3 then continue
    print, ii
endfor
END
