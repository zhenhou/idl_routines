;;;
; NAME: script_0607
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) plot_sv_err
; 2) compare run_04 to run_05 stuff
;
; MODIFICATION HISTORY:
;  06/07/2012: (KTS) Created
;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot sv error bars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_sv_err

;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120608_174259_kweight.sav'
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120608_220412_kweight.sav'
;restore, '/home/kstory/lps12/end2end/sv_err_05_kweight.sav' ; gets sv_err

v1 = indgen(50)+8
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

window, 1, xsize=700, ysize=500
plot, l[v1], sqrt(cov_all_nobeam[v1,v1]), /ylog, xtitle='ell', ytitle='dDl [K^2]', title='Simulated SV errors 05'
oplot, l[v1], sqrt(cov_all_nobeam_e2e[v1,v1]), color=!red                                                    
oplot, l[v1], sqrt(cov_all_nobeam_data[v1,v1]), color=!skyblue
oplot, l[v1], sqrt(cov_all_nobeam_mc[v1,v1]), color=!lavender
oplot, l[v1], sv_err[v1], color=!green
oplot, l[v1], sqrt(cov_all_nobeam[v1,v1])
legend, ['Sim_total_err', 'original_err', 'data_err', 'mc_err', 'sv_err'], linestyle=[0,0,0,0,0], $
  colors=[!white, !red, !skyblue, !lavender, !green], pos=[2000, 5.e-11]

window, 2
plot, l[v1], sqrt(cov_all_nobeam[v1,v1]) / sqrt(cov_all_nobeam_e2e[v1,v1]) - 1., xtitle='ell', ytitle='(err_SVsim - err_old) / err_old', title='Fraction change in error bars'

stop
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot sv error bars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO comp_0405
; run_04
restore, '/data/kstory/projects/lps12/end2end/sv_err_04_kweight.sav' ; get dl_all
restore, '/data/kstory/projects/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav'
cov_all_nobeam_04 = cov_all_nobeam
cov_all_nobeam_tot_04 = cov_all_nobeam_tot
cov_all_nobeam_data_04 = cov_all_nobeam_data
cov_all_nobeam_mc_04 = cov_all_nobeam_mc
sv_err_04 = sv_err

; run_05
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120608_174259_kweight.sav'
cov_all_nobeam_05 = cov_all_nobeam
cov_all_nobeam_tot_05 = cov_all_nobeam_tot
cov_all_nobeam_data_05 = cov_all_nobeam_data
cov_all_nobeam_mc_05 = cov_all_nobeam_mc
sv_err_05 = sv_err
stop
END


PRO comp_tf, idx
f = lps12_fieldstruct()

restore, '/home/kstory/lps12/end2end/end_'+f[idx].name+'_04_kweight.sav'
tf_04 = transfer
restore, '/home/kstory/lps12/end2end/end_'+f[idx].name+'_05_kweight.sav'
tf_05 = transfer

plot, ellkern, (tf_04 / tf_05) - 1., yr=[-0.02, 0.02], xr=[0,4500], title=f[idx].name, ytitle='(TF04 - TF05)/TF05)'

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
err=tvread(/png,/nodialog,filename=fdir+'tf_'+strtrim(string(idx),2)+'_0608')
END


PRO comp_err
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120608_220412_kweight.sav'
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120608_231258_kweight.sav'

df=dblarr(4,58)
de=dblarr(4,58)
dd=dblarr(4,58)
dm=dblarr(4,58)

for i=0, 57 do begin
    df[*,i] = sqrt(cov_y_nobeam_final[*,i,i])
    de[*,i]=sqrt(cov_y_nobeam_e2e[*,i,i])
    dd[*,i]=sqrt(cov_y_nobeam_data[*,i,i])
    dm[*,i]=sqrt(cov_y_nobeam_mc[*,i,i])
endfor

v1 = indgen(50) + 8
plot, l, de[0,*], /ylog, yr=[10e-13, 10e-9]
stop
END
