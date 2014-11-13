;;;
; NAME: pipeline_test1
; PURPOSE:
;   Check that the power spectrum of the sims re-produce the input spectrum
;
; NOTES:
; 1) first pipeline check.
;
; MODIFICATION HISTORY:
;  06/06/2012: (KTS) Created
;  06/08/2012: (KTS) use 'wf2use' to choose window function.
;;;

PRO pipeline_test1, run

; setup
nbin=58
end_dir = '/home/kstory/lps12/end2end/'

; Get the bandpower window function, wf_all_sim
case run of
    '04' : begin
        restore, '/home/kstory/lps12/end2end/wf_all_sim_run_04.sav'
        wf2use = wf_all_sim
        ;restore, '/data/kstory/projects/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav'
    end
    '05' : begin
        ;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120608_174259_kweight.sav'
        ;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120608_193944_kweight.sav'
        ;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120608_194306_kweight.sav'
        restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_022240_kweight.sav'
        wf2use = wf_all_sim
    end
    '07' : begin
        restore, '/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_144542_kweight.sav'
        wf2use = wf_all_sim
    end
    '09' : begin
        restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120828_170101_kweight.sav'
        wf2use = wf_all_sim
    end
    else : begin
        print, 'un-recognized run.  Return -1'
        RETURN
    endelse
endcase

; Get calculated sim-spec
restore, end_dir+'cov_sv_extra_'+run+'_kweight.sav' ; get dl_all
l = this_banddef-25
mndl = fltarr(nbin)
for ii=0, nbin-1 do mndl[ii] = mean(this_dl_all[ii,*])
mndl *= 1d12

; get theory
readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_uK2

;dl_th = interpol(dl_uK2, l_vec, l)
;convolve the theory spectrum with the bandpower window function
dl_uK2_2 = dl_uK2[50:3998]
l_vec_2  = l_vec[50:3998]

dl_th_2 = dl_uK2_2 # wf2use

; plot
v = indgen(48)+8

window, 0
plot, l[v], dl_th_2[v], /ylog, title='sim and theory spectrum, run_'+run, xtitle='ell', ytitle='Dl [uK^2]', thick=2,chars=1.8,charthick=1
oplot, l[v], dl_th_2[v], color=!red,thick=2
oplot, l[v], mndl[v],thick=2
legend, ['Dl_sim', 'Dl_th'], colors=[!white, !red], pos=[2000, 5000], linestyle=[0,0],thick=2,charthick=1

window, 2, xsize=800, ysize=800
!p.multi=[0,1,2]
;plot, l[v], (mndl/dl_th_2)[v]-1., ytitle='(sim-dl_th) / dl_th', xtitle='ell', title='Pipeline Test1, run_'+run
plot, l[v], ((mndl-dl_th_2)/dl_th_2)[v], ytitle='(sim-dl_th) / dl_th', xtitle='ell', title='Pipeline Test1, run_'+run, thick=2,chars=1.8,charthick=1
oplot, l[v], ((mndl-dl_th_2)/dl_th_2)[v], color=!red,thick=2
plot, l[v], (mndl/dl_th_2)[v]-1., ytitle='(sim-dl_th) / dl_th', xtitle='ell', title='Pipeline Test1, run_'+run+', zoomed out', yr=[-0.02, 0.02], thick=2,chars=1.8,charthick=1
oplot, l[v], ((mndl-dl_th_2)/dl_th_2)[v], color=!red,thick=2
!p.multi=0

window, 3, xsize=800, ysize=800
!p.multi=[0,1,2]
;plot, l[v], (mndl/dl_th_2)[v]-1., ytitle='(sim-dl_th) / dl_th', xtitle='ell', title='Pipeline Test1, run_'+run
plot, l[v], (mndl-dl_th_2)[v], ytitle='sim-dl_th', xtitle='ell', title='Pipeline Test1, run_'+run, thick=2,chars=1.8,charthick=1
oplot, l[v], (mndl-dl_th_2)[v], color=!red,thick=2
;plot, l[v], (mndl)[v]-1., ytitle='(sim-dl_th) / dl_th', xtitle='ell', title='Pipeline Test1, run_'+run+', zoomed out', yr=[-0.02, 0.02], thick=2,chars=1.8,charthick=1
;oplot, l[v], (mndl-dl_th_2)[v], color=!red,thick=2
!p.multi=0

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
stop
END
