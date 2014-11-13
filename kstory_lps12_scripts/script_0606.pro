;;;
; NAME: script_0606
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Plotting azrms 10% cuts
; 2) make_expected_dls_run_05
;
; MODIFICATION HISTORY:
;  06/06/2012: (KTS) Created
;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot azrms 10% cuts
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_azrms


files = ['/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav', $
         '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_145116_kweight_azrms_90.sav', $
         '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_145300_kweight_randcut_90_seed1.sav', $
         '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_145423_kweight_randcut_90_seed2.sav', $
         '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_150600_kweight_randcut_90_seed3.sav']
nfiles = 5
nbins = 58

; make arrays
dls   = dblarr(nfiles, nbins)
diags = dblarr(nfiles, nbins)

for ii=0, nfiles-1 do begin
    restore, files[ii]
    dls[ii,*] = dl_all
    diags[ii,*] = diag_nobeam
endfor


v1 = indgen(50)+8
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
color_arr = [!white, !red, !skyblue, !lavender, !green]

window, 1, xsize=700, ysize=500
plot, l[v1], (dls[0,*])[v1], /ylog, xtitle='ell', ytitle='dDl [K^2]', title='Simulated SV errors'
for ii=0, nfiles-1 do begin
    oplot, l[v1], (dls[ii,*])[v1], color=color_arr[ii]
endfor
legend, ['full_set', 'azrms_90', 'rand_1', 'rand_2', 'rand_3'], linestyle=[0,0,0,0,0], $
  colors=[!white, !red, !skyblue, !lavender, !green], pos=[2000, 5.e-9]
err = tvread(/png,/nodialog,filename=fdir+'azrms_ps_0606')

window, 2
plot, l[v1], ((dls[0,*] - dls[1,*])/dls[0,*])[v1], $
  xtitle='ell', ytitle='(dl_full - dl_new) / dl_full', title='Fraction change in dls'
for jj=1, nfiles-1 do begin
    oplot, l[v1], ((dls[0,*] - dls[jj,*])/dls[0,*])[v1], color=color_arr[jj]
endfor
legend, ['azrms_90', 'rand_1', 'rand_2', 'rand_3'], linestyle=[0,0,0,0], $
  colors=[!red, !skyblue, !lavender, !green], pos=[600, 0.025]
err = tvread(/png,/nodialog,filename=fdir+'azrms_ps_frac_0606')

stop
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make expected_dls sav files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO make_expected_dls_run_05
jackdir = '/home/kstory/lps12/jacks/sims/'
run='05'

nbin=5
expected_dls=fltarr(nbin,20)

stub = 'lr'
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
save, expected_dls, filename=savfile

stub = '12'
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
save, expected_dls, filename=savfile

stub = 'azrms'
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
save, expected_dls, filename=savfile

stub = 'azrms_95'
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
save, expected_dls, filename=savfile

stub = 'azrms_90'
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
save, expected_dls, filename=savfile

stub = 'tweight'
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
save, expected_dls, filename=savfile

stub = 'moon'
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
save, expected_dls, filename=savfile

stub = 'sun'
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
save, expected_dls, filename=savfile

END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate pixel window function for sims, run_04
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO mk_area
run='04'
; setup
f = lps12_fieldstruct()
end_dir = '/home/kstory/lps12/end2end/'
nbin = 58
nsim = 100
;calib = 0.76
reso = 1.0

; Get list of fields
idx_list = [0,1,indgen(17)+3]
nfields = n_elements(idx_list)

; Make arrays to store output
area        = dblarr(nfields)

; retrieve the information        
for ilist=0, nfields-1 do begin
    idx = idx_list[ilist]
    print, 'process field ',idx, ', '+f[idx].name
    restore, end_dir+'end_'+f[idx].name+'_'+run+'_kweight.sav'

    ; get the spectra
    area[ilist] = total(mask_padded)*(reso/60.)^2.
endfor

save, area, 'tmp.sav'
END


; 
PRO sss
restore, 'tmp1.sav'

; get field weights
w0 = dblarr(nfields,nbin) ; weight array
for i=0,nfields-1 do w0[i,*] = area[i]
for i=0,nbin-1 do w0[*,i] /= total(w0[*,i])

; combine
for i=0,nbin-1 do begin
    tmp = dblarr(nwf)
    for j=0,nfields-1 do begin
        tmp += w0[j,i]*wfs_sim[*,i,j]
    endfor
    wf_all_sim[*,i] = tmp
endfor
; save
savfile = end_dir+'wf_all_sim_run_'+run+'.sav'
print, 'Saving WF in file: ' + savfile
save, wf_all_sim, l_wf_sim, filename=savfile

stop
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Re-run end2end run_05
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO end2end_run_05_0606
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    if idx eq 0 then idx++
    if idx eq 1 then idx++
    print, 'run end2end_05, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_05, idx, /resume, /use_kweight, out_dir='/data18/kstory/lps12/end2end/'

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Jacks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;--------------------------
; Jackknives for tonight
PRO night_jacks_1
; azrms
; for idx=0, 19 do begin
;     print, 'jack azrms, field ', idx
;     t0 = systime(0, /seconds)
;     lps12_jack, idx, 'azrms'
;     print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
; endfor

; azrms_95
for idx=0, 19 do begin
    print, 'jack azrms_95, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_95', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; azrms_90
for idx=0, 19 do begin
    print, 'jack azrms_90, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_90', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; lr
for idx=0, 19 do begin
    print, 'jack lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

END



;;; night 2
PRO night_jacks_2
; 12
for idx=18, 19 do begin ; only last 2 files left to run
    print, 'jack 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; tweight
for idx=0, 19 do begin
    print, 'jack tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; moon
for idx=0, 19 do begin
    print, 'jack moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; sun
idx_list = [4, 6, 11, 17, 18, 19]
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', savdir='/data18/kstory/lps12/jacks/'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

;;;--------------------------
; Sim Jackknives for tonight
PRO night_sim_jacks

; ; coadd azrms
; for idx=0, 19 do begin
;     print, 'coadd jack sims azrms, field ', idx
;     coadd_lps12_jack_sims, idx, 'azrms'
; endfor
; ; azrms sim jack
; for idx=0, 19 do begin
;     print, 'jack sims azrms, field ', idx
;     t0 = systime(0, /seconds)
;     lps12_jack, idx, 'azrms', /sim
;     print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
; endfor


; coadd azrms_95
for idx=0, 19 do begin
    print, 'coadd jack sims azrms_95, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms_95', out_dir='/data18/kstory/lps12/sims/'
endfor
; azrms_95 sim jack
for idx=0, 19 do begin
    print, 'jack sims azrms_95, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_95', /sim, sim_dir='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd azrms_90
for idx=0, 19 do begin
    print, 'coadd jack sims azrms_90, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms_90', out_dir='/data18/kstory/lps12/sims/'
endfor
; azrms_90 sim jack
for idx=0, 19 do begin
    print, 'jack sims azrms_90, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms_90', /sim, sim_dir='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd moon
for idx=0, 19 do begin
    print, 'coadd jack sims moon, field ', idx
    coadd_lps12_jack_sims, idx, 'moon', out_dir='/data18/kstory/lps12/sims/'
endfor
; jack moon
for idx=0, 19 do begin
    print, 'jack sims moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', /sim, sim_dir='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor


; coadd sun
idx_list = [4, 6, 11, 17, 18, 19]

for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'coadd jack sims sun, field ', idx
    coadd_lps12_jack_sims, idx, 'sun', out_dir='/data18/kstory/lps12/sims/'
endfor
; jack sun
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sims sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', /sim, sim_dir='/data18/kstory/lps12/sims/', savdir='/data18/kstory/lps12/jacks/sims/'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

END
;;;--------------------------

