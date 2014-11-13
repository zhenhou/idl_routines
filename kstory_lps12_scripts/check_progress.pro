;;;
; NAME: check_progress.pro
; PURPOSE:
;   Check processing progress
;
; INPUTS: all optional
;
; NOTES:
;   hard-coded run = '03'
;
; MODIFICATION HISTORY:
;  05/08/2012: (KTS) Created
;;;


;...................................................................
; Check the progress of processing
PRO check_progress, $
                    field_idx=field_idx, $   ; check everything on one field
                    end2end=end2end, $       ; check first end2end run
                    end2end_kw=end2end_kw, $ ; check end2end with kweight
                    kern=kern, $             ; check coupling kernel
                    tf_prep=tf_prep, $       ; check tf_prep file
                    transfer=transfer, $     ; check twod_tf
                    kweight=kweight, $       ; check twod_kweight
                    sims=sims, $             ; check coadded sims file
                    noise=noise              ; check noise_psd

bdir = '/home/kstory/lps12/'
f = lps12_fieldstruct()
run='03'

if n_elements(field_idx) ne 0 then begin
    restore, bdir+'scripts/tf_prep03.sav' ; for tf_prep
    ii = field_idx
    print, 'Check field ', ii, ', '+f[ii].name
    print, 'end2end:  ', file_test(bdir+'end2end/end_'+f[ii].name+'_03.sav')
    print, 'kern:     ', file_test(bdir+'masks/masks_50mJy/kern_'+f[ii].name+'.sav')
    print, 'tf_prep:  ',  maxmin(mode_couplings[ii,*])
    ex=execute('help, s'+strtrim(string(ii),2))
    print, 'transfer: ', file_test(bdir+'twod_tfs/tf_'+f[ii].name+'.sav')
    print, 'kweight:  ', file_test(bdir+'twod_kweights/weight_2d_'+f[ii].name+'.sav')
    print, 'noise:    ', file_test(bdir+'noise_psds/noise_psd_'+f[ii].name+'.sav')
    print, 'sims:     ', file_test(bdir+'sims/coaddsim_'+f[ii].name+'.dat')
endif


if keyword_set(end2end) then begin
    for ii=0, 19 do print, ii, file_test(bdir+'end2end/end_'+f[ii].name+'_03.sav')
endif

if keyword_set(end2end_kw) then begin
    for ii=0, 19 do print, ii, file_test(bdir+'end2end/end_'+f[ii].name+'_03_kweight.sav')
endif

if keyword_set(kern) then begin
    for ii=0, 19 do print, ii, file_test(bdir+'masks/masks_50mJy/kern_'+f[ii].name+'.sav')
endif

if keyword_set(tf_prep) then begin
    restore, bdir+'scripts/tf_prep03.sav'
    for ii=0, 19 do begin
        print, ii, maxmin(mode_couplings[ii,*])
        ex=execute('help, s'+strtrim(string(ii),2))
    endfor
endif

if keyword_set(transfer) then begin
    for ii=0, 19 do print, ii, file_test(bdir+'twod_tfs/tf_'+f[ii].name+'.sav')
endif

if keyword_set(kweight) then begin
    for ii=0, 19 do print, ii, file_test(bdir+'twod_kweights/weight_2d_'+f[ii].name+'.sav')
endif

if keyword_set(noise) then begin
    for ii=0, 19 do print, ii, file_test(bdir+'noise_psds/noise_psd_'+f[ii].name+'.sav')
endif

if keyword_set(sims) then begin
    for ii=0, 19 do print, ii, file_test(bdir+'sims/coaddsim_'+f[ii].name+'.dat')
endif

END

