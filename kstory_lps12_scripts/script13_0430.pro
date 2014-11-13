;;;
; NAME: script13_0430.pro
; PURPOSE:
;   make ptsrc masks
;
; NOTES:
;  1) Check point-source masks
;;;

;;;;;;;;;;;;;;;;;;;
;
; Check ptsrc masks
;
;;;;;;;;;;;;;;;;;;;
PRO check_pt
f=lps12_fieldstruct()
list = [1,2,7,9,12,15,18]
n = n_elements(list)

for i=0, n-1 do begin
    idx = list[i]
    print, "*************", idx, ', '+f[idx].name
    savdir='/data/kstory/projects/lps12/masks/0428_ptscor/'
    restore, savdir+'ptsrc_mask_'+f[idx].name+'_50mJy_proj5.sav'
    pt1 = ptsrc_mask

    restore,'/home/kstory/lps12/masks/ptsrc/ptsrc_mask_'+f[idx].name+'_50mJy_proj5.sav'
    pt0 = ptsrc_mask

    tv_spt_map, pt1-pt0, /norms, winnum=2
    stop
endfor
END



;;;;;;;;;;;;;;;;;;;
;
; Make apod masks
;
;;;;;;;;;;;;;;;;;;;
PRO make_apod_ptscor
f=lps12_fieldstruct()
list = [1,2,7,9,12,15,18]
n = n_elements(list)

ptdir = '/data/kstory/projects/lps12/masks/0428_ptscor/'
apoddir = '/data/kstory/projects/lps12/masks/apod/'

for i=0, n-1 do begin
    idx = list[i]
    MAKE_MASK_LPS12, field_idx=idx, $
      apod_dir='/data/kstory/projects/lps12/masks/apod/', $
      ptsrc_dir='/data/kstory/projects/lps12/masks/0428_ptscor/', $
      outdir = '/data/kstory/projects/lps12/masks/0428_ptscor/'

endfor
print, "end of make_apod_ptscor"
stop
END    

; Check apod_ptscor
PRO check_apod_ptscor
f=lps12_fieldstruct()
list = [1,2,7,9,12,15,18]
n = n_elements(list)
ptscor_dir = '/data/kstory/projects/lps12/masks/0428_ptscor/'
orig_dir   = '/data/kstory/projects/lps12/masks/masks_50mJy/'

for i=0, n-1 do begin
    idx = list[i]
    print, '*** ',i,', '+f[idx].name
    restore, orig_dir+'mask_'+f[idx].name+'.sav'
    m0 = mask

    restore, ptscor_dir+'mask_'+f[idx].name+'.sav'
    m1 = mask

    tv_spt_map, m0 - m1, winnum=2
    stop
end
stop
END




;;;;;;;;;;;;;;;;;;;
;
; run end2end
;
;;;;;;;;;;;;;;;;;;;
PRO run_e2e
list = [7,9,18]
n = n_elements(list)
for i=0, n-1 do begin
    idx = list[i]
    
    print, 'run try_end2end_ptscor, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_ptscor, idx, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

END

