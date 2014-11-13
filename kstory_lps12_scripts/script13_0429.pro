;;;
; NAME: script13_0429.pro
; PURPOSE:
;   make ptsrc masks
;
; NOTES:
;;;

;;;;;;;;;;;;;;;;;;;
; Make point source masks
;;;;;;;;;;;;;;;;;;;
PRO mk_ptsrc_mask
list = [1,2,7,9,12,15,18]
n = n_elements(list)

for i=0, n-1 do begin
    idx = list[i]
    print, "*************", idx
    make_ptsrc_mask_lps12, idx, config_dir='/home/kstory/lps12/ptsrc_lists/20130428_corrected/', savdir='/data/kstory/projects/lps12/masks/0428_ptscor/'
endfor
;make_ptsrc_mask_lps12, 1, config_dir='/home/kstory/lps12/ptsrc_lists/20130428_corrected/', savdir='/data/kstory/projects/lps12/masks/0428_ptscor/', /stopit



END

