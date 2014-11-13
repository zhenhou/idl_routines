;;;
; NAME: script_0423
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) make apod masks for 2008, 2009 fields
;
; MODIFICATION HISTORY:
;  04/23/2012: (KTS) Created
;;;

;...................................................................
; Make apod masks
PRO make_apod_masks

f = lps12_fieldstruct()
for ii=0, 19 do begin
    fst = f[ii]

    ; get old runlist name
    if fst.lead_trail then begin
        rx = '/home/kstory/lps12/runlists/lt/runlist_lt_'+fst.name+'.txt'
    endif else begin
        rx = '/home/kstory/lps12/runlists/runlist_xspec_lps12_'+fst.name+'.txt'
    endelse

    ; get new runilst name
    bdir = '/data/kstory/projects/lps12/runlists/xspec2/'
    rx2 = bdir+'runlist_xspec2_lps12_'+fst.name+'.txt'

    print, 'old: ', rx
    print, 'new: ', rx2
    print, ' '
    spawn, 'ln -s '+rx+' '+rx2
endfor

END


;...................................................................
; Check progress of current idf job
PRO sav2txt_lt_runlist
         stopit=stopit
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()
idx_list = [3,4,5]
nfields = n_elements(idx_list)

; dirs
bdir = '/home/kstory/lps12/runlists/lt/'

; Write out runlists into .txt file
;for ii=0, nfields-1 do begin

END



