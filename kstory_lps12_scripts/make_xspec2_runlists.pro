;;;
; NAME: make_xspec2_runlists
; PURPOSE:
;   Make xspec2 runlists, which have coadded dates by xspec map
;
; NOTES:
; 1) output: /data/kstory/projects/lps12/runlists/xspec2/runlist_xspec2_lps12_[field].txt
; 2) sav2txt_lt_runlists, takes runlists from
;    /home/kstory/lps12/runlists/lt/list_[field_name]_150.sav' and turns them into 
;    /home/kstory/lps12/runlists/lt/runlist_lt_[field_name].txt

;
; MODIFICATION HISTORY:
;  04/20/2012: (KTS) Created
;;;

;...................................................................
; Make runlists for CR to coadd sims
PRO make_xspec2_runlists

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
; convert .sav LT runlists into .txt files
PRO sav2txt_lt_runlist
         stopit=stopit
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()
idx_list = [3,4,5]
nfields = n_elements(idx_list)

; dirs
bdir = '/home/kstory/lps12/runlists/lt/'

; Write out runlists into .txt file
for ii=0, nfields-1 do begin

    idx = idx_list[ii]
    fst = f[idx]
    list_txt = bdir+'runlist_lt_'+fst.name+'.txt' ; DEBUGGING

    ; get the .sav runlist
    list_sav = bdir+'list_'+fst.name+'_150.sav'
    restore, list_sav
    ndates = n_elements(lead_list)

    lead_dates = extract_date_from_filename(lead_list)
    trail_dates = extract_date_from_filename(trail_list)

    get_lun, lun1
    openw, lun1, list_txt
    for jj=0, ndates-1 do begin
        printf, lun1, lead_dates[jj], ' ', trail_dates[jj]
    endfor
    close, lun1
    free_lun,lun1
endfor

stop
END



