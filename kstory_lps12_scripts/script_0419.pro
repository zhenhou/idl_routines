;;;
; NAME: script_0419
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) check_idf_progress, checks how many idfs have been made for each field.
; 2) plot_0419, make plots of coadded map, dmap for ra5h30dec-45
; 3) make_mask_10, make apodization mask for ra5h30dec-45
; 4) make_kern_10, make apodization mask for ra5h30dec-45
;
; MODIFICATION HISTORY:
;  04/19/2012: (KTS) Created
;;;

;...................................................................
; Check progress of current idf job
PRO check_idf_progress, $
         stopit=stopit
compile_opt IDL2, HIDDEN

;--------------------
; setup
;--------------------
f = lps12_fieldstruct()


for ii=0, 19 do begin
;ii=11
;ii=14
    if ii eq 9 then ii =10

    fst = f[ii]
    idf_dir = '/home/rkeisler/lowellfits/'+fst.name+'/'

    rdates = get_lps12_runlist(ii, /obs_dates)
    nr = n_elements(rdates)

    files = file_search(idf_dir+'field_scan_150_*.fits')
    nf = n_elements(files)

    ; compare against runlist
    tmp = intarr(nr)
    for jj=0, nr-1 do begin
        tmp[jj] = file_test(idf_dir+'field_scan_150_'+rdates[jj]+'.fits')
    endfor

    whg = where(tmp ne 0, nwhg, complement=whb)

    gfiles = nwhg gt 0 ? files[ where(tmp ne 0)] : 0
    bfiles = whb[0] ne -1 ? files[ where(tmp eq 0)] : 0

    print, strtrim(string(ii),2), ' ' + fst.name + ': nf / nr = ', strtrim(string(nf),2) + '/'+ strtrim(string(nr),2)
    print, '          : ncomplete / nr = ', strtrim(string(nwhg),2) , '/', strtrim(string(nr),2)
    print, ' '

endfor
stop

END

;...................................................................
; Check progress of current idf job
PRO plot_0419, $
               stopit=stopit
compile_opt IDL2, HIDDEN

;--------------------
; setup
;--------------------
f = lps12_fieldstruct()
idx = 10
fst = f[idx]
d = krf('/home/kstory/lps12/maps/20120420/coadds/coadd_ra5h30dec-45_50mJy.fits')
;d = krf('/home/kstory/lps12/maps/20120305/coadds/coadd_ra5h30dec-45_50mJy.fits')

; overall plot
tv_spt_map, d.map.map[300:600, 450:700], winnum=4
tv_spt_map, d.dmap.map[300:600, 450:700], winnum=3

; zoomed in plot
tv_spt_map, d.dmap.map[380:400, 500:510], winnum=7
tv_spt_map, d.map.map[380:400, 500:510], winnum=8

; ptsrc cut
window, 9 & plot, d.map.map[380:400, 508], title='MAP, ptsrc cut'
window, 10 & plot, d.dmap.map[380:400, 508], title='DMAP, ptsrc cut'

; dipole #1
print, max(d.dmap.map[350:370, 572]) / max(d.map.map[350:370, 572])

; dipole #2
print, max(d.dmap.map[395:415, 560]) / max(d.map.map[395:415, 560])

stop

END


;...................................................................
; Make mask for ra5h30dec-45
PRO make_mask_10
compile_opt IDL2, HIDDEN
make_apod_mask_lps12, 10, coadd_dir = '/data/kstory/projects/lps12/maps/20120420/coadds/'
make_ptsrc_mask_lps12, 10
make_mask_lps12, field_idx=10
END

;...................................................................
; Make mask for ra5h30dec-45
PRO make_kern_10
compile_opt IDL2, HIDDEN
make_coupling_kernel, 10
END

