;;;
; NAME: calculate_total_area_lps12.pro
; PURPOSE:
;   Calculate the total effective area of this analysis
;
; NOTES:
;
; MODIFICATION HISTORY:
;  05/23/2012: (KTS) Created
;;;

;...................................................................
; Check overlap in apod masks between fields
pro calculate_total_area_lps12
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
f = lps12_fieldstruct()
reso_arcmin = 1.0

idx_list = [0,1,indgen(17)+3]
nfield = n_elements(idx_list)

; make big map
radec0_big = [20, -65]
npixels_big = [floor(960*7L), floor(960*4L)]
bigmap = fltarr(npixels_big[0], npixels_big[1])

for ii=0, nfield-1 do begin
    idx = idx_list[ii]
    print, "loop ",idx, ii

    ; get apod mask and field info
    mask = get_lps12_mask(idx)
    ;restore, '/home/kstory/lps12/masks/0522_run03/masks_50mJy/mask_'+f[idx].name+'.sav'
    radec0  = [ f[idx].ra0, f[idx].dec0 ]
    npx0 = (get_lps12_fieldinfo(idx)).npix 

    ; convert to ra,dec
    pix2ang_proj5, npx0, radec0, reso_arcmin, ra0, dec0


    ;------------------------
    ; Inset into big map
    ;------------------------

    npix_sm = long(npx0[0]) * long(npx0[1])
    for j=0, npix_sm-1 do begin
        if (mask[j] ne 0) then begin
            ang2pix_proj5, ra0[j], dec0[j], npixels_big, radec0_big, 1.0, ipix
            bigmap[ipix] += mask[j]
        endif
    endfor

endfor
;;;


;------------------------
; Plots
;------------------------

; Plot 1: apod masks
;tv_spt_map, bigmap[1000:*,300:*], /norms, scale=0.6, /forcesize, winnum=1, title='5h-7h, dec-65 to -50'
tv_spt_map, bigmap, /norms, winnum=1, title='masks for all maps'

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
;err= tvread(/png, filename=figname, /nodialog)

print, "Total area in deg^2 = ", total(bigmap) *(1/60.)^2.

stop
END

