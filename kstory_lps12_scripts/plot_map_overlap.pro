;;;
; NAME: plot_map_overlap.pro
; PURPOSE:
;   Plot all apod masks in 2500 deg^2 footprint to check overlap
;
; INPUTS:
;   apod,                     plot apod masks
;   mask,                     plot final masks
;
; NOTES:
; 1) must specify either apod or mask keyword
;
; MODIFICATION HISTORY:
;  03/13/2012: (KTS) Created
;  03/26/2012: (KTS) use get_lps12_fieldinfo instead of get_lps12_map_npix
;;;


;...................................................................
; make final apod*ptsrc_mask mask for one field
pro plot_map_overlap, apod=apod, mask=mask
compile_opt IDL2, HIDDEN

; output plot name
out_tag = '0313'

; Read command line args.  Note you must specify either apod_dir or mask_dir
if ( ~keyword_set(apod) and ~keyword_set(mask) ) then begin
    print, "PLOT_MAP_OVERLAP: you must specify either apod or mask.  Exit."
    return
endif

; plot apod masks
if keyword_set(apod) then begin
    mask_dir = '/home/kstory/lps12/masks/apod_0312/'

    ; define the file name
    pre  = 'apod_'
    post = '_60_0.0500_30.sav'
    
    my_title = 'Apodization masks'  ; plot title
endif

if keyword_set(mask) then begin
    mask_dir = '/data/kstory/projects/lps12/masks/masks_50mJy/'

    ; define the file name
    pre  = 'mask_'
    post = '.sav'
    
    my_title = 'Full masks'
endif
    

;------------------------
; Setup
;------------------------
field_arr = lps12_fieldstruct()
reso_arcmin = 1.0


;------------------------
; Make big map
;------------------------

radec_big = [20, -60]
npixels_big = [floor(960*8), floor(960*4)]
bigmap = fltarr(npixels_big[0], npixels_big[1])

pix2ang_proj5, npixels_big, radec_big, reso_arcmin, ra_big, dec_big ; for debugging

; Loop over all fields
for ii=0, 19 do begin

    if ii eq 2 then ii=3 ; skip ra23h30dec-55_2010

    fst = field_arr[ii] & field_name = fst.name & npixels = (get_lps12_fieldinfo(ii)).npix & radec0 = [ fst.ra0, fst.dec0 ] ; Get the apod masks
    print, "PLOT_MAP_OVERLAP: process field ", field_name

    restore, mask_dir + pre + field_name + post
    if keyword_set(apod) then mask = apod
    if keyword_set(mask) then mask = mask
    
    pix2ang_proj5, npixels, radec0, reso_arcmin, ra, dec
    
    ; insert mask into bigmap
    print, "loop ", ii
    npix_sm = long(npixels[0]) * long(npixels[1])
    for jj=0, npix_sm-1 do begin
        if (mask[jj] ne 0) then begin
            ang2pix_proj5, ra[jj], dec[jj], npixels_big, radec_big, 1.0, ipix
            bigmap[ipix] += mask[jj]
        endif
    endfor
endfor

;------------------------
; Plots
;------------------------

; Plot 1: show overlap greater than 1
ss = bigmap*0
wh = where(bigmap gt 1)
ss[wh] = 1
fname='/data/kstory/projects/lps12/figs/ss_'+out_tag+'.ps'
print, 'Save file: ', fname
tv_spt_map, ss, /norms, scale=1.0, /forcesize, winnum=1, title='2500 deg^2, show where above 1', psfile=fname

; Plot 2: apod masks
;save plot
fname='/data/kstory/projects/lps12/figs/bigmap_'+out_tag+'.ps'
print, 'Save file: ', fname
tv_spt_map, bigmap, /norms, scale=1.0, /forcesize, winnum=2, title='2500 deg^2 masks', psfile=fname

stop
end

