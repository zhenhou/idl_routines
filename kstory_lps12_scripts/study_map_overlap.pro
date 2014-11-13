;;;
; NAME: study_map_overlap
; PURPOSE:
;   Study map overlap to make apod masks
;
; NOTES:
;
; MODIFICATION HISTORY:
;  03/08/2012: (KTS) Created
;  05/23/2012: (KTS) re-use to check map-overlap
;;;

;...................................................................
; Check overlap in apod masks between fields
pro study_map_overlap
compile_opt IDL2, HIDDEN

;restore, '/home/kstory/lps12/scripts/apod_test1_0305.sav'

;------------------------
; Setup
;------------------------
f = lps12_fieldstruct()
reso_arcmin = 1.0

ii = 14
field_name1 = f[ii].name & npixels1 = (get_lps12_fieldinfo(ii)).npix & radec0_1  = [ f[ii].ra0, f[ii].dec0 ]
apod1 = 1 & make_apod_mask_lps12, ii, overlap=[0.75, 0.75, 0.48, 0.48], ret_apod=apod1
;apod1 = get_lps12_mask(ii)

ii = 15
field_name2 = f[ii].name & npixels2 = (get_lps12_fieldinfo(ii)).npix & radec0_2  = [ f[ii].ra0, f[ii].dec0 ]
apod2 = get_lps12_mask(ii)

ii = 5
field_name3 = f[ii].name & npixels3 = (get_lps12_fieldinfo(ii)).npix & radec0_3  = [ f[ii].ra0, f[ii].dec0 ]
apod3 = get_lps12_mask(ii)

ii = 12
field_name4 = f[ii].name & npixels4 = (get_lps12_fieldinfo(ii)).npix & radec0_4  = [ f[ii].ra0, f[ii].dec0 ]
apod4 = get_lps12_mask(ii)

ii = 3
field_name5 = f[ii].name & npixels5 = (get_lps12_fieldinfo(ii)).npix & radec0_5  = [ f[ii].ra0, f[ii].dec0 ]
apod5 = get_lps12_mask(ii)

; convert to ra,dec
pix2ang_proj5, npixels1, radec0_1, reso_arcmin, ra1, dec1
pix2ang_proj5, npixels2, radec0_2, reso_arcmin, ra2, dec2
pix2ang_proj5, npixels3, radec0_3, reso_arcmin, ra3, dec3
pix2ang_proj5, npixels4, radec0_4, reso_arcmin, ra4, dec4
pix2ang_proj5, npixels5, radec0_5, reso_arcmin, ra5, dec5


;------------------------
; Inset into big map
;------------------------

radec0 = [330, -52]
npixels = [floor(960*5), floor(960*2.5)]
bigmap = fltarr(npixels[0], npixels[1])

print, "loop 1"
npix_sm = long(npixels1[0]) * long(npixels1[1])
for ii=0, npix_sm-1 do begin
    if (apod1[ii] ne 0) then begin
        ang2pix_proj5, ra1[ii], dec1[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod1[ii]
    endif
endfor

print, "loop 2"
npix_sm = long(npixels2[0]) * long(npixels2[1])
for ii=0, npix_sm-1 do begin
    if (apod2[ii] ne 0) then begin
        ang2pix_proj5, ra2[ii], dec2[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod2[ii]
    endif
endfor

print, "loop 3"
npix_sm = long(npixels3[0]) * long(npixels3[1])
for ii=0, npix_sm-1 do begin
    if (apod3[ii] ne 0) then begin
        ang2pix_proj5, ra3[ii], dec3[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod3[ii]
    endif
endfor

print, "loop 4"
npix_sm = long(npixels4[0]) * long(npixels4[1])
for ii=0, npix_sm-1 do begin
    if (apod4[ii] ne 0) then begin
        ang2pix_proj5, ra4[ii], dec4[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod4[ii]
    endif
endfor

print, "loop 5"
npix_sm = long(npixels5[0]) * long(npixels5[1])
for ii=0, npix_sm-1 do begin
    if (apod5[ii] ne 0) then begin
        ang2pix_proj5, ra5[ii], dec5[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod5[ii]
    endif
endfor


;------------------------
; Plots
;------------------------

; Plot 1: show overlap greater than 1
ss = bigmap*0
wh = where(bigmap gt 1)
ss[wh] = 1
tv_spt_map, ss[1000:*,300:*], /norms, scale=0.6, /forcesize, winnum=2, title='show where above 1'

; Plot 2: show where equal to 0
; ss0 = bigmap*0
; wh = where(bigmap le 0)
; ss0[wh] = 1
; tv_spt_map, ss0[*,150:*], /norms, scale=0.6, /forcesize, winnum=2, title='show where equal to 0'

; Plot 3: apod masks
tv_spt_map, bigmap[1000:*,300:*], /norms, scale=0.6, /forcesize, winnum=1, title='5h-7h, dec-65 to -50'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_0.5overlap_0306'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_Nooverlap_0306'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_1deg_0306'
;err= tvread(/png, filename=figname, /nodialog)

; save output for easy recycling
;print, 'save current apod masks into apod_test1_0305.sav'
;save, apod1, apod2, apod3, apod4, apod5, bigmap, ss, filename='/home/kstory/lps12/scripts/apod_test1_0305.sav'
stop
end

