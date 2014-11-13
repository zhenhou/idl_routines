;;;
; NAME: script_0306
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) plot_overlap: plot the overlap between apod masks
;
; MODIFICATION HISTORY:
;  03/02/2012: (KTS) Created
;;;



;...................................................................
; Check overlap in apod masks between
; ra23h30dec-55, ra22h30dec-55
pro plot_overlap
compile_opt IDL2, HIDDEN

field_arr = lps12_fieldstruct()

reso_arcmin = 1.0

;restore, '/home/kstory/lps12/scripts/apod_test1_0305.sav'

ii = 0 ; ra5h30dec-55
;ii = 1 ; ra23h30dec-55
fst1 = field_arr[ii]
field_name1 = fst1.name
 npixels1 = get_lps12_map_npix(ii)
radec0_1  = [ fst1.ra0, fst1.dec0 ]
apod1 = 1
;make_apod_mask_0306, ii, overlap=[0.65, 0.65, 0.48, 0.48], ret_apod=apod1 ; BEST GUESS
;make_apod_mask_0306, ii, overlap=[0., 0., 0., 0.], ret_apod=apod1
;make_apod_mask_0306, ii, overlap=[-1, -1, -1, -1], ret_apod=apod1
restore, '/home/kstory/lps12/masks/apod_0305/apod_'+field_name1+'_60_0.0500_40.sav'
apod1 = apod

ii = 11 ; ra6h30dec-55
;ii = 14 ; ra22h30dec-55
fst2 = field_arr[ii]
field_name2 = fst2.name
npixels2 = get_lps12_map_npix(ii)
radec0_2  = [ fst2.ra0, fst2.dec0 ]
apod2 = 1
;make_apod_mask_0306, ii, overlap=[-1, 0., 0., 0.], ret_apod=apod2
;make_apod_mask_0306, ii, overlap=[-1, -1, -1, -1], ret_apod=apod2
restore, '/home/kstory/lps12/masks/apod_0305/apod_ra6h30dec-55_60_0.0500_40.sav'
apod2 = apod


ii = 16 ; ra6hdec-62.5
fst3 = field_arr[ii]
field_name3 = fst3.name
npixels3 = get_lps12_map_npix(ii)
radec0_3  = [ fst3.ra0, fst3.dec0 ]
apod3 = 1
;make_apod_mask_0306, ii, overlap=[-1, 0., -1, 0.], ret_apod=apod3
;make_apod_mask_0306, ii, overlap=[-1, -1, -1, -1], ret_apod=apod3
restore, '/home/kstory/lps12/masks/apod_0305/apod_ra6hdec-62.5_60_0.0500_40.sav'
apod3 = apod

; convert to ra,dec
pix2ang_proj5, npixels1, radec0_1, reso_arcmin, ra1, dec1
pix2ang_proj5, npixels2, radec0_2, reso_arcmin, ra2, dec2
pix2ang_proj5, npixels3, radec0_3, reso_arcmin, ra3, dec3


npixels = [960*1.7, 1200]
bigmap = fltarr(npixels[0], npixels[1])
;radec0 = [345, -55]
radec0 = [90, -57]

print, "loop 1"
npix_sm = 960L*960L
for ii=0, npix_sm-1 do begin
    if (apod1[ii] ne 0) then begin
        ang2pix_proj5, ra1[ii], dec1[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod1[ii]
    endif

    if (apod2[ii] ne 0) then begin
        ang2pix_proj5, ra2[ii], dec2[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod2[ii]
    endif
endfor

; ra6h30dec-62.5
print, "loop 2"
npix_sm = long(npixels3[0]) * long(npixels3[1])
for ii=0, npix_sm-1 do begin
    if (apod3[ii] ne 0) then begin
        ang2pix_proj5, ra3[ii], dec3[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod3[ii]
    endif
endfor


; cut out pixels that equal 2:
;bigmap[ where(bigmap eq 2) ] = 1

ss = bigmap*0
wh = where(bigmap gt 1)
ss[wh] = 1
;tv_spt_map, ss[*,150:*], /norms, scale=0.6, /forcesize, winnum=2, title='show where above 1'

ss0 = bigmap*0
wh = where(bigmap le 0)
ss0[wh] = 1
tv_spt_map, ss0[*,150:*], /norms, scale=0.6, /forcesize, winnum=2, title='show where equal to 0'

tv_spt_map, bigmap[*,150:*], /norms, scale=0.6, /forcesize, winnum=1, title='5h-7h, dec-65 to -50'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_0.5overlap_0306'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_Nooverlap_0306'
figname = '/home/kstory/lps12/masks/figs_0305/apod_1deg_0306'
;err= tvread(/png, filename=figname, /nodialog)

save, apod1, apod2, apod3, bigmap, ss, filename='apod_test1_0305.sav'
stop
end



pro test_overlap

restore, '/home/kstory/lps12/scripts/apod_test1_0305.sav'

npixels = [960*1.7, 1200]
bigmap = fltarr(npixels[0], npixels[1])
;radec0 = [345, -55]
radec0 = [90, -60]

field_arr = lps12_fieldstruct()
ii = 16 ; ra6hdec-62.5
fst3 = field_arr[ii]
field_name3 = fst3.name
npixels3 = get_lps12_map_npix(ii)
radec0_3  = [ fst3.ra0, fst3.dec0 ]
pix2ang_proj5, npixels3, radec0_3, 1.0, ra3, dec3
npixels3 = get_lps12_map_npix(16)

npix_sm = long(npixels3[0]) * long(npixels3[1])
for ii=0L, npix_sm-1 do begin
    if (apod3[ii] ne 0) then begin
        ang2pix_proj5, ra3[ii], dec3[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod3[ii]
    endif
endfor


; cut out pixels that equal 2:
;bigmap[ where(bigmap eq 2) ] = 1

ss = bigmap*0
wh = where(bigmap gt 1)
ss[wh] = 1
tv_spt_map, ss, /norms, scale=0.75, /forcesize, winnum=2

tv_spt_map, bigmap, /norms, scale=0.75, /forcesize, winnum=1

stop
end
