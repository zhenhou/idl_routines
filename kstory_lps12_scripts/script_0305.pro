;;;
; NAME: script_0305
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) test_nmiss_all: plot nmiss arrays for individual obs
; 2) check_nmiss: plot final nmiss (combined all obs) for all fields
; 3) make_apod_0305: make apodization masks, probably not used
; 4) test_ptsrc: plot ptsrc masks to test coverage
; 5) plot_ptsrc_ra5h: ptsrc masks to test coverage
; 6) plot_overlap: moved to script_0306
;
; MODIFICATION HISTORY:
;  03/05/2012: (KTS) Created
;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Make nmiss arrays
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;...................................................................
; Plot nmiss arrays for individual observations
pro test_nmiss_all, field_idx, plot=plot
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
outdir = '/data/kstory/projects/lps12/masks/nmiss_sav/'

;------------------------
; Get field information
;------------------------
;field_idx = 1
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

;------------------------
; get the map files
;------------------------
maps = get_lps12_runlist(field_idx, /xspec_maps)
nmaps = n_elements(maps)

;------------------------
; Get nmiss for this observation
;------------------------
info = expand_fits_struct(read_spt_fits(maps[0])) ; get map shape
sss = size(info.map.map)
nmiss_all = fltarr(nmaps, sss[1], sss[2])

print, "MAKE_NMISS_ARRAY, nmaps = ", nmaps
for ii=0, nmaps-1 do begin
;for ii=0, 9 do begin

    ; make some noise
    if (ii mod 50 eq 2) then print, "processing " + strtrim(string(ii),2) + "/" + strtrim(string(nmaps),2)

    d = expand_fits_struct(read_spt_fits(maps[ii]))
    weight = d.weight.map
    weight /=max(weight)
    nmiss     = 0*weight

    ; find pixels that have zero obs
    wh_0 = where(weight eq 0, nwh, complement=wh_1)

    nmiss[wh_0] = 1
    nmiss_all[ii,*,*] = nmiss

    if keyword_set(plot) then begin
        ;tv_spt_map, nmiss[*,200:*], /norms, scale=0.75, /forcesize, winnum=1
        tv_spt_map, nmiss, /norms, scale=0.75, /forcesize, winnum=1
        com = ''
        read, com, prompt='mapnum = ' + strtrim(string(ii))
    endif

endfor

; Save nmiss
outname = outdir+'nmiss_all_'+field_name+'.sav'
print, "Save file here: ", outname
;save, nmiss_all, filename=outname
end


;...................................................................
; Check nmiss savs
pro check_nmiss
compile_opt IDL2, HIDDEN

field_arr_ = lps12_fieldstruct()

for ii=0, 19 do begin
    if ii eq 3 then ii=4

    fst = field_arr_[ii]
    restore, '/data/kstory/projects/lps12/masks/nmiss_sav/nmiss_'+fst.name+'.sav'

    print, "nmiss for field " + fst.name
    tv_spt_map, reform(nmiss), /norms, scale=0.75, /forcesize, winnum=2
    com = ''
    read, com, prompt = 'mapnum = ' + strtrim(string(ii),2)
    if (com eq 5) then return
endfor

end


;...................................................................
; Make all apod masks
pro make_apod_0305
compile_opt IDL2, HIDDEN

for ii=0, 19 do begin
    if (ii eq 3) then ii=4 ; skip ra21hdec-60 for now

    make_apod_mask_0305, ii
endfor
end



;...................................................................
; Check ptsrc and apod masks
pro test_ptsrc;, field_idx
compile_opt IDL2, HIDDEN

for ii=0, 19 do begin

if ii eq 3 then ii=4
field_arr = lps12_fieldstruct()
fst = field_arr[ii]
field_name = fst.name

; get ptsrc mask
restore, '/home/kstory/lps12/masks/ptsrc_0305/ptsrc_mask_'+field_name+'_50mJy_proj5.sav'

; get apod mask
restore, '/home/kstory/lps12/masks/apod_0305/apod_'+field_name+'_60_0.0500_40.sav'

mask = apod*ptsrc_mask

tv_spt_map, mask, /norms, scale=0.75, /forcesize, winnum=1

figname = '/home/kstory/lps12/masks/figs_0305/fullMask_0305_'+field_name
err= tvread(/png, filename=figname, /nodialog)

endfor

end


;...................................................................
; Check ptsrc and apod masks
pro plot_ptsrc_ra5h
compile_opt IDL2, HIDDEN

ii = 0
field_arr = lps12_fieldstruct()
fst = field_arr[ii]
field_name = fst.name

; get ptsrc mask
restore, '/home/kstory/lps12/masks/ptsrc_0305/ptsrc_mask_'+field_name+'_6p4mJy_proj5.sav'

; get apod mask
restore, '/home/kstory/lps12/masks/apod_0305/apod_'+field_name+'_60_0.0500_40.sav'

mask = apod*ptsrc_mask
tv_spt_map, mask, /norms, scale=0.75, /forcesize, winnum=1
figname = '/home/kstory/lps12/masks/figs_0305/fullMask_6p4mJy_0305_'+field_name
;err= tvread(/png, filename=figname, /nodialog)
stop
end


;...................................................................
; Check overlap in apod masks between
; ra23h30dec-55, ra22h30dec-55
pro plot_overlap
compile_opt IDL2, HIDDEN

field_arr = lps12_fieldstruct()

reso_arcmin = 1.0

; ra23h30dec-55
ii = 1
fst1 = field_arr[ii]
field_name1 = fst1.name
npixels1 = get_lps12_map_npix(ii)
radec0_1  = [ fst1.ra0, fst1.dec0 ]
;restore, '/home/kstory/lps12/masks/ptsrc_0305/ptsrc_mask_'+field_name1+'_6p4mJy_proj5.sav'
;pts1 = ptsrc_mask
restore, '/home/kstory/lps12/masks/apod_0305/apod_'+field_name1+'_60_0.0500_40.sav'
apod1 = apod

; ra22h30dec-55
ii = 14
fst2 = field_arr[ii]
field_name2 = fst2.name
npixels2 = get_lps12_map_npix(ii)
radec0_2  = [ fst2.ra0, fst2.dec0 ]
;restore, '/home/kstory/lps12/masks/ptsrc_0305/ptsrc_mask_'+field_name2+'_6p4mJy_proj5.sav'
;pts2 = ptsrc_mask
restore, '/home/kstory/lps12/masks/apod_0305/apod_'+field_name2+'_60_0.0500_40.sav'
apod2 = apod

; convert to ra,dec
pix2ang_proj5, npixels1, radec0_1, reso_arcmin, ra1, dec1
pix2ang_proj5, npixels2, radec0_2, reso_arcmin, ra2, dec2


npixels = [960*2, 960]
bigmap = fltarr(npixels[0], npixels[1])
radec0 = [345, -55]

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

tv_spt_map, bigmap, /norms, scale=0.75, /forcesize, winnum=1
;figname = '/home/kstory/lps12/masks/figs_0305/blah'
;err= tvread(/png, filename=figname, /nodialog)
stop
end

