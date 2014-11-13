;;;
; NAME: script_0220
; PURPOSE:
;   make point source masks
;
; INPUTS:
;   field_idx,           index in lps12_fieldstruct()
;   cut_type,            either '50mJy' or '5snr', specific to my
;                        naming scheme
;
; NOTES:
; 1) Runs on config files made with make_ptsrc_config_0220
;
; MODIFICATION HISTORY:
;  02/20/2012: (KTS) Created from /home/rkeisler/ps09/make_ptsrc_mask.pro
;;;

; Make point source masks for all fields
pro run_0220
for ii=0, 19 do begin
    mk_ptsrc_mask_0220, ii, cut_type='50mJy'
    mk_ptsrc_mask_0220, ii, cut_type='5snr'
endfor
end

;pro mk_ptsrc_mask_0220, field_idx, rdisc=rdisc, rtaper=rtaper, cut_type=cut_type
pro mk_ptsrc_mask_0220, field_idx, cut_type=cut_type, stopit=stopit

; Get the field information
field_arr = lps12_fieldstruct()

field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name

; deal with special case of naming
if(field_name eq 'ra5h30dec-55_2008') then begin
    coadd_field_name = 'ra5h30dec-55'
endif else begin
    coadd_field_name = field_name
endelse


band = 150

; Parse command line arguments
;if n_elements(rdisc) eq 0 then rdisc=5.
;if n_elements(rtaper) eq 0 then rtaper=5.
if n_elements(cut_type) eq 0 then cut_type='50mJy'

;------------------------
; Setup
ptsrc_file='/home/kstory/lps12/ptsrc_lists/ptsrc_config_'+cut_type+'_'+field_name+'_fixed_*.txt'
coadd_dir = '/data17/rkeisler/beams_2010_2011/all_coadds/'

; Read the ptsrc config file
readcol,ptsrc_file,ind_src,ra_src,dec_src,maskrad_src

coaddname = coadd_dir + 'coadd_'+coadd_field_name+'.fits'
d=read_spt_fits(coaddname)
map = d.map.map
;proj = d.mapinfo.projection
proj = 5
nx = n_elements(map[*,0])
ny = n_elements(map[0,*])
savename='/home/kstory/lps12/masks/ptsrc_mask_'+field_name+'_'+cut_type+'_proj'+strtrim(floor(proj),2)+'.sav'

radec0 = [d.mapinfo.ra0, d.mapinfo.dec0]
reso = d.mapinfo.reso_arcmin
pix2ang_anyproj, [nx,ny], radec0, reso, ra_map, dec_map, proj=proj

nsrc = n_elements(ra_src)
print,'MAKE_PTSRC_MASK: field=',field_name
print,'MAKE_PTSRC_MASK: masking ',nsrc,' sources.'
print, 'MAKE_PTSRC_MASK: output file name =', savename

ptsrc_mask = fltarr(nx,ny) + 1.
for i=0,nsrc-1 do begin
    this_ra = ra_src[i]
    this_dec = dec_src[i]
    this_rdisc  = maskrad_src[i]*60.
    this_rtaper = this_rdisc
    dr = sphdist(this_ra, this_dec, ra_map, dec_map, /degrees)*60.

    tmp = disc_taper(dr, r_disc=this_rdisc, sigma_taper=this_rtaper)
    ptsrc_mask *= tmp
endfor

save,ptsrc_mask,filename=savename

if keyword_set(stopit) then stop

end



