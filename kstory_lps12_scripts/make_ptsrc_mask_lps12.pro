;;;
; NAME: make_ptsrc_mask_lps12
; PURPOSE:
;   Make point source masks
;
; INPUTS
;   field_idx,         lps12_fieldstruct index
;   cut_type,          string, either '50mJy' or '6p4mJy', matches naming scheme of ptsrc_config files
;   config_dir,        directory for ptsrc_config files
;   savdir,            directory where output masks are saved
; 
; NOTES:
; 1) expects a ptsrc config file
;
; MODIFICATION HISTORY:
;  03/05/2012: (KTS) Created
;  03/21/2012: (KTS) finalize script
;  03/26/2012: (KTS) use get_lps12_fieldinfo instead of get_lps12_map_npix
;;;


;...................................................................
; make masks for all fields
PRO run_all_masks
for ii=0, 19 do begin
    make_ptsrc_mask_lps12, ii
endfor
END

;...................................................................
; Make point source mask for one field
PRO make_ptsrc_mask_lps12, field_idx, $
                           cut_type=cut_type, $
                           config_dir = config_dir, $
                           savdir = savdir, $
                           stopit=stopit
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------

if ~keyword_set(config_dir) then config_dir = '/home/kstory/lps12/ptsrc_lists/'
if ~keyword_set(savdir) then savdir = '/home/kstory/lps12/masks/ptsrc/'

; Get the field information
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

band = 150

; Parse command line arguments
if n_elements(cut_type) eq 0 then cut_type='50mJy'


ptsrc_file = config_dir + 'ptsrc_config_'+cut_type+'_'+field_name+'_*.txt'

; Read the ptsrc config file
print, 'MAKE_PTSRC_MASK: reading in ptsrc_config file='+ptsrc_file
readcol,ptsrc_file,ind_src,ra_src,dec_src,maskrad_src

proj = 5
reso = 1.0

npix = (get_lps12_fieldinfo(field_idx)).npix
nx = npix[0]
ny = npix[1]
savename=savdir + 'ptsrc_mask_'+field_name+'_'+cut_type+'_proj'+strtrim(floor(proj),2)+'.sav'

radec0 = [fst.ra0, fst.dec0]
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

END




