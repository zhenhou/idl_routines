;;;
; NAME: make_ptsrc_config_0305
; PURPOSE:
;   make point-source config files
;
; NOTES:
; 1) with expanded coverage, add_rk_ptsrc_to_config is no longer necessary
;
; MODIFICATION HISTORY:
;  02/17/2012: (KTS) Created
;  02/20/2012: (KTS) added procedure add_rk_ptsrc_to_config
;  03/01/2012: (KTS) call add_rk_ptsrc_to_config automatically
;  03/05/2012: (KTS) Add 5 deg on all sides to be sure we have covered
;                    all observed area.
;  03/05/2012: (KTS) Remove add_rk_ptsrc_to_config
;;;

;--------------------------------
; function to call make_ptsrc_config_single_field
; on all fields
pro make_ptsrc_config_0305

for ii=0, 20 do begin
    make_ptsrc_config_single_field, ii
endfor
end


;--------------------------------
; Make the ptsrc_config file for a single field
; indexed by lps12_fieldstruct
;
pro make_ptsrc_config_single_field, field_idx

; get the field properties
field_arr = lps12_fieldstruct()

field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name

tmp = read_spt_fields()
wh = where( strcmp(tmp.name, field_dir_name) eq 1, nwh)
if(nwh eq 0) then begin
    print, "field name miss-match"
    return
endif

fstruct = tmp[wh]

;------------
; Set the input values to compile_source_mask_list_from_tcfiles
;------------
radec0 = [fstruct.ra0, fstruct.dec0]
pixel_mask = 0.
ra_min  = radec0[0] - ( (fstruct.dx / 2.) + 5.)
ra_max  = radec0[0] + ( (fstruct.dx / 2.) + 5.)
dec_min = radec0[1] - ( (fstruct.dy / 2.) + 5.)
dec_max = radec0[1] + ( (fstruct.dy / 2.) + 5.)
proj = 5
mytime = '20120305';date_toolkit(systime(/julian), 'file')


;------------
; Mask 1: mask all sources above 50mJy
maskrad_arcmin_50mJy  = 5
maskthresh_50mJy      = 50
thresh4maskrad_50mJy  = 50
configfile_out_50mJy  = '/home/kstory/lps12/ptsrc_lists/ptsrc_config_50mJy_'+field_name+'_'+mytime+'.txt'

;------------
; Mask 2: mask all sources with flux in [6.4mJy, 50mJy] with 2 arcmin mask,
;         and all sources above 50mJy with 5 arcmin mask

maskrad_arcmin_6p4mJy  = [5, 2]
maskthresh_6p4mJy      = [50, 6.4]
thresh4maskrad_6p4mJy  = [50, 6.4]
configfile_out_6p4mJy  = '/home/kstory/lps12/ptsrc_lists/ptsrc_config_6p4mJy_'+field_name+'_'+mytime+'.txt'


;;; print out arguments:
print, '******** field ', field_name, ', ', field_idx
; print, '               ', fstruct.name
; print, 'radec0 = ', radec0
; print, 'pixel_mask = ', pixel_mask
; print, 'maskthresh_50mJy = ', maskthresh_50mJy
; print, 'ra_min = ', ra_min
; print, 'ra_max = ', ra_max
; print, 'dec_min = ', dec_min
; print, 'dec_max = ', dec_max
; print, 'maskrad_arcmin = ', maskrad_arcmin
; print, 'proj = ', proj
print, ' - configfile_out_50mJy = ',configfile_out_50mJy
print, ' - configfile_out_6p4mJy  = ',configfile_out_6p4mJy
print, ''
;stop

;----------------------
; PS cut at 50 mJy
print, '--- Make config for 50 mJy cut ---'

compile_source_mask_list_from_tcfiles, radec0, pixel_mask, maskthresh_50mJy, output_struct, $
  ra_min=ra_min, ra_max=ra_max, $
  dec_min=dec_min, dec_max=dec_max, $
  maskrad_arcmin=maskrad_arcmin_50mJy, $
  bands4thresh=150, $
  thresh4maskrad=thresh4maskrad_50mJy, $
  proj=proj, $
  /flux_select, $
  configfile_out=configfile_out_50mJy

;----------------------
; PS cut at 6.4 mJy, 50 mJy
print, '--- Make config for split cut, [6.4-50, 50+] ---'
compile_source_mask_list_from_tcfiles, radec0, pixel_mask, maskthresh_6p4mJy, output_struct, $
  ra_min=ra_min, ra_max=ra_max, $
  dec_min=dec_min, dec_max=dec_max, $
  maskrad_arcmin=maskrad_arcmin_6p4mJy, $
  bands4thresh=[150, 150], $
  thresh4maskrad=thresh4maskrad_6p4mJy, $
  proj=proj, $
  /flux_select, $
  configfile_out=configfile_out_6p4mJy

end
