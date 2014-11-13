;;;
; NAME: script_0217
; PURPOSE:
;   make point-source config files
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/17/2012: (KTS) Created
;;;

;--------------------------------
; function to call make_ptsrc_config_single_field
; on all fields
pro run_0217

for ii=0, 19 do begin
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

; Set the input values to compile_source_mask_list_from_tcfiles
radec0 = [fstruct.ra0, fstruct.dec0]
;radec0 = [82.978, -55.203] ;; from CR

pixel_mask = 0.
maskthresh_50mJy = 50 & flux_select_50mJy = 1
maskthresh_5snr = 5 & flux_select_5snr = 0 
ra_min  = radec0[0] - (fstruct.dx / 2.)
ra_max  = radec0[0] + (fstruct.dx / 2.)
dec_min = radec0[1] - (fstruct.dy / 2.)
dec_max = radec0[1] + (fstruct.dy / 2.)

mytime = date_toolkit(systime(/julian), 'file')
configfile_out_50mJy = '/home/kstory/lps12/ptsrc_lists/ptsrc_config_50mJy_'+field_name+'_'+mytime+'.txt'
configfile_out_5snr = '/home/kstory/lps12/ptsrc_lists/ptsrc_config_5snr_'+field_name+'_'+mytime+'.txt'


;;; print out arguments:
print, '******** field ', field_name, ', ', field_idx
print, '               ', fstruct.name
print, 'radec0 = ', radec0
print, 'pixel_mask = ', pixel_mask
print, 'maskthresh_50mJy = ', maskthresh_50mJy
print, 'flux_select_50mJy = ', flux_select_50mJy
print, 'maskthresh_5snr = ', maskthresh_5snr
print, 'flux_select_5snr = ', flux_select_5snr
print, 'ra_min = ', ra_min
print, 'ra_max = ', ra_max
print, 'dec_min = ', dec_min
print, 'dec_max = ', dec_max
print, ' - configfile_out_50mJy = ',configfile_out_50mJy
print, ' - configfile_out_5snr  = ',configfile_out_5snr
print, ''
;stop

;----------------------
; PS cut at 50 mJy
compile_source_mask_list_from_tcfiles, radec0, pixel_mask, maskthresh_50mJy, output_struct, $
                                      ra_min=ra_min, ra_max=ra_max, $
                                      dec_min=dec_min, dec_max=dec_max, $
                                      flux_select=flux_select_50mJy, $
                                      configfile_out=configfile_out_50mJy

;----------------------
; PS cut at 5 SNR
compile_source_mask_list_from_tcfiles, radec0, pixel_mask, maskthresh_5snr, output_struct, $
                                      ra_min=ra_min, ra_max=ra_max, $
                                      dec_min=dec_min, dec_max=dec_max, $
                                      flux_select=flux_select_5snr, $
                                      configfile_out=configfile_out_5snr

                                      

end

