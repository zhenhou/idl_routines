;;;
; NAME: make_coadds_0820
; PURPOSE:
;   TEMP file to make coadds from runlists
;
; MODIFICATION HISTORY:
;  08/20/2011: (KTS) Created
;;;

;...................................................................
; Main function, call this
;
pro make_coadds_ra21hdec60_0820
compile_opt IDL2, HIDDEN

test_fields = [2] ; ra21hdec-60
cut_list = [0, 6, 5, 8, 10, 9, 7, 11]

for ii=0, n_elements(test_fields)-1 do begin
    for cc=0, n_elements(cut_list)-1 do begin
        make_coadds_for_field, test_fields[ii], cut_list[cc]
    endfor
endfor

end

;...................................................................
; Make coadds for one field
;
pro make_coadds_for_field, field_idx, cut_idx
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'
field_name = all_fields[field_idx]
field_dir = all_map_dirs[field_idx]

; Get cut set names
cutSetNames = cut_set_names()

restore, '/data/kstory/projects/lps12/runlists/ra21hdec-60/runlist_wtrms_default2_'+field_name+'.sav'
output_coadd_fname = 'coadd_wtrms_'+cutSetNames[cut_idx]+'_'+field_name+'.fits'

; Make some noise
print, "Making coadd for field ", field_name, ", cut ", cutSetNames[cut_idx]

;; If field is Lead-Trail, you already have DMAPS
if (lead_trail[field_idx]) then begin
    coadd_fits_maps,dates=dates150, dir=field_dir, $
      band=150,fileout='/data/kstory/projects/lps12/runlists/coadds/'+output_coadd_fname
endif else begin
    ; Not Lead-Trail
    print, 'Not Lead-Trail, IMPLEMENT THIS LATER'
endelse

    

end

