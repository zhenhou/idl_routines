;;;
; NAME: make_wtrms_coadds
; PURPOSE:
;   Make coadds from runlists
;
; INPUTS:
;    field_name,     name of field
;
; NOTES:
;   requires .sav files to be made with make_wtrms_runlists.pro
;
; MODIFICATION HISTORY:
;  08/20/2011: (KTS) Created
;  08/22/2011: (KTS) Make field_name an input
;;;

;...................................................................
; Main function, make coadds for wtrms cuts
;
pro make_wtrms_coadds, field_name, cut_idx_list
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'
field_idx = where(all_fields eq field_name)
;test_fields = [2] ; ra21hdec-60
;cut_list = [0, 6, 5, 8, 10, 9, 7, 11]

for cc=0, n_elements(cut_idx_list)-1 do begin
    make_coadds_single_field, field_idx, cut_idx_list[cc]
endfor

end

;...................................................................
; Make coadds for one field
;
pro make_coadds_single_field, field_idx, cut_idx
compile_opt IDL2, HIDDEN

restore, 'lps12_fieldnames.sav'
field_name = all_fields[field_idx]
field_dir = all_map_dirs[field_idx]

; Get cut set names
cutSetNames = cut_set_names()

runlist_name = '/data/kstory/projects/lps12/runlists/'+field_name+'/runlist_wtrms_'+cutSetNames[cut_idx]+'_'+field_name+'.sav'

print, "Reading runlist: " + runlist_name
restore, runlist_name

; Scramble and split the dates, then make a coadd of each half.
ndates = n_elements(dates150)
dates_scr = scramble(dates150, seed=17)
dates_first_half = dates_scr[0:ndates/2-1]
dates_second_half = dates_scr[ndates/2:ndates-1]

output_coadd_fname1 = 'coadd_wtrms_'+cutSetNames[cut_idx]+'_'+field_name+'_first_half.fits'
output_coadd_fname2 = 'coadd_wtrms_'+cutSetNames[cut_idx]+'_'+field_name+'_second_half.fits'

; Make some noise
print, "Making coadd for field " + field_name + ", cut " + cutSetNames[cut_idx]
print, "   output coadd: "+output_coadd_fname1
print, "                 "+output_coadd_fname2

coadd_fits_maps,dates=dates_first_half, dir=field_dir, $
  band=150,fileout='/data/kstory/projects/lps12/runlists/coadds/'+field_name+'/'+output_coadd_fname1
coadd_fits_maps,dates=dates_second_half, dir=field_dir, $
  band=150,fileout='/data/kstory/projects/lps12/runlists/coadds/'+field_name+'/'+output_coadd_fname2

end

; ;If field is Lead-Trail, you already have DMAPS
; if (lead_trail[field_idx]) then begin
;     ; Make some noise
;     output_coadd_fname = 'coadd_wtrms_'+cutSetNames[cut_idx]+'_'+field_name+'.fits'
;     print, "Making coadd for field " + field_name + ", cut " + cutSetNames[cut_idx]
;     print, "   output coadd: "+output_coadd_fname

;     coadd_fits_maps,dates=dates150, dir=field_dir, $
;       band=150,fileout='/data/kstory/projects/lps12/runlists/coadds/'+field_name+'/'+output_coadd_fname

