;;;
; NAME: script_0427
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) run end2end on 10, 11, 16
;
; MODIFICATION HISTORY:
;  04/27/2012: (KTS) Created
;;;

;...................................................................
; run end2end
PRO end_6101116 ; done
try_end2end, 6, run='02', /resume
try_end2end, 10, run='02', /resume
try_end2end, 11, run='02', /resume
try_end2end, 16, run='02', /resume
END

PRO make_noise2
make_noise_psd_lps12, 10
make_noise_psd_lps12, 11
make_noise_psd_lps12, 16
END

; Process a field
PRO process_field, idx
coadd_maps_0420, idx
make_mask_lps12, field_idx=idx
make_coupling_kernel, idx
;make_noise_psd_lps12, idx
;try_end2end, idx ??
;make_twod_tfs, idx ??
;make_twod_kweights_lps12, idx ??
END

PRO process_0427
print, "********* process field 10 ****************"
process_field, 10

print, "********* process field 11 ****************"
process_field, 11

print, "********* process field 16 ****************"
process_field, 16
END

PRO maps_8
; make the maps
map_script = '/data/kstory/projects/lps12/maps/20120420/make_maps_ra2h30dec-50_150.txt'
log_file   = '/data/kstory/projects/lps12/maps/20120420/run0_ra2h30dec-50.out'
spawn, '/home/kstory/sptcdevel/mapping/mpibatch.x '+map_script+' > & '+log_file

END

;...................................................................
; process field 7
PRO process_7

; make the maps
map_script = '/data/kstory/projects/lps12/maps/20120420/make_maps_ra0h50dec-50_150.txt'
log_file   = '/data/kstory/projects/lps12/maps/20120420/run0_ra0h50dec-50.out'
;spawn, 'mpirun -n 4 /home/kstory/sptcdevel/mapping/mpibatch.x '+map_script+' & > '+log_file
spawn, '/home/kstory/sptcdevel/mapping/mpibatch.x '+map_script+' > & '+log_file

; coadd the maps
coadd_maps_0420, 7

; make mask
make_mask_lps12, field_idx=7

; make coupling kernel
make_coupling_kernel, 7

END

;...................................................................
; Check noise psds
PRO check_noise

restore, '/home/kstory/lps12/noise_psds/noise_psd_ra5h30dec-55_2008.sav'





END

