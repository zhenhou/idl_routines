;;;
; NAME: merge_and_coadd_lt_maps_ra23h30dec-55_2008.pro
; PURPOSE:
;   Merge lead-trail pairs, and coadd multiple pairs in the
;   ra23h30dec-55_2008 field.  The result will be the maps that we
;   intend to use for xspec
;
; INPUTS: None - all hard-coded
;
; OUTPUTS:
;   maps in '/data/kstory/projects/lps12/maps/20120224/ra23h30dec-55_2008_lt/'
;
; NOTES:
; 1) Current map directory is /data/kstory/projects/lps12/maps/20120224/
; 2) Resulting coadded maps are stored in /data/kstory/projects/lps12/maps/20120224/*field*_lt/
;
; MODIFICATION HISTORY:
;  03/01/2012: (KTS) Created
;  03/07/2012: (KTS) Change output dir to 20120305
;  03/10/2012: (KTS) Change name to '4Xpair_'
;;;


;...................................................................
; Merge and coadd maps
pro merge_and_coadd_now
compile_opt IDL2, HIDDEN

; get the runlist
restore, '/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.sav'
nmaps = n_elements(final_maps[*,0])

map_dir = '/data/kstory/projects/lps12/maps/20120305/'

print, "MERGE_AND_COADD_LT_MAPS_RA23H30DEC-55_2008, make coadds."
for ii=0, nmaps-1 do begin

    ; Get a set of indicies to use in lead_set, trail_set
    print, "*** Coadding ", ii, " / ", nmaps

    ; name the output file with the date of the first lead obs in this coadd
    coadd_name = map_dir+'ra23h30dec-55_2008_lt/4Xpair_map_ra23h30dec-55_150_'+keep_date[ii,0]+'.fits'

    filelist = reform(final_maps[ii,*])
    coadd_fits_maps, filelist=filelist, fileout=coadd_name
endfor

end


