;;;
; NAME: script_0302
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  03/02/2012: (KTS) Created
;;;


;...................................................................
; Make nmiss arrays for all fields
pro sss
compile_opt IDL2, HIDDEN

maps = get_lps12_runlist(1, /obs_maps)
nmaps = n_elements(maps)
min_dec = fltarr(nmaps)

radec0 = [352.5, -55.0]
pix2ang_proj5, [960, 960], radec0, 1.0, ra, dec

print, "Nmaps = ", nmaps
for ii=0, nmaps-1 do begin
    d = expand_fits_struct(read_spt_fits(maps[ii]))

    wh = where(d.map.map ne 0, nwh)
    min_dec[ii] = min( abs(dec[wh]) )
    print, "ii = ", ii, ", ", min_dec[ii]
endfor

save, min_dec, filename='min_dec_0302.sav'
stop
end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Make nmiss arrays
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;...................................................................
; Make nmiss arrays for all fields
pro make_nmsiss_0302;, field_idx
compile_opt IDL2, HIDDEN

for ii=4, 19 do begin
    make_nmiss_array, ii
endfor
end

;...................................................................
; Make a test nmiss array for each obs of ra23h30dec-55_2008
pro make_nmiss_all, field_idx, plot=plot
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
        tv_spt_map, nmiss[*,200:*], /norms, scale=1, /forcesize
        com = ''
        read, com, prompt='mapnum = ' + strtrim(string(ii))
    endif

endfor

; Save nmiss
outname = outdir+'nmiss_all_'+field_name+'.sav'
print, "Save file here: ", outname
save, nmiss_all, filename=outname
end


;...................................................................
; Study the result
pro study_nmiss_all
compile_opt IDL2, HIDDEN

print, 'Restore sav file'
restore, '/data/kstory/projects/lps12/masks/nmiss_sav/nmiss_all_ra23h30dec-55_2008.sav'

sss = size(nmiss_all)
nmaps = sss[1]

print, 'Loop over files'
for ii=0, nmaps-1 do begin
    tv_spt_map, reform(nmiss_all[ii,*,300:*]), /norms, scale=1, /forcesize, winnum=2
    com = ''
    read, com, prompt = 'mapnum = ' + strtrim(string(ii),2)
    if (com eq 5) then return
endfor

end


;...................................................................
; new strategy for merge_and_coadd_lt
pro mac
compile_opt IDL2, HIDDEN

; get lead_list
restore, '/data/kstory/projects/lps12/runlists/lt/list_ra23h30dec-55_2008_150.sav'
lead_date  = extract_date_from_filename(lead_list)
trail_date = extract_date_from_filename(trail_list)
nobs = n_elements(lead_list)
group_list = intarr(nobs)

; get jitter group list
js=read_jitterlist_config(fieldname='ra23h30dec-55')
x = js.date

for ii=0, nobs-1 do begin
    wh = where(x eq lead_date[ii], nwh)
    group_list[ii] = js[wh].group
endfor

; Make arrays
keep = intarr(1,2)
keep_date  = strarr(1,2)
keep_group = intarr(1,2)

residual_g1 = intarr(1)
residual_g2 = intarr(1)
    
;------------------------
; Make good pairs first
;------------------------
for ii=0, 7 do begin

    ; Groups to be matched
    g1 = ii
    g2 = ii+8

    ; Find where the obs are in each group
    wh1 = where(group_list eq g1, nwh1)
    wh2 = where(group_list eq g2, nwh2)

    ; Find which one is bigger
    if nwh1 le nwh2 then npair = nwh1 else npair = nwh2

    tmp = strarr(npair, 2)

    ; loop over the number of good pairs in this group matching
    for pp=0, npair-1 do begin
        tmp[pp,0] = wh1[pp]
        tmp[pp,1] = wh2[pp]
    endfor

    ; Array of indicies in lead_list
    keep = [keep, tmp]
    ;keep_date = [keep_date, lead_date[keep]]
    ;keep_group = [keep_date, group_list[keep]]

    ; Figure out the residual, if any
    if npair lt nwh1 then begin
        residual_g1 = [ residual_g1, wh1[npair:*] ]
    endif
    if npair lt nwh2 then begin
        residual_g2 = [ residual_g2, wh2[npair:*] ]
    endif
endfor

; Remove first elements
keep       = keep[1:*,*]

residual_g1 = residual_g1[1:*,*]
residual_g2 = residual_g2[1:*,*]


;------------------------
; Deal with the residuals
;------------------------
lo = [residual_g1, residual_g2]
nlo = n_elements(lo[*,0])
nlo2 = floor(nlo / 2)
l1 = lo[0:nlo2-1]
l2 = lo[nlo2:nlo2*2-1]

tmp = [ [l1], [l2] ]
keep = [keep, tmp]

keep_date  = lead_date[keep]
keep_group = group_list[keep]


; Add trail info
final_dates = [ [lead_date[keep]], [trail_date[keep]] ]


;------------------------
; sav and write out runlist
;------------------------
save, keep, keep_date, keep_group, final_dates, filename='/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.sav'

new_runlist = '/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.txt'

nxspec = n_elements(final_dates[*,0])

get_lun, lun1
openw, lun1, new_runlist
for ii=0, nxspec-1 do begin
    tmp = final_dates[ii,0] + ' ' + final_dates[ii,1] + ' ' + final_dates[ii,2] + ' ' + final_dates[ii,3]
    printf, lun1, tmp
endfor
close, lun1
free_lun,lun1

stop

end


;...................................................................
; new strategy for merge_and_coadd_lt
pro check_mac
compile_opt IDL2, HIDDEN

; get the runlist
restore, '/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.sav'
nmaps = n_elements(final_maps[*,0])

tmpx = reform(final_dates[*,0])
tmpy = reform(final_dates[*,1])
ndates = n_elements(tmpx)
xx = fltarr(ndates)
yy = fltarr(ndates)

for ii=0, ndates-1 do begin
    xx[ii] = date_toolkit(tmpx[ii], 'mjd')
    yy[ii] = date_toolkit(tmpy[ii], 'mjd')
endfor

xx = xx - min(xx)
yy = yy - min(yy)

plot, xx, yy, psym=4
stop

end
