;;;
; NAME: make_runlist_xspec_ra23h30dec-55_2008
; PURPOSE:
;   Make runlists of maps that we plan to cross in the power spectrum
;   analysis for the ra23h30dec-55_2008 field
;
;   The algorithm is as follows:
;   1) Only consider obs that have matched LT pairs.
;   2) Sort everything by lead_list, which is sorted by date.
;   3) Split obs up by jitter group.  Make 'good pairs,' which are
;     defined as obs that have jitter groups differing by 8 (there are
;     16 total).  Make as many good pairs as possible
;   4) Order the remaining 38 observations by jitter number, and match
;     the first half to the second half (0-19, 1-20, 2-21, ... 18-37)
;
; INPUTS: None
;
; OUTPUTS:
;   - sav file in '/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.sav'
;   - txt file in '/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.txt'
;   - xspec runlist at '/data/kstory/projects/lps12/runlists/runlist_xspec_lps12_ra23h30dec-55_2008.txt'
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  03/02/2012: (KTS) Created
;;;

;...................................................................
;
pro make_runlist_xspec
compile_opt IDL2, HIDDEN

; get lead - trail pairs
restore, '/data/kstory/projects/lps12/runlists/lt/list_ra23h30dec-55_2008_150.sav'
lead_date  = extract_date_from_filename(lead_list)
trail_date = extract_date_from_filename(trail_list)
nobs = n_elements(lead_list)
group_list = intarr(nobs) ; jitter group

; get jitter group list
js=read_jitterlist_config(fieldname='ra23h30dec-55')
x = js.date

for ii=0, nobs-1 do begin
    wh = where(x eq lead_date[ii], nwh)
    group_list[ii] = js[wh].group
endfor

; Make output arrays
keep = intarr(1,2) ; indicies in lead_list
keep_date  = strarr(1,2) ; dates in lead_list
keep_group = intarr(1,2) ; jitter group number

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
final_maps  = [ [lead_list[keep]], [trail_list[keep]] ]

;------------------------
; sav and write out runlist with all 4 dates
;------------------------
save, keep, keep_date, keep_group, final_dates, final_maps, filename='/data/kstory/projects/lps12/runlists/lt/4Xpair_ra23h30dec-55_2008.sav'

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


;------------------------
; Write out xspec runlist
;------------------------
xspec_runlist = '/data/kstory/projects/lps12/runlists/runlist_xspec_lps12_ra23h30dec-55_2008.txt'
get_lun, lun1
openw, lun1, xspec_runlist
for ii=0, nxspec-1 do begin
    tmp = final_dates[ii,0]
    printf, lun1, tmp
endfor
close, lun1
free_lun,lun1


stop

end
