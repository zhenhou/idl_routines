;;;
; NAME: make_randcut_runlists
; PURPOSE:
;   Make runlists with a random percentage of observations cut.
;
; NOTES:
;   1) This operates on PREAZ runlists.
;   2) Default is to keep 95% of the data
;   3) To make the numbers match the current runlists, nkeep must be even.
;
; MODIFICATION HISTORY:
;  06/15/2012: (KTS) Copied from script_0531.pro
;;;

; Make runlists with random 10% data cuts.
PRO make_randcut_runlists, seed=seed, keep_pct=keep_pct, make_goodfiles=make_goodfiles

; setup
if n_elements(keep_pct) eq 0 then keep_pct = 95
skeep = strtrim(string(keep_pct), 2)
if (n_elements(seed) eq 0) then seed = 1
sseed = strtrim(string(seed), 2)

; setup
f = lps12_fieldstruct()
runlist_dir  = '/home/kstory/lps12/runlists/'

; loop over fields
for idx=0, 19 do begin
    
    new_runlist = runlist_dir+'runlist_randcut_'+skeep+'_seed'+sseed+'_'+f[idx].name+'.txt'
    
    ; Read the first column of the preaz runlist
    preaz_runlist = runlist_dir+'runlist_preaz_'+f[idx].name+'.txt'
    readcol,/silent,preaz_runlist,dates,format='a'
    n_dates = n_elements(dates)

    ; Cut a random 10%
    seed2use = seed
    rand_idx = randperm(n_dates, SEED=seed2use)
    rand_dates = dates[rand_idx]

    ; Calculate the number of obs to keep.
    neven = n_dates - (n_dates mod 2) ; even number, because jack setdefs have to match up.
    nkeep = floor((keep_pct/100.) * neven)
    nkeep = nkeep - (nkeep mod 2)     ; even number, because jack setdefs have to match up.
    new_dates = rand_dates[0:nkeep]
    new_dates = new_dates[sort(new_dates)]
    
    ; write new runlist
    print, 'write out new runlist: ' + new_runlist
    get_lun, lun1
    openw, lun1, new_runlist
    for ii=0, nkeep-1 do begin
        printf, lun1, new_dates[ii]
    endfor
    close, lun1
    free_lun,lun1

    ; Also make jack goodfiles?
    if keyword_set(make_goodfiles) then begin
        print, 'MAKE_RANDCUT_RUNLISTS: also make jackgoodfiles.'
        seed2use = seed
        make_azrms_jack_randcut_list, idx, seed2use, keep_pct
    endif

endfor


END


