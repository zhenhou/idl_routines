;;;
; NAME: script_0616
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) sim jacks
;
; MODIFICATION HISTORY:
;  06/16/2012: (KTS) Created
;;;


; Check runlist_randcut_95_seed1
PRO sss
f = lps12_fieldstruct()
rdir = '/home/kstory/lps12/runlists/'
rdir2 = '/home/kstory/lps12/runlists_0612/'
for i=0, 19 do begin
    print, i, ', Field ' + f[i].name
    spawn, 'less '+rdir+'runlist_randcut_95_seed1_'+f[i].name+'.txt | wc'
    spawn, 'less '+rdir+'runlist_lps12_'+f[i].name+'.txt | wc'
    ;spawn, 'less '+rdir2+'runlist_randcut_95_seed1_'+f[i].name+'.txt | wc'
endfor
END

; check randcut_95_seed1_az goodfiles
PRO sss2
f = lps12_fieldstruct()
rdir = '/home/kstory/lps12/runlists/'
gdir = '/home/kstory/lps12/jacks/jackgoodfiles/'

for i=0, 19 do begin

    rlist = rdir+'runlist_randcut_95_seed1_'+f[i].name+'.txt'
    readcol,/silent,rlist,dates_randcut,format='a'

    ; read good-files
    glist1 = gdir+'goodfiles_'+f[i].name+'_randcut_95_seed1_az_lowrms.txt'
    glist2 = gdir+'goodfiles_'+f[i].name+'_randcut_95_seed1_az_highrms.txt'
    readcol,/silent,glist1,dates_g1,format='a'
    readcol,/silent,glist1,dates_g2,format='a'
    dates_g = [dates_g1, dates_g2]
    ng = n_elements(dates_g)
    
    whmatch = intarr(ng)
    for j=0, ng-1 do begin
        whmatch[j] = where(dates_randcut eq dates_g[j])
    endfor
    
    wh = where(whmatch eq -1, nwh)

    ; print results
    print, i, 'field '+f[i].name
    print, nwh
    print, ng, n_elements(dates_randcut)
endfor

END



;;;;;;;;;;;;;;;;;;;;;;
; azmrs jack
;;;;;;;;;;;;;;;;;;;;;;

PRO jack_azrms
; azrms ; Note, run this as run_06, but run_07 is linked to 06.
for idx=0, 19 do begin
    print, 'jack azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms', run='06'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

;;;;;;;;;;;;;;;;;;;;;;
; randcut azmrs jack
;;;;;;;;;;;;;;;;;;;;;;

PRO jack_randcut_95_seed1_azrms
for idx=0, 19 do begin
    print, 'jack randcut_95_seed1_azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'randcut_95_seed1_azrms', run='06'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_randcut_95_seed2_azrms
for idx=0, 19 do begin
    print, 'jack randcut_95_seed2_azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'randcut_95_seed2_azrms', run='06'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_randcut_95_seed3_azrms
for idx=0, 19 do begin
    print, 'jack randcut_95_seed3_azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'randcut_95_seed3_azrms', run='06'
    print, 'That jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END



;;;;;;;;;;;;;;;;;;;;;;
; Jack_sims, per type
;;;;;;;;;;;;;;;;;;;;;;

PRO jack_sims_azrms
; coadd azrms
for idx=0, 19 do begin
     print, 'coadd jack sims azrms, field ', idx
    coadd_lps12_jack_sims, idx, 'azrms', run='07'
endfor
; jack azrms
for idx=0, 19 do begin
    print, 'jack sims azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'azrms', run='07', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_lr
; coadd lr
for idx=0, 19 do begin
     print, 'coadd jack sims lr, field ', idx
    coadd_lps12_jack_sims, idx, 'lr', run='07'
endfor
; jack lr
for idx=0, 19 do begin
    print, 'jack sims lr, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'lr', run='07', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_12
; coadd 12
for idx=0, 19 do begin
     print, 'coadd jack sims 12, field ', idx
    coadd_lps12_jack_sims, idx, '12', run='07'
endfor
; jack 12
for idx=0, 19 do begin
    print, 'jack sims 12, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, '12', run='07', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_tweight
; coadd tweight
for idx=0, 19 do begin
     print, 'coadd jack sims tweight, field ', idx
    coadd_lps12_jack_sims, idx, 'tweight', run='07'
endfor
; jack tweight
for idx=0, 19 do begin
    print, 'jack sims tweight, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'tweight', run='07', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_moon
; coadd moon
for idx=0, 19 do begin
     print, 'coadd jack sims moon, field ', idx
    coadd_lps12_jack_sims, idx, 'moon', run='07'
endfor
; jack moon
for idx=0, 19 do begin
    print, 'jack sims moon, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'moon', run='07', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


PRO jack_sims_sun
; coadd sun
idx_list = [4, 6, 11, 17, 18, 19]

for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'coadd jack sims sun, field ', idx
    coadd_lps12_jack_sims, idx, 'sun', run='07'
endfor
; jack sun
for ii=0, n_elements(idx_list)-1 do begin
    idx = idx_list[ii]
    print, 'jack sims sun, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'sun', /sim, run='07'
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


PRO jack_sims_randcut_95_seed1_azrms
; coadd randcut_95_seed1_azrms
for idx=0, 19 do begin
     print, 'coadd jack sims randcut_95_seed1_azrms, field ', idx
    coadd_lps12_jack_sims, idx, 'randcut_95_seed1_azrms', run='07'
endfor
; jack randcut_95_seed1_azrms
for idx=0, 19 do begin
    print, 'jack sims randcut_95_seed1_azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'randcut_95_seed1_azrms', run='07', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_randcut_95_seed2_azrms
; coadd randcut_95_seed2_azrms
for idx=0, 19 do begin
     print, 'coadd jack sims randcut_95_seed2_azrms, field ', idx
    coadd_lps12_jack_sims, idx, 'randcut_95_seed2_azrms', run='07'
endfor
; jack randcut_95_seed2_azrms
for idx=0, 19 do begin
    print, 'jack sims randcut_95_seed2_azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'randcut_95_seed2_azrms', run='07', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO jack_sims_randcut_95_seed3_azrms
; coadd randcut_95_seed3_azrms
for idx=0, 19 do begin
     print, 'coadd jack sims randcut_95_seed3_azrms, field ', idx
    coadd_lps12_jack_sims, idx, 'randcut_95_seed3_azrms', run='07'
endfor
; jack randcut_95_seed3_azrms
for idx=0, 19 do begin
    print, 'jack sims randcut_95_seed3_azrms, field ', idx
    t0 = systime(0, /seconds)
    lps12_jack, idx, 'randcut_95_seed3_azrms', run='07', /sim
    print, 'That sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END








;;;;;;;;;;;;;;;;;;;;;;
; Processing scripts
;;;;;;;;;;;;;;;;;;;;;;
PRO jacks1
;jack_azrms
jack_sims_azrms
jack_sims_lr
END

PRO jacks2
jack_sims_12
jack_sims_tweight
jack_sims_moon
jack_sims_sun
END

PRO jacks3
jack_randcut_95_seed1_azrms
jack_sims_randcut_95_seed1_azrms
END

PRO jacks4
jack_randcut_95_seed2_azrms
jack_sims_randcut_95_seed2_azrms
END

PRO jacks5
jack_randcut_95_seed3_azrms
jack_sims_randcut_95_seed3_azrms
END

; PRO jack_sims_1
;jack_sims_azrms
; jack_sims_lr
; jack_sims_12
; jack_sims_tweight
; jack_sims_moon
; jack_sims_sun
;jack_sims_randcut_95_seed1_azrms
;jack_sims_randcut_95_seed2_azrms
;jack_sims_randcut_95_seed3_azrms
; END

; PRO other_jacks
;jack_azrms
;jack_randcut_95_seed1_azrms
;jack_randcut_95_seed2_azrms
;jack_randcut_95_seed3_azrms
; END


PRO coadd_sims_run_07
coadd_all_lps12_sims
coadd_all_lps12_sims, /jacklr
END


;;;;;;;;;;;;;;;;;;;;;;
; Make azrms jackgoodfiles
;;;;;;;;;;;;;;;;;;;;;;

PRO make_azrms_goodfiles
f = lps12_fieldstruct()

for idx=0, 19 do begin
    print, 'Make lists for field ', idx, ': ' + f[idx].name
    make_azrms_jack_list, idx
    make_azrms_jack_list, idx, /preaz
endfor
END


;;;;;;;;;;;;;;;;;;;;;;
; End2End, run_07
;;;;;;;;;;;;;;;;;;;;;;

PRO end2end_run_07
tbegin = systime(0, /seconds)
for idx=0, 19 do begin
    print, 'run end2end_07, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_07, idx, /resume, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
tend = systime(0, /seconds)
print, 'The whole thing took ', (tend - tbegin)/60., ' minutes.'
END
