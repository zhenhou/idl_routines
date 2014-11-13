;;;
; NAME: script13_0318
;
; NOTES:
;  1) re-process 1e5 sim
;;;


PRO save_plots;,tag=tag
;tag = '0318_1e4_realistic' & yr=[-20,20]
;tag = '0318_1e4_sn1' & yr=[.5, 1.5]
tag = '0317_1e5_realistic' & yr=[-20,20]

print, tag
;restore, 'sav_files/sim_sigbias_map_'+tag+'_stub.sav'
restore, 'sav_files/sim_sigbias_map_'+tag+'.sav'
s = sim_st
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

plot_sim, sim_st, sim2, tag=tag,pn=pn
err = tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+tag+'_tsplot')
stop

plot_bias, sim_st, sim0,sim1,sim2,s0,s1,s2,whmask,tag=tag,yr=yr
err = tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+tag+'_tbias')
stop

plot_pbias,sim_st,s0,s1,s2,tag=tag
wset, 5
err = tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+tag+'_presid')
wset, 0
err = tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+tag+'_pbias')
stop

END





PRO scratch, npix=npix

info = get_lps12_fieldinfo(6)
fname = 'ra4h10dec-50'
restore, '/home/kstory/lps12/masks/apod/apod_'+fname+'_60_0.0500_30.sav'
mask = apod
mask_row = mask[*,info.npix[1]/2]
; expand or collapse mask if necessary
if (npix gt info.npix[0]) then begin ; make mask bigger
    mm = fltarr(npix) + 1.
    mm[0:info.npix[0]/2-1] = mask_row[0:info.npix[0]/2-1]
    mm[npix-info.npix[0]/2:npix-1] = mask_row[info.npix[0]/2:info.npix[0]-1]
    mask_row = mm
endif
if (npix lt info.npix[0]) then begin ; make mask smaller    
    wh = where(mask_row lt 0.9, nwh)
    if (info.npix[0] - nwh ) lt 200 then begin
        print, '*** too few pixels: ', npix, ' .  Stopping'
        stop
    endif
    mm = fltarr(npix)
    mm[0:npix/2-1] = mask_row[0:npix/2-1]
    mm[npix/2:npix-1] = mask_row[info.npix[0]-npix/2:info.npix[0]-1]
    mask_row = mm
endif


stop
END
