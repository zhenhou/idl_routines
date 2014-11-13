;;;
; NAME: script13_0419.pro
; PURPOSE:
;   Make a fields plot
;
; NOTES:
;;;

;;;;;;;;;;;;;;;;;;;
; Make fields overplot
;;;;;;;;;;;;;;;;;;;
PRO plots_s12_fields, type, huge=huge
;;;
; Make a plot of the spt fields for observation
; Copied from TC: spt:/home/tcrawfor/spt_analysis/temp/mkplot_spt_obsreg_2010version_batch
; 09/23/2010 (KTS): Created
;;;

type = 'coadd90_icon'

if keyword_set(huge) then begin
    pxsize = 1600
    pysize = 1600
endif
pix2ang_proj5,[6144,3072],[22.5,-62.3],1.,rapix,decpix
apod_temp=fltarr(6144,3072) 
apod_temp[where(rapix ge 20.*15. and decpix ge -65. and decpix le -40.)]=1.
apod_temp[where(rapix le 7.*15. and decpix ge -65. and decpix le -40.)]=1
apod_tempsm=apod_temp & smooth_map,apod_tempsm,1.,30.

; Set default plot options
use_boarder_outline = 1
cb = 0 ; black,     ;cb = 255 ; white
uthick=5
lthick=1.5


case type of
    'coadd90_icon': begin
        sc1=read_spt_fits('/data17/rkeisler/weakly_filtered/super_coadd_90_old_runlist.fits')
        mapsm=sc1.map.map*apod_tempsm
        smooth_map,mapsm,1.,8.
        
        nside = 8192
        filename = strcompress('spt2500deg_'+string(fix(nside))+'.eps',/rem)
        
        fac_rebin = max([round(float(nside)/2048.),1])
        dusttemp = proj_to_healp(mapsm,5,1.,[22.5,-62.3],nside,rebin=fac_rebin[0])

        masktemp = proj_to_healp(apod_temp,5,1.,[22.5,-62.3],nside,rebin=fac_rebin[0])
        wh0 = where(masktemp lt 0.02,n0)
        if n0 gt 0 then dusttemp[wh0] = dusttemp[0]
        whbad = where(dusttemp eq dusttemp[0],nbad)
        if nbad gt 0 then dusttemp[whbad] = 1d30
        
        myrot = [202.5,240,180]
        ;orthview,/online,dusttemp*1e6,grat=[15,10],colt=39,min=-400.,max=400.,eul_mat=emtemp,rot=myrot,/half,title=' ',/nobar
        orthview,/online,dusttemp*1e6,colt=39,min=-400.,max=400.,eul_mat=emtemp,rot=myrot,/half,title=' ',/nobar,pxsize=1400;,pysize=pysize
        ;f_annotate = '/home/kstory/lps12/scripts/plotting/annotate_fields_lps12_coadd90_icon.dat'
        plot_field_boundaries = 0
    endcase

endcase

stop

;legstr=['4000 sq. deg.','Published fields (Lueker 09','Finished fields','Future fields (2010, 2011)']
;legstr=['2008 fields','2009 fields','2010 fields','2011 fields']
legstr=['           ','           ','           ','           ']
loadct,39

;DEVICE, SET_FONT = 'Helvetica*Bold' 
;DEVICE, SET_FONT = 'Palatino*Bold*Italic' 

; color names, 190=yellow, 140=green, 250=red, 40=purple, 80=blue
;color_full=190
color_2008=140
color_2009=250
color_2010=80
color_2011=190
;legend,legstr,line=0,psym=0,thick=2,color=[color_2008,color_2009,color_2010,color_2011],position=[-1.01,1.25],charsize=1.9,charthick=1.5,font=-1

; Plot the full 4000 field
; ra_full=[20.,7.] & dec_full=[-65.,-30.]
; oplot_rectangle,ra_full,dec_full,emtemp,proj='ORTH',/half,color=color_full,thick=2


if plot_field_boundaries then begin
    
; Plot black outlines
    if use_boarder_outline then begin
;2008
    oplot_rectangle,[5.,6.],[-60.,-50.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[23.,0.],[-60.,-50.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
;2009
    oplot_rectangle,[20.,22.],[-65.,-55.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[2.,5.],[-65.,-55.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[20.,22.],[-55.,-45.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
;2010
    oplot_rectangle,[0.,2.],[-65.,-55.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[0.,25.]/15.,[-55.,-45.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[25.,50.]/15.,[-55.,-45.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[50.,75.]/15.,[-55.,-45.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[5.,6.],[-50.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[6.,7.],[-60.,-50.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
;2011
    oplot_rectangle,[20.,22.],[-45.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[22.,0.],[-65.,-60.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[22.,23.],[-60.,-50.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[22.,0.],[-50.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[0.,2.],[-45.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[2.,5.],[-45.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[5.,7.],[-65.,-60.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop
    oplot_rectangle,[6.,7.],[-50.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=cb,thick=uthick,/crop

    endif

; Plot colored outlines
;2008
    oplot_rectangle,[5.,6.],[-60.,-50.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2008,thick=lthick,/crop
    oplot_rectangle,[23.,0.],[-60.,-50.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2008,thick=lthick,/crop
;2009
    oplot_rectangle,[20.,22.],[-65.,-55.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2009,thick=lthick,/crop
    oplot_rectangle,[2.,5.],[-65.,-55.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2009,thick=lthick,/crop
    oplot_rectangle,[20.,22.],[-55.,-45.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2009,thick=lthick,/crop
;2010
    oplot_rectangle,[0.,2.],[-65.,-55.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
    oplot_rectangle,[0.,25.]/15.,[-55.,-45.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
    oplot_rectangle,[25.,50.]/15.,[-55.,-45.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
    oplot_rectangle,[50.,75.]/15.,[-55.,-45.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
    oplot_rectangle,[5.,6.],[-50.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
    oplot_rectangle,[6.,7.],[-60.,-50.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2010,thick=lthick,/crop
;2011
    oplot_rectangle,[20.,22.],[-45.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
    oplot_rectangle,[22.,0.],[-65.,-60.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
    oplot_rectangle,[22.,23.],[-60.,-50.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
    oplot_rectangle,[22.,0.],[-50.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
    oplot_rectangle,[0.,2.],[-45.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
    oplot_rectangle,[2.,5.],[-45.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
    oplot_rectangle,[5.,7.],[-65.,-60.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
    oplot_rectangle,[6.,7.],[-50.,-40.],emtemp,rot=myrot,proj='ORTH',/half,color=color_2011,thick=lthick,/crop
endif

;annotate,load_file='/home/kstory//projects/field_obs/preview_2010_08/annotate_preview_0923.dat'
;annotate,load_file='/home/tcrawfor/SPTData/annotate_obsreg.dat'
;annotate,load_file=f_annotate

stop
END



;;;;;;;;;;;;;;;;;;;
; 
; Plot ss0809 / ss1011
; 
;;;;;;;;;;;;;;;;;;;
PRO plot_year
restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120830_154612_kweight_0809.sav'
s08 = {dl:dl_all, er:diag_nobeam}

restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120830_210134_kweight_1011.sav'
s10 = {dl:dl_all, er:diag_nobeam}

restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_ApJ2012.sav'
sall = {dl:dl_all, er:diag_nobeam}


vec=indgen(47)+9

!p.multi=[0,2,1]
window, 0, xsize=800,ysize=300
plot, l[vec],(s08.dl/sall.dl)[vec], yr=[0.85,1.15], $
  title='Bandpower Comparison',xtitle='!12l!X!N',ytitle='Ratio of D!D!12l!X!'
oplot, l[vec],(s08.dl/sall.dl)[vec], color=!red
oplot, l[vec],(s10.dl/sall.dl)[vec], color=!blue
legend,['s0809/sall','s1011/sall'],colors=[!red,!blue],pos=[510,1.19],linestyle=[0,0]

plot, l[vec],(s08.dl/sall.dl)[vec] / mean((s08.dl/sall.dl)[vec]), yr=[0.85,1.15],$
    title='Re-scaled in power',xtitle='!12l!X!N',ytitle='Ratio of D!D!12l!X!'
oplot, l[vec],(s08.dl/sall.dl)[vec] / mean((s08.dl/sall.dl)[vec]), color=!red
oplot, l[vec],(s10.dl/sall.dl)[vec] / mean((s10.dl/sall.dl)[vec]), color=!blue

!p.multi=0
;err=tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/pow_ratio_year')

stop
END

