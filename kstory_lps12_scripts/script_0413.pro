;;;
; NAME: script_0413
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) plot_data; plot coadds, looking at point sources and power leaking into DMAP
; 2) plot_tcs; make plots for my webpage

; MODIFICATION HISTORY:
;  04/13/2012: (KTS) Created
;;;

;...................................................................
; 
PRO plot_data, idx
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()

coadd_dir = '/data/kstory/projects/lps12/maps/20120305/coadds/'
file = coadd_dir + 'coadd_'+f[idx].name+'_50mJy.fits'
d = krf(file)

; plot
mnmx_x = [400, 900]
mnmx_y = [300, 700]
;pix_mn = 100 & pix_mx = 639
while 1 do begin
    tv_spt_map, d.map.map[mnmx_x[0]:mnmx_x[1], mnmx_y[0]:mnmx_y[1]], title='map', scale=1.5, /forcesize, winnum=6
    fig_name = fig_dir+'map_'+name+'_0413'
    ;err = tvread(/png, filename=fig_name, /nodialog)
    pause

    tv_spt_map, d.dmap.map[mnmx_x[0]:mnmx_x[1], mnmx_y[0]:mnmx_y[1]], scale=1.5, /forcesize, winnum=6 
    fig_name = fig_dir+'map_'+name+'_0413'
    ;err = tvread(/png, filename=fig_name, /nodialog)
    pause
endwhile


stop
END

;...................................................................
; Time-constant plots
PRO plot_tcs
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()

fig_dir   = '/data/kstory/projects/lps12/figs/'
coadd_dir = '/data/kstory/projects/lps12/maps/20120305/coadds/'

; figure name stub
;stub = '_norm'
stub = '_zoom'

;--------------------
; plot good field
;--------------------
idx = 3
file = coadd_dir + 'coadd_'+f[idx].name+'_50mJy.fits'
d = krf(file)
name = f[idx].name

; get new map
;mnmx_x = [400, 900] & mnmx_y = [300, 700] ;; scale = 1
mnmx_x = [600, 900] & mnmx_y = [500, 700]  ;; scale = 2
mymap  = d.map.map[mnmx_x[0]:mnmx_x[1], mnmx_y[0]:mnmx_y[1]]
mydmap = d.dmap.map[mnmx_x[0]:mnmx_x[1], mnmx_y[0]:mnmx_y[1]]

; normalize by T_source
mx = max(mymap)
gsrange = [-0.1*mx, 0.1*mx]

;tv_spt_map, mymap, title='map, '+name, minval=gsrange[0], maxval=gsrange[1], scale=1, /forcesize, winnum=6
tv_spt_map, mymap, title='map, '+name, scale=2, /forcesize, winnum=6
;tv_spt_map, mymap, title='map, '+name, scale=1, /forcesize, winnum=6
fig_name = fig_dir+'map_'+name+stub+'_0413'
pause
err = tvread(/png, filename=fig_name, /nodialog)

;tv_spt_map, mydmap, title='dmap, '+name, scale=1, minval=gsrange[0], maxval=gsrange[1], /forcesize, winnum=6 
tv_spt_map, mydmap, title='dmap, '+name, scale=2, /forcesize, winnum=6 
;tv_spt_map, mydmap, title='dmap, '+name, scale=1, /forcesize, winnum=6 
fig_name = fig_dir+'dmap_'+name+stub+'_0413'
pause
err = tvread(/png, filename=fig_name, /nodialog)



;--------------------
; plot bad field
;--------------------
idx = 7
file = coadd_dir + 'coadd_'+f[idx].name+'_50mJy.fits'
d = krf(file)

; get new map
mnmx_x = [400, 900] & mnmx_y = [300, 700] ;; scale=1
mnmx_x = [600, 900] & mnmx_y = [500, 700] ;; scale=2
mymap  = d.map.map[mnmx_x[0]:mnmx_x[1], mnmx_y[0]:mnmx_y[1]]
mydmap = d.dmap.map[mnmx_x[0]:mnmx_x[1], mnmx_y[0]:mnmx_y[1]]

; normalize by T_source
mx = max(mymap)
gsrange = [-0.1*mx, 0.1*mx]

name = f[idx].name
;tv_spt_map, mymap, title='map, '+name, scale=1, minval=gsrange[0], maxval=gsrange[1], /forcesize, winnum=6
tv_spt_map, mymap, title='map, '+name, scale=2, /forcesize, winnum=6
;tv_spt_map, mymap, title='map, '+name, scale=1, /forcesize, winnum=6
fig_name = fig_dir+'map_'+name+stub+'_0413'
pause
err = tvread(/png, filename=fig_name, /nodialog)

;tv_spt_map, mydmap, title='dmap, '+name, scale=1, minval=gsrange[0], maxval=gsrange[1], /forcesize, winnum=6 
tv_spt_map, mydmap, title='dmap, '+name, scale=2, /forcesize, winnum=6 
;tv_spt_map, mydmap, title='dmap, '+name, scale=1, /forcesize, winnum=6 
fig_name = fig_dir+'dmap_'+name+stub+'_0413'
pause
err = tvread(/png, filename=fig_name, /nodialog)

;stop
END

;...................................................................
; Time-constant plots
PRO plot_ffts, idx
compile_opt IDL2, HIDDEN

f = lps12_fieldstruct()

coadd_dir = '/data/kstory/projects/lps12/maps/20120305/coadds/'
file = coadd_dir + 'coadd_'+f[idx].name+'_50mJy.fits'
d = krf(file)

mymap  = d.map.map
mydmap = d.dmap.map

; fft's
myfft  = fft(mymap)
mydfft = fft(mydmap)

; plot
stop

END

;...................................................................
; Time-constant plots
PRO sim_filter, idx
compile_opt IDL2, HIDDEN

sim_map = [0,0]

stop
END
