;;;
; NAME: plots_s12_alens
; PURPOSE:
;   Make alens plot
;
; NOTES:
; 1) Test SV averaging
;
; MODIFICATION HISTORY:
;  08/09/2012: (KTS) Created
;;;


;;;;;;;;;;;;;;;;;;;
; WMAP only, for testing
PRO p1_alens_w

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

; S12
dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/'
files_s12 = file_search(dir_s12+'c24_lcdm_alens_camb_w7s12*.txt')
pname_s12 = dir_s12+'c24_lcdm_alens_camb_w7s12.paramnames'

; WMAP
dir_w = '/data23/hou/lps12/paramfits/chains_final/lcdm_alens_camb_w7/chains/'
files_w = file_search(dir_w+'lcdm_alens_camb_w7*.txt')
pname_w = dir_w+'lcdm_alens_camb_w7.paramnames'

; K11
dir_k11   = '/data17/rkeisler/ps09/chains_20110202/'
files_k11 = file_search(dir_k11+'chain_k10_alens_1.txt')
pname_k11 = dir_k11+'chain_k10_alens.paramnames'


;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1700 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.
thickspt = 1
lfont = 0
ls = [0,2]

;*******************************
; colors
color_this_work = pscolors[8] ; blue
color_wmap = pscolors[2] ; red
color_k11 = pscolors[5] ; red

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[-2, 5]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, -3.0, 0.45, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 1.8, -0.16, charsize=1.5, charthick=3, font=0, '!3A!DL!X!N'

oplot,[0,0],[0,0],color=3

; S12
plot_like1dname,files_s12,pname_s12,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=ls[0], color=color_this_work, /oplot
xyouts,2,0.9,'SPT+WMAP',charsize=csize, color=color_this_work,font=lfont

; WMAP
i=1
plot_like1dname,files_w,pname_w,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=ls[i], color=color_wmap, /oplot
xyouts,2,0.8,'WMAP',charsize=csize, color=color_wmap,font=lfont

; K11
plot_like1dname,files_k11,pname_k11,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=3, color=color_k11, /oplot
xyouts,2,0.7,'K11',charsize=csize, color=color_K11,font=lfont

; vertical line
oplot, [1.,1.], [0,1], linestyle=1, thick=3

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/alens_w_0810.pdf'
stop
END



;;;;;;;;;;;;;;;;;;;
; WMAP only, for testing
PRO p1_alens2_w

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

; S12
dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/'
files_s12 = file_search(dir_s12+'c24_lcdm_alens_camb_w7s12*.txt')
pname_s12 = dir_s12+'c24_lcdm_alens_camb_w7s12.paramnames'

; WMAP
dir_w = '/data23/hou/lps12/paramfits/chains_final/lcdm_alens_camb_w7/chains/'
files_w = file_search(dir_w+'lcdm_alens_camb_w7*.txt')
pname_w = dir_w+'lcdm_alens_camb_w7.paramnames'

; K11
dir_k11   = '/data17/rkeisler/ps09/chains_20110202/'
files_k11 = file_search(dir_k11+'chain_k10_alens_1.txt')
pname_k11 = dir_k11+'chain_k10_alens.paramnames'


;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1700 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.
thickspt = 1
lfont = 0
ls = [0,2]

;*******************************
; colors
color_this_work = pscolors[8] ; blue
color_wmap = pscolors[2] ; red
color_k11 = pscolors[5] ; red

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[-0.2, 2.2]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, -0.55, 0.45, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 1., -0.16, charsize=1.5, charthick=3, font=0, '!3A!DL!X!N'

oplot,[0,0],[0,0],color=3

; S12
plot_like1dname,files_s12,pname_s12,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=ls[0], color=color_this_work, /oplot
xyouts,1.32,0.55,'SPT+WMAP',charsize=csize, color=color_this_work,font=lfont

; WMAP
i=1
plot_like1dname,files_w,pname_w,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=ls[i], color=color_wmap, /oplot
xyouts,1.32,0.45,'WMAP',charsize=csize, color=color_wmap,font=lfont

; K11
plot_like1dname,files_k11,pname_k11,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=3, color=color_k11, /oplot
xyouts,1.32,0.35,'K11',charsize=csize, color=color_K11,font=lfont

; vertical line
oplot, [1.,1.], [0,1], linestyle=1, thick=3

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/alens2_w_0810.pdf'
stop
END





;;;;;;;;;;;;;;;;;
; Figure: alens
;;;;;;;;;;;;;;;;;
PRO p1_alens
; Setup from CR idl_startup
WINDOW, /PIXMAP & WDELETE
DEVICE, BYPASS_TRANSLATION=0
device, retain=2
device, true_color=24
device, decomp=1
!p.charsize=1.
!p.color = !black
!p.background = !white
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]


;----------------------
; Get the data

; S12
dir = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/'
files = file_search(dir+'c24_lcdm_alens_camb_w7s12*.txt')
pname=dir+'c24_lcdm_alens_camb_w7s12.paramnames'
;n=n_elements(files)

; WMAP
;;;;


;*******************************
; colors
color_this_work = 60 ; blue
color_wmap = !red

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait
loadct,39,ncolors=12

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1700 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.
thickspt = 1
lfont = 0
ls = [0,2]



;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[-1, 5]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, -2.1, 0.45, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 1.8, -0.16, charsize=1.5, charthick=3, font=0, '!3A!DL!X!N'

oplot,[0,0],[0,0],color=3

; S12
plot_like1dname,files,pname,'alens',subsamp=3,thick=7,linestyle=ls[0], color=color_this_work, /oplot
xyouts,2,0.8,'SPT+WMAP',charsize=csize*1.3, color=color_this_work,font=lfont

; WMAP
i=1 ;& plot_like1dname,files[i],pname[i],'alens',subsamp=3,thick=4,linestyle=ls[i], color=cc[i], /oplot
xyouts,2,0.67,'WMAP',charsize=csize*1.3, color=color_wmap,font=lfont

; vertical line
oplot, [1.,1.], [0,1], linestyle=2, thick=5

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/ns_r_0815.pdf'
stop
END



