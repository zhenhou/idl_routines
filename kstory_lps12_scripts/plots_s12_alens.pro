;;;
; NAME: plots_s12_alens
; PURPOSE:
;   Make alens plot
;
; NOTES:
; 1) Use p1_alens_w, others are obsolete
;
; MODIFICATION HISTORY:
;  08/09/2012: (KTS) Created
;  09/10/2012: (KTS) Update with 0828 chains
;;;

;;;;;;;;;;;;;;;;;
; Figure: alens
;;;;;;;;;;;;;;;;;
PRO p1_alens
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

; S12+WMAP
dir_cmb = cdir+'c54_lcdm_alens_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c54_lcdm_alens_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c54_lcdm_alens_camb_w7s12.paramnames'

; WMAP
dir_w = '/data23/hou/lps12/paramfits/chains_0717/lcdm_alens_camb_w7/chains/'
files_w = file_search(dir_w+'lcdm_alens_camb_w7*.txt')
pname_w = dir_w+'lcdm_alens_camb_w7.paramnames'

; S12-only
dir_s12   = cdir+'c52_lcdm_alens_camb_s12tau/chains/'
files_s12 = file_search(dir_s12+'c52_lcdm_alens_camb_s12tau*.txt')
pname_s12 = dir_s12+'c52_lcdm_alens_camb_s12tau.paramnames'

; K11
dir_k11   = '/data17/rkeisler/ps09/chains_20110202/'
files_k11 = file_search(dir_k11+'chain_k10_alens_1.txt')
pname_k11 = dir_k11+'chain_k10_alens.paramnames'


;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.0;1.3
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/alens'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1.48;1.8 ; xyouts label
ytxt=0.68
dytxt = 0.065
bigspt = 1.
thickspt = 1
lfont = 0
ls = [0,2,3]

;*******************************
; colors
color_cmb   = pscolors[0] ; black
color_s12   = pscolors[8] ; blue
color_wmap  = pscolors[2] ; !red

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
xyouts, -0.45, 0.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, '!3Likelihood'
xyouts, 0.95, -0.16, charsize=1.5, charthick=3, font=0, '!3A!DL!X!N'

oplot,[0,0],[0,0],color=3
;subsamp = 5

; CMB
plot_like1dname,files_cmb,pname_cmb,'alens',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[0], color=color_cmb, /oplot

; WMAP
plot_like1dname,files_w,pname_w,'alens',subsamp=5,nskip=1000,$
  thick=7,linestyle=ls[1], color=color_wmap, /oplot

; S12
plot_like1dname,files_s12,pname_s12,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=3, color=color_s12, /oplot

;xyouts,xtxt,ytxt,'SPT+WMAP7',charsize=csize, color=color_cmb,font=lfont
;oplot,[xtxt-0.2, xtxt-0.005],[ytxt+0.02,ytxt+0.02],thick=7,linestyle=ls[0],color=color_cmb
xyouts,xtxt,ytxt+dytxt-0.02,'SPT+',charsize=csize, color=color_cmb,font=lfont
xyouts,xtxt-0.03,ytxt,'WMAP7',charsize=csize, color=color_cmb,font=lfont
oplot,[xtxt-0.2, xtxt-0.005],[ytxt+dytxt,ytxt+dytxt],thick=7,linestyle=ls[0],color=color_cmb
xyouts,xtxt,ytxt-dytxt,'WMAP7',charsize=csize, color=color_wmap,font=lfont
oplot,[xtxt-0.2, xtxt-0.005],[ytxt-dytxt+0.02,ytxt-dytxt+0.02],thick=7,linestyle=ls[1],color=color_wmap
xyouts,xtxt,ytxt-2*dytxt,'SPT',charsize=csize, color=color_S12,font=lfont
oplot,[xtxt-0.2, xtxt-0.005],[ytxt-2*dytxt+0.02,ytxt-2*dytxt+0.02],thick=7,linestyle=ls[2],color=color_s12

; ; K11
; plot_like1dname,files_k11,pname_k11,'alens',subsamp=3,nskip=1000,$
;   thick=7,linestyle=3, color=color_s12, /oplot
; xyouts,xtxt,0.35,'K11',charsize=csize, color=color_s12,font=lfont
; oplot,[xtxt-0.3, xtxt-0.001],[0.35+0.02,0.35+0.02],thick=7,linestyle=ls[2],color=color_s12


; vertical line
oplot, [1.,1.], [0,1], linestyle=1, thick=3

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
;spawn, 'cp figs/tmp.pdf figs/alens.pdf'
print, 'output: ', filename_stub+'.pdf'
stop
END

