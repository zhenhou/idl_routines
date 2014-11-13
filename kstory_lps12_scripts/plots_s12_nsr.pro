;;;
; NAME: plots_s12_nsr.pro
; PURPOSE:
;   Make 2D contour plot of ns v.s. r
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  08/17/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;
; Use this one
PRO plots_s12_nsr

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

; WMAP
dir_w = '/data23/hou/lps12/paramfits/chains_final/lcdm_r_camb_w7/chains/'
files_w = file_search(dir_w+'lcdm_r_camb_w7_*.txt')
pname_w = dir_w+'lcdm_r_camb_w7.paramnames'

; CMB
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c26_lcdm_r_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c26_lcdm_r_camb_w7s12_*.txt')
pname_cmb = dir_cmb+'c26_lcdm_r_camb_w7s12.paramnames'

; CMB+BAO
dir_ext = '/data23/hou/lps12/paramfits/chains_final/c28_lcdm_r_camb_w7s12_BAO/chains/'
files_ext = file_search(dir_ext+'c28_lcdm_r_camb_w7s12_BAO_*.txt')
pname_ext = dir_ext+'c28_lcdm_r_camb_w7s12_BAO.paramnames'


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
xtxt=0.99 ; xyouts label
ytxt=0.45
dytxt = 0.05
bigspt = 1.
thickspt = 1
lfont = 0
ls = [0,2]

;*******************************
; colors
;colors = [!red,!green,!blue,!yellow,!black,!white,!skyblue,!purple]
cc = [pscolors[2], pscolors[5], pscolors[8], pscolors[6], pscolors[0], pscolors[1], pscolors[7], pscolors[10]]

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

nskip=5000
subsamp=3

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[0.90, 1.07]
!y.window=[y0,y0+ddy]
!y.crange=[0, 0.8]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ylog=0,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, 0.98, -0.1, alignment=0.5, charsize=1.5, font=0, charthick=3, '!3n!Ds!X!N'
xyouts, 0.87, 0.4, charsize=1.5, charthick=3, font=0, '!3r!X!N'

oplot,[0,0],[0,0],color=3

;oplot,[0,1]
plot_like2dname,files_w,pname_w,'ns','r',/overplot,icolors=[cc[0],cc[3]],sigma=[1,2];,smooth=0.15;,priorname=['czero_ksz'],priormin=[0],priormax=[100]
plot_like2dname,files_cmb,pname_cmb,'ns','r',/overplot,icolors=[cc[2],cc[1]],sigma=[1,2];,smooth=0.15;,priorname=['czero_ksz'],priormin=[0],priormax=[100]
plot_like2dname,files_ext,pname_ext,'ns','r',/overplot,icolors=[cc[7],cc[6]],sigma=[1,2];,smooth=0.15;,priorname=['czero_ksz'],priormin=[0],priormax=[100]
; plot_like1dname,files_cmb,pname_cmb,'ns',subsamp=subsamp,nskip=nskip,scale=scale,$
;   yl=0,thick=7,linestyle=3, color=color_cmb, /oplot

xyouts,0.91,0.7,'WMAP',charsize=csize*1.3, color=cc[0],font=lfont
xyouts,0.91,0.64,'CMB',charsize=csize*1.3, color=cc[2],font=lfont
xyouts,0.91,0.58,'CMB+BAO',charsize=csize*1.3, color=cc[7],font=lfont

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/nsr_0817.pdf'
stop
END

