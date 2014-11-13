;;;
; NAME: plots_s12_ns2.pro
; PURPOSE:
;   Make ns, ns+models plot
;
; NOTES:
;
; MODIFICATION HISTORY:
;  08/22/2012: (KTS) created
;  09/10/2012: (KTS) Update with 0828 chains
;  10/04/2012: (KTS) 
;;;


;;;;;;;;;;;;;;;;;;;
; 
; 2d ns plot
;
;;;;;;;;;;;;;;;;;;;

PRO ns2

cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;;; LCDM
; WMAP
dir_w = cdir+'c1_lcdm_pico_w7/chains/'
files_w = file_search(dir_w+'c1_lcdm_pico_w7*.txt')
pname_w = dir_w+'c1_lcdm_pico_w7.paramnames'

; cmb
dir_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c2_lcdm_pico_w7s12*.txt')
pname_cmb = dir_cmb+'c2_lcdm_pico_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = cdir+'c4_lcdm_pico_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c4_lcdm_pico_w7s12_H0*.txt')
pname_cmb_h0 = dir_cmb_h0+'c4_lcdm_pico_w7s12_H0.paramnames'

; CMB+BAO 
dir_cmb_bao = cdir+'c3_lcdm_pico_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c3_lcdm_pico_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c3_lcdm_pico_w7s12_BAO.paramnames'

;----------------------
; external datasets

;;; cmb-only
dir_lcdm_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'
files_lcdm_cmb = file_search(dir_lcdm_cmb+'c2_lcdm_pico_w7s12*.txt')
pname_lcdm_cmb = dir_lcdm_cmb+'c2_lcdm_pico_w7s12.paramnames'

; neff_cmb
type = 'c12_lcdm_neff_pico_w7s12'
dir_neff_cmb = cdir+type+'/chains/'
files_neff_cmb = file_search(dir_neff_cmb+type+'_*.txt')
pname_neff_cmb = dir_neff_cmb+type+'.paramnames'

; mnu_cmb
type = 'c19_lcdm_mnu_pico_w7s12'
dir_mnu_cmb = cdir+type+'/chains/'
files_mnu_cmb = file_search(dir_mnu_cmb+type+'_*.txt')
pname_mnu_cmb = dir_mnu_cmb+type+'.paramnames'


;;; cmb_h0
; lcdm_h0
type = 'c4_lcdm_pico_w7s12_H0'
dir_lcdm_h0 = cdir+type+'/chains/'
files_lcdm_h0 = file_search(dir_lcdm_h0+type+'_*.txt')
pname_lcdm_h0 = dir_lcdm_h0+type+'.paramnames'

; neff_h0
type = 'c14_lcdm_neff_pico_w7s12_H0'
dir_neff_h0 = cdir+type+'/chains/'
files_neff_h0 = file_search(dir_neff_h0+type+'_*.txt')
pname_neff_h0 = dir_neff_h0+type+'.paramnames'

; mnu_h0
type = 'c21_lcdm_mnu_pico_w7s12_H0'
dir_mnu_h0 = cdir+type+'/chains/'
files_mnu_h0 = file_search(dir_mnu_h0+type+'_*.txt')
pname_mnu_h0 = dir_mnu_h0+type+'.paramnames'


;;; cmb_bao
; lcdm_bao
type = 'c3_lcdm_pico_w7s12_BAO'
dir_lcdm_bao = cdir+type+'/chains/'
files_lcdm_bao = file_search(dir_lcdm_bao+type+'_*.txt')
pname_lcdm_bao = dir_lcdm_bao+type+'.paramnames'

; neff_bao
type = 'c13_lcdm_neff_pico_w7s12_BAO'
dir_neff_bao = cdir+type+'/chains/'
files_neff_bao = file_search(dir_neff_bao+type+'_*.txt')
pname_neff_bao = dir_neff_bao+type+'.paramnames'

; mnu_bao
type = 'c20_lcdm_mnu_pico_w7s12_BAO'
dir_mnu_bao = cdir+type+'/chains/'
files_mnu_bao = file_search(dir_mnu_bao+type+'_*.txt')
pname_mnu_bao = dir_mnu_bao+type+'.paramnames'


;;; cmb_h0_bao
; lcdm_h0_bao
type = 'c5_lcdm_pico_w7s12_BAO_H0'
dir_lcdm_h0_bao = cdir+type+'/chains/'
files_lcdm_h0_bao = file_search(dir_lcdm_h0_bao+type+'_*.txt')
pname_lcdm_h0_bao = dir_lcdm_h0_bao+type+'.paramnames'

; neff_h0_bao
type = 'c15_lcdm_neff_pico_w7s12_BAO_H0'
dir_neff_h0_bao = cdir+type+'/chains/'
files_neff_h0_bao = file_search(dir_neff_h0_bao+type+'_*.txt')
pname_neff_h0_bao = dir_neff_h0_bao+type+'.paramnames'

; mnu_h0_bao
type = 'c22_lcdm_mnu_pico_w7s12_BAO_H0'
dir_mnu_h0_bao = cdir+type+'/chains/'
files_mnu_h0_bao = file_search(dir_mnu_h0_bao+type+'_*.txt')
pname_mnu_h0_bao = dir_mnu_h0_bao+type+'.paramnames'


;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 11
csize = 1.3
yb = 4.5

filename_stub = '~kstory/lps12/scripts/figs/ns2'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.10
ddx=.38
dx_off = .5
xoff=.1
y0=.16
ddy=.77

; Define the output plot and font sizes
!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot

;-----------------
; Plot 1
; colors; use setup_ps_plotting
color_w = pscolors[2]           ; red
color_cmb = pscolors[0]         ; black
color_cmb_h0 = pscolors[3]      ; orange
color_cmb_bao = pscolors[0]     ; black
color_cmb_h0_bao = pscolors[8] ; blue
; 5 - forest

ls = [2, 0, 4, 1, 3] ; w, cmb, cmb_h0, cmb_bao, cmb_h0_bao

thick=3 ; points
xtxt=0.993 ; xyouts label
ytxt=0.93
dytxt = 0.08
bigspt = 1.2
thickspt = 1
lfont = 0

!x.window=[x0,x0+ddx]
!y.window=[y0,y0+ddy]
!x.crange=[0.92, 1.05]          ; LCDM
!y.crange=[0, 1]

param = 'ns'
subsamp = 5
nskip=1000

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=5
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=5

axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=5 ,yl=0
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=5, yl=0

oplot, [1,1],[0,0.1], linestyle=0, thick=10

plot_like1dname,files_w,pname_w,param,subsamp=subsamp,nskip=nskip,scale=scale,$
  yl=1,thick=7,linestyle=ls[0], color=color_w, /oplot
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,$
  yl=1,thick=7,linestyle=ls[1], color=color_cmb, /oplot
plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,$
  yl=1,thick=7,linestyle=ls[2], color=color_cmb_h0, /oplot
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[3], color=color_cmb_bao, /oplot
plot_like1dname,files_lcdm_h0_bao,pname_lcdm_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[4], color=color_cmb_h0_bao, /oplot

xtitle='!3n!Ds!X!N'
ytitle='!3Likelihood!X!N'
xyouts, 0.904, 0.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, ytitle
xyouts, 0.975, -0.16, charsize=1.5, charthick=3, font=0, xtitle

xyouts,xtxt,ytxt,     'WMAP7',chars=bigspt,charthick=thickspt,font=lfont,color=color_w
oplot, [xtxt-0.01, xtxt-0.001], [ytxt+0.04,ytxt+0.04], thick=7,linestyle=ls[0],color=color_w
xyouts,xtxt,ytxt-dytxt*1, 'SPT+WMAP7',chars=bigspt, charthick=thickspt,font=lfont,color=color_cmb
oplot, [xtxt-0.01, xtxt-0.001], [(ytxt-dytxt*1)+0.04,(ytxt-dytxt*1)+0.04], thick=7,linestyle=ls[1],color=color_cmb
xyouts,xtxt,ytxt-dytxt*2,'SPT+WMAP7+H!D0!X!N',chars=bigspt, charthick=thickspt,font=lfont,color=color_cmb_h0
oplot, [xtxt-0.01, xtxt-0.001], [(ytxt-dytxt*2)+0.04,(ytxt-dytxt*2)+0.04], thick=7,linestyle=ls[2],color=color_cmb_h0
xyouts,xtxt,ytxt-dytxt*3,'SPT+WMAP7+BAO',chars=bigspt, charthick=thickspt,font=lfont,color=color_cmb_bao
oplot, [xtxt-0.01, xtxt-0.001], [(ytxt-dytxt*3)+0.04,(ytxt-dytxt*3)+0.04], thick=7,linestyle=ls[3],color=color_cmb_bao
i=4 & xyouts,xtxt,ytxt-dytxt*i,'SPT+WMAP7+',chars=bigspt, charthick=thickspt,font=lfont,color=color_cmb_h0_bao
xyouts,xtxt+.0165,ytxt-dytxt*i-0.06,'H!D0!N+BAO',chars=bigspt, charthick=thickspt,font=lfont,color=color_cmb_h0_bao
oplot, [xtxt-0.01, xtxt-0.001], [(ytxt-dytxt*i)+0.04,(ytxt-dytxt*i)+0.04], thick=7,linestyle=ls[4],color=color_cmb_h0_bao


;-----------------
; Plot 2
; colors
cc = [ pscolors[0], $ ; black
       pscolors[8], $ ; blue
       pscolors[2], $ ; red
       pscolors[10],$ ; purple
       pscolors[3], $   ; orange
       pscolors[5]]   ; forest

; c_h0  = pscolors[8] ; blue
; c_bao = pscolors[2] ; red

;ls = [0,3,2]
ls = [0,2,3]

bigspt = 1.
thickspt = 1
lfont = -1
lthick = 5
lsize=1.7

xtxt=1.00 ; xyouts label

i=1
!x.window=[x0+i*(ddx+xoff),x0+(i+1)*ddx+i*xoff]
!x.crange=[0.92, 1.05]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=5
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=5
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=5
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=5

; axes labels
xyouts, 0.904, 0.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, ytitle
xyouts, 0.975, -0.16, charsize=1.5, charthick=3, font=0, xtitle

oplot,[0,0],[0,0],color=3
oplot, [1,1],[0,0.1], linestyle=0, thick=10
;subsamp = 5

; ;;; H0
; plot_like1dname,files_lcdm_h0,pname_lcdm_h0,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[0], color=cc[0], ylog=0,/oplot
; plot_like1dname,files_mnu_h0,pname_mnu_h0,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[1], color=cc[1], /oplot
; plot_like1dname,files_neff_h0,pname_neff_h0,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[2], color=cc[2], /oplot
 
; ;;; BAO
; plot_like1dname,files_lcdm_bao,pname_lcdm_bao,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[0], color=c_bao, /oplot
; plot_like1dname,files_mnu_bao,pname_mnu_bao,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[1], color=c_bao, /oplot
; plot_like1dname,files_neff_bao,pname_neff_bao,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[2], color=cc[4], /oplot

;if 0 then begin ;DEB
;;; H0_BAO
plot_like1dname,files_lcdm_h0_bao,pname_lcdm_h0_bao,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[0], color=cc[0], /oplot
plot_like1dname,files_mnu_h0_bao,pname_mnu_h0_bao,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[1], color=cc[5], /oplot
plot_like1dname,files_neff_h0_bao,pname_neff_h0_bao,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[2], color=cc[3], /oplot
;endif ; DEB

; ;;; CMB-only
; plot_like1dname,files_lcdm_cmb,pname_lcdm_cmb,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[0], color=cc[0], /oplot
; plot_like1dname,files_mnu_cmb,pname_mnu_cmb,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[1], color=cc[1], /oplot
; plot_like1dname,files_neff_cmb,pname_neff_cmb,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[2], color=cc[2], /oplot

xyouts,xtxt-0.01,0.92,'SPT+WMAP7+',charsize=1.4,color=cc[0],font=0,charthick=lthick
xyouts,xtxt+0.009,0.86,'H!D0!X!N+BAO',charsize=1.4,color=cc[0],font=0,charthick=lthick

; lcdm
i=0 & xyouts,xtxt,0.74-i*0.08,'!4K!X!NCDM',charsize=lsize,color=cc[0],font=lfont,charthick=lthick
oplot,[xtxt-0.01, xtxt-0.001],[0.74-i*0.08+0.02,0.74-i*0.08+0.02],thick=7,linestyle=ls[0],color=cc[0]

; mnu
i=1 & xyouts,xtxt,0.74-i*0.08,'!4K!X!NCDM+!4R!Xm!D!4m!X!N',charsize=lsize,color=cc[5],font=lfont,charthick=lthick
oplot,[xtxt-0.01, xtxt-0.001],[0.74-i*0.08+0.02,0.74-i*0.08+0.02],thick=7,linestyle=ls[1],color=cc[5]

; neff
i=2 & xyouts,xtxt,0.74-i*0.08,'!4K!X!NCDM+!3N!Deff!X!N',charsize=lsize,color=cc[3],font=lfont,charthick=lthick
oplot,[xtxt-0.01, xtxt-0.001],[0.74-i*0.08+0.02,0.74-i*0.08+0.02],thick=7,linestyle=ls[2],color=cc[3]

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
;spawn, 'cp '+filename_stub+'.pdf figs/.'
print, 'output: ', filename_stub+'.pdf'
stop
END
