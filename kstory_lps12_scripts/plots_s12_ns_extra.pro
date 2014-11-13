;;;
; NAME: plots_s12_ns_extra.pro
; PURPOSE:
;   Make ns + other models plot
;
; NOTES:
;
; MODIFICATION HISTORY:
;  08/22/2012: (KTS) created
;  09/10/2012: (KTS) Update with 0828 chains
;;;


;;;;;;;;;;;;;;;;;;;
; 
; 2d ns plot
;
;;;;;;;;;;;;;;;;;;;

PRO ns_plot_2d, nsigma=nsigma
if n_elements(nsigma) eq 0 then nsigma = 2

cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data, CMB+BAO+H0

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

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 11
csize = 1.3
yb = 4.5

filename_stub = '~kstory/lps12/scripts/figs/ns_extra_tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.10
ddx=.38
dx_off = .5
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1.00 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.
thickspt = 1
lfont = -1
lthick = 5
lsize=1.7

;*******************************
; colors
cc = [ pscolors[0], $ ; black
       pscolors[8], $ ; blue
       pscolors[2], $ ; red
       pscolors[10],$ ; magenta
       pscolors[3]]   ; orange

; c_h0  = pscolors[8] ; blue
; c_bao = pscolors[2] ; red

ls = [0,3,2]

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[0.92, 1.05]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=5
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=5
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=5
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=5

; axes labels
xyouts, 0.904, 0.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 0.975, -0.16, charsize=1.5, charthick=3, font=0, '!3n!Ds!X!N'

oplot,[0,0],[0,0],color=3
;subsamp = 5

;;; H0
plot_like1dname,files_lcdm_h0,pname_lcdm_h0,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[0], color=cc[0], ylog=0,/oplot
plot_like1dname,files_mnu_h0,pname_mnu_h0,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[1], color=cc[1], /oplot
plot_like1dname,files_neff_h0,pname_neff_h0,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[2], color=cc[2], /oplot
 
; ;;; BAO
; plot_like1dname,files_lcdm_bao,pname_lcdm_bao,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[0], color=c_bao, /oplot
; plot_like1dname,files_mnu_bao,pname_mnu_bao,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[1], color=c_bao, /oplot
; plot_like1dname,files_neff_bao,pname_neff_bao,'ns',subsamp=8,nskip=1000,$
;   thick=7,linestyle=ls[2], color=c_bao, /oplot


xyouts,xtxt-0.01,0.92,'SPT+WMAP7+H!D0!X!N',charsize=1.4,color=cc[0],font=lfont,charthick=lthick

; lcdm
xyouts,xtxt,0.82,'!4K!X!NCDM',charsize=lsize,color=cc[0],font=lfont,charthick=lthick
oplot,[xtxt-0.01, xtxt-0.001],[0.82+0.02,0.82+0.02],thick=7,linestyle=ls[0],color=cc[0]

; mnu
xyouts,xtxt,0.74,'!4K!X!NCDM+!4R!Xm!D!4m!X!N',charsize=lsize,color=cc[1],font=lfont,charthick=lthick
oplot,[xtxt-0.01, xtxt-0.001],[0.74+0.02,0.74+0.02],thick=7,linestyle=ls[1],color=cc[1]

; neff
xyouts,xtxt,0.66,'!4K!X!NCDM+!3N!Deff!X!N',charsize=lsize,color=cc[2],font=lfont,charthick=lthick
oplot,[xtxt-0.01, xtxt-0.001],[0.66+0.02,0.66+0.02],thick=7,linestyle=ls[2],color=cc[2]


;-----------------
; Plot 2
!x.window=[dx_off+x0,dx_off+x0+ddx]
!y.window=[y0,y0+ddy]
if nsigma eq 3 then begin
    !x.crange=[0.92, 1.03]
    !y.crange=[2.0, 5.0]
endif else begin
    !x.crange=[0.93, 1.01]
    !y.crange=[2.5, 4.5]
endelse


axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=5
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=5
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=5
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=5

oplot,[0,0],[0,0],color=3

if nsigma eq 3 then begin
    xyouts, 0.905, 3.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, '!3N!Deff!X!N' ; axis label
    xyouts, 0.975, 1.6, charsize=1.5, charthick=3, font=0, '!3n!Ds!X!N'
    
    plot_like2dname,files_neff_h0,pname_neff_h0,'ns','Neff',/overplot,sigma=[1,2,3],icolors=[pscolors[2],pscolors[3],pscolors[6]],smooth=0.002,subsamp=2,yl=0;,priorname=['czero_ksz'],priormin=[0],priormax=[100]
    xyouts,0.93,4.7,'!4K!X!NCDM+!3N!Deff!X!N',charsize=lsize, color=cc[2],font=lfont,charthick=lthick

endif 

if nsigma eq 2 then begin
    xyouts, 0.919, 3.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, '!3N!Deff!X!N' ; axis label
    xyouts, 0.975, 2.18, charsize=1.5, charthick=3, font=0, '!3n!Ds!X!N'
    
    plot_like2dname,files_neff_h0,pname_neff_h0,'ns','Neff',/overplot,sigma=[1,2],icolors=[pscolors[2],pscolors[6]],smooth=0.002,subsamp=2,yl=0 ;,priorname=['czero_ksz'],priormin=[0],priormax=[100]
    xyouts,0.935,4.3,'!4K!X!NCDM+!3N!Deff!X!N',charsize=lsize, color=cc[2],font=lfont,charthick=lthick
endif

oplot, [1,1], !y.crange,linestyle=2,thick=7
;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp '+filename_stub+'.pdf figs/ns_extra_2d_0924.pdf'
stop
END





;;;;;;;;;;;;;;;;;;;
;
; Obsolete
;
;;;;;;;;;;;;;;;;;;;

PRO ns_plot
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data, CMB+BAO+H0

; nmu, sum_nu, omegaK, w

; cmb+bao+h0
dir_cmb_bao_h0 = cdir+'/c5_lcdm_pico_w7s12_BAO_H0/chains/'
files_cmb_bao_h0 = file_search(dir_cmb_bao_h0+'c5_lcdm_pico_w7s12_BAO_H0*.txt')
pname_cmb_bao_h0 = dir_cmb_bao_h0+'c5_lcdm_pico_w7s12_BAO_H0.paramnames'

; neff
type = 'c15_lcdm_neff_pico_w7s12_BAO_H0'
dir_neff = cdir+type+'/chains/'
files_neff = file_search(dir_neff+type+'_*.txt')
pname_neff = dir_neff+type+'.paramnames'

; mnu
type = 'c22_lcdm_mnu_camb_w7s12_BAO_H0'
dir_mnu = cdir+type+'/post_chains/'
files_mnu = file_search(dir_mnu+type+'_*.txt')
pname_mnu = dir_mnu+type+'.paramnames'

; ; omk
; type = ''
; dir_omk = cdir+type+'/post_chains/'
; files_omk = file_search(dir_omk+type+'_*.txt')
; pname_omk = dir_omk+type+'.paramnames'

; ; w
; type = '';'c30_lcdm_w_camb_w7s12_BAO'
; dir_w = cdir+type+'/post_chains/'
; files_w = file_search(dir_w+type+'_*.txt')
; pname_w = dir_w+type+'.paramnames'

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
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
lfont = -1
ls = [0,1,2,3,4]

;*******************************
; colors
cc = [ pscolors[8], $ ; blue
       pscolors[2], $ ; blue
       pscolors[5], $ ; forest
       pscolors[10],$ ; magenta
       pscolors[3]]   ; orange


!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[0.92, 1.02]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=5
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=5
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=5
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=5

; axes labels
xyouts, 0.904, 0.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 0.975, -0.16, charsize=1.5, charthick=3, font=0, '!3n!Ds!X!N'

oplot,[0,0],[0,0],color=3
;subsamp = 5

; LCDM
plot_like1dname,files_cmb,pname_cmb,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[0], color=cc[0], /oplot
xyouts,0.985,0.95,'!4X!X!NCDM',charsize=csize, color=cc[0],font=lfont

; neff
plot_like1dname,files_neff,pname_neff,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[1], color=cc[1], /oplot
xyouts,0.985,0.9,'!4X!X!NCDM+neff',charsize=csize, color=cc[1],font=lfont

; mnu
plot_like1dname,files_mnu,pname_mnu,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[2], color=cc[2], /oplot
xyouts,0.985,0.85,'!4X!X!NCDM+mnu',charsize=csize, color=cc[2],font=lfont

; omk
plot_like1dname,files_omk,pname_omk,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[3], color=cc[3], /oplot
xyouts,0.985,0.8,'!4X!X!NCDM+omk',charsize=csize, color=cc[3],font=lfont

; w
plot_like1dname,files_w,pname_w,'ns',subsamp=8,nskip=1000,$
  thick=7,linestyle=ls[4], color=cc[4], /oplot
xyouts,0.985,0.75,'!4X!X!NCDM+w',charsize=csize, color=cc[4],font=lfont


; vertical line
;oplot, [1.,1.], [0,1], linestyle=1, thick=3

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/ns_extra_0822.pdf'
stop
END




;;;;;;;;;
; CMB+BAO

; ; cmb
; dir_cmb = cdir+'c4_lcdm_camb_w7s12/chains/'
; files_cmb = file_search(dir_cmb+'c4_lcdm_camb_w7s12*.txt')
; pname_cmb = dir_cmb+'c4_lcdm_camb_w7s12.paramnames'

; ; neff
; type = 'c12_lcdm_neff_camb_w7s12_BAO'
; dir_neff = cdir+type+'/chains/'
; files_neff = file_search(dir_neff+type+'_*.txt')
; pname_neff = dir_neff+type+'.paramnames'

; ; mnu
; type = 'c55_lcdm_mnu_pico_w7s12_BAO'
; dir_mnu = cdir+type+'/post_chains/'
; files_mnu = file_search(dir_mnu+type+'_*.txt')
; pname_mnu = dir_mnu+type+'.paramnames'

; ; omk
; type = 'c25_lcdm_omk_camb_w7s12_BAO'
; dir_omk = cdir+type+'/post_chains/'
; files_omk = file_search(dir_omk+type+'_*.txt')
; pname_omk = dir_omk+type+'.paramnames'

; ; w
; type = 'c30_lcdm_w_camb_w7s12_BAO'
; dir_w = cdir+type+'/post_chains/'
; files_w = file_search(dir_w+type+'_*.txt')
; pname_w = dir_w+type+'.paramnames'

