;;;
; NAME: script12_0917
;
; NOTES:
;  1) transparency plots
;;;

;; Transparent plots
PRO trans1

cdir = '/data23/hou/lps12/paramfits/chains_0828/' & 

dir_w7 = cdir+'c1_lcdm_pico_w7/chains/' &  files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7*.txt') & pname_w7 = dir_w7+'c1_lcdm_pico_w7.paramnames'
dir_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'  & files_cmb = file_search(dir_cmb+'c2_lcdm_pico_w7s12*.txt') & pname_cmb = dir_cmb+'c2_lcdm_pico_w7s12.paramnames'

;;;;--------------------
x0=.20
ddx=.75
y0=.16
ddy=.77

csize = 1.2

!x.window=[x0,x0+ddx]
!x.crange=[60,80]
!y.window=[y0,y0+ddy]
!y.crange=[1/15., 1/12.];!y.crange=[12, 15]

xb = 550
yb = 550


; p1
window,2,xsize=xb,ysize=yb
oplot,[0,0],[0,0]
scale2=1. / (154.66/150.82)
plot_like2dname, files_w7,pname_w7, 'H0*','Dvp57ors',sigma=[1,2,3],icolors=[!black,!green,!blue], $
  scale2=scale,power2=(-1.0),$
  smooth=0.0002,/overplot
  ;yr=!y.crange,xr=!x.crange
;oplot, [60,80], [1/15., 1/12.]
p1 = TVREAD(TRUE=3)
;stop

; p2
window,2,xsize=xb,ysize=yb
oplot,[0,0],[0,0]
plot_like2dname, files_cmb,pname_cmb, 'H0*','Dvp57ors',sigma=[1,2,3],icolors=[!red,!orange,!yellow], $
  scale2=scale,power2=(-1.0),$
  smooth=0.0002,/overplot
; x2=[60, 75] & y2=[0.065, 0.08]
; POLYFILL, [x2[0], x2[1], x2[1], x2[0], x2[0]], [y2[0],  y2[0], y2[1], y2[1], y2[0]], /DATA, COLOR=FSC_COLOR('deep pink')
p2 = TVREAD(TRUE=3)
;stop

; p3
; window,2,xsize=xb,ysize=yb
; oplot,[0,0],[0,0]
; x2=[65, 80] & y2=[0.07, 0.085]
; POLYFILL, [x2[0], x2[1], x2[1], x2[0], x2[0]], [y2[0],  y2[0], y2[1], y2[1], y2[0]], /DATA, COLOR=FSC_COLOR('blue')
; p3 = TVREAD(TRUE=3)

window,3,xsize=xb,ysize=yb
TV, p1*0.5 + p2*0.5, TRUE=3
axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

stop

END

;; Transparent plots
PRO trans
;setup_ps_plotting
cdir = '/data23/hou/lps12/paramfits/chains_0828/' & 

dir_w7 = cdir+'c1_lcdm_pico_w7/chains/' &  files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7*.txt') & pname_w7 = dir_w7+'c1_lcdm_pico_w7.paramnames'
dir_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'  & files_cmb = file_search(dir_cmb+'c2_lcdm_pico_w7s12*.txt') & pname_cmb = dir_cmb+'c2_lcdm_pico_w7s12.paramnames'

x0=.20
ddx=.75
y0=.16
ddy=.77

csize = 1.0

!x.window=[x0,x0+ddx]
!x.crange=[60,80]
!y.window=[y0,y0+ddy]
!y.crange=[12, 15]

wset, 1
plot_like2dname, files_w7,pname_w7, 'H0*','Dvp57ors',sigma=[1,2,3],colors=[!black,!green,!blue], yr=[12,16],xr=[60,85],/oplot
fig1 = tvread(true=3)

wset, 2
plot_like2dname, files_cmb,pname_cmb, 'H0*','Dvp57ors',sigma=[1,2,3],colors=[!red,!orange,!yellow],yr=[12,16],xr=[60,85],/oplot
fig2 = tvread(true=3)


setup_ps_plotting
filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

alpha = 0.5
tv, fig1*alpha + fig2*alpha, true=3

ps_close

END



