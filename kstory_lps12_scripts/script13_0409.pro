;;;
; NAME: script13_0409
;
; NOTES:
;  1) plot s12 and k11 lcdm constraints
;;;


PRO plot_k11

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

cdir = '/data23/kstory/lps12/chains/'
chains  = ['c1_lcdm_pico_w7_newKp','c27_lcdm_camb_s12tau_newKp','c2_lcdm_pico_w7s12_newKp']



; S12-only
dir_s12   = cdir+'c27_lcdm_camb_s12tau_newKp/chains/'
files_s12 = file_search(dir_s12+'c27_lcdm_camb_s12tau_newKp*.txt')
pname_s12 = dir_s12+'c27_lcdm_camb_s12tau_newKp.paramnames'

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7_newKp/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7_newKp*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7_newKp.paramnames'

; S12+wmap
dir_w7s12 = cdir+'c2_lcdm_pico_w7s12_newKp/chains/'
files_w7s12 = file_search(dir_w7s12+'c2_lcdm_pico_w7s12_newKp*.txt')
pname_w7s12 = dir_w7s12+'c2_lcdm_pico_w7s12_newKp.paramnames'

; K11+wmap
file_w7k11 = '/home/rkeisler/ps09/cons_chain_baseline.sav'
pname_w7k11 = '/home/rkeisler/ps09/paramnames_alens.txt'

restore, file_w7k11


;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 8
csize = 0.8
yb = xb * 0.5

filename_stub = '~kstory/lps12/scripts/figs/lcdm_8panel'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/landscape

;Define the axis size relative to the window size
x0=.095
ddx=.22
y0=.12
ddy=.35

; Define the output plot and font sizes
axlCthick = 3
axlCsize  = 0.9;0.8
ytxt = [-0.23, .5] ; y position of [x,y]-axis label
;bigspt = 1.4
;thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_w7s12 = pscolors[0] ; black
color_s12   = pscolors[8] ; blue
color_w7    = pscolors[2] ; !red
color_k11   = pscolors[10] ; !purple

ls      = [3,2,0] ; [SPT-only, WMAP-only, SPT+WMAP]
mycolor = [color_s12, color_w7, color_w7s12] ; [SPT-only, WMAP-only, SPT+WMAP]

ls      = [3,2,0] ; [SPT-only, WMAP-only, SPT+WMAP]
mycolor = [color_s12, color_w7, color_w7s12] ; [SPT-only, WMAP-only, SPT+WMAP]


!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

; axis properties
!y.crange=[0,1]
ytitle='!X!3Likelihood!X!N'

; TMP
!x.crange=[0,9.5]
xtmp=indgen(10)
ytmp=0.1*xtmp

i=0
!x.window=[x0+i*ddx,x0+(i+1)*ddx]
!y.window=[(y0+ddy)+y0,(y0+ddy)+y0+ddy]

; -------------
; Plot
            !x.crange=[1.8, 2.92]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xyouts, 1.61, ytxt[1], alignment=0.5, orientation=90, charsize=0.9, font=-1, charthick=axlCthick, ytitle
            xtitle='100!4X!3!Db!Nh!U2!X!N'
            xyouts, 2.21, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            ; K11
            scale=100
            subsamp=3 
            pp = echain[2,nskip:*]*scale
            sig = stddev(pp)
            binsize=sig/subsamp
            h=double(histogram(pp,binsize=binsize,omin=omn,omax=omx))
            h/=max(h)
            bins = (findgen(n_elements(h))+.5)*binsize +omn

            oplot,[0,0],[0,0],color=pscolors[0]
;             plot_like1dname,files_s12,pname_s12,'omegabh2',subsamp=3,nskip=1000,scale=100,$
;               thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            oplot, bins,h,thick=4,linestyle=ls[1],color=color_k11

;             plot_like1dname,files_w7,pname_w7,'omegabh2',subsamp=3,nskip=1000,scale=100,$
;               thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot

            oplot, [2.6,2.9], [0.9+0.02,0.9+0.02], thick=4,linestyle=ls[2],color=mycolor[2]
            oplot, [2.6,2.9], [0.8+0.02,0.8+0.02], thick=4,linestyle=ls[1],color=color_k11
            ;oplot, [2.6,2.9], [0.7+0.02,0.7+0.02], thick=4,linestyle=ls[0],color=mycolor[0]

;------------------


ps_close

print, 'output file: figs/lcdm_8panel.pdf'
spawn,'epstopdf '+filename_stub+'.ps'


stop
END

