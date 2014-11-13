;;;
; NAME: plots_s12_lrange.pro
; PURPOSE:
;   Make LCDM plot
;
; MODIFICATION HISTORY:
;  08/12/2012: (KTS) Created
;  09/10/2012: (KTS) Use 0828 chains
;  04/19/2013: (KTS) Copied from plot_s12_lcdm
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 8-panel LCDM plot
PRO plot_lrange_8panel, zoom=zoom

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

cdir = '/data23/kstory/lps12/chains/'
;cdir = '/data23/hou/lps12/paramfits/chains_0828/'
;chains  = ['c1_lcdm_pico_w7_newKp','c27_lcdm_camb_s12tau_newKp','c2_lcdm_pico_w7s12_newKp']

;;; SPT+WMAP7, lmax = 1500
dir_lmx15   = cdir+'c2_lcdm_pico_w7s12_lmax1500_newKp/chains/'
files_lmx15 = file_search(dir_lmx15+'c2_lcdm_pico_w7s12_lmax1500_newKp*.txt')
pname_lmx15 = dir_lmx15+'c2_lcdm_pico_w7s12_lmax1500_newKp.paramnames'

;;; SPT+WMAP7, 1500 < l < 3000
dir_lmx30 = cdir+'c2_lcdm_pico_w7s12_1500ell3000_newKp/chains/'
files_lmx30 = file_search(dir_lmx30+'c2_lcdm_pico_w7s12_1500ell3000_newKp*.txt')
pname_lmx30 = dir_lmx30+'c2_lcdm_pico_w7s12_1500ell3000_newKp.paramnames'


; S12-only
; dir_s12   = cdir+'c27_lcdm_camb_s12tau_newKp/chains/'
; files_s12 = file_search(dir_s12+'c27_lcdm_camb_s12tau_newKp*.txt')
; pname_s12 = dir_s12+'c27_lcdm_camb_s12tau_newKp.paramnames'

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7_newKp/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7_newKp*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7_newKp.paramnames'

; S12+wmap
dir_w7s12 = cdir+'c2_lcdm_pico_w7s12_newKp/chains/'
files_w7s12 = file_search(dir_w7s12+'c2_lcdm_pico_w7s12_newKp*.txt')
pname_w7s12 = dir_w7s12+'c2_lcdm_pico_w7s12_newKp.paramnames'


;;; As(k0=0.002)
;cdir = '/data23/hou/lps12/paramfits/chains_0828/'
; ; S12-only
; dir_s12   = cdir+'c27_lcdm_camb_s12tau/chains/'
; files_s12 = file_search(dir_s12+'c27_lcdm_camb_s12tau*.txt')
; pname_s12 = dir_s12+'c27_lcdm_camb_s12tau.paramnames'

; ; WMAP-only
; dir_w7 = cdir+'c1_lcdm_pico_w7/chains/'
; files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7*.txt')
; pname_w7 = dir_w7+'c1_lcdm_pico_w7.paramnames'

; ; S12+wmap
; dir_w7s12 = cdir+'c2_lcdm_pico_w7s12/chains/'
; files_w7s12 = file_search(dir_w7s12+'c2_lcdm_pico_w7s12*.txt')
; pname_w7s12 = dir_w7s12+'c2_lcdm_pico_w7s12.paramnames'

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 8
csize = 0.8
yb = xb * 0.5

filename_stub = '~kstory/lps12/scripts/figs/lcdm_lrange'
if keyword_set(zoom) then filename_stub = '~kstory/lps12/scripts/figs/lcdm_lrange_zoom'
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
color_lmx15   = pscolors[5] ; forest
color_lmx30    = pscolors[10] ; purple
color_w7    = pscolors[2] ; red

ls      = [3,4,0,2] ; [SPT-only, WMAP-only, SPT+WMAP]
mycolor = [color_lmx15, color_lmx30, color_w7s12, color_w7] ; [SPT-only, WMAP-only, SPT+WMAP]


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

;;; Top Row
for i=0,3 do begin
    !x.window=[x0+i*ddx,x0+(i+1)*ddx]
    !y.window=[(y0+ddy)+y0,(y0+ddy)+y0+ddy]

    case i of
        0: begin ; omegahb plot
            !x.crange=[1.8, 2.92]
            if keyword_set(zoom) then !x.crange=[2.0, 2.44]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xyouts, 1.61, ytxt[1], alignment=0.5, orientation=90, charsize=0.9, font=-1, charthick=axlCthick, ytitle
            xtitle='100!4X!3!Db!Nh!U2!X!N'
            xyouts, 2.18, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_w7,pname_w7,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[3], color=mycolor[3], /oplot
            plot_like1dname,files_lmx15,pname_lmx15,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_lmx30,pname_lmx30,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot

            if ~keyword_set(zoom) then begin
                oplot, [2.6,2.9], [0.9+0.02,0.9+0.02], thick=4,linestyle=ls[2],color=mycolor[2]
                oplot, [2.6,2.9], [0.8+0.02,0.8+0.02], thick=4,linestyle=ls[1],color=mycolor[1]
                oplot, [2.6,2.9], [0.7+0.02,0.7+0.02], thick=4,linestyle=ls[0],color=mycolor[0]
                oplot, [2.6,2.9], [0.6+0.02,0.6+0.02], thick=4,linestyle=ls[3],color=mycolor[3]
            endif 
            
        endcase
        1: begin ; omegadmh2
            !x.crange=[0.075, 0.139]
            if keyword_set(zoom) then !x.crange=[0.09, 0.135]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.02
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!4X!3!Dc!Nh!U2!X!N'
            xyouts, 0.105, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_w7,pname_w7,'omegadmh2',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[3], color=mycolor[3], /oplot
            plot_like1dname,files_lmx15,pname_lmx15,'omegadmh2',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_lmx30,pname_lmx30,'omegadmh2',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'omegadmh2',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
            if ~keyword_set(zoom) then begin
                xyouts, 0.076, 0.9, ' 650< l <3000', charsize=0.7, charthick=2, font=0, color=mycolor[2]
                xyouts, 0.076, 0.8, '1500< l <3000', charsize=0.7, charthick=2, font=0, color=mycolor[1]
                xyouts, 0.076, 0.7, ' 650< l <1500', charsize=0.7, charthick=2, font=0, color=mycolor[0]
                xyouts, 0.076, 0.6, ' WMAP7-only', charsize=0.7, charthick=2, font=0, color=mycolor[3]
            endif

            if keyword_set(zoom) then begin
                oplot, [0.122,0.134], [0.9+0.02,0.9+0.02], thick=4,linestyle=ls[2],color=mycolor[2]
                oplot, [0.122,0.134], [0.8+0.02,0.8+0.02], thick=4,linestyle=ls[1],color=mycolor[1]
                oplot, [0.122,0.134], [0.7+0.02,0.7+0.02], thick=4,linestyle=ls[0],color=mycolor[0]
                oplot, [0.122,0.134], [0.6+0.02,0.6+0.02], thick=4,linestyle=ls[3],color=mycolor[3]
            endif

        endcase
        2: begin ; theta_s
            !x.crange=[1.028, 1.053]
            if keyword_set(zoom) then !x.crange=[1.028, 1.048]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.01
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!3 100 !4h!D!3s!X!N'
            xyouts, 1.039, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_w7,pname_w7,'theta_s',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[3], color=mycolor[3], /oplot
            plot_like1dname,files_lmx15,pname_lmx15,'theta_s',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_lmx30,pname_lmx30,'theta_s',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'theta_s',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
            if keyword_set(zoom) then begin
                xyouts, 1.0285, 0.9, ' 650< l <3000', charsize=0.7, charthick=2, font=0, color=mycolor[2]
                xyouts, 1.0285, 0.8, '1500< l <3000', charsize=0.7, charthick=2, font=0, color=mycolor[1]
                xyouts, 1.0285, 0.7, ' 650< l <1500', charsize=0.7, charthick=2, font=0, color=mycolor[0]
                xyouts, 1.0285, 0.6, ' WMAP7-only', charsize=0.7, charthick=2, font=0, color=mycolor[3]
            endif

        endcase
        3: begin ; omegal
            !x.crange=[0.56, 0.89]
            if keyword_set(zoom) then !x.crange=[0.61,0.85]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.1 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!N!4X!DK!X!N'
            xyouts, 0.73, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_w7,pname_w7,'omegal*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[3], color=mycolor[3], /oplot
            plot_like1dname,files_lmx15,pname_lmx15,'omegal*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_lmx30,pname_lmx30,'omegal*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'omegal*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        else: begin
            print, 'i = ', i
        endcase
    endcase
endfor

;;; Bottom Row
for i=0,3 do begin
    !x.window=[x0+i*ddx,x0+(i+1)*ddx]
    !y.window=[y0,y0+ddy]

    case i of
        0: begin ; tau
            !x.crange=[0.041, 0.149]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xyouts, 0.022, ytxt[1], alignment=0.5, orientation=90, charsize=0.9, font=-1, charthick=axlCthick, ytitle
            xtitle='!4s!X!N'
            xyouts, 0.095, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            ;;   Do not plot s12, since we use the WMAP prior on tau
            plot_like1dname,files_w7,pname_w7,'tau',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[3], color=mycolor[3], /oplot
            plot_like1dname,files_lmx15,pname_lmx15,'tau',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_lmx30,pname_lmx30,'tau',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'tau',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        1: begin ; ns
            !x.crange=[0.81, 1.05]
            if keyword_set(zoom) then !x.crange=[0.91,1.02]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!3n!Ds!X!N'
            xyouts, 0.96, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_w7,pname_w7,'ns',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[3], color=mycolor[3], /oplot
            plot_like1dname,files_lmx15,pname_lmx15,'ns',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_lmx30,pname_lmx30,'ns',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'ns',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        2: begin ; logA
            !x.crange=[1.81, 2.49]
            ;!x.crange=[2., 3.]
            if keyword_set(zoom) then !x.crange=[1.81, 2.49]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            ;xtitle='!3log[10!U10!N!4D!3!U2!DR!N]!X'
            xtitle='!310!U9!N!4D!3!U2!DR!N!X'
            ;xyouts, 2.38, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle
            xyouts, 2.1, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            scale=1.

            oplot,[0,0],[0,0],color=pscolors[0]
;             plot_like1dname,files_lmx15,pname_lmx15,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
;               thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
              thick=4,linestyle=ls[3], color=mycolor[3], /oplot
            plot_like1dname,files_lmx15,pname_lmx15,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_lmx30,pname_lmx30,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        3: begin ; H0
            !x.crange=[61, 89.9]
            if keyword_set(zoom) then !x.crange=[61,83]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!N!3H!D0!X!N'
            xyouts, 72, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_w7,pname_w7,'H0*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[3], color=mycolor[3], /oplot
            plot_like1dname,files_lmx15,pname_lmx15,'H0*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_lmx30,pname_lmx30,'H0*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,'H0*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        else: begin
            !x.crange=[0,9.5]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            oplot,[0,0],[0,0],color=pscolors[0]
            oplot, xtmp, ytmp
        endcase
    endcase
            
endfor
;*******************************
ps_close

print, 'output file: '+filename_stub+'.pdf'
spawn,'epstopdf '+filename_stub+'.ps'

stop
END


