;;;
; NAME: script12_0910
; PURPOSE:
;   Make WMAP-only comparison plot
;
; MODIFICATION HISTORY:
;  09/10/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 8-pannel LCDM plot
PRO plot_lcdm_8pannel
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7.paramnames'

; WMAP-0717
dir_old_w7 = '/data23/hou/lps12/paramfits/chains_0717/lcdm_camb_w7/chains/'
files_old_w7 = file_search(dir_old_w7+'lcdm_camb_w7*.txt')
pname_old_w7 = dir_old_w7+'lcdm_camb_w7.paramnames'

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 12
csize = 1.2
yb = xb * 0.5

filename_stub = '~kstory/lps12/scripts/figs/tmp_lcdm'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.095
ddx=.22
y0=.12
ddy=.35

; Define the output plot and font sizes
axlCthick = 4
axlCsize  = 1.5
ytxt = [-0.25, .45] ; y position of [x,y]-axis label
;bigspt = 1.4
;thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_s12   = pscolors[0] ; black
color_w7    = pscolors[2] ; !red
color_w7s12 = pscolors[8] ; blue

ls      = [2,3,0] ; [SPT-only, WMAP-only, SPT+WMAP]
mycolor = [color_s12, color_w7, color_w7s12] ; [SPT-only, WMAP-only, SPT+WMAP]


!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

; axis properties
!y.crange=[0,1]
ytitle='!3Likelihood!X!N'

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
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xyouts, 1.6, ytxt[1], alignment=0.5, orientation=90, charsize=axlCsize, font=-1, charthick=axlCthick, ytitle
            xtitle='!4X!3!Db!Nh!U2!X!N'
            xyouts, 2.3, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_old_w7,pname_old_w7,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=7,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=7,linestyle=ls[1], color=mycolor[1], /oplot

            oplot, [2.6,2.9], [0.8+0.02,0.8+0.02], thick=7,linestyle=ls[1],color=mycolor[1]
            oplot, [2.6,2.9], [0.7+0.02,0.7+0.02], thick=7,linestyle=ls[0],color=mycolor[0]
            
        endcase
        1: begin ; omegadmh2
            !x.crange=[0.075, 0.139]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.02
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!4X!3!Dc!Nh!U2!X!N'
            xyouts, 0.105, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_old_w7,pname_old_w7,'omegadmh2',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'omegadmh2',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[1], color=mycolor[1], /oplot
            
            xyouts, 0.076, 0.8, 'WMAP_0828', charsize=1.0, charthick=2, font=0, color=mycolor[1]
            xyouts, 0.076, 0.7, 'WMAP_0717', charsize=1.0, charthick=2, font=0, color=mycolor[0]

        endcase
        2: begin ; theta_s
            !x.crange=[1.028, 1.053]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.01
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!3 100 !4h!D!3s!X!N'
            xyouts, 1.039, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_old_w7,pname_old_w7,'theta_s',subsamp=3,nskip=1000,scale=100,$
              thick=7,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'theta_s',subsamp=3,nskip=1000,scale=100,$
              thick=7,linestyle=ls[1], color=mycolor[1], /oplot
            
        endcase
        3: begin ; omegal
            !x.crange=[0.56, 0.89]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.1 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!U*!N!4X!DK!X!N'
            xyouts, 0.73, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_old_w7,pname_old_w7,'omegal*',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'omegal*',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[1], color=mycolor[1], /oplot

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

            xyouts, 0.022, ytxt[1], alignment=0.5, orientation=90, charsize=axlCsize, font=-1, charthick=axlCthick, ytitle
            xtitle='!4s!X!N'
            xyouts, 0.095, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            ;;   Do not plot s12, since we use the WMAP prior on tau
            plot_like1dname,files_old_w7,pname_old_w7,'tau',subsamp=3,nskip=1000,$
             thick=7,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'tau',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[1], color=mycolor[1], /oplot
            
        endcase
        1: begin ; ns
            !x.crange=[0.81, 1.05]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!3n!Ds!X!N'
            xyouts, 0.93, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_old_w7,pname_old_w7,'ns',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'ns',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[1], color=mycolor[1], /oplot
            
        endcase
        2: begin ; As
            ;!x.crange=[2.81, 3.9]
            !x.crange=[1.5, 4.2]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            ;xtitle='!3log[10!U10!N!4D!3!U2!DR!N]!X'
            xtitle='!310!U9!N!4D!3!U2!DR!N!X'
            xyouts, 2.5, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            scale=1.

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_old_w7,pname_old_w7,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
              thick=7,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
              thick=7,linestyle=ls[1], color=mycolor[1], /oplot
            
        endcase
        3: begin ; H0
            !x.crange=[61, 90]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!U*!N!3H!D0!X!N'
            xyouts, 75, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,files_old_w7,pname_old_w7,'H0*',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files_w7,pname_w7,'H0*',subsamp=3,nskip=1000,$
              thick=7,linestyle=ls[1], color=mycolor[1], /oplot
            
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

print, 'output file: figs/lcdm_8pannel.pdf'
spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp_lcdm.pdf figs/lcdm_w7_8pannel.pdf'

stop
END



