;;;
; NAME: plots_s12_nslike
; PURPOSE:
;   Make ns likelihood and CDF plot
;
; MODIFICATION HISTORY:
;  08/15/2012: (KTS) Created
;  09/10/2012: (KTS) Update with 0828 chains
;;;

;;;;;;;;;;;;;;;;;;;
; 
PRO plot_ns
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

param = 'ns'
subsamp = 5
nskip=2000

;----------------------
; Get the data

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

; ; CMB+BAO+H0
; dir_cmb_bao_h0 = cdir+'/c5_lcdm_pico_w7s12_BAO_H0/chains/'
; files_cmb_bao_h0 = file_search(dir_cmb_bao_h0+'c5_lcdm_pico_w7s12_BAO_H0*.txt')
; pname_cmb_bao_h0 = dir_cmb_bao_h0+'c5_lcdm_pico_w7s12_BAO_H0.paramnames'

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 10
csize = 1.2
yb = 4                          ;xb * 0.45

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.095
ddx=.39
xoff=.1
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3                         ; points
pthick=5                        ; points
xtxt=0.99                       ; xyouts label [LCDM]
;xtxt=1.01 ; xyouts label [LCDM+r]
ytxt=0.93
dytxt = 0.08
bigspt = 1.2
thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_w = pscolors[2]           ; red
color_cmb = pscolors[0]         ; black
color_cmb_h0 = pscolors[3]      ; orange
color_cmb_bao = pscolors[5]     ; forest

ls = [2, 0, 3, 1] ; w, cmb, cmb_h0, cmb_bao

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]


for i=0,1 do begin
    !x.window=[x0+i*(ddx+xoff),x0+(i+1)*ddx+i*xoff]
    !x.crange=[0.92, 1.05]      ; LCDM
    !y.window=[y0,y0+ddy]

    axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
    axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
    xtitle='!3n!Ds!X!N'

    ; Plot #1
    if i eq 0 then begin
        !y.crange=[0, 1]
        axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 ,yl=0
        axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2, yl=0

        plot_like1dname,files_w,pname_w,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[0], color=color_w, /oplot
        plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[1], color=color_cmb, /oplot
        plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[2], color=color_cmb_h0, /oplot
        plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[3], color=color_cmb_bao, /oplot

        ytitle='!3Likelihood!X!N'
        xyouts, 0.903, 0.5, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
        xyouts, 0.98, -0.15, charsize=1.4, charthick=3, xtitle

        xyouts,xtxt,ytxt,     'WMAP7',chars=bigspt,charthick=thickspt,font=lfont,color=color_w
        oplot, [xtxt-0.01, xtxt-0.001], [ytxt+0.04,ytxt+0.04], thick=7,linestyle=ls[0],color=color_w
        xyouts,xtxt,ytxt-dytxt*1, 'SPT+WMAP7',chars=bigspt, charthick=thickspt,font=lfont,color=color_cmb
        oplot, [xtxt-0.01, xtxt-0.001], [(ytxt-dytxt*1)+0.04,(ytxt-dytxt*1)+0.04], thick=7,linestyle=ls[1],color=color_cmb
        xyouts,xtxt,ytxt-dytxt*2,'SPT+WMAP7+H!D0!X!N',chars=bigspt, charthick=thickspt,font=lfont,color=color_cmb_h0
        oplot, [xtxt-0.01, xtxt-0.001], [(ytxt-dytxt*2)+0.04,(ytxt-dytxt*2)+0.04], thick=7,linestyle=ls[2],color=color_cmb_h0
        xyouts,xtxt,ytxt-dytxt*3,'SPT+WMAP7+BAO',chars=bigspt, charthick=thickspt,font=lfont,color=color_cmb_bao
        oplot, [xtxt-0.01, xtxt-0.001], [(ytxt-dytxt*3)+0.04,(ytxt-dytxt*3)+0.04], thick=7,linestyle=ls[3],color=color_cmb_bao

    ; Plot #2
    endif else begin
        !y.crange=[-4.5,0]
        axis,yaxis=0,font=0,ystyle=1,/yl,/save,ycharsize=csize,ythick=2,ytickname=['10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','10!U0!N']
        axis,yaxis=1,ytickname=empty,ystyle=1,/yl,ythick=2

        oplot,[0,0],[0,0],color=pscolors[0]

        plot_like1dname,files_w,pname_w,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[0],color=color_w,/cdf,/oplot
        plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[1],color=color_cmb,/cdf,/oplot
        plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[2],color=color_cmb_h0,/cdf,/oplot
        plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,$
          yl=1,thick=7,linestyle=ls[3],color=color_cmb_bao,/cdf,/oplot


        ytitle='!3P(>n!Ds!N)!X!N'
        xyouts, 0.903, 5e-3, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle ; LCDM
;        xyouts, 0.90, 5e-3, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle ; LCDM+r
        xyouts, 0.98, 0.6e-5, charsize=1.4, charthick=3, font=-1, xtitle
        
    endelse

endfor
;*******************************
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/nslike_0921.pdf'

stop
END



