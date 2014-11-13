;;;
; NAME: plots_s12_dl_all.pro
; PURPOSE:
;   Make bandpower comparison plot
;
; NOTES:
; 1) Test SV averaging
;
; MODIFICATION HISTORY:
;  08/09/2012: (KTS) Created
;  09/12/2012: (KTS) Update to use 0828 bandpowers
;;;

;;;;;;;;;;;;;;;;;;
; Procedure to make bandpower comparison plot.
; other options in script_0801.pro
PRO plot_dl_all_1

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120828_170101_kweight.sav'
uk2mk = 1d-6

; load theory
;restore, '/home/kstory/lps12/scripts/plotting/bestfit_WMAP7_lps12.sav'
readcol, '/home/kstory/lps12/best_fit/lcdm_w7s12_ML_dl.txt', l_total, dl_total
best_cal = 0.98773074 ; from plotting/bestfit_WMAP7_lps12.sav
cl4_fit = dl_total *(l_total^4.) / (l_total*(l_total+1)) * uk2mk

; readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
; cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
; dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)
; cl4_th = dl_th *(l_vec^4.) / (l_vec*(l_vec+1)) * uk2mk
; cl4_fit = cl4_th
; l_total = l_vec


; S12
restore,file
dl_all_lps12 = dl_all*(1d12) * best_cal
diag_nobeam_lps12 = diag_nobeam*(1d12)
l_lps12 = l

istart=9
istop=55
dl_all_lps12 = dl_all_lps12[istart:istop]
diag_nobeam_lps12 = diag_nobeam_lps12[istart:istop]
l_lps12 = l_lps12[istart:istop]
cl4_lps12 = dl_all_lps12 *(l_lps12^4.) / (l_lps12*(l_lps12+1)) * uk2mk
dcl4_lps12 = diag_nobeam_lps12 *(l_lps12^4.) / (l_lps12*(l_lps12+1)) * uk2mk

; load  K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
cal_rk = 1;0.76
dl_all_k11      = dl_all*1d12*cal_rk^2. * 0.97
diag_nobeam_k11 = diag_nobeam*1d12*cal_rk^2.
l_k11 = l

dl_all_k11 = dl_all_k11[istart:istop]
diag_nobeam_k11 = diag_nobeam_k11[istart:istop]
l_k11 = l_k11[istart:istop]
cl4_k11 = dl_all_k11 *(l_k11^4.) / (l_k11*(l_k11+1)) * uk2mk
dcl4_k11 = diag_nobeam_k11 *(l_k11^4.) / (l_k11*(l_k11+1)) * uk2mk

sf=1.

; load act (das et al)
readcol,pdir+'act_150_use.txt',l_act,dl_act,ddl_act
l_act   = l_act[0:39]
dl_act  = dl_act[0:39]
ddl_act = ddl_act[0:39]
cl4_act = dl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk
dcl4_act = ddl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk

; load acbar
readcol,pdir+'cl_dc1_v3.dat',iband,dl_acbar,delta_dl_acbar,log_normal_offset
ledge=[[350,550],[550,650],[650,730],[730,790],[790,850],[850,910],[910,970],[970,1030],[1030,1090],[1090,1150],[1150,1210],[1210,1270],[1270,1330],[1330,1390],[1390,1450],[1450,1510],[1510,1570],[1570,1650],[1650,1750],[1750,1850],[1850,1950],[1950,2100],[2100,2300],[2300,2500],[2500,3000]]
llo=reform(ledge[0,*])
lhi=reform(ledge[1,*])
l_acbar=0.5*(llo+lhi)
ddl_acbar = delta_dl_acbar
cl4_acbar = dl_acbar *(l_acbar^4.) / (l_acbar*(l_acbar+1)) * uk2mk
dcl4_acbar = ddl_acbar *(l_acbar^4.) / (l_acbar*(l_acbar+1)) * uk2mk


;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 10
csize = 1.2
yb = 4;xb * 0.45

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.095
ddx=.39
xoff=.1
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3  ; points
pthick=3 ; points
athick=5 ; axes
xtxt=1700 ; xyouts label
ytxt=2700
dytxt = 0.63
bigspt = 1.4
thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_th = pscolors[0]        ; black
color_act = pscolors[3]       ; orange
color_acbar = pscolors[9]     ; purple
color_k11 = pscolors[5]       ; forest
color_this_work = pscolors[8] ; blue

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]


for i=0,1 do begin
    !x.window=[x0+i*(ddx+xoff),x0+(i+1)*ddx+i*xoff]
    ;!x.window=[x0+i*ddx,x0+(i+1)*ddx]
    !x.crange=[450, 3200]
    !y.window=[y0,y0+ddy]

    axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=athick 
    axis,xaxis=1,xtickname=empty,xstyle=1,xthick=athick
    
    if i eq 0 then begin
        !y.crange=[ALOG10(30),ALOG10(4000)]
        axis,yaxis=0,font=0,ystyle=1,/yl,/save,ycharsize=csize,ythick=athick 
        axis,yaxis=1,ytickname=empty,ystyle=1,/yl,ythick=athick
    endif else begin
        !y.crange=[0, 2000]
        axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=athick ,yl=0
        axis,yaxis=1,ytickname=empty,ystyle=1,ythick=athick, yl=0
    endelse


    ; Plot #1
    if i eq 0 then begin
        oplot,[0,0],[0,0],color=pscolors[0]

        oplot,l_lps12/sf,dl_all_lps12,ps=3,color=color_this_work
        errplot,l_lps12/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=4,color=color_this_work

        xtitle='!12l!X!N'

        ;ytitle='!12l(l+1)C!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'
        ytitle='!12D!Dl!X!N (!7l!6K!3!U2!X!N)'
        xyouts, 0, 300, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
        xyouts, 1800, 13, charsize=1.4, charthick=3, font=-1, xtitle

        xyouts,xtxt,ytxt*dytxt^2.,'SPT',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
        
    ; Plot #2
    endif else begin
        oplot,[0,0],[0,0],color=3
        ;oplot,l_vec[450:3200],cl4_th[450:3200],color=color_th, thick=thick

        wh = where(l_total ge 450 and l_total lt 3200)
        oplot,l_total[wh],cl4_fit[wh],color=color_th, thick=thick

        errplot,l_acbar,cl4_acbar-dcl4_acbar,cl4_acbar+dcl4_acbar,thick=pthick,color=color_acbar
        errplot,l_act,cl4_act-dcl4_act,cl4_act+dcl4_act,thick=pthick,color=color_act
        errplot,l_k11,(cl4_k11-dcl4_k11),(cl4_k11+dcl4_k11),thick=pthick,color=color_k11
        errplot,l_lps12,(cl4_lps12-dcl4_lps12),(cl4_lps12+dcl4_lps12),thick=pthick+2,color=color_this_work

        ytitle='!12l!U4!N!12C!Dl!3!N/2!7p!X!N (!6mK!3!U2!X!N)'
        xyouts, 50, 1000, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
        xyouts, 1800, -330, charsize=1.4, charthick=3, xtitle

        xyouts,xtxt+180,1800,       'ACBAR',chars=bigspt,charthick=thickspt,font=lfont,color=color_acbar
        xyouts,xtxt+180,1650,         'ACT',chars=bigspt,charthick=thickspt,font=lfont,color=color_act
        xyouts,xtxt+180,1500,   'SPT, K11',chars=bigspt, charthick=thickspt,font=lfont,color=color_k11
        xyouts,xtxt+180,1350,'SPT, this work',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
    endelse

endfor
;*******************************
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/dl_all.pdf'
print, 'output file: figs/dl_all.pdf'

stop
END





