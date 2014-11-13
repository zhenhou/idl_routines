;;;
; NAME: script_0801
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Analyze svnv test, run_08_lr
; 2) Pipeline test2
; 3) testing dl_all plots
;;;

;;;;;;;;;;;;;;;;;;;;;;;
;;; 1) Analyze svnv test
;;;;;;;;;;;;;;;;;;;;;;;
PRO svnv
;f = lps12_fieldstruct()
;!p.font = 1
edir='/home/kstory/lps12/end2end/'
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

restore, edir+'run_08/combined_spectrum_20120801_151258_kweight_lr.sav'
cov_n = cov_all_data
cov_sv = cov_sv

nl = 58
diag_n  = dblarr(nl)
diag_sv = dblarr(nl)

for i=0, nl-1 do begin
    diag_n[i]  = sqrt(cov_n[i,i]) * 1d12
    diag_sv[i] = sqrt(cov_sv[i,i]) * 1d12
endfor

; plots
vec=indgen(47) + 9

window, 1
plot, l[vec], diag_sv[vec], /yl, xtitle='ell', ytitle='diag_cov [uK^2]', title='Sample v.s. Noise Variance', charsize=2, thick=2
oplot, l[vec], diag_sv[vec], color=!blue, thick=2
oplot, l[vec], diag_n[vec], color=!red, thick=2
legend, ['diag_SV', 'diag_noise'], linestyle=0, thick=2, chars=2, colors=[!blue, !red], position=[1700, 60]
xyouts, 1700, 10, 'crossover at ell ~ 2900', chars=2

err = tvread(/png,/nodialog,filename=fdir+'svnv_0801')

window, 2
plot, l[vec], (diag_sv / diag_n)[vec]


stop

END


;;;;;;;;;;;;;;;;;;;;;;;
;;; 2) Pipeline test2
;;;;;;;;;;;;;;;;;;;;;;;
PRO pipe2
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
file = edir+'end_ra3h30dec-60_08p2_kweight.sav'

; get pipe2
restore, file
;calib = 0.825
calib = 0.76
dl_p2 = spectrum*(1d12)^2.* (calib)
cov_p2 = cov

; get final combined
restore, edir+'run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all = dl_all*1d12
cov_all = cov
wf2use = wf_all_sim

;;; get dl_th
; readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_uK2
; dl_uK2_2 = dl_uK2[50:3998]
; l_vec_2  = l_vec[50:3998]
; dl_th_2 = dl_uK2_2 # wf2use
theoryspectrumfiles=[ ['/home/cr/cmb_models/wmap7_lcdm_lensedCls_extended.dat',$
                       '/home/cr/paramfits/cosmomc.s10/data/dl_ksz_sehgal.txt',$
                       '/home/cr/cmb_models/foreground_sim09_150.txt']]

;----------------------------
;Get the theory spectrum

; Plotting
vec=indgen(47) + 9

window, 1
plot, l[vec], dl_all[vec], xr=[500,3000], yr=[10, 5000], /yl
oplot, l[vec], dl_all[vec]
oplot, l[vec], dl_p2[vec], color=!red, linestyle=3

window, 2
plot, l[vec], (dl_th_2/dl_p2)[vec]
stop
END




;;;;;;;;;;;;;;;;;;
; Procedure to make over-plot
;;;;;;;;;;;;;;;;;;
PRO pl4
;----------------------
; Get the data
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

; wmap
;readcol,pdir+'wmap7cmb_bestfit_dl_k11.txt',lth_cmb,dlth_cmb

uk2mk = 1d-6

; S12
restore,file
dl_all_lps12 = dl_all*(1d12)
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
dl_all_k11      = dl_all*1d12*cal_rk^2.
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
l_act   = l_act[3:39]
dl_act  = dl_act[3:39]
ddl_act = ddl_act[3:39]
cl4_act = dl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk
dcl4_act = ddl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk

; load theory
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)
cl4_th = dl_th *(l_vec^4.) / (l_vec*(l_vec+1)) * uk2mk

;------------------
; PLOT
vec=indgen(47) + 9
color_act = 250 ;!red
color_k11 = 105 ; !cyan
color_this_work = 60 ;!blue
thick = 2

window, 0
plot, l_lps12, dl_all_lps12,ps=3, /yl, xr=[500,3000]
oplot,l_vec,dl_th,color=!red
oplot,l_lps12,dl_all_lps12,ps=3,color=color_this_work
errplot,l_lps12,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work

window, 1
plot, l_lps12, cl4_lps12,ps=3, xr=[500,3000]
oplot,l_vec,cl4_th,color=!green, thick=thick
errplot,l_act,cl4_act-dcl4_act,cl4_act+dcl4_act,color=color_act,thick=thick
errplot,l_k11,(cl4_k11-dcl4_k11),(cl4_k11+dcl4_k11),thick=thick,color=color_k11
errplot,l_lps12,(cl4_lps12-dcl4_lps12),(cl4_lps12+dcl4_lps12),thick=thick,color=color_this_work


stop

END



;;;;;;;;;;;;;;;;;;
; l^2, l^2
PRO plot_dl_all_0

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

; wmap
;readcol,pdir+'wmap7cmb_bestfit_dl_k11.txt',lth_cmb,dlth_cmb

; S12
restore,file
dl_all_lps12 = dl_all*(1d12)
diag_nobeam_lps12 = diag_nobeam*(1d12)
l_lps12 = l

istart=9
istop=55
dl_all_lps12 = dl_all_lps12[istart:istop]
diag_nobeam_lps12 = diag_nobeam_lps12[istart:istop]
l_lps12 = l_lps12[istart:istop]

color_this_work = 60 ; blue

; load  K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
cal_rk = 1;0.76
dl_all_k11      = dl_all*1d12*cal_rk^2.
diag_nobeam_k11 = diag_nobeam*1d12*cal_rk^2.
l_k11 = l

dl_all_k11 = dl_all_k11[istart:istop]
diag_nobeam_k11 = diag_nobeam_k11[istart:istop]
l_k11 = l_k11[istart:istop]

sf=1.

; load act (das et al)
readcol,pdir+'act_150_use.txt',l_act,dl_act,ddl_act
l_act   = l_act[3:39]
dl_act  = dl_act[3:39]
ddl_act = ddl_act[3:39]


;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 8
csize = xb/7.
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.45

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.11
ddx=.44
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1700 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.4
thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_th = pscolors[0]        ; black
color_act = pscolors[2]       ; !red
color_k11 = pscolors[5]       ;!cyan
color_this_work = pscolors[8] ; blue

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]


for i=0,1 do begin
    !x.window=[x0+i*ddx,x0+(i+1)*ddx]
    !x.crange=[200, 3400]
    !y.window=[y0,y0+ddy]
    !y.crange=[ALOG10(30),ALOG10(4000)]

    axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
    axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
    
    if i eq 0 then begin
        axis,yaxis=0,font=0,ystyle=1,/yl,/save,ycharsize=csize,ythick=2 
    endif else $
      axis,yaxis=0,ytickname=empty,ystyle=1,/yl,/save,ythick=2

    axis,yaxis=1,ytickname=empty,ystyle=1,/yl,ythick=2

    ; Plot #1
    if i eq 0 then begin
        oplot,[0,0],[0,0],color=3

        oplot,l_lps12/sf,dl_all_lps12,ps=3,color=color_this_work
        errplot,l_lps12/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work

        ;xyouts,-320,200,alignment=0.5,orientation=90,'!12l(l+1)C!D!12l!3!N/2!7p!X  !N',charsize=1.2
        ;xyouts,-320,930,alignment=0.5,orientation=90,'(!7l!XK!U2!N)',charsize=1.2

        xtitle='!12l!X!N'
        ;ytitle='!12l(l+1)C!D!12l!3!N/2!7p!X!N (!7l!X!NK!U2!X!N)'
        ytitle='!12l(l+1)C!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'

        xyouts, -320, 300, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
        xyouts, 1800, 13, charsize=1.4, charthick=3, font=-1, xtitle

        xyouts,xtxt,ytxt*dytxt^2.,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
        
    ; Plot #2
    endif else begin
        oplot,[0,0],[0,0],color=3
        errplot,l_act,dl_act-ddl_act,dl_act+ddl_act,color=color_act,thick=thick
        errplot,l_k11/sf,(dl_all_k11-diag_nobeam_k11),(dl_all_k11+diag_nobeam_k11),thick=thick,color=color_k11
        errplot,l_lps12/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work

        xyouts, 1800, 13, charsize=1.4, charthick=3, xtitle

        xyouts,xtxt,ytxt,         'ACT',chars=bigspt,charthick=thickspt,font=lfont,color=color_act
        xyouts,xtxt,ytxt*dytxt,   'SPT, K11',chars=bigspt, charthick=thickspt,font=lfont,color=color_k11
        xyouts,xtxt,ytxt*dytxt^2.,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
    endelse

endfor
;*******************************
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/dl_all0.pdf'

stop
END



;;;;;;;;;
; l^2, l^4
;;;;;;;;;
PRO plot_dl_all_1

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
uk2mk = 1d-6

; S12
restore,file
dl_all_lps12 = dl_all*(1d12)
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
dl_all_k11      = dl_all*1d12*cal_rk^2.
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
l_act   = l_act[3:39]
dl_act  = dl_act[3:39]
ddl_act = ddl_act[3:39]
cl4_act = dl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk
dcl4_act = ddl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk

; load theory
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)
cl4_th = dl_th *(l_vec^4.) / (l_vec*(l_vec+1)) * uk2mk


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
pthick=5 ; points
xtxt=1700 ; xyouts label
ytxt=2700
dytxt = 0.63
bigspt = 1.4
thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_th = pscolors[0]        ; black
color_act = pscolors[2]       ; !red
color_k11 = pscolors[5]       ;!cyan
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

    axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
    axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
    
    if i eq 0 then begin
        !y.crange=[ALOG10(30),ALOG10(4000)]
        axis,yaxis=0,font=0,ystyle=1,/yl,/save,ycharsize=csize,ythick=2 
        axis,yaxis=1,ytickname=empty,ystyle=1,/yl,ythick=2
    endif else begin
        !y.crange=[0, 2000]
        axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 ,yl=0
        axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2, yl=0
    endelse


    ; Plot #1
    if i eq 0 then begin
        oplot,[0,0],[0,0],color=pscolors[0]

        oplot,l_lps12/sf,dl_all_lps12,ps=3,color=color_this_work
        errplot,l_lps12/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work

        xtitle='!12l!X!N'

        ytitle='!12l(l+1)C!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'
        xyouts, 0, 300, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
        xyouts, 1800, 13, charsize=1.4, charthick=3, font=-1, xtitle

        xyouts,xtxt,ytxt*dytxt^2.,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
        
    ; Plot #2
    endif else begin
        oplot,[0,0],[0,0],color=3
        oplot,l_vec[450:3200],cl4_th[450:3200],color=color_th, thick=thick
        errplot,l_act,cl4_act-dcl4_act,cl4_act+dcl4_act,thick=pthick,color=color_act
        errplot,l_k11,(cl4_k11-dcl4_k11),(cl4_k11+dcl4_k11),thick=pthick,color=color_k11
        errplot,l_lps12,(cl4_lps12-dcl4_lps12),(cl4_lps12+dcl4_lps12),thick=pthick,color=color_this_work

        ytitle='!12l!U4!X!NC!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'
        xyouts, 50, 1000, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
        xyouts, 1800, -330, charsize=1.4, charthick=3, xtitle

        xyouts,xtxt,1800,         'ACT',chars=bigspt,charthick=thickspt,font=lfont,color=color_act
        xyouts,xtxt,1650,   'SPT, K11',chars=bigspt, charthick=thickspt,font=lfont,color=color_k11
        xyouts,xtxt,1500,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
    endelse

endfor
;*******************************
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/dl_all1.pdf'

stop
END



;;;;;;;;;
; l^4, l^4
;;;;;;;;;
PRO plot_dl_all_2

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
uk2mk = 1d-6

; S12
restore,file
dl_all_lps12 = dl_all*(1d12)
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
dl_all_k11      = dl_all*1d12*cal_rk^2.
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
l_act   = l_act[3:39]
dl_act  = dl_act[3:39]
ddl_act = ddl_act[3:39]
cl4_act = dl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk
dcl4_act = ddl_act *(l_act^4.) / (l_act*(l_act+1)) * uk2mk

; load theory
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)
cl4_th = dl_th *(l_vec^4.) / (l_vec*(l_vec+1)) * uk2mk

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 10
csize = 1.2
yb = 4;xb * 0.45

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.11
ddx=.44
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3  ; Lines, text
pthick=5 ; points
xtxt=1700 ; xyouts label
ytxt=2700
dytxt = 0.63
bigspt = 1.4
thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_th = pscolors[0]        ; black
color_act = pscolors[2]       ; !red
color_k11 = pscolors[5]       ;!cyan
color_this_work = pscolors[8] ; blue

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]


for i=0,1 do begin
    !x.window=[x0+i*ddx,x0+(i+1)*ddx]
    !x.crange=[450, 3200]
    !y.window=[y0,y0+ddy]
    !y.crange=[0, 2000]

    axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
    axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
    
    if i eq 0 then begin
        axis,yaxis=0,font=0,ystyle=1,yl=0,/save,ycharsize=csize,ythick=2 
    endif else begin
        axis,yaxis=0,font=0,ytickname=empty,ystyle=1,/save,ycharsize=csize,ythick=2 ,yl=0
    endelse

    axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2, yl=0

    ; Plot #1
    if i eq 0 then begin
        oplot,[0,0],[0,0],color=3

;         oplot,l_lps12/sf,dl_all_lps12,ps=3,color=color_this_work
;         errplot,l_lps12/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work
        oplot,l_vec[450:3200],cl4_th[450:3200],color=color_th, thick=thick
        ;oplot,cl4_lps12,dcl4_lps12,ps=3,color=color_this_work
        errplot,l_lps12,(cl4_lps12-dcl4_lps12),(cl4_lps12+dcl4_lps12),thick=pthick,color=color_this_work

        xtitle='!12l!X!N'

        ytitle='!12l!U4!X!NC!Dl!3!N/2!7p!X!N (!7l!6K!3!U2!X!N)'
        xyouts, 50, 1000, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
        xyouts, 1800, -330, charsize=1.4, charthick=3, font=-1, xtitle

        xyouts,xtxt,1500.,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
        
    ; Plot #2
    endif else begin
        oplot,[0,0],[0,0],color=3
        oplot,l_vec[450:3200],cl4_th[450:3200],color=color_th, thick=thick
        errplot,l_act,cl4_act-dcl4_act,cl4_act+dcl4_act,color=color_act,thick=pthick
        errplot,l_k11,(cl4_k11-dcl4_k11),(cl4_k11+dcl4_k11),thick=pthick,color=color_k11
        errplot,l_lps12,(cl4_lps12-dcl4_lps12),(cl4_lps12+dcl4_lps12),thick=pthick,color=color_this_work

        ;xyouts, -50, 1000, alignment=0.5, orientation=90, charsize=1.4, font=-1, charthick=3, ytitle
        xyouts, 1800, -330, charsize=1.4, charthick=3, xtitle

        xyouts,xtxt,1800,         'ACT',chars=bigspt,charthick=thickspt,font=lfont,color=color_act
        xyouts,xtxt,1650,   'SPT, K11',chars=bigspt, charthick=thickspt,font=lfont,color=color_k11
        xyouts,xtxt,1500,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
    endelse

endfor
;*******************************
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/dl_all2.pdf'

stop
END


;;;;;;;;;;;;;;;;;;
; Procedure to make over-plot that is in K11
PRO plot_dl_all_k11
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
xmargin=[8.5,2]
ymargin=[4,2]

readcol,pdir+'wmap7cmb_bestfit_dl_k11.txt',lth_cmb,dlth_cmb


;size=3.9
size=6.
;chars=0.9
chars=1.5
lchars=chars*1.6
lfont=1
filename = 'figs/dl_all.eps'
openplotps,filen=filename,xs=size*1.8,ys=0.8*size,/eps
!p.multi = [0,2,1]
;window, 1, xsize=1000, ysize=400
setcolors2,/sys
thick=4

xtitle='!12l!X!N/1000'
xtitle='!12l!X!N'
ytitle='!12l(l+1)C!D!12l!3!N/2!7p!X'+textoidl(' [\muK^2]')
sf=1.
xr=[200,3400]/sf
xxx = 200. & xr = [650-xxx,3000+xxx]
yr=[30,4e3]
lines_neg = 2
nsmooth=5

; Make the plot
plot,[0],[0],/nodata,xr=xr,yr=yr,/xst,/yst,xtitle=xtitle,ytitle=ytitle, $
  chars=chars,/yl,xtickint=1000/sf,xmargin=xmargin,ymargin=ymargin

; Get the data
restore,file

dl_all_lps12 = dl_all*(1d12)
diag_nobeam_lps12 = diag_nobeam*(1d12)


istart=9
istop=55
dl_all_lps12 = dl_all_lps12[istart:istop]
diag_nobeam_lps12 = diag_nobeam_lps12[istart:istop]
l = l[istart:istop]

color_this_work = !black
;color_this_work = !blue
oplot,l/sf,dl_all_lps12,ps=3,color=color_this_work
errplot,l/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work
;oplot,lth_cmb/sf,dlth_cmb,lines=2

xtxt=2060
ytxt=2000*0.9
dytxt = 0.63
bigspt = 1.05
xyouts,xtxt,ytxt,'SPT,',chars=lchars*bigspt,font=lfont,color=color_this_work
xyouts,xtxt,ytxt*dytxt,'Full Survey',chars=lchars*bigspt,font=lfont,color=color_this_work

;;;;;;;;;;;;;;;
; RIGHT PANEL
;;;;;;;;;;;;;;;
color_quad = !darkgreen
color_acbar = !blue
color_act = !darkorange
color_shiro = !orange


; load quad
readcol,pdir+'quad_lowell.txt',l_quad1,dl_quad1,delta_dl_quad1,dl_lcdm_quad1
readcol,pdir+'quad_highell.txt',l_quad2,dl_quad2,delta_dl_quad2,dl_lcdm_quad2
l_quad = [l_quad1, l_quad2]
dl_quad = [dl_quad1, dl_quad2]
ddl_quad = [delta_dl_quad1, delta_dl_quad2]

; load acbar
readcol,pdir+'cl_dc1_v3.dat',iband,dl_acbar,delta_dl_acbar,log_normal_offset
ledge=[[350,550],[550,650],[650,730],[730,790],[790,850],[850,910],[910,970],[970,1030],[1030,1090],[1090,1150],[1150,1210],[1210,1270],[1270,1330],[1330,1390],[1390,1450],[1450,1510],[1510,1570],[1570,1650],[1650,1750],[1750,1850],[1850,1950],[1950,2100],[2100,2300],[2300,2500],[2500,3000]]
llo=reform(ledge[0,*])
lhi=reform(ledge[1,*])
l_acbar=0.5*(llo+lhi)
ddl_acbar = delta_dl_acbar

; load act (das et al)
readcol,pdir+'act_150_use.txt',l_act,dl_act,ddl_act

; load spt shirokoff et al
readcol,pdir+'shiro.txt',llo_shiro,lhi_shiro,l_shiro,dl150,ddl150,dl150220,ddl150220,dl220,ddl220
dl_shiro = dl150
ddl_shiro = ddl150

; load  K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
dl_all_k11      = dl_all*1d12
diag_nobeam_k11 = diag_nobeam*1d12

xtitle='!12l!X!N/1000'
xtitle='!12l!X!N'
ytitle='!12l(l+1)C!D!12l!3!N/2!7p!X'+textoidl(' [\muK^2]')
sf=1.
xr=[200,3400]/sf
xxx = 200. & xr = [650-xxx,3000+xxx]
yr=[30,4e3]
lines_neg = 2
nsmooth=5

; Plot #2
plot,[0],[0],/nodata,xr=xr,yr=yr,/xst,/yst,xtitle=xtitle,ytitle=ytitle, $
  chars=chars,/yl,xtickint=1000/sf,xmargin=xmargin,ymargin=ymargin

restore,file

dl_all_lps12 = dl_all*(1d12)
diag_nobeam_lps12  = diag_nobeam(1d12)

istart=9
istop=55
dl_all_lps12 = dl_all_lps12[istart:istop]
diag_nobeam_lps12 = diag_nobeam_lps12[istart:istop]
l = l[istart:istop]

errplot,l_acbar,dl_acbar-ddl_acbar,dl_acbar+ddl_acbar,color=color_acbar,thick=thick
errplot,l_quad,dl_quad-ddl_quad,dl_quad+ddl_quad,color=color_quad,thick=thick
errplot,l_act,dl_act-ddl_act,dl_act+ddl_act,color=color_act,thick=thick
errplot,l_shiro,dl_shiro-ddl_shiro,dl_shiro+ddl_shiro,color=color_shiro,thick=thick

xyouts,xtxt,ytxt*dytxt^0.,'ACBAR',chars=lchars,font=lfont,color=color_acbar
xyouts,xtxt,ytxt*dytxt^1.,'QUaD',chars=lchars,font=lfont,color=color_quad
xyouts,xtxt,ytxt*dytxt^2.,'ACT',chars=lchars,font=lfont,color=color_act
xyouts,xtxt,ytxt*dytxt^3.,'SPT, S10',chars=lchars,font=lfont,color=color_shiro

closeps
spawn,'convert '+filename+' figs/dl_all_k11.pdf'

stop
END




;;;;;
PRO ptest
WINDOW, /PIXMAP & WDELETE
DEVICE, BYPASS_TRANSLATION=0
device, retain=2
device, true_color=24
device, decomp=0

     names = [ 'black',  $
               'red',    $
               'orange', $
               'green',  $
               'forest', $
               'yellow', $
               'cyan',   $
               'blue',   $
               'magenta',$
               'purple', $
               'gray',   $
               'white']

     r_plot = 255b*[0, 1,  1, 0,  0, 1, 0, 0, 1,  0,  0, 1]$
       +byte([0, 0,  0, 0, 35, 0, 0, 0, 0,153,127, 0])
     g_plot = 255b*[0, 0,  0, 1,  0, 1, 1, 0, 0,  0,  0, 1]$
       +byte([0, 0,127, 0,142, 0, 0, 0, 0, 50,127, 0])
     b_plot = 255b*[0, 0,  0, 0,  0, 0, 1, 1, 1,  0,  0, 1]$
       +byte([0, 0,  0, 0, 35, 0, 0, 0, 0,205,127, 0])

val = color_quan(r_plot, g_plot, b_plot, rr, gg, bb)
tvlct, rr, gg, bb

i=4
filename_stub = '~kstory/lps12/scripts/figs/tmp' & ps_open,filename_stub,/color,xsize=5,ysize=5,/inches,/portrait
plot, indgen(5), color=val[4]
ps_close

stop
END

PRO ptest2
SET_PLOT, 'PS'

;Set the PostScript device to *8* bits per color, not 24:

filename_stub = '~kstory/lps12/scripts/figs/tmp.ps'
DEVICE, FILE=filename_stub, /COLOR, BITS=8


r=[70, 255, 0] & g=[70, 255, 255] & b=[70, 0, 0]
val = color_quan(r, g, b, rr, gg, bb)
tvlct, rr, gg, bb
;Plot, Findgen(11), Color=val[1], Background=val[0], /NoData
Oplot, Findgen(11), Color=val[2]

;plot, indgen(5), color=[[[35]], [[142]], [[35]]], TRUE=3
;TV, 
ps_close
END

PRO ptest3
;!path = !path+':'+expand_path('+/home/kstory/idl/setup/')
;.run '/home/kstory/idl/setup/setup_ps_plotting.pro'
setup_ps_plotting, pscolors=pscolors
filename_stub = '~kstory/lps12/scripts/figs/tmp' & ps_open,filename_stub,/color,xsize=5,ysize=5,/inches,/portrait
plot, indgen(5), color=pscolors[5]
ps_close
END
