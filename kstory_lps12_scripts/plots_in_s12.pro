;;;
; NAME: plots_in_s12.pr
; PURPOSE:
;   Make the plots that are in the paper.
;
; NOTES:
; 1) Test SV averaging
;
; MODIFICATION HISTORY:
;  07/20/2012: (KTS) Created
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



;;;;;;;;;;;;;;;;;
; Figure: alens
;;;;;;;;;;;;;;;;;
PRO p1_alens
; Setup from CR idl_startup
WINDOW, /PIXMAP & WDELETE
DEVICE, BYPASS_TRANSLATION=0
device, retain=2
device, true_color=24
device, decomp=1
!p.charsize=1.
!p.color = !black
!p.background = !white
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]


;----------------------
; Get the data

; S12
dir = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/'
files = file_search(dir+'c24_lcdm_alens_camb_w7s12*.txt')
pname=dir+'c24_lcdm_alens_camb_w7s12.paramnames'
;n=n_elements(files)

; WMAP
;;;;


;*******************************
; colors
color_this_work = 60 ; blue
color_wmap = !red

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait
loadct,39,ncolors=12

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
lfont = 0
ls = [0,2]



;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[-1, 5]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, -2.1, 0.45, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 1.8, -0.16, charsize=1.5, charthick=3, font=0, '!3A!DL!X!N'

oplot,[0,0],[0,0],color=3

; S12
plot_like1dname,files,pname,'alens',subsamp=3,thick=7,linestyle=ls[0], color=color_this_work, /oplot
xyouts,2,0.8,'SPT+WMAP',charsize=csize*1.3, color=color_this_work,font=lfont

; WMAP
i=1 ;& plot_like1dname,files[i],pname[i],'alens',subsamp=3,thick=4,linestyle=ls[i], color=cc[i], /oplot
xyouts,2,0.67,'WMAP',charsize=csize*1.3, color=color_wmap,font=lfont

; vertical line
oplot, [1.,1.], [0,1], linestyle=2, thick=5

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/alens_0726.pdf'
stop
END



;;;;;;;;;;;;;;;;;
; Figure: r_ns
; Now this figure is being made by Brent
;;;;;;;;;;;;;;;;;
PRO p2_r_ns

dir = '/data23/hou/lps12/paramfits/chains_final/c35_lcdm_r_camb_w7s12_BAO_SNe/chains/'
files = file_search(dir+'c35_lcdm_r_camb_w7s12_BAO_SNe*.txt')
pname=dir+'c35_lcdm_r_camb_w7s12_BAO_SNe.paramnames'


cc=[!blue, !red]
ls=[0,2]
xb = 7
csize = xb/7.
yb = 2.5
csize=2
lfont=1
!p.multi=[0,2,1]

window, 0, xsize=900,ysize=400

plot_like1dname,files,pname,'r',subsamp=3,thick=4,linestyle=ls[0], xr=[0,0.5],yr=[0,1.],/xstyle,xtitle='r'
i=0 & plot_like1dname,files,pname,'r',subsamp=3,thick=4,linestyle=ls[i], color=cc[i], /oplot
xyouts,0.3,0.85,'SPT+WMAP',charsize=csize*1.3, color=cc[0],font=lfont
xyouts,0.3,0.75,'WMAP',charsize=csize*1.3, color=cc[1],font=lfont

plot_like2dname,files,pname,'ns','r',xtitle='n_s',ytitle=textoidl('r');,icolors=[3,8,10] ;,smooth=0.15,priorname=['czero_ksz'],priormin=[0],priormax=[100]
xyouts,0.94,0.35,'SPT+WMAP',charsize=csize*1.3, color=cc[0],font=lfont


fdir='/home/kstory/lps12/scripts/figs/'
ff=fdir+'r_ns_0724' & fname=ff
;err = tvread(/png,/nodialog,filename=fname)
spawn, 'convert '+fname+' '+ff+'.pdf '
stop

END

;;;;;;;;;;;;;;;;;
; Table 2, bandpowers
;;;;;;;;;;;;;;;;;
PRO tab_bandpower
restore, '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

; calculate l_effective
leff = fltarr(58)
l2 = fltarr(58)
for i=0, 57 do leff[i] = total(wf_all[*,i]*l_wf) ;/ total(wf_all[*,i])

istart=9
istop=55
dl_all_lps12 = dl_all[istart:istop] * 1d12
diag_nobeam_lps12 = diag_nobeam[istart:istop] * 1d12
l = l[istart:istop]
leff = leff[istart:istop]

nbin = n_elements(l)
n2bin = nbin/2
for i=0, n2bin-1 do begin
    i2 = i+n2bin+1
    print, strtrim(string(round(l[i])-24),2), ' - ', strtrim(string(round(l[i])+25),2), ' & ', round(leff[i]), ' & ', dl_all_lps12[i], ' & ', diag_nobeam_lps12[i], ' & ', $
      strtrim(string(round(l[i2])-24),2), ' - ', strtrim(string(round(l[i2])+25),2), ' & ', round(leff[i2]), ' & ', dl_all_lps12[i2], ' & ', diag_nobeam_lps12[i2], ' \\', $
      format = '(A4, A3, A5, A3, i6, A3, f6.1, A3, f4.1, A3,   A6, A3, A5, A3, i6, A3, f6.1, A3, f4.1, A3)'

endfor
i=n2bin
print, strtrim(string(round(l[i])-24),2), ' - ', strtrim(string(round(l[i])+25),2), ' & ', round(leff[i]), ' & ', dl_all_lps12[i], ' & ', sigfig(diag_nobeam_lps12[i], 3), ' & & & \\', $
  format = '(A4, A3, A5, A3, i6, A3, f6.1, A3, f4.1, A9)'

stop
END
