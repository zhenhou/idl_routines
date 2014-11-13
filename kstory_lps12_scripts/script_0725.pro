;;;
; NAME: script_0725
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Beam test
; 2) calculate leff for the table in S12
;;;

PRO beam_test

;;-------------
; Make the beam function
beamfiles_orig = get_lps12_beams(2009, l_beam, bl_beam)

; Test changing the beam by 5% at ell=3000
fl = -2.128*10^(-5.) *(l_beam - 650) + 1
fn = 1/fl^2.
;;-------------

restore, '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
dl_all_08 = dl_all

restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_08_kweight.sav'
dl_08 = spectrum * (combined_calib)^2.

restore, '/home/kstory/lps12/end2end/end_ra3h30dec-60_08b1_kweight.sav'
combined_calib = 0.825
dl_08b1 = spectrum * (combined_calib)^2.
;diag_08b1 = 

istart = 9
iend = 54
whplot = indgen(iend-istart) + istart

plot, l[whplot], (dl_08b1 / dl_08)[whplot], $
  title='Beam test, ra6h30dec-60', xtitle=textoidl(' '), ytitle='dl_08b1 / dl_orig'
oplot, l[whplot], (dl_08b1 / dl_all_08)[whplot], color=!blue
oplot, l[whplot], (dl_08b1 / dl_08)[whplot], color=!red
oplot, l_beam, fn, color=!green, linestyle=2

stop
END


;----------------------------
; calculating leff for the table
PRO leff
restore, '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

leffa = fltarr(58) ; debugging
leff_nonmz = fltarr(58) ; debugging
leffa_nonmz = fltarr(58) ; debugging

; calculate l_effective
leff = fltarr(58)
l2 = fltarr(58)
;for i=0, 57 do begin
for i=5, 57 do begin
    leff[i] = total(wf_all[*,i]*l_wf) / total(wf_all[*,i])
    leff_nonmz[i] = total(wf_all[*,i]*l_wf)

    wh = where(l_wf gt l[i]-25 and l_wf le l[i]+25)
    leffa[i] = total(wf_all[wh,i]*l_wf[wh]) / total(wf_all[wh,i])

    leffa_nonmz[i] = total(wf_all[wh,i]*l_wf[wh])
endfor

print, "leff ", leff[9:12]
print, "leff_nonmz ", leff_nonmz[9:12]
print, "leffa ", leffa[9:12]
print, "leffa_nonmz ", leffa_nonmz[9:12]
stop


END


;-------------------------------
; Calculate where SPT becomes more sensitive than WMAP
PRO sensitive
readcol,'wmap_binned_tt_spectrum_7yr_v4p1.txt',l_wmap,llo,lhi,dl_wmap,ddl_wmap,ddl_noise,ddl_sv
restore, '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'

ddl = dblarr(58)
for i=0, 57 do begin
    ddl[i] = sqrt(cov_all[i,i])*1d12
    print, l[i], l_wmap[27+i], ddl_wmap[27+i], ddl[i]
endfor
END



;-------------------------------
; test program
PRO ss

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

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 6.5
csize = xb/7.
yb = xb *.88/.9*(.44/.9)

ps_open,'~kstory/lps12/scripts/figs/tmp',/color,xsize=xb,ysize=yb,/inches,/portrait
loadct,39,ncolors=12

;Define the axis size relative to the window size
x0=.11
ddx=.8
y0=.07
ddy=.86

!x.window=[x0, x0+ddx] ;*xb
!x.crange=[0,20]
!y.window=[y0,y0+ddy]           ;*yb
!y.crange=[0,400]

; Make the axes
axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

xyouts,max(!x.crange)*0.9,100,alignment=0.5,'D!S!D3000!R!N!UtSZ  !N',font=0,charsize=1.2
xyouts,max(!x.crange)*0.93,100,alignment=0.5,'(!7l!XK!U2!N)',charsize=1.2
xyouts,5,300,alignment=0.5,orientation=90,'D!S!D3000!R!N!UkSZ  !N',font=0,charsize=1.2
xyouts,5,330,alignment=0.5,orientation=90,'(!7l!XK!U2!N)',charsize=1.2


; Now make the plot
x = findgen(20)
y = x^2.

oplot,[0,0],[0,0],color=3
oplot, x, y
ps_close

stop
END


;-------------------------------
; Make plot for s12
PRO plot_dl_all
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


;----------------------
; Get the data
pdir = '/home/kstory/lps12/scripts/plotting/'
file = '/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
xmargin=[8.5,2]
ymargin=[4,2]

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

; colors
color_this_work = 60 ; blue
color_k11 = !orange
color_act = !green

;*******************************

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 8
csize = xb/7.
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.45

ps_open,'~kstory/lps12/scripts/figs/tmp',/color,xsize=xb,ysize=yb,/inches,/portrait
loadct,39,ncolors=12

;Define the axis size relative to the window size
x0=.11
ddx=.44
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1900 ; xyouts label
ytxt=2000*0.9
dytxt = 0.63
bigspt = 1.6
thickspt = 5
lfont = 1

; Axes lables

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
        ytitle='!12l(l+1)C!D!12l!3!N/2!7p!X'+textoidl(' (\muK^2)')
        xyouts, -320, 300, alignment=0.5, orientation=90, ytitle, charsize=1.2 ; ytitle
        xyouts, 1800, 13, xtitle, charsize=1.2  ; xtitle

        xyouts,xtxt,ytxt,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
        
    ; Plot #2
    endif else begin
        oplot,[0,0],[0,0],color=3

        errplot,l_act,dl_act-ddl_act,dl_act+ddl_act,color=color_act,thick=thick
        errplot,l_k11/sf,(dl_all_k11-diag_nobeam_k11),(dl_all_k11+diag_nobeam_k11),thick=thick,color=color_k11
        errplot,l_lps12/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work

        xtitle='!12l!X!N'
        xyouts, 1800, 13, xtitle, charsize=1.2  ; xtitle

        xyouts,xtxt,ytxt,'ACT',chars=bigspt,charthick=thickspt,font=lfont,color=color_act
        xyouts,xtxt,ytxt*dytxt,'SPT, K11',chars=bigspt, charthick=thickspt,font=lfont,color=color_k11
        xyouts,xtxt,ytxt*dytxt*dytxt,'SPT, Full Survey',chars=bigspt, charthick=thickspt,font=lfont,color=color_this_work
    endelse

endfor
;*******************************

ps_close

stop
END


;;; Keep this for later
PRO spt_other_dl_plot
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
filename = 'figs/tmp.eps'
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
l_lps12 = l

istart=9
istop=55
dl_all_lps12 = dl_all_lps12[istart:istop]
diag_nobeam_lps12 = diag_nobeam_lps12[istart:istop]
l_lps12 = l_lps12[istart:istop]

color_this_work = !black
;color_this_work = !blue
oplot,l_lps12/sf,dl_all_lps12,ps=3,color=color_this_work
errplot,l_lps12/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work
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
;color_quad = !darkgreen
;color_acbar = !blue
color_act = !darkorange
color_k11 = !orange


; ; load quad
; readcol,pdir+'quad_lowell.txt',l_quad1,dl_quad1,delta_dl_quad1,dl_lcdm_quad1
; readcol,pdir+'quad_highell.txt',l_quad2,dl_quad2,delta_dl_quad2,dl_lcdm_quad2
; l_quad = [l_quad1, l_quad2]
; dl_quad = [dl_quad1, dl_quad2]
; ddl_quad = [delta_dl_quad1, delta_dl_quad2]

; ; load acbar
; readcol,pdir+'cl_dc1_v3.dat',iband,dl_acbar,delta_dl_acbar,log_normal_offset
; ledge=[[350,550],[550,650],[650,730],[730,790],[790,850],[850,910],[910,970],[970,1030],[1030,1090],[1090,1150],[1150,1210],[1210,1270],[1270,1330],[1330,1390],[1390,1450],[1450,1510],[1510,1570],[1570,1650],[1650,1750],[1750,1850],[1850,1950],[1950,2100],[2100,2300],[2300,2500],[2500,3000]]
; llo=reform(ledge[0,*])
; lhi=reform(ledge[1,*])
; l_acbar=0.5*(llo+lhi)
; ddl_acbar = delta_dl_acbar

; ; load act (das et al)
; readcol,pdir+'act_150_use.txt',l_act,dl_act,ddl_act

; ; load spt shirokoff et al
; readcol,pdir+'shiro.txt',llo_shiro,lhi_shiro,l_shiro,dl150,ddl150,dl150220,ddl150220,dl220,ddl220
; dl_shiro = dl150
; ddl_shiro = ddl150

; load  K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
cal_rk = 1.;0.76
dl_all_k11      = dl_all*1d12*cal_rk^2.
diag_nobeam_k11 = diag_nobeam*1d12*cal_rk^2.

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
diag_nobeam_lps12  = diag_nobeam*(1d12)

istart=9
istop=55
dl_all_lps12 = dl_all_lps12[istart:istop]
diag_nobeam_lps12 = diag_nobeam_lps12[istart:istop]
dl_all_k11 = dl_all_k11[istart:istop]
diag_nobeam_k11 = diag_nobeam_k11[istart:istop]
l_k11 = l_lps12

;errplot,l_acbar,dl_acbar-ddl_acbar,dl_acbar+ddl_acbar,color=color_acbar,thick=thick
;errplot,l_quad,dl_quad-ddl_quad,dl_quad+ddl_quad,color=color_quad,thick=thick
;errplot,l_act,dl_act-ddl_act,dl_act+ddl_act,color=color_act,thick=thick
;errplot,l_shiro,dl_shiro-ddl_shiro,dl_shiro+ddl_shiro,color=color_shiro,thick=thick
errplot,l_k11/sf,(dl_all_k11-diag_nobeam_k11),(dl_all_k11+diag_nobeam_k11),thick=thick,color=color_k11
errplot,l_lps12/sf,(dl_all_lps12-diag_nobeam_lps12),(dl_all_lps12+diag_nobeam_lps12),thick=thick,color=color_this_work

xyouts,xtxt,ytxt*dytxt^0.,'SPT, S12',chars=lchars,font=lfont,color=color_this_work
xyouts,xtxt,ytxt*dytxt^1.,'SPT, K11',chars=lchars,font=lfont,color=color_k11
;xyouts,xtxt,ytxt*dytxt^0.,'ACBAR',chars=lchars,font=lfont,color=color_acbar
;xyouts,xtxt,ytxt*dytxt^1.,'QUaD',chars=lchars,font=lfont,color=color_quad
;xyouts,xtxt,ytxt*dytxt^2.,'ACT',chars=lchars,font=lfont,color=color_act
;xyouts,xtxt,ytxt*dytxt^3.,'SPT, S10',chars=lchars,font=lfont,color=color_shiro

closeps
spawn,'convert '+filename+' figs/tmp.pdf'

stop
END
