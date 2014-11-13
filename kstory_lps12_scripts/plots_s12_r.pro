;;;
; NAME: plots_s12_r
; PURPOSE:
;   Make r plot
;
; NOTES:
;
; MODIFICATION HISTORY:
;  09/27/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;
; Figure: r
;;;;;;;;;;;;;;;;;
PRO plot_r
; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors


; Get the chains
cdir = ['/data23/hou/lps12/paramfits/chains_0828/', '/data23/hou/lps12/paramfits/chains_0717/']

chains  = ['lcdm_r_camb_w7','c56_lcdm_r_camb_w7s12','c58_lcdm_r_camb_w7s12_H0','c57_lcdm_r_camb_w7s12_BAO','c59_lcdm_r_camb_w7s12_BAO_H0','c73_lcdm_r_camb_w7_BAO']
nchains = n_elements(chains)
dirs    = cdir[0]+chains+'/chains/'
dirs[0] = cdir[1]+chains[0]+'/chains/'
dirs[5] = cdir[1]+chains[5]+'/post_chains/'
nfiles = n_elements(file_search(dirs[1]+chains[1]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
ndsets = n_elements(chains)

;indicies: 
; 0 - w7
; 1 - w7s12
; 2 - w7s12_H0
; 3 - w7s12_BAO
; 4 - w7s12_BAO_H0
; 5 - w7_BAO

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.2
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/r_1d'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=0.3;1.8 ; xyouts label
ytxt=0.9
dytxt = 0.065
bigspt = 1.
thickspt = 1
lfont = 0
ls = [2, 0, 3, 1, 4, 4] ; w, cmb, cmb_h0, cmb_bao, cmb_h0_bao, w7_bao

;*******************************
; colors
cc = [pscolors[2],$ ; red
      pscolors[0],$ ; black
      pscolors[3],$ ; orange
      pscolors[0],$ ; forest
      pscolors[10],$ ; purple
      pscolors[3]] ; purple

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[0., 0.6]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, -0.08, 0.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 0.3, -0.16, charsize=1.5, charthick=3, font=0, 'r'

oplot,[0,0],[0,0],color=3
;subsamp = 5

set = [0,1,5,3]
for j=0, n_elements(set)-1 do begin
    ii = set[j]
    plot_like1dname,files[ii].ff,pnames[ii],'r',subsamp=subsamp,nskip=1000,$
      thick=7,linestyle=ls[ii], color=cc[ii], /oplot
endfor
;plot_like1dname,files[5].ff,pnames[5],'r',subsamp=subsamp,nskip=1000,thick=7,linestyle=ls[5], color=cc[5], psym=2,/oplot


xyouts,xtxt,ytxt,'WMAP7',charsize=csize,font=lfont,color=cc[0]
oplot,[xtxt-0.1, xtxt-0.005],[ytxt+0.02,ytxt+0.02],thick=7,linestyle=ls[0],color=cc[0]

xyouts,xtxt,ytxt-1*dytxt,'WMAP7+BAO',charsize=csize,font=lfont,color=cc[5]
oplot,[xtxt-0.1, xtxt-0.005],[ytxt-1*dytxt+0.02,ytxt-1*dytxt+0.02],thick=7,linestyle=ls[5],color=cc[5]

pp=1
xyouts,xtxt,ytxt-(1+pp)*dytxt,'SPT+WMAP7',charsize=csize,font=lfont,color=cc[1]
oplot,[xtxt-0.1, xtxt-0.005],[ytxt-(1+pp)*dytxt+0.02,ytxt-(1+pp)*dytxt+0.02],thick=7,linestyle=ls[1],color=cc[1]

; xyouts,xtxt,ytxt-(2+pp)*dytxt,'SPT+WMAP7+H!D0!N',charsize=csize,font=lfont,color=cc[2]
; oplot,[xtxt-0.1, xtxt-0.005],[ytxt-(2+pp)*dytxt+0.02,ytxt-(2+pp)*dytxt+0.02],thick=7,linestyle=ls[2],color=cc[2]

pp=0
xyouts,xtxt,ytxt-(3+pp)*dytxt,'SPT+WMAP7+BAO',charsize=csize,font=lfont,color=cc[3]
oplot,[xtxt-0.1, xtxt-0.005],[ytxt-(3+pp)*dytxt+0.02,ytxt-(3+pp)*dytxt+0.02],thick=7,linestyle=ls[3],color=cc[3]

; xyouts,xtxt,ytxt-4*dytxt,'WMAP7+BAO',charsize=csize,font=lfont,color=cc[4]
; ;xyouts,xtxt+0.085,ytxt-5*dytxt,'H!D0!N+BAO',charsize=csize,font=lfont,color=cc[4]
; oplot,[xtxt-0.1, xtxt-0.005],[ytxt-4*dytxt+0.02,ytxt-4*dytxt+0.02],thick=7,linestyle=ls[4],color=cc[4]



; vertical line
oplot, [1.,1.], [0,1], linestyle=1, thick=3

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp '+filename_stub+'.pdf '+filename_stub+'_0928.pdf'
stop
END

