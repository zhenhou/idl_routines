;;;
; NOTES:
;  1) check degeneracy between r and tau
;;;

;--------------------
; degeneracy between r and tau
;--------------------
PRO r_tau
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

set = [0,1]
cc = [[!purple,!green],[!blue,!orange]]
; ii=0
; plot_like2dname,files[ii].ff,pnames[ii],'r','tau',sigma=[1,2],icolors=[cc[0,ii],cc[1,ii] ],smooth=0.01,subsamp=4,yl=0

; ii=1
; plot_like2dname,files[ii].ff,pnames[ii],'r','tau',/overplot,sigma=[1,2],icolors=[cc[0,ii],cc[1,ii] ],smooth=0.01,subsamp=4,yl=0

ii=1
plot_like2dname,files[ii].ff,pnames[ii],'r','tau',sigma=[1,2],icolors=[cc[0,1],cc[1,1] ],smooth=0.01,subsamp=4,yl=0

ii=4
plot_like2dname,files[ii].ff,pnames[ii],'r','tau',/overplot,sigma=[1,2],icolors=[cc[0,0],cc[1,0] ],smooth=0.01,subsamp=4,yl=0




stop
END




;--------------------
; degeneracy between r and tau
;--------------------
PRO r_tau_2
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

filename_stub = '~kstory/lps12/scripts/figs/r_tau'

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
!y.crange=[0, 0.12]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, -0.08, 0.5, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Tau'
xyouts, 0.3, -0.16, charsize=1.5, charthick=3, font=0, 'r'

oplot,[0,0],[0,0],color=3
;subsamp = 5

set = [0,1]
cc = [[5,10],[8,3]]
for j=0, n_elements(set)-1 do begin
    ii = set[j]
    plot_like2dname,files[ii].ff,pnames[ii],'r','tau',/overplot,sigma=[1,2],$
      icolors=[pscolors[cc[0,j]],pscolors[cc[1,j]] ],$
      smooth=0.01,subsamp=4,yl=0 ;,priorname=['czero_ksz'],priormin=[0],priormax=[100]
endfor

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
;spawn, 'cp '+filename_stub+'.pdf '+filename_stub+'_0928.pdf'

stop
END
