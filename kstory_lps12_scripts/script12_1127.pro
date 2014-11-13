;;;
; NAME: scripts12_1127.pro
;
; NOTES:
; 1) Make a txt file with tabulated plotting data for the A_L plot in S12
;
;;;

;;;;;;;;;;;;;;;;;
; Figure: alens
;;;;;;;;;;;;;;;;;
PRO p1_alens
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; Set up the environment for ps plotting, get colors
;setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

; S12+WMAP
dir_cmb = cdir+'c54_lcdm_alens_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c54_lcdm_alens_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c54_lcdm_alens_camb_w7s12.paramnames'

; WMAP
dir_w = '/data23/hou/lps12/paramfits/chains_0717/lcdm_alens_camb_w7/chains/'
files_w = file_search(dir_w+'lcdm_alens_camb_w7*.txt')
pname_w = dir_w+'lcdm_alens_camb_w7.paramnames'

; S12-only
dir_s12   = cdir+'c52_lcdm_alens_camb_s12tau/chains/'
files_s12 = file_search(dir_s12+'c52_lcdm_alens_camb_s12tau*.txt')
pname_s12 = dir_s12+'c52_lcdm_alens_camb_s12tau.paramnames'

; K11
dir_k11   = '/data17/rkeisler/ps09/chains_20110202/'
files_k11 = file_search(dir_k11+'chain_k10_alens_1.txt')
pname_k11 = dir_k11+'chain_k10_alens.paramnames'


; CMB
plot_like1dname,files_cmb,pname_cmb,'alens',subsamp=8,nskip=1000,/stopit

; WMAP
plot_like1dname,files_w,pname_w,'alens',subsamp=5,nskip=1000,$
  thick=7,linestyle=ls[1], color=color_wmap, /oplot, /stopit

; S12
plot_like1dname,files_s12,pname_s12,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=3, color=color_s12, /oplot, /stopit

;xyouts,xtxt,ytxt,'SPT+WMAP7',charsize=csize, color=color_cmb,font=lfont
;oplot,[xtxt-0.2, xtxt-0.005],[ytxt+0.02,ytxt+0.02],thick=7,linestyle=ls[0],color=color_cmb
xyouts,xtxt,ytxt+dytxt-0.02,'SPT+',charsize=csize, color=color_cmb,font=lfont
xyouts,xtxt-0.03,ytxt,'WMAP7',charsize=csize, color=color_cmb,font=lfont
oplot,[xtxt-0.2, xtxt-0.005],[ytxt+dytxt,ytxt+dytxt],thick=7,linestyle=ls[0],color=color_cmb
xyouts,xtxt,ytxt-dytxt,'WMAP7',charsize=csize, color=color_wmap,font=lfont
oplot,[xtxt-0.2, xtxt-0.005],[ytxt-dytxt+0.02,ytxt-dytxt+0.02],thick=7,linestyle=ls[1],color=color_wmap
xyouts,xtxt,ytxt-2*dytxt,'SPT',charsize=csize, color=color_S12,font=lfont
oplot,[xtxt-0.2, xtxt-0.005],[ytxt-2*dytxt+0.02,ytxt-2*dytxt+0.02],thick=7,linestyle=ls[2],color=color_s12

; ; K11
; plot_like1dname,files_k11,pname_k11,'alens',subsamp=3,nskip=1000,$
;   thick=7,linestyle=3, color=color_s12, /oplot
; xyouts,xtxt,0.35,'K11',charsize=csize, color=color_s12,font=lfont
; oplot,[xtxt-0.3, xtxt-0.001],[0.35+0.02,0.35+0.02],thick=7,linestyle=ls[2],color=color_s12


; vertical line
oplot, [1.,1.], [0,1], linestyle=1, thick=3

;*******************************

; Save the plot
;ps_close

spawn,'epstopdf '+filename_stub+'.ps'
;spawn, 'cp figs/tmp.pdf figs/alens.pdf'
print, 'output: ', filename_stub+'.pdf'
stop
END

