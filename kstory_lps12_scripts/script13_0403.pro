;;;
; NAME: script13_0321
;
; NOTES:
;  1) plot s12 and k11 lcdm constraints
;;;


PRO plot_k11

; Set up the environment for ps plotting, get colors
;setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

cdir = '/data23/kstory/lps12/chains/'
chains  = ['c1_lcdm_pico_w7_newKp','c27_lcdm_camb_s12tau_newKp','c2_lcdm_pico_w7s12_newKp']



; S12-only
dir_s12   = cdir+'c27_lcdm_camb_s12tau_newKp/chains/'
files_s12 = file_search(dir_s12+'c27_lcdm_camb_s12tau_newKp*.txt')
pname_s12 = dir_s12+'c27_lcdm_camb_s12tau_newKp.paramnames'

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7_newKp/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7_newKp*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7_newKp.paramnames'

; S12+wmap
dir_w7s12 = cdir+'c2_lcdm_pico_w7s12_newKp/chains/'
files_w7s12 = file_search(dir_w7s12+'c2_lcdm_pico_w7s12_newKp*.txt')
pname_w7s12 = dir_w7s12+'c2_lcdm_pico_w7s12_newKp.paramnames'

; K11+wmap
file_w7k11 = '/home/rkeisler/ps09/cons_chain_baseline.sav'
pname_w7k11 = '/home/rkeisler/ps09/paramnames_alens.txt'

restore, file_w7k11


;*******************************
; colors; use setup_ps_plotting
color_w7s12 = !black;pscolors[0] ; black
color_s12   = !blue;pscolors[8] ; blue
color_w7    = !red;pscolors[2] ; !red
color_k11   = !purple

ls      = [3,2,0] ; [SPT-only, WMAP-only, SPT+WMAP]
mycolor = [color_s12, color_w7, color_w7s12] ; [SPT-only, WMAP-only, SPT+WMAP]



; omegahb plot
nskip=1000
subsamp=3
scale=100

pp = echain[2,nskip:*]*scale
sig = stddev(pp)
binsize=sig/subsamp
h=double(histogram(pp,binsize=binsize,omin=omn,omax=omx))
h/=max(h)
bins = (findgen(n_elements(h))+.5)*binsize +omn

!x.crange=[1.8, 2.92]
plot_like1dname,files_s12,pname_s12,'omegabh2',subsamp=3,nskip=1000,scale=100,$
  thick=4,linestyle=ls[0], color=mycolor[0], /oplot
; plot_like1dname,files_w7,pname_w7,'omegabh2',subsamp=3,nskip=1000,scale=100,$
;   thick=4,linestyle=ls[1], color=mycolor[1], /oplot
; plot_like1dname,files_w7s12,pname_w7s12,'omegabh2',subsamp=3,nskip=1000,scale=100,$
;   thick=4,linestyle=ls[2], color=mycolor[2], /oplot

oplot, bins,h,thick=4,linestyle=4,color=color_k11


stop
END

