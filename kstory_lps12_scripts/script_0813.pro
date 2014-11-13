;;;
; NAME: script_0813
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) constraints on sigma8 and H0
;;;


;;;;;;;;;;;;;;;;;;;;;;
; Alens^0.6
;;;;;;;;;;;;;;;;;;;;;;
PRO h0

; S12-only
dir_s12   = '/data23/hou/lps12/paramfits/chains_final/c1_lcdm_camb_s12tau/chains/'
files_s12 = file_search(dir_s12+'c1_lcdm_camb_s12tau*.txt')
pname_s12 = dir_s12+'c1_lcdm_camb_s12tau.paramnames'

; WMAP-only
dir_w7 = '/data23/hou/lps12/paramfits/chains_final/lcdm_camb_w7/chains/'
files_w7 = file_search(dir_w7+'lcdm_camb_w7*.txt')
pname_w7 = dir_w7+'lcdm_camb_w7.paramnames'

; S12+wmap
dir_w7s12 = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
files_w7s12 = file_search(dir_w7s12+'c4_lcdm_camb_w7s12*.txt')
pname_w7s12 = dir_w7s12+'c4_lcdm_camb_w7s12.paramnames'

subsamp=100000
print, '************ S12 *****************'
plot_like1dname,files_s12,pname_s12,'H0*',subsamp=subsamp,nskip=1000,$
  thick=7,linestyle=ls[0], color=mycolor[0], /oplot
print, '************ WMAP7 *****************'
plot_like1dname,files_w7,pname_w7,'H0*',subsamp=subsamp,nskip=1000,$
  thick=7,linestyle=ls[1], color=mycolor[1], /oplot
print, '************ WMAP7+S12 *****************'
plot_like1dname,files_w7s12,pname_w7s12,'H0*',subsamp=subsamp,nskip=1000,$
  thick=7,linestyle=ls[2], color=mycolor[2], /oplot
END
