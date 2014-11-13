;;;
; NOTES:
;  1) get omega_c from lcdm+yp
;;;

;--------------------
; degeneracy between r and tau
;--------------------
PRO ss
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

;----------------------
; Get the data

; S12+WMAP
dir_cmb = cdir+'c11_lcdm_yp_pico_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c11_lcdm_yp_pico_w7s12*.txt')
pname_cmb = dir_cmb+'c11_lcdm_yp_pico_w7s12.paramnames'

; WMAP
dir_w = cdir+'c28_lcdm_yp_pico_w7/chains/'
files_w = file_search(dir_w+'c28_lcdm_yp_pico_w7*.txt')
pname_w = dir_w+'c28_lcdm_yp_pico_w7.paramnames'

; S12-only
dir_s12   = cdir+'c52_lcdm_alens_camb_s12tau/chains/'
files_s12 = file_search(dir_s12+'c52_lcdm_alens_camb_s12tau*.txt')
pname_s12 = dir_s12+'c52_lcdm_alens_camb_s12tau.paramnames'

; K11
dir_k11   = '/data17/rkeisler/ps09/chains_20110202/'
files_k11 = file_search(dir_k11+'chain_k10_alens_1.txt')
pname_k11 = dir_k11+'chain_k10_alens.paramnames'

plot_like1dname,files_cmb,pname_cmb,'ommh2',subsamp=1000,nskip=1000
plot_like1dname,files_w,pname_w,'ommh2',subsamp=1000,nskip=1000,color=!red,/oplot



stop
END
