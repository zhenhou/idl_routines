;;;
; NAME: script_0814
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) constraints on sigma8 and H0
;;;


;;;;;;;;;;;;;;;;;;;;;;
; H0
;;;;;;;;;;;;;;;;;;;;;;
PRO h0

name = 'H0*'
subsamp=100000
nskip = 2000
scale=1.

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

print, '************ S12 *****************'
plot_like1dname,files_s12,pname_s12,name,subsamp=subsamp,nskip=nskip,scale=scale
print, '' & print, '************ WMAP7 *****************'
plot_like1dname,files_w7,pname_w7,name,subsamp=subsamp,nskip=nskip,scale=scale
print, '' & print, '************ WMAP7+S12 *****************'
plot_like1dname,files_w7s12,pname_w7s12,name,subsamp=subsamp,nskip=nskip,scale=scale
stop
END


;;;;;;;;;;;;;;;;;;;;;;
; sigma8
;;;;;;;;;;;;;;;;;;;;;;
PRO sigma8

name = 'sigma8*'
subsamp=100000
nskip = 2000
scale=1.

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

print, '************ S12 *****************'
plot_like1dname,files_s12,pname_s12,name,subsamp=subsamp,nskip=nskip,scale=scale
print, '' & print, '************ WMAP7 *****************'
plot_like1dname,files_w7,pname_w7,name,subsamp=subsamp,nskip=nskip,scale=scale
print, '' & print, '************ WMAP7+S12 *****************'
plot_like1dname,files_w7s12,pname_w7s12,name,subsamp=subsamp,nskip=nskip,scale=scale
stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Power-law [LCDM]
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm, param=param
;name = 'ns'
subsamp=100000
nskip = 2000
scale=1.

; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c4_lcdm_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c4_lcdm_camb_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = '/data23/kstory/lps12/chains/c64_lcdm_camb_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c64_lcdm_camb_w7s12_H0*.txt')
pname_cmb_h0 = dir_cmb_h0+'c64_lcdm_camb_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c11_lcdm_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c11_lcdm_camb_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c11_lcdm_camb_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = '/data23/hou/lps12/paramfits/chains_final/c15_lcdm_camb_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c15_lcdm_camb_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c15_lcdm_camb_w7s12_BAO_H0.paramnames'

print, '' & print, '************ CMB *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
print, '' & print, '************ CMB+H0 *****************'
plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
print, '' & print, '************ CMB+H0+BAO *****************'
plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
print, 'CMB+H0: ', results_cmb_h0.median, results_cmb_h0.err_plus, results_cmb_h0.err_minus, results_cmb_h0.avg_err
print, 'CMB+BAO: ', results_cmb_bao.median, results_cmb_bao.err_plus, results_cmb_bao.err_minus, results_cmb_bao.avg_err
print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END


;;;;;;;;;;;;;;;;;;;;;;
; LCDM parameters
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_table3, param=param, scale=scale
;name = 'ns'
subsamp=100000
nskip = 2000
if n_elements(scale) eq 0 then scale=1.

; S12-only
dir_spt   = '/data23/hou/lps12/paramfits/chains_final/c1_lcdm_camb_s12tau/chains/'
files_spt = file_search(dir_spt+'c1_lcdm_camb_s12tau*.txt')
pname_spt = dir_spt+'c1_lcdm_camb_s12tau.paramnames'

; WMAP-only
dir_w7 = '/data23/hou/lps12/paramfits/chains_final/lcdm_camb_w7/chains/'
files_w7 = file_search(dir_w7+'lcdm_camb_w7*.txt')
pname_w7 = dir_w7+'lcdm_camb_w7.paramnames'

; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c4_lcdm_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c4_lcdm_camb_w7s12.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = '/data23/hou/lps12/paramfits/chains_final/c15_lcdm_camb_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c15_lcdm_camb_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c15_lcdm_camb_w7s12_BAO_H0.paramnames'

print, '' & print, '************ SPT *****************'
plot_like1dname,files_spt,pname_spt,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_spt
print, '' & print, '************ WMAP *****************'
plot_like1dname,files_w7,pname_w7,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_w7
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
print, '' & print, '************ CMB+H0+BAO *****************'
plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'param = ', param
print, 'WMAP: ', results_w7.median, results_w7.err_plus, results_w7.err_minus, results_w7.avg_err
print, 'SPT: ', results_spt.median, results_spt.err_plus, results_spt.err_minus, results_spt.avg_err
print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END



;;;;;;;;;;;;;;;;;;;;;;
; LCDM parameters
;;;;;;;;;;;;;;;;;;;;;;
PRO ml_lcdm_table3
; S12-only
dir_spt   = '/data23/hou/lps12/paramfits/chains_final/c1_lcdm_camb_s12tau/chains/'
files_spt = file_search(dir_spt+'c1_lcdm_camb_s12tau*.txt')
pname_spt = dir_spt+'c1_lcdm_camb_s12tau.paramnames'

; WMAP-only
dir_w7 = '/data23/hou/lps12/paramfits/chains_final/lcdm_camb_w7/chains/'
files_w7 = file_search(dir_w7+'lcdm_camb_w7*.txt')
pname_w7 = dir_w7+'lcdm_camb_w7.paramnames'

; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c4_lcdm_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c4_lcdm_camb_w7s12.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = '/data23/hou/lps12/paramfits/chains_final/c15_lcdm_camb_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c15_lcdm_camb_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c15_lcdm_camb_w7s12_BAO_H0.paramnames'

print, '' & print, '************ WMAP *****************'
print_ml_chain_all, files_w7
print, '' & print, '************ SPT *****************'
print_ml_chain_all, files_spt
print, '' & print, '************ CMB *****************'
print_ml_chain_all, files_cmb
print, '' & print, '************ CMB+H0+BAO *****************'
print_ml_chain_all, files_cmb_h0_bao
END



;;;;;;;;;;;;;;;;;;;;;;
; Tensor [LCDM+r]
;   param = 'ns' or 'r'
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_r, param=param
; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c26_lcdm_r_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c26_lcdm_r_camb_w7s12_*.txt')
pname_cmb = dir_cmb+'c26_lcdm_r_camb_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = '/data23/hou/lps12/paramfits/chains_final/c54_lcdm_r_camb_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c54_lcdm_r_camb_w7s12_H0_*.txt')
pname_cmb_h0 = dir_cmb_h0+'c54_lcdm_r_camb_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c28_lcdm_r_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO_*.txt')
pname_cmb_bao = dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = '/data23/hou/lps12/paramfits/chains_final/c35_lcdm_r_camb_w7s12_BAO_SNe/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c35_lcdm_r_camb_w7s12_BAO_SNe_*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c35_lcdm_r_camb_w7s12_BAO_SNe.paramnames'

;param = 'ns'
subsamp=100000
nskip = 2000
scale=1.

print, '' & print, '************ CMB *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
print, '' & print, '************ CMB+H0 *****************'
plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
print, '' & print, '************ CMB+H0+BAO *****************'
plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
print, 'CMB+H0: ', results_cmb_h0.median, results_cmb_h0.err_plus, results_cmb_h0.err_minus, results_cmb_h0.avg_err
print, 'CMB+BAO: ', results_cmb_bao.median, results_cmb_bao.err_plus, results_cmb_bao.err_minus, results_cmb_bao.avg_err
print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END


;;;;;;;;;;;;;;;;;;;;;;
; r, LCDM+r
;;;;;;;;;;;;;;;;;;;;;;
; PRO r_lcdm_r
; name = 'r'
; subsamp=100000
; nskip = 2000
; scale=1.

; ; cmb
; dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c26_lcdm_r_camb_w7s12/chains/'
; files_cmb = file_search(dir_cmb+'c26_lcdm_r_camb_w7s12_*.txt')
; pname_cmb = dir_cmb+'c26_lcdm_r_camb_w7s12.paramnames'

; ; CMB+H0
; dir_cmb_h0 = '/data23/hou/lps12/paramfits/chains_final/c54_lcdm_r_camb_w7s12_H0/chains/'
; files_cmb_h0 = file_search(dir_cmb_h0+'c54_lcdm_r_camb_w7s12_H0_*.txt')
; pname_cmb_h0 = dir_cmb_h0+'c54_lcdm_r_camb_w7s12_H0.paramnames'

; ; CMB+BAO
; dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c28_lcdm_r_camb_w7s12_BAO/chains/'
; files_cmb_bao = file_search(dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO_*.txt')
; pname_cmb_bao = dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO.paramnames'

; ; CMB+H0+BAO
; dir_cmb_h0_bao = '/data23/hou/lps12/paramfits/chains_final/c35_lcdm_r_camb_w7s12_BAO_SNe/chains/'
; files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c35_lcdm_r_camb_w7s12_BAO_SNe_*.txt')
; pname_cmb_h0_bao = dir_cmb_h0_bao+'c35_lcdm_r_camb_w7s12_BAO_SNe.paramnames'

; print, '' & print, '************ CMB *****************'
; plot_like1dname,files_cmb,pname_cmb,name,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
; print, '' & print, '************ CMB+H0 *****************'
; plot_like1dname,files_cmb_h0,pname_cmb_h0,name,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
; print, '' & print, '************ CMB+BAO *****************'
; plot_like1dname,files_cmb_bao,pname_cmb_bao,name,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
; print, '' & print, '************ CMB+H0+BAO *****************'
; plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,name,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao
; stop
; END




;;;;;;;;;;;;;;;;;;;;;;
; Running [LCDM+nrun]
;   param = 'nrun' or 'ns'
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_nrun, param=param
; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c10_lcdm_nrun_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c10_lcdm_nrun_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c10_lcdm_nrun_camb_w7s12.paramnames'

; CMB+H0
; CMB+BAO
dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c13_lcdm_nrun_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c13_lcdm_nrun_camb_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c13_lcdm_nrun_camb_w7s12_BAO.paramnames'
; ; CMB+H0+BAO


subsamp=100000
nskip = 2000
scale=1.

print, '' & print, '************ CMB *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
; print, '' & print, '************ CMB+H0 *****************'
; plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
; print, '' & print, '************ CMB+H0+BAO *****************'
; plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
;print, 'CMB+H0: ', results_cmb_h0.median, results_cmb_h0.err_plus, results_cmb_h0.err_minus, results_cmb_h0.avg_err
print, 'CMB+BAO: ', results_cmb_bao.median, results_cmb_bao.err_plus, results_cmb_bao.err_minus, results_cmb_bao.avg_err
;print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Running + tensor [LCDM+nrun+r
;   param = 'ns', 'r', or 'nrun'
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_nrun_r, param=param
; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c27_lcdm_nrun_r_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c27_lcdm_nrun_r_camb_w7s12_*.txt')
pname_cmb = dir_cmb+'c27_lcdm_nrun_r_camb_w7s12.paramnames'

; CMB+H0

; CMB+BAO
dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c29_lcdm_nrun_r_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c29_lcdm_nrun_r_camb_w7s12_BAO_*.txt')
pname_cmb_bao = dir_cmb_bao+'c29_lcdm_nrun_r_camb_w7s12_BAO.paramnames'

; CMB+H0+BAO

subsamp=100000
nskip = 2000
scale=1.

print, '' & print, '************ CMB *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
; print, '' & print, '************ CMB+H0 *****************'
; plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
; print, '' & print, '************ CMB+H0+BAO *****************'
; plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
;print, 'CMB+H0: ', results_cmb_h0.median, results_cmb_h0.err_plus, results_cmb_h0.err_minus, results_cmb_h0.avg_err
print, 'CMB+BAO: ', results_cmb_bao.median, results_cmb_bao.err_plus, results_cmb_bao.err_minus, results_cmb_bao.avg_err
;print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END








;;;;;;;;;;;;;;;;;;;;;;
; ML: LCDM
;;;;;;;;;;;;;;;;;;;;;;
PRO ml_lcdm
; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c4_lcdm_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c4_lcdm_camb_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = '/data23/hou/lps12/paramfits/chains_final/c54_lcdm_r_camb_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c54_lcdm_r_camb_w7s12_H0_*.txt')
pname_cmb_h0 = dir_cmb_h0+'c54_lcdm_r_camb_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c11_lcdm_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c11_lcdm_camb_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c11_lcdm_camb_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = '/data23/hou/lps12/paramfits/chains_final/c15_lcdm_camb_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c15_lcdm_camb_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c15_lcdm_camb_w7s12_BAO_H0.paramnames'

; print, '' & print, '************ CMB *****************'
; print_ml_chain_all, files_cmb
print, '' & print, '************ CMB+H0 *****************'
print_ml_chain_all, files_cmb_h0
; print, '' & print, '************ CMB+BAO *****************'
; print_ml_chain_all, files_cmb_bao
; print, '' & print, '************ CMB+H0+BAO *****************'
; print_ml_chain_all, files_cmb_h0_bao

END


;;;;;;;;;;;;;;;;;;;;;;
; ML: LCDM+r
;;;;;;;;;;;;;;;;;;;;;;
PRO ml_lcdm_r
; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c26_lcdm_r_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c26_lcdm_r_camb_w7s12_*.txt')
pname_cmb = dir_cmb+'c26_lcdm_r_camb_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = '/data23/hou/lps12/paramfits/chains_final/c54_lcdm_r_camb_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c54_lcdm_r_camb_w7s12_H0_*.txt')
pname_cmb_h0 = dir_cmb_h0+'c54_lcdm_r_camb_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c28_lcdm_r_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO_*.txt')
pname_cmb_bao = dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = '/data23/hou/lps12/paramfits/chains_final/c35_lcdm_r_camb_w7s12_BAO_SNe/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c35_lcdm_r_camb_w7s12_BAO_SNe_*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c35_lcdm_r_camb_w7s12_BAO_SNe.paramnames'

print, '' & print, '************ CMB *****************'
print_ml_chain_all, files_cmb
print, '' & print, '************ CMB+H0 *****************'
print_ml_chain_all, files_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
print_ml_chain_all, files_cmb_bao
print, '' & print, '************ CMB+H0+BAO *****************'
print_ml_chain_all, files_cmb_h0_bao

END


;;;;;;;;;;;;;;;;;;;;;;
; ML: LCDM+nrun
;;;;;;;;;;;;;;;;;;;;;;
PRO ml_lcdm_nrun

; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c10_lcdm_nrun_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c10_lcdm_nrun_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c10_lcdm_nrun_camb_w7s12.paramnames'

; CMB+H0
; CMB+BAO
dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c13_lcdm_nrun_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c13_lcdm_nrun_camb_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c13_lcdm_nrun_camb_w7s12_BAO.paramnames'
; ; CMB+H0+BAO

print, '' & print, '************ CMB *****************'
print_ml_chain_all, files_cmb
;print, '' & print, '************ CMB+H0 *****************'
;print_ml_chain_all, files_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
print_ml_chain_all, files_cmb_bao
;print, '' & print, '************ CMB+H0+BAO *****************'
;print_ml_chain_all, files_cmb_h0_bao
END
