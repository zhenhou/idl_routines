;;;
; NAME: script_0808
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) alens^0.6
;;;


;;;;;;;;;;;;;;;;;;;;;;
; Alens^0.6
;;;;;;;;;;;;;;;;;;;;;;
PRO test_alens
dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/' 
files_s12 = [file_search(dir_s12+'c24_lcdm_alens_camb_w7s12*.txt')]  
pname_s12 = dir_s12+'c24_lcdm_alens_camb_w7s12.paramnames'

;plot_like1d,files_s12,14,nskip=5000,power=0.6,results=results
plot_like1dname,files_s12,pname_s12,'alens',power=0.6,nskip=5000,subsamp=10000,results=results

print,' '
print,results.median,' +/- ',results.err

stop
END
