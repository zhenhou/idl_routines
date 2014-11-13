;;;
; NAME: script13_0501.pro
;
; NOTES:
;  1) Check ptscor 2500d spectrum
;;;

;;;;;;;;;;;;;;;;;;;
;
; 
;
;;;;;;;;;;;;;;;;;;;
PRO get_k11_ml
readcol, '/home/rkeisler/ps09/paramnames_alens.txt', idx,name,c,format='I,A,A'
restore, '/home/rkeisler/ps09/cons_chain_baseline.sav'

paramset = ['omegabh2','omegadmh2','theta','tau','ns','logA','czero_tsz','czero_dg_po','czero_dg_cl']
n = n_elements(paramset)

; get ML index
like = echain[1,*]
ind = (where(like eq min(like)))[0]

for i=0, n-1 do begin
    wh = where(name eq paramset[i])
    ii = idx[wh]
    print, paramset[i], echain[ii,ind]
endfor
stop

END
