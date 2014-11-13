;;;
; NAME: test.pro
; PURPOSE:
;   test script, keep nothing important here!
;;;


;...................................................................
; the program
pro test, dosave=dosave, dosave_y=dosave_y
compile_opt IDL2, HIDDEN

restore, '/home/rkeisler/ps09/1.37/end_ra21hdec-50_1.37_kweight.sav'
kmask_tmp=weight_2d
nbig = 4320

kmask_tmp = shift(kmask_tmp, 1080, 1080)
small = dblarr(1080,1080)

nn = 1080
for i=0, nn-1 do begin
    small[i,*] = kmask_tmp[(nn/2)+i, (nn/2):(3*nn/2)-1 ]
endfor

; now expand to 4320x4320
kmask_rk = congrid(small, nbig, nbig)
kmask_rk = shift(kmask_rk, nbig/2, nbig/2)

kmask_ks = get_lps12_kweight(5)


stop
end
