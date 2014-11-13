;;;
; Daily script
;
; Notes:
;  1) Make template for Poisson Point Sounces
;;;




;;;;;;;;;;
;
; Make template for Poisson Point Sounces
;
;;;;;;;;;;
PRO mk_dl_poisson_point_source

nell = 10000

dl_3000 = 19.3
ell = lindgen(10000) + 1

dl_ps = dl_3000 * ell*(ell+1) / (3000*3001.)

print, "dl_3000 = ", dl_ps[where(ell eq 3000)]

script_fname = '/home/kstory/lps12/webpage/dl_poisson_point_source.txt'
print, 'write output to: ', script_fname

get_lun, lun1
openw,lun1,script_fname
for i=0, n_elements(ell)-1 do begin
    printf,lun1, ell[i], dl_ps[i]
endfor
close, lun1
free_lun,lun1

END

