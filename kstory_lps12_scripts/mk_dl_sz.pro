;;;
; Script for making SZ template for the s12 data products website
;;;

;;;;;;;;;;
;
; Make template for SZ sources
;
;;;;;;;;;;
PRO mk_dl_sz

nell = 10000
ell = lindgen(10000) + 1

readcol, '/home/hou/paramfits/cosmomc.lps12/scripts/data/dl_shaw_tsz_s10_153ghz.txt', ell, dl_sz, format='i,d'

wh = (where(ell eq 3000))[0 ]
dl_sz /= dl_sz[wh]

script_fname = '/home/kstory/lps12/webpage/dl_sz.txt'
print, 'write output to: ', script_fname

get_lun, lun1
openw,lun1,script_fname
for i=0, n_elements(ell)-1 do begin
    printf,lun1, ell[i], dl_sz[i]
endfor
close, lun1
free_lun,lun1

stop
END

