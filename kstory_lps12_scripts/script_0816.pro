;;;
; NAME: script_0816
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) n_s constraints
; 2) make_ini, for calculating best-fit chisq
;;;

PRO ns

;;; LCDM
; cmb
dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c4_lcdm_camb_w7s12*.txt')
pname_cmb = dir_cmb+'c4_lcdm_camb_w7s12.paramnames'

;;; LCDM+r
; ; cmb
; dir_cmb = '/data23/hou/lps12/paramfits/chains_final/c26_lcdm_r_camb_w7s12/chains/'
; files_cmb = file_search(dir_cmb+'c26_lcdm_r_camb_w7s12_*.txt')
; pname_cmb = dir_cmb+'c26_lcdm_r_camb_w7s12.paramnames'

; ; CMB+BAO
; dir_cmb_bao = '/data23/hou/lps12/paramfits/chains_final/c28_lcdm_r_camb_w7s12_BAO/chains/'
; files_cmb_bao = file_search(dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO_*.txt')
; pname_cmb_bao = dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO.paramnames'

subsamp=10000
nskip=5000
stale=1
plot_like1dname,files_cmb,pname_cmb,'ns',subsamp=subsamp,nskip=nskip,scale=scale,/cdf,/stopit

;; Run the following on the command line when the program stops
;wh = where(bins ge 1.0)
;print, ff[wh[0]]
;print, gauss_cvf(ff[wh[0]])

stop
END



;;;
; Make ini file for best-fit LCDM model
;;;
PRO make_ML_ini, chain

dir = '/data23/hou/lps12/paramfits/chains_final/'+chain+'/chains/'
files = file_search(dir+chain+'*.txt')
pname = dir+chain+'.paramnames'

out_dir = '/home/kstory/paramfits/cosmomc.lps12/scripts/'
fout = out_dir + 'params_'+chain+'_ML.ini'

readcol,pname,names,format='a'

;------------
; Get the ML values of all params
a = read_ascii(files[0])
a=a.field01[*,*]

nparam = n_elements(a[*,0])
nf = n_elements(files)
final = fltarr(nparam, 1)

; read all files
for i=1, nf-1 do begin
    print, 'Read file ', i
    b = read_ascii(files[i])
    b = b.field01[*,*]
    final = [ [final], [b] ]
endfor
final = final[*,1:*]

like = final[1,*]
ind = (where(like eq min(like)))[0]
ml = final[*,ind]
;print,ml[2:17]
;print,ml[0:1]
print,ml

; Make ml range
ml_min = ml - abs(ml/100.)
ml_max = ml + abs(ml/100.)


;------------
; write output file
n_param = n_elements(names)

get_lun, lun1
openw,lun1,fout

for i=0, n_param-1 do begin
    if (names[i] eq 'A_Pois') then begin
        printf,lun1,'param['+names[i]+'] = ' + $
          strtrim(string(ml[i+2]),2) + ' ' + $
          strtrim(string(ml_min[i+2]),2) + ' ' + $
          strtrim(string(ml_max[i+2]),2) + ' ' + $
          '0.1 0.1'
    endif else begin        
        printf,lun1,'param['+names[i]+'] = ' + $
          strtrim(string(ml[i+2]),2) + ' ' + $
          strtrim(string(ml_min[i+2]),2) + ' ' + $
          strtrim(string(ml_max[i+2]),2) + ' ' + $
          '0.000 0.000'
    endelse
endfor

close, lun1
free_lun,lun1

stop
END




;;;
; make full chain
PRO make_cons_chain



END
