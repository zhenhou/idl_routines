;;;
; NAME: make_ML_ini.pro
; PURPOSE:
;   Make ini files for ML chains
;
; NOTES:
;
; Modification History
; 08/17/2012: (KTS) Created
; 09/12/2012: (KTS) Use 0828 chains
;;;

;;;
; Make ini file 
;  Example, chain = 'c4_lcdm_camb_w7s12'
;;;
PRO make_ML_ini, chain, chain_dir=chain_dir,keep_derived=keep_derived, stopit=stopit,fout=fout
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

dir = cdir+chain+'/chains/'
if n_elements(chain_dir) ne 0 then dir = chain_dir
files = file_search(dir+chain+'*.txt')
pname = dir+chain+'.paramnames'

out_dir = '/home/kstory/paramfits/cosmomc.lps12/scripts/'
if n_elements(fout) eq 0 then fout = out_dir + 'params_'+chain+'_ML.ini'

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
print,ml

; Make ml range
ml_min = ml - abs(ml/100.)
ml_max = ml + abs(ml/100.)


;------------
; write output file
n_param = n_elements(names)

print, "write to output file: ", fout
get_lun, lun1
openw,lun1,fout

for i=0, n_param-1 do begin
    ; Check for derived parameters
    if ~keyword_set(keep_derived) then $
      if strlen( strsplit(names[i],'*',/extract)) ne strlen(names[i]) then break

    print, names[i]
    if (names[i] eq 'A_Pois') then begin
        printf,lun1,'param['+names[i]+'] = ' + $
          strtrim(string(ml[i+2]),2) + ' ' + $
          strtrim(string(ml_min[i+2]),2) + ' ' + $
          strtrim(string(ml_max[i+2]),2) + ' ' + $
          '0.0000001 0.0000001'
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

if keyword_set(stopit) then stop
END

