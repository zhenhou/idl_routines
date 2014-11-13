;;;
; NAME: script_0524
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Write Dls_theory.txt
;
; MODIFICATION HISTORY:
;  05/24/2012: (KTS) Created
;;;

; write Dls_theory.txt
PRO write_Dls_theory
bdir = '/home/kstory/lps12/cls_theory/'
readcol,bdir+'Cls_theory.txt',l_vec,cl_uK2
outfile = bdir+'Dls_theory.txt'

dl_uK2 = cl_uK2 * ( l_vec*(l_vec+1)/(2*!pi) )

nn = n_elements(l_vec)

; Write out new file
get_lun, lun1
openw, lun1, outfile
for jj=0, nn-1 do begin
    ss = strtrim(string(l_vec[jj]),2) + ' ' + strtrim(string(dl_uK2[jj]),2)
    printf, lun1, ss
endfor
close, lun1
free_lun,lun1

stop
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; end2end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO end2end_0to9
for idx=0, 9 do begin
    ; skip fields that have finished
    if idx eq 0 then idx++

    print, 'run end2end_04, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_04, idx, run='04', /resume, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END

PRO end2end_10to19
for idx=10, 19 do begin
    ; skip fields that have finished
    if idx eq 10 then idx++
    if idx eq 11 then idx++

    print, 'run end2end_04, field ', idx
    t0 = systime(0, /seconds)

    try_end2end_04, idx, run='04', /resume, /use_kweight

    print, 'That end2end run took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Jacks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO azrms_jack
for idx=0, 19 do begin
    print, 'azrms jack, field ', idx
    t0 = systime(0, /seconds)

    lps12_jack, idx, 'azrms'

    print, 'That azrms jack took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor
END


;-----------------------------------------
; run lots of stuff, part 1
PRO process_0524_1

; azrms sim jack
for idx=0, 19 do begin
    print, 'azrms jack sim, field ', idx
    t0 = systime(0, /seconds)

    lps12_jack, idx, 'azrms', /sim

    print, 'That azrms jack sim took ', (systime(0, /seconds) - t0)/60., ' minutes.'
endfor

; azrms_95 jack, sim
for idx=0, 19 do begin
    print, 'azrms_95 jack, field ', idx
    t0 = systime(0, /seconds)

    ; jack
    lps12_jack, idx, 'azrms_95'
    t1 = (systime(0, /seconds))
    print, 'That azrms_95 jack took ', (t1 - t0)/60., ' minutes.'

    ; sim
    lps12_jack, idx, 'azrms_95', /sim

    print, 'That azrms_95 jack sim took ', (systime(0, /seconds) - t1)/60., ' minutes.'
endfor

; azrms_90 jack, sim
for idx=0, 19 do begin
    print, 'azrms_90 jack, field ', idx
    t0 = systime(0, /seconds)

    ; jack
    lps12_jack, idx, 'azrms_90'
    t1 = (systime(0, /seconds))
    print, 'That azrms_90 jack took ', (t1 - t0)/60., ' minutes.'

    ; sim
    lps12_jack, idx, 'azrms_90', /sim

    print, 'That azrms_90 jack sim took ', (systime(0, /seconds) - t1)/60., ' minutes.'
endfor

; azrms_70 jack, sim
for idx=0, 19 do begin
    print, 'azrms_70 jack, field ', idx
    t0 = systime(0, /seconds)

    ; jack
    lps12_jack, idx, 'azrms_70'
    t1 = (systime(0, /seconds))
    print, 'That azrms_70 jack took ', (t1 - t0)/60., ' minutes.'

    ; sim
    lps12_jack, idx, 'azrms_70', /sim

    print, 'That azrms_70 jack sim took ', (systime(0, /seconds) - t1)/60., ' minutes.'
endfor

END




;;;;;;;;;;;;;;;;;;;;
;-----------------------------------------
; run lots of stuff, part 2
PRO process_0524_2

; lr jack, sim
for idx=0, 19 do begin
    print, 'lr jack, field ', idx
    t0 = systime(0, /seconds)

    ; jack
    lps12_jack, idx, 'lr'
    t1 = (systime(0, /seconds))
    print, 'That lr jack took ', (t1 - t0)/60., ' minutes.'

    ; sim
    lps12_jack, idx, 'lr', /sim

    print, 'That lr jack sim took ', (systime(0, /seconds) - t1)/60., ' minutes.'
endfor

; 12 jack, sim
for idx=0, 19 do begin
    print, '12 jack, field ', idx
    t0 = systime(0, /seconds)

    ; jack
    lps12_jack, idx, '12'
    t1 = (systime(0, /seconds))
    print, 'That 12 jack took ', (t1 - t0)/60., ' minutes.'

    ; sim
    lps12_jack, idx, '12', /sim

    print, 'That 12 jack sim took ', (systime(0, /seconds) - t1)/60., ' minutes.'
endfor

; tweight jack, sim
for idx=0, 19 do begin
    print, 'tweight jack, field ', idx
    t0 = systime(0, /seconds)

    ; jack
    lps12_jack, idx, 'tweight'
    t1 = (systime(0, /seconds))
    print, 'That tweight jack took ', (t1 - t0)/60., ' minutes.'

    ; sim
    lps12_jack, idx, 'tweight', /sim

    print, 'That tweight jack sim took ', (systime(0, /seconds) - t1)/60., ' minutes.'
endfor

; moon jack, sim
for idx=0, 19 do begin
    print, 'moon jack, field ', idx
    t0 = systime(0, /seconds)

    ; jack
    lps12_jack, idx, 'moon'
    t1 = (systime(0, /seconds))
    print, 'That moon jack took ', (t1 - t0)/60., ' minutes.'

    ; sim
    lps12_jack, idx, 'moon', /sim

    print, 'That moon jack sim took ', (systime(0, /seconds) - t1)/60., ' minutes.'
endfor

END
