;;;
; NAME: script_0209
; PURPOSE:
;   Make IDL sav files for 1hz bolos
;
; NOTES:
; 1) Now using the lpf-ds IDF's from Ryan, Jan3 2012
;
; MODIFICATION HISTORY:
;  02/09/2012: (KTS) Created
;;;

pro get_1hz_0209_0to5
compile_opt IDL2, HIDDEN
field_idx_list = [0,1,2,3,4,5]

field_arr_ = lps12_fieldstruct()
for ii=0, n_elements(field_idx_list) -1 do begin
    print, '*********** get_1hz_bolo_info for field, ',field_arr_[field_idx_list[ii]].name
    get_1hz_bolo_info, ii
endfor
end

; plot stuff
pro myplot
restore, '/data/kstory/projects/lps12/1hzBolo/near_one_hz_ra5h30dec-55.sav'

;-----------------------
; fitamp
;-----------------------
whneg = where(fitamp lt 0, nwh)
if nwh gt 0 then fitamp[whneg] = 0

; Plot all observations in one hist
wset, 1
mx = 2e-4
nbins = 100
bins = indgen(nbins)*mx/nbins + 0.
hist = histogram(fitamp, min=0, max=mx, nbins=nbins)

plot, bins[1:*], hist[1:*], ytitle='Nbolos',xtitle='fitamp of 1hz line'


; Plot individual observations
if 0 then begin
for i=0, nobs-1 do begin
    x = fitamp[i,*]

    bins = indgen(50)*( max(x) - min(x) ) / 50. + min(x)
    hist = histogram(x, nbins=50)

    plot, bins, hist
    pause
endfor
endif


;-----------------------
; fitfreq
;-----------------------
wset, 2
whneg = where(fitfreq lt 0, nwh)
if nwh gt 0 then fitfreq[whneg] = 0

; Plot all observations in one hist
mx = 1.3
mn = 0.7
nbins = 100
bins = indgen(nbins)*(mx-mn)/nbins + mn
hist = histogram(fitfreq, min=mn, max=mx, nbins=nbins)

plot, bins[1:*], hist[1:*], ytitle='Nbolos',xtitle='fitted location of 1hz line'


;-----------------------
; Find bolos to cut
;-----------------------
p_cutoff = 8.e-5
nbad_cutoff = 27


nbolos = n_elements(fitamp[0,*])
nbad = intarr(nbolos)

count = 0
for i=0, nobs-1 do begin

    for ibolo = 0, nbolos-1 do begin
        if fitamp[i, ibolo] gt p_cutoff then begin
            nbad[ibolo]++
            print, 'ibolo =', ibolo, ', nbad[ibolo] = ', nbad[ibolo]
            count++
        endif

    endfor
endfor

wset, 3
plot, nbad, xtitle='BoloIndex', ytitle='Nobs with high 1hz power'
oplot, [0,nbolos], [nbad_cutoff, nbad_cutoff], lines=2, color=!red


badidx = where(nbad gt nbad_cutoff)
badbolo = boloind[badidx]

; From RK
bad1hzbolos = [22, 58, 66, 93,99, 116, 127, 326, 341, 344, 355, 356, 362, 370, 412, 426, 435, 472, 502, 543, 639]

print, "======== Bad Bolometers ==========="
print, "KTS: ", badbolo
print, "RK:  ", bad1hzbolos

stop

end

