;;;
; NAME: find_1hz_bolo_cuts
; PURPOSE:
;   Find bolometer indicies of bolometers that are systematically
;   affected by 1hz noise.
;
; OUTPUT:
;   Print bad bolos to the screen, make plots
;
; NOTES:
; 1) Uses .sav files that were made with near_one_hz.pro
;
; MODIFICATION HISTORY:
;  02/10/2012: (KTS) Created
;;;

;*****************************************************************
; Analyze 1hz bolo cuts for a single field
pro get_1field, field_idx, save_plots=save_plots
compile_opt IDL2, HIDDEN
field_arr = lps12_fieldstruct()
field_name = field_arr[field_idx].name

restore, '/data/kstory/projects/lps12/1hzBolo/near_one_hz_'+field_name+'.sav'

if (intersect(field_idx, [0,1]) ne -1) then year=2008
if (intersect(field_idx, [3,4,5]) ne -1) then year=2009
if (intersect(field_idx, [2,6,7,8,9,10]) ne -1) then year=2010
if (intersect(field_idx, indgen(9)+11) ne -1) then year=2011

print, "year = ", year

if keyword_set(save_plots) then begin
    find_1hz_bolo_cuts, fitamp, fitfreq, field_name, year, boloind, /save_plots
endif else begin
    find_1hz_bolo_cuts, fitamp, fitfreq, field_name, year, boloind
endelse
end
    




;*****************************************************************
; Analyze all 2008 fields simultaneously
pro get_2008, save_plots=save_plots
compile_opt IDL2, HIDDEN
field_arr = lps12_fieldstruct()

field_list = [0,1] ; List of 2009 fields
nbolos = 394 ; from sav files
nfield = n_elements(field_list)

tmpamp = fltarr(nfield*100, nbolos)
tmpfreq = fltarr(nfield*100, nbolos)

for i=0, nfield-1 do begin
    field_idx = field_list[i]

    field_name = field_arr[field_idx].name
    restore, '/data/kstory/projects/lps12/1hzBolo/near_one_hz_'+field_name+'.sav'
    tmpamp[i*100+0:i*100+99,*] = fitamp
    tmpfreq[i*100+0:i*100+99,*] = fitfreq
endfor

fitamp = tmpamp
fitfreq = tmpfreq

; Call find_1hz_bolo_cuts
if keyword_set(save_plots) then begin
    find_1hz_bolo_cuts, fitamp, fitfreq, 'all_2008', 2008.1, boloind, /save_plots
endif else begin
    find_1hz_bolo_cuts, fitamp, fitfreq, 'all_2008', 2008.1, boloind
endelse

end



;*****************************************************************
; Analyze all 2009 fields simultaneously
pro get_2009, save_plots=save_plots
compile_opt IDL2, HIDDEN
field_arr = lps12_fieldstruct()

field_list = [3,4,5] ; List of 2009 fields
nbolos = 514 ; from sav files
nfield = n_elements(field_list)

tmpamp = fltarr(nfield*100, nbolos)
tmpfreq = fltarr(nfield*100, nbolos)

for i=0, nfield-1 do begin
    field_idx = field_list[i]

    field_name = field_arr[field_idx].name
    restore, '/data/kstory/projects/lps12/1hzBolo/near_one_hz_'+field_name+'.sav'
    tmpamp[i*100+0:i*100+99,*] = fitamp
    tmpfreq[i*100+0:i*100+99,*] = fitfreq
endfor

fitamp = tmpamp
fitfreq = tmpfreq

; Call find_1hz_bolo_cuts
if keyword_set(save_plots) then begin
    find_1hz_bolo_cuts, fitamp, fitfreq, 'all_2009', 2009.1, boloind, /save_plots
endif else begin
    find_1hz_bolo_cuts, fitamp, fitfreq, 'all_2009', 2009.1, boloind
endelse

end


;*****************************************************************
; Analyze all 2010 fields simultaneously
pro get_2010, save_plots=save_plots
compile_opt IDL2, HIDDEN
field_arr = lps12_fieldstruct()

field_list = [2,6,7,8,9,10] ; List of 2010 fields
nbolos = 514 ; from sav files
nfield = n_elements(field_list)

tmpamp = fltarr(nfield*100, nbolos)
tmpfreq = fltarr(nfield*100, nbolos)

for i=0, nfield-1 do begin
    field_idx = field_list[i]

    field_name = field_arr[field_idx].name
    restore, '/data/kstory/projects/lps12/1hzBolo/near_one_hz_'+field_name+'.sav'
    tmpamp[i*100+0:i*100+99,*] = fitamp
    tmpfreq[i*100+0:i*100+99,*] = fitfreq
endfor

fitamp = tmpamp
fitfreq = tmpfreq

; Call find_1hz_bolo_cuts
if keyword_set(save_plots) then begin
    find_1hz_bolo_cuts, fitamp, fitfreq, 'all_2010', 2010.1, boloind, /save_plots
endif else begin
    find_1hz_bolo_cuts, fitamp, fitfreq, 'all_2010', 2010.1, boloind
endelse

end


;*****************************************************************
; Analyze all 2011 fields simultaneously
pro get_2011, save_plots=save_plots
compile_opt IDL2, HIDDEN
field_arr = lps12_fieldstruct()

field_list = indgen(9)+11 ; List of 2010 fields
nbolos = 514 ; from sav files
nfield = n_elements(field_list)

tmpamp = fltarr(nfield*100, nbolos)
tmpfreq = fltarr(nfield*100, nbolos)

for i=0, nfield-1 do begin
    field_idx = field_list[i]

    field_name = field_arr[field_idx].name
    restore, '/data/kstory/projects/lps12/1hzBolo/near_one_hz_'+field_name+'.sav'
    tmpamp[i*100+0:i*100+99,*] = fitamp
    tmpfreq[i*100+0:i*100+99,*] = fitfreq
endfor

fitamp = tmpamp
fitfreq = tmpfreq

; Call find_1hz_bolo_cuts
if keyword_set(save_plots) then begin
    find_1hz_bolo_cuts, fitamp, fitfreq, 'all_2011', 2011.1, boloind, /save_plots
endif else begin
    find_1hz_bolo_cuts, fitamp, fitfreq, 'all_2011', 2011.1, boloind
endelse

end



; plot 1hz bolo stuff
pro find_1hz_bolo_cuts, fitamp, fitfreq, field_name, year, boloind, save_plots=save_plots
compile_opt IDL2, HIDDEN

;-----------------------
; Cut Values
;-----------------------
; Cut 1: cut bolos if fitamp > p_cutoff
; Cut 2: cut bolos if nobs with (fitamp > p_cutoff) is >= nbad_obs_cutoff
; Cut 3: cut bolos if average fitfreq not in [ 1 +- freq_cutoff]

case year of 
    ;;; 2008 cuts
    2008: begin
        p_cutoff = 7.e-5 
        nbad_obs_cutoff = 28
        freq_cutoff = 0.05
        bad1hzbolos = [22, 58, 66, 93,99, 116, 127, 326, 341, 344, 355, 356, 362, 370, 412, 426, 435, 472, 502, 543, 639] ; From RK
    end

    ;;; All of 2009
    2008.1: begin
        p_cutoff = 7.e-5 
        nbad_obs_cutoff = 28*2
        freq_cutoff = 0.05
        bad1hzbolos = [22, 58, 66, 93,99, 116, 127, 326, 341, 344, 355, 356, 362, 370, 412, 426, 435, 472, 502, 543, 639] ; From RK
    end

    ;;; 2009 cuts
    2009: begin
        p_cutoff = 5.e-5 
        nbad_obs_cutoff = 20
        freq_cutoff = 0.05
        bad1hzbolos = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768] ; From RK
    end

    ;;; all of 2009
    2009.1: begin
        p_cutoff = 5.e-5 
        nbad_obs_cutoff = 60
        freq_cutoff = 0.05
        bad1hzbolos = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768] ; From RK
    end

    ;;; 2010 cuts
    2010: begin
        p_cutoff = 5.e-5 
        nbad_obs_cutoff = 30
        freq_cutoff = 0.05
        bad1hzbolos = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768] ; From RK, 2009
    end

    ;;; all of 2010
    2010.1: begin
        p_cutoff = 5.e-5 
        nbad_obs_cutoff = 150;30*6
        freq_cutoff = 0.05
        bad1hzbolos = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768] ; From RK, 2009
    end

    ;;; 2011 cuts
    2011: begin
        p_cutoff = 5.e-5 
        nbad_obs_cutoff = 30
        freq_cutoff = 0.05
        bad1hzbolos = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768] ; From RK, 2009
    end

    ;;; all of 2011
    2011.1: begin
        p_cutoff = 5.e-5 
        nbad_obs_cutoff = 150;30*6
        freq_cutoff = 0.05
        bad1hzbolos = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768] ; From RK, 2009
    end

endcase


; Remove fitamp, fitfreq < 0
whneg = where(fitamp lt 0, nwh)
if nwh gt 0 then fitamp[whneg] = 0

whneg = where(fitfreq lt 0, nwh)
if nwh gt 0 then fitfreq[whneg] = 0

;-----------------------
; Find bolos to cut
;-----------------------

nbolos = n_elements(fitamp[0,*])
nbad = intarr(nbolos)
nobs = n_elements(fitamp[*,0])

; Cut 1
for i=0, nobs-1 do begin
    for ibolo = 0, nbolos-1 do begin
        if fitamp[i, ibolo] gt p_cutoff then begin
            nbad[ibolo]++
        endif

    endfor
endfor

; Cut 2
badidx = where(nbad ge nbad_obs_cutoff, complement=goodidx)
badbolo2 = boloind[badidx]

; Cut 3
failfreq12 = fitfreq[*,badidx] ; frequency arrays of bolos that fail Cut 1 and Cut 2
passfreq12 = fitfreq[*,goodidx]
n_bad12 = n_elements(failfreq12[0,*])
mvec = fltarr(n_bad12)
for i=0, n_bad12 - 1 do mvec[i] = median(failfreq12[*,i]) ; find mean freq for each bad bolo
whbad3 = where( abs(mvec - 1) lt freq_cutoff, nwh)
badbolo3 = boloind[ badidx[whbad3] ]

;-----------------------
; Print results
;-----------------------

print, "======== Bad Bolometers:  =========== field: ", field_name, ": ", p_cutoff, freq_cutoff, nbad_obs_cutoff
print, "init ", badbolo2
print, "KTS: ", badbolo3
print, "RK:  ", bad1hzbolos





;-----------------------
; Plot fitamp
;-----------------------
; Plot all observations in one hist
wset, 1
mx = 2e-4
nbins = 100
bins = indgen(nbins)*mx/nbins + 0.
hist = histogram(fitamp, min=0, max=mx, nbins=nbins)

plot, bins[1:*], hist[1:*], ytitle='Nbolos',xtitle='fitamp of 1hz line', $
  title='Fitamp of 1hz line ' + field_name
oplot, [p_cutoff, p_cutoff], [0,max(hist)], lines=2, color=!red


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
; Plot fitfreq
;-----------------------
wset, 2
; Plot all observations in one hist
mx = 1.3
mn = 0.7
nbins = 100
bins = indgen(nbins)*(mx-mn)/nbins + mn
hist = histogram(fitfreq, min=mn, max=mx, nbins=nbins)

plot, bins[1:*], hist[1:*], ytitle='Nbolos',xtitle='fitted location of 1hz line', $
  title='Fitfreq of 1hz line ' + field_name


;-----------------------
; Plot cut 2
;-----------------------
wset, 3
plot, nbad, xtitle='BoloIndex', ytitle='Nobs with high 1hz power', $
  title='Bolos to cut for high 1hz Power ' + field_name
oplot, [0,nbolos], [nbad_obs_cutoff, nbad_obs_cutoff], lines=2, color=!red


;-----------------------
; Plot cut 3
;-----------------------
wset, 4
mx = 1.3 & mn = 0.7 & nbins = 100
bins = indgen(nbins)*(mx-mn)/nbins + mn
hist = histogram(failfreq12, min=mn, max=mx, nbins=nbins)
plot, bins, hist, ytitle='Nbolos',xtitle='fitted location of 1hz line', $
  title='fitfreq of Bolos to cut ' + field_name

; plot the remainder
wset, 6
hist = histogram(passfreq12, min=mn, max=mx, nbins=nbins)
plot, bins, hist, ytitle='Nbolos',xtitle='fitted location of 1hz line', $
  title='fitfreq of good Bolos ' + field_name


wset, 5
mx = 1.3 & mn = 0.7 & nbins = 100
bins = indgen(nbins)*(mx-mn)/nbins + mn
hist = histogram(mvec, min=mn, max=mx, nbins=nbins)
plot, bins, hist, ytitle='Nbolos',xtitle='<fitted location of 1hz line>_obs', $
  title='<fitfreq>_obs of Bolos to cut ' + field_name
oplot, [1-freq_cutoff, 1-freq_cutoff], [0, max(hist)], lines=2, color=!red
oplot, [1+freq_cutoff, 1+freq_cutoff], [0, max(hist)], lines=2, color=!red


if keyword_set(save_plots) then begin
    figs_dir = '/data/kstory/projects/lps12/1hzBolo/figs/'

    wset, 1
    err= tvread(/png, filename=figs_dir+'fitamp_'+field_name, /nodialog)

    wset, 2
    err= tvread(/png, filename=figs_dir+'fitfreq_'+field_name, /nodialog)

    wset, 3
    err= tvread(/png, filename=figs_dir+'amp_cut_'+field_name, /nodialog)

    wset, 4
    err= tvread(/png, filename=figs_dir+'fitfreq_cut_'+field_name, /nodialog)

    wset, 5
    err= tvread(/png, filename=figs_dir+'fitfreq_avg_cut_'+field_name, /nodialog)
endif    

stop

end



;*****************************************************************
; Check if flagged bolometers live on the same squids
; This uses the sets of flagged bolometers obtained above.
pro check_squids

; 2008
bb08 = [22,  58,  66,  99, 116, 127, 129, 137, 326, 341, 355, 356, 362, 370, 426, 472, 502, 529, 543, 639]
rk08 = [22,  58,  66,  93,  99, 116, 127, 326, 341, 344, 355, 356, 362, 370, 412, 426, 435, 472, 502, 543, 639]

ac08 = read_array_config(starttime='10-Mar-2008:06:24:22', endtime='10-Mar-2008:16:54:26')

sq08 = ac08.squid_addr[bb08]
rk_sq08 = ac08.squid_addr[rk08]

uniq_sq08 = uniqval( union(sq08, rk_sq08) )

print, '=== 2008 ==='
tab = STRING(9B)
print, 'squid'+tab+tab+' kst'+tab+'     rk'
for ii=0, n_elements(uniq_sq08)-1 do begin
    wh = where(sq08 eq uniq_sq08[ii], nwh)
    wh = where(rk_sq08 eq uniq_sq08[ii], rk_nwh)
    print, uniq_sq08[ii], nwh, rk_nwh
endfor

; 2009
bb09 = [58, 158, 166, 171, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 220, 235, 268, 362, 389, 407, 465, 644, 768]
rk09 = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768]

ac09 = read_array_config(starttime='25-Mar-2009:11:18:49', endtime='25-Mar-2009:12:21:38')

sq09 = ac09.squid_addr[bb09]
rk_sq09 = ac09.squid_addr[rk09]

uniq_sq09 = uniqval( union(sq09, rk_sq09) )

print, '=== 2009 ==='
tab = STRING(9B)
print, 'squid'+tab+tab+' kts'+tab+'     rk'
for ii=0, n_elements(uniq_sq09)-1 do begin
    wh = where(sq09 eq uniq_sq09[ii], nwh)
    wh = where(rk_sq09 eq uniq_sq09[ii], rk_nwh)
    print, uniq_sq09[ii], nwh, rk_nwh
endfor
stop
end
