;;;
; NAME: plots_0807
; PURPOSE:
;   stack histograms of runlist cuts for medwt, rms, wtrms
;
; MODIFICATION HISTORY:
;  08/07/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro plots_0807, save_plots=save_plots, stopit=stopit
compile_opt IDL2, HIDDEN

; Field names
restore, 'lps12_fieldnames.sav'
colors = get_color_array()

window, 1, xsize=600, ysize=500
first_plot=1

;!p.multi = [0,3,2]
for ii=0, n_elements(all_fields)-1 do begin
;for ii=0, 5 do begin
    field_name = all_fields[ii]
    if (field_name eq 'ra23hdec-45') then continue

    print, field_name
    restore, '/data/kstory/projects/lps12/runlists/runlist_cuts_wtrms_'+field_name+'.sav'

    ; Make histograms
    ;;; medwts
    ;data = medwts & hist_title='Medwts' & sav_name='medwts'
    ;data = rms_1am_midmap & hist_title='RMS' & sav_name='rms'
    data = medwts*rms_1am_midmap^2 & hist_title='Medwts * RMS^2' & & sav_name='wtrms'

    wset, 1
    if (first_plot eq 1) then begin
        oplot_hist, data, colors[ii], title=hist_title,/no_oplot
        first_plot = 0          ; Set first plot to zero
    endif else begin
        oplot_hist, data, colors[ii]
    endelse

    if keyword_set(stopit) then stop
endfor
;!p.multi=0

; Save plots
wset, 1
if keyword_set(save_plots) then err=tvread(/png, filename='figs/wtrms_hist_all_fields_'+sav_name+'_0807', /nodialog)


end

;;; over-Plot histograms
pro oplot_hist, data, hist_color, title=title, no_oplot=no_oplot
hmin=min(data) & hmax=max(data)
hh=  histogram(data, MIN=hmin, MAX=hmax, nbins=50)
hh_norm = DOUBLE(hh) / max(hh)
xvec = hmin + findgen(50) * (hmax-hmin) / 50.
low_medwt = median(data) / 1e6 & high_medwt = median(data) * 2.

if keyword_set(no_oplot) then begin
    plot, xvec, hh_norm, title=title
    oplot, xvec, hh_norm, color=hist_color
endif else begin
    oplot, xvec, hh_norm, color=hist_color
endelse
;oplot, [low_medwt, low_medwt], [0,1], color=!red
;oplot, [high_medwt, high_medwt], [0,1], color=!red
end

;;; over-Plot histograms, align by peak value
pro oplot_hist, data, hist_color, title=title, no_oplot=no_oplot
hmin=min(data) & hmax=max(data)
hh=  histogram(data, MIN=hmin, MAX=hmax, nbins=50)
hh_norm = DOUBLE(hh) / max(hh)
xvec = hmin + findgen(50) * (hmax-hmin) / 50.
low_medwt = median(data) / 1e6 & high_medwt = median(data) * 2.

if keyword_set(no_oplot) then begin
    plot, xvec, hh_norm, title=title
    oplot, xvec, hh_norm, color=hist_color
endif else begin
    oplot, xvec, hh_norm, color=hist_color
endelse
;oplot, [low_medwt, low_medwt], [0,1], color=!red
;oplot, [high_medwt, high_medwt], [0,1], color=!red
end




function get_color_array
device, decomposed=1
device, true_color=24

cArr = [255L, $ 
256L * (255), $
256L * (256L * 255), $
3000010, $
38650L, $
16581006, $
111617L, $
!blue+255L*256L, $
6782758L, $
10564482, $
7310724, $
7036682, $
2499467, $
7331449, $
16480249L, $
874399L, $
1195345L, $
961733, $
10136312L, $
13799838L, $
10866611L, $
15127481L]

return, cArr
end

