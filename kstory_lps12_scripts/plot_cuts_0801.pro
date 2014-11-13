;;;
; NAME: plot_cuts_0801
; PURPOSE:
;   Plot runlist cuts
;
; MODIFICATION HISTORY:
;  08/01/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro plot_cuts_0801, save_plots=save_plots, stopit=stopit
compile_opt IDL2, HIDDEN

; Field names
restore, 'lps12_fieldnames.sav'

window, 1, xsize=1200, ysize=900

!p.multi = [0,3,2]
;for ii=0, n_elements(all_fields)-1 do begin
for ii=0, 3 do begin
    field_name = fields_2010[ii]
    restore, '/data/kstory/projects/lps12/runlists/runlist_cuts_wtrms_'+field_name+'.sav'

    print, field_name

    wset, 1
    plot, medwts, title='medtws'
    plot, rms_1am_midmap, title=field_name+': rms'
    plot, medwts*(rms_1am_midmap^2.), title='medwts*rms^2'

    ; Make histograms
    data = medwts & kplot_hist, data, 'medtws'
    data = rms_1am_midmap & kplot_hist, data, field_name+': rms'
    data = medwts*rms_1am_midmap^2 & kplot_hist, data, 'medwts*rms^2'

    if keyword_set(save_plots) then err=tvread(/png, filename='figs/wtrms_cuts_'+field_name, /nodialog)

    if keyword_set(stopit) then stop
endfor
!p.multi=0

end

pro kplot_hist, data, ktitle
hmin=min(data) & hmax=max(data)
hh=  histogram(data, MIN=hmin, MAX=hmax, nbins=50)
hh_norm = DOUBLE(hh) / max(hh)
xvec = hmin + findgen(50) * (hmax-hmin) / 50.
low_medwt = median(data) / 1e6 & high_medwt = median(data) * 2.

plot, xvec, hh_norm, title=ktitle
;oplot, [low_medwt, low_medwt], [0,1], color=!red
;oplot, [high_medwt, high_medwt], [0,1], color=!red
end
