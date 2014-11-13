;;;
; NAME: plots_0820
; PURPOSE:
;   Iteratively plot the runlist cuts for all fields
;
; MODIFICATION HISTORY:
;  08/20/2011: (KTS) Created
;;;

;...................................................................
; Main function
;
pro plots_0820

restore, 'lps12_fieldnames.sav'

for ii=0, n_elements(all_fields)-1 do begin
    field_name = all_fields[ii]
    print, 'Field name = ', field_name
    restore, '/data/kstory/projects/lps12/runlists/keep_0820/runlist_cuts_wtrms_'+field_name+'.sav'
    
    ;data = medwts & mytitle='medwts '+field_name
    ;data = rms_1am_midmap & mytitle='rms '+field_name
    ;data = medwts*rms_1am_midmap^2 & mytitle='medwts * rms^2 '+field_name
    xmean = mean(data)
    wset, 1 & plot, data
    
    wset, 2
    kplot_hist, data, mytitle
    oplot, [1.5*xmean, 1.5*xmean], [0,1], color=!red
    oplot, [xmean/3., xmean/3.], [0,1], color=!green
    oplot, [xmean/2., xmean/2.], [0,1], color=!green
    oplot, [xmean/1.5, xmean/1.5], [0,1], color=!green
    oplot, [xmean/1.3, xmean/1.3], [0,1], color=!green
    
    wh = where(data lt xmean/3., nwh) & print, nwh
    wh = where(data lt xmean/2., nwh) & print, nwh
    wh = where(data lt xmean/1.5, nwh) & print, nwh
    wh = where(data lt xmean/1.3, nwh) & print, nwh
    
    stop
endfor
end
