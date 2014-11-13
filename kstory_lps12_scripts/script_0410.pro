;;;
; NAME: script_0410
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) lr jacks, over-plotting
;
; MODIFICATION HISTORY:
;  04/10/2012: (KTS) Created
;;;

;...................................................................
; 
PRO plot_lr
compile_opt IDL2, HIDDEN

restore, '/home/kstory/lps12/scripts/sav_files/lr_0410.sav'
i0=1 & i1=5
stub='lr'
f = lps12_fieldstruct()
legend_str = strarr(20)

window, 2, xsize=1000, ysize=800
color_arr = kscolor_array()
yr = minmax(dls*1e12)
plot,banddef[i0:i1],mndl*1e12,$
  title='combined '+stub+' jack, '+strtrim(string(chisq_tot),2)+' chisq',$
  xr=[min(banddef[i0:i1])-100,max(banddef[i0:i1])+100],/xst,xtitle='!12l!X!N',$
  ytitle='Dl ('+textoidl('\mu K^2')+')', yr=yr;,/yst

; make array of colors
for ii=0, 19 do begin
    fdl  = dls[*,ii]*1e12
    fddl = sqrt(diagcovs[*,ii])*1e12

    oplot, banddef[i0:i1], fdl, color=color_arr[ii], linestyle=3
    errplot,banddef[i0:i1],(fdl-fddl),(fdl+fddl), color=color_arr[ii], linestyle=3

    legend_str[ii] = f[ii].name
endfor

; make legend

legend, legend_str, linestyle=3, colors=color_arr[0:19], position=[2300, 7.5]

stop
END

