;;;
; NAME: script_0411
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) lr jacks, over-plotting
;
; MODIFICATION HISTORY:
;  04/11/2012: (KTS) Created
;;;

;...................................................................
; 
PRO plot_lr
compile_opt IDL2, HIDDEN

restore, '/home/kstory/lps12/scripts/sav_files/lr_0410.sav'
i0=1 & i1=5
stub='lr'
field_arr_ = lps12_fieldstruct()
legend_str = strarr(20)

; make plots

; plot all jack ontop of each other
window, 2, xsize=1000, ysize=800
legend_str = strarr(20)
color_arr = kscolor_array()
yr = minmax(dls*1e12)
plot,banddef[i0:i1],mndl*1e12,$
  title='combined '+stub+' jack, '+strtrim(string(chisq_tot),2)+' chisq',$
  xr=[min(banddef[i0:i1])-100,max(banddef[i0:i1])+100],/xst,xtitle='!12l!X!N',$
  ytitle='Dl ('+textoidl('\mu K^2')+')', yr=yr;,/yst

; make array of colors
for ii=0, 19 do begin
    fdl  = dls[*,ii]*1e12 &   fddl = sqrt(diagcovs[*,ii])*1e12
    oplot, banddef[i0:i1], fdl, color=color_arr[ii], linestyle=3
    errplot,banddef[i0:i1],(fdl-fddl),(fdl+fddl), color=color_arr[ii], linestyle=3
    legend_str[ii] = field_arr_[ii].name
endfor
legend, legend_str, linestyle=3, colors=color_arr[0:19], position=[2300, ((yr[1] - yr[0])*0.95 + yr[0])]

oplot,banddef[i0:i1],mndl*1e12,thick=3

; save plots
if 0 then begin
    plot_name = '/home/kstory/lps12/figs/'+stub+'_jack_0411'
    err= tvread(/png, filename=plot_name, /nodialog)
endif


stop
END

