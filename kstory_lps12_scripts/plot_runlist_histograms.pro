;;;
; NAME: plot_runlist_histograms
; PURPOSE:
;   Iteratively plot the runlist cuts for all fields
;
; MODIFICATION HISTORY:
;  08/20/2011: (KTS) Created
;  08/30/2011: (KTS) Add total_weight as a plotting option
;  09/13/2011: (KTS) Modify to work with lps12_fieldstruct
;;;

;...................................................................
; Main function
;
pro plot_runlist_histograms, save_plots=save_plots
compile_opt IDL2, HIDDEN

;restore, 'lps12_fieldnames.sav'
field_arr_ = lps12_fieldstruct()

;field_list = all_fields[[1,5,11]]
field_list = [11]
;field_list = indgen(15)

for ii=0, n_elements(field_list)-1 do begin

    idx = field_list[ii]
    field_name = field_arr_[idx].name

    print, 'Field name = ', field_name, ' idx = ',idx
    ;restore, '/data/kstory/projects/lps12/runlists/default_sav/runlist_cuts_wtrms_'+field_name+'.sav'
    ;restore, '/data/kstory/projects/lps12/runlists/tweight_sav/runlist_tweight_'+field_name+'.sav'
    restore, '/data/kstory/projects/lps12/runlists/sav_files/runlist_lps12_'+field_name+'.sav'

    ;wset, 1 & plot, data
    
    for jj=0, 3 do begin
        if (jj eq 0) then begin 
            data = medwts & type = 'medwts' & mytitle=type+' '+field_name & fx=1.25 & fn=2.
        endif else if (jj eq 1) then begin
            data = rms_1am_midmap & type = 'rms' & mytitle=type+' '+field_name & fx=2. & fn=3.
        endif else if (jj eq 2) then begin
            data = medwts*rms_1am_midmap^2 & type = 'wtrms' & mytitle=type+' '+field_name & fx=2. & fn=2.
        endif else if (jj eq 3) then begin
            data = tweight & type = 'tweight' & mytitle=type+' '+field_name & fx=1.2 & fn=1.5
        endif
        xmedian = median(data)
        

        window, jj
        kplot_hist, data[where(data lt 3*xmedian)], mytitle, nbins=100, xmax=xmedian*3.
;         oplot, [1.25*xmedian, 1.25*xmedian], [0,1], color=!red
;         oplot, [1.3*xmedian, 1.3*xmedian], [0,1], color=!red
;         oplot, [1.5*xmedian, 1.5*xmedian], [0,1], color=!red
;         oplot, [2.*xmedian, 2.*xmedian], [0,1], color=!red
;         oplot, [xmedian/3., xmedian/3.], [0,1], color=!green
;         oplot, [xmedian/2., xmedian/2.], [0,1], color=!green
;         oplot, [xmedian/1.5, xmedian/1.5], [0,1], color=!green
;         oplot, [xmedian/1.3, xmedian/1.3], [0,1], color=!green
;         wh = where(data gt xmedian*2., nwh) & print, nwh
;         wh = where(data gt xmedian*1.5, nwh) & print, nwh
;         wh = where(data gt xmedian*1.3, nwh) & print, nwh
;         wh = where(data gt xmedian*1.2, nwh) & print, nwh
;         wh = where(data gt xmedian*1.1, nwh) & print, nwh
    

        ;;; Final cut set:
        oplot, [fx*xmedian, fx*xmedian], [0,1], color=!red
        oplot, [xmedian/fn, xmedian/fn], [0,1], color=!green

        print, '  cut type: ', type
        wh = where(data lt xmedian/fn, nwh) & print, '            Number cut by LOW cut (',fn,'): ', nwh
        wh = where(data gt xmedian*fx, nwh) & print, '            Number cut by HIGH cut(',fx,'): ', nwh
    
        if keyword_set(save_plots) then err = tvread(/png,filename=type+'_cut_hist_'+field_name, /nodialog) 
    
    endfor
    stop
endfor
end
