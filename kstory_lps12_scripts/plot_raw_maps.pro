;;;
; NAME: plot_raw_maps.pro
; PURPOSE:
;   Make plots of all raw maps
;
; CALLING SEQUENCE: plot_raw_maps
;
; INPUTS:
;   field_name,  name of the field
;   dates,       a string array of dates in the format
;                '10-Aug-2008:10:00:00' to use for selecting which files to coadd.
;   scale,       the scale for tv_spt_map, usually 0.17
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
; NOTES:
;   1) 
;
; MODIFICATION HISTORY:
;   07/22/2011: (KTS) Created
;   08/31/2011: (KTS) Modified to run on updated runlists
;   09/20/2011: (KTS) Modified to use lps12_fieldstruct
;;;

;...................................................................
; Call plot_raw_maps
;
pro call_plot_raw_maps, field_list
compile_opt IDL2, HIDDEN

field_arr_ = lps12_fieldstruct()
;restore, 'lps12_fieldnames.sav'
;field_list = [1,2,6,12]
;field_list = [6]
;field_list = setdifference(lindgen(n_elements(all_fields)), exclude_list)

for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name
    print, "************* make raw maps, " + name
    plot_raw_maps, idx, scale=0.17
    make_summary_pdf, [idx]
endfor
end

;...................................................................
; Main function
;
pro plot_raw_maps, idx, scale=scale
compile_opt IDL2, HIDDEN

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
loadct, 0
device, decomposed=1
device, true_color=24
field_arr_ = lps12_fieldstruct()

if n_elements(scale) eq 0 then scale=-1

name = field_arr_[idx].name
map_dirs = field_arr_[idx].autoprocessed_map_dirs

; Get dates
runlist = '/data/kstory/projects/lps12/runlists/runlist_lps12_'+name+'.txt'
readcol, runlist, dates, path_to_maps, format='a,a',delim=' '

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Get the filenames of all maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

list = ['']
for md=0, n_elements(map_dirs) - 1 do begin
    ls = 'ls '+ map_dirs[md] + '/' + '*' + '_150_*.fits'
    spawn, ls, tmp
    list = [list, tmp]
endfor
list = list[1:*]
listdates=convert_filelist_to_mjds(list)
inputdates=date_string_to_mjd(dates)

count=0
for i=0,n_elements(inputdates)-1 do begin
    wh=where(listdates eq inputdates[i],cnt)
    if cnt eq 1 then begin
        if i eq 0 then begin
            whgood=wh
        endif else begin
            whgood=[whgood,wh]
        endelse
        ++count
    endif else begin
        print, 'Did not find file for date: ',dates[i]
    endelse
    
endfor

if count ne n_elements(dates) then begin
    print, 'did not find all of the requested dates'
endif

if count eq 0 then begin
    print, 'found none of requested dates.  quitting coadd_fits_maps!'
    return
endif

list=list[whgood]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot all maps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print, '   Plotting raw maps: number of maps = ', n_elements(list)
for ii=0, n_elements(list)-1 do begin
;for ii=0, 3 do begin
    data = spt_expand_smallmap(read_spt_fits( list[ii] ))
    tmp = strsplit(list[ii], '_',/extract)
    n = n_elements(tmp)
    plot_title = name+'_'+tmp[n-3]+'_'+tmp[n-2]+'_'+( strsplit(tmp[n-1],'.',/extract))[0]
    if (ii mod 50 eq 2 ) then print, '  plot_title for number [' + strtrim(string(ii),2) + ']: ' + plot_title

    if (scale eq -1) then begin
        tv_spt_map, data.map.map, title=plot_title+': MAP', winnum=1
    endif else begin
        tv_spt_map, data.map.map, scale=scale, /forcesize, title=plot_title+': MAP', winnum=1
    endelse
    fname='/data/kstory/projects/lps12/raw_maps/'+name+'/'+plot_title+'_map'
    err=tvread(/png, filename=fname, /nodialog) 

    if (scale eq -1) then begin
        tv_spt_map, data.weight.map, min=0., max=max(data.weight.map), title=plot_title+': WEIGHT', winnum=1
    endif else begin
        tv_spt_map, data.weight.map, min=0., max=max(data.weight.map), scale=scale, /forcesize, title=plot_title+': WEIGHT', winnum=1
    endelse
    fname='/data/kstory/projects/lps12/raw_maps/'+name+'/'+plot_title+'_weight'
    err=tvread(/png, filename=fname, /nodialog) 
endfor
end


;...................................................................
; Make a summary pdf of all raw maps
;
pro make_summary_pdf, field_list
compile_opt IDL2, HIDDEN

field_arr_ = lps12_fieldstruct()

for ii=0, n_elements(field_list)-1 do begin
    idx = field_list[ii]
    name = field_arr_[idx].name

    base_dir = '/data/kstory/projects/lps12/raw_maps/'
    pdf_name = base_dir + 'all_raw_maps_' + name + '.pdf'
    
    command = 'convert ' + base_dir + name + '/*.png ' + pdf_name

    print, '************* Execute: $ ' + command
    spawn, command
    
endfor
end

