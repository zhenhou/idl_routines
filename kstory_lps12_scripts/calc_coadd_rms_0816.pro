;;;
; NAME: calc_coadd_rms_0816
; PURPOSE:
;   Calculate the RMS of coadds for runlist wtrms cut tests
;
; MODIFICATION HISTORY:
;  08/16/2011: (KTS) Created
;;;

;...................................................................
; Mean
;
pro calc_coadd_rms_0816 
compile_opt IDL2, HIDDEN

;fields = ['ra23h30dec-55','ra21hdec-60','ra4h10dec-50','ra23hdec-62.5']
fields = ['ra21hdec-60']
cut_vals=['1.5', '2.0', '3.0']
cut_vals_num=[1.5, 2., 3.]
;field_list = [ 1, 2, 5, 11]

;window, 1, xsize=1200, ysize=1200
;!p.multi=[0,4,3]

for ii=0, n_elements(fields) - 1 do begin
    print, "Starting field", fields[ii]

    ; HIGH
    type = 'medwt'
    mrms=fltarr(3) & mmean=fltarr(3)
    for cv=0, 2 do begin ; cv = cut value
        coadd_fname = 'runlist_cuts_coadd_'+type+'_high'+cut_vals[cv]+'_low2.0_'+fields[ii]+'.fits'
        output = calc_coadd_rms_and_mean(coadd_fname)
        mrms[cv] = output[0] & mmean[cv] = output[1]
        stop
    endfor
    ;plot, cut_vals_num, mrms, title=type + '_high: RMS'
    ;plot, cut_vals_num, mmean, title=type + '_high: MEAN'
    
    print, "RMS, MEAN"
    print, "High: ", mrms, mmean


    ; LOW
    mrms=fltarr(3) & mmean=fltarr(3)
    for cv=0, 2 do begin ; cv = cut value
        coadd_fname = 'runlist_cuts_coadd_'+type+'_high2.0_low'+cut_vals[cv]+'_'+fields[ii]+'.fits'
        output = calc_coadd_rms_and_mean(coadd_fname)
        mrms[cv] = output[0] & mmean[cv] = output[1]
    endfor
    ;plot, cut_vals_num, mrms, title=type + '_low: RMS'
    ;plot, cut_vals_num, mmean, title=type + '_low: MEAN'
    print, "Low:  ", mrms, mmean
    
endfor
stop

;!p.multi=0
end

;...................................................................
; Calculate the dmap rms and mean
;
function calc_coadd_rms_and_mean, coadd_fname
compile_opt IDL2, HIDDEN

print, "  *** Restore File"
data = read_spt_fits('/data/kstory/projects/lps12/runlists/coadds/'+coadd_fname)

print, "  *** Calculate rms and mean"
dmap = data.dmap.map
npix = [(size(dmap))[1], (size(dmap))[2] ]
meshgrid, npix[0], npix[1], xx, yy

xmax = round(npix[0]*(2./3)) & xmin = round(npix[0]*(1./3))
ymax = round(npix[1]*(2./3)) & ymin = round(npix[1]*(1./3))

wh = where( (xx lt xmax) and (xx gt xmin) and (yy lt xmax) and (yy lt xmin) , nwh)
cmap = dmap[wh]

return, [rms(cmap), mean(cmap)]
;stop
end

;...................................................................
; Calculate the dmap rms and mean
;
pro make_coadd_dmap, runlist_fname, dmap
compile_opt IDL2, HIDDEN

end
