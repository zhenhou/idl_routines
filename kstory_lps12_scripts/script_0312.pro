;;;
; NAME: script_0312
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) make_all_apodmasks_0312, exactly what it says
; 2) make_final_masks_0312, combine apod and ptsrc masks into final mask
; 3) check_apod_
; 4) plot_fullmap_overlap, plot apod masks in 2500 deg^2 footprint
; 5) pro run_lr_jack_0312
;
; MODIFICATION HISTORY:
;  03/12/2012: (KTS) Created
;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; apod masks
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;...................................................................
; Check overlap in apod masks between fields
pro make_all_apodmasks_0312
compile_opt IDL2, HIDDEN

for ii=0, 19 do begin
    make_apod_mask_0312, ii
endfor

end



;...................................................................
; make final apod*ptsrc_mask masks
pro make_final_masks_0312
compile_opt IDL2, HIDDEN

outdir = '/home/kstory/lps12/masks/'

field_arr = lps12_fieldstruct()

for ii=0, 19 do begin
    fst = field_arr[ii] & field_name = fst.name

    print, "make final mask for field ", field_name
    ; get the apod mask
    restore, '/home/kstory/lps12/masks/apod_0312/apod_'+field_name+'_60_0.0500_30.sav'

    ; get the ptsrc mask
    restore, '/home/kstory/lps12/masks/ptsrc_0305/ptsrc_mask_'+field_name+'_50mJy_proj5.sav'
    final_mask = apod*ptsrc_mask
    save, final_mask, filename=outdir+'masks_50mJy/mask_'+field_name+'.sav'
endfor

end


;...................................................................
; check difference between exp(-nmiss) and weight apod masks
pro check_apod_diff
compile_opt IDL2, HIDDEN

field_arr = lps12_fieldstruct()

for ii=0, 19 do begin
    fst = field_arr[ii] & field_name = fst.name

    print, "make final mask for field ", field_name
    ; get the apod mask
    restore, '/home/kstory/lps12/masks/apod_0312/apod_'+field_name+'_60_0.0500_30.sav'
    apod_w = apod

    restore, '/home/kstory/lps12/masks/apod_nmiss_0312/apod_'+field_name+'_60_0.0500_30.sav'
    apod_m = apod

    print, "apod_w = ", total(apod_w)
    print, "apod_m = ", total(apod_m)
endfor
end
    




;...................................................................
; Check overlap in apod masks between fields
pro plot_fullmap_overlap
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr = lps12_fieldstruct()
reso_arcmin = 1.0


;------------------------
; Make big map
;------------------------

radec_big = [20, -60]
npixels_big = [floor(960*8), floor(960*4)]
bigmap = fltarr(npixels_big[0], npixels_big[1])

pix2ang_proj5, npixels_big, radec_big, reso_arcmin, ra_big, dec_big ; for debugging

for ii=0, 19 do begin

    if ii eq 2 then ii=3 ; skip ra23h30dec-55_2010

    fst = field_arr[ii] & field_name = fst.name & npixels = get_lps12_map_npix(ii) & radec0 = [ fst.ra0, fst.dec0 ] ; Get the apod masks
    print, "PLOT_FULLMAP_OVERLAP: process field ", field_name
    restore, '/home/kstory/lps12/masks/apod_0312/apod_'+field_name+'_60_0.0500_30.sav'

    pix2ang_proj5, npixels, radec0, reso_arcmin, ra, dec
    
    ; insert apod into bigmap
    print, "loop ", ii
    npix_sm = long(npixels[0]) * long(npixels[1])
    for jj=0, npix_sm-1 do begin
        if (apod[jj] ne 0) then begin
            ang2pix_proj5, ra[jj], dec[jj], npixels_big, radec_big, 1.0, ipix
            bigmap[ipix] += apod[jj]
        endif
    endfor

endfor


;------------------------
; Plots
;------------------------

; Plot 1: show overlap greater than 1
ss = bigmap*0
wh = where(bigmap gt 1)
ss[wh] = 1
fname='/data/kstory/projects/lps12/figs/ss_0312.ps'
tv_spt_map, ss, /norms, scale=1.0, /forcesize, winnum=1, title='2500 deg^2, show where above 1', psfile=fname

; Plot 2: apod masks
;save plot
fname='/data/kstory/projects/lps12/figs/bigmap_0312.ps'
tv_spt_map, bigmap, /norms, scale=1.0, /forcesize, winnum=2, title='2500 deg^2 apod masks', psfile=fname

stop
end


; test ps plotting
pro test_ps
restore, '/home/kstory/lps12/masks/apod_0312/apod_ra21hdec-60_60_0.0500_30.sav'
;set_plot, 'ps'
;device, filename='/data/kstory/projects/lps12/figs/test_0312.ps'
fname='/data/kstory/projects/lps12/figs/test_0312.ps'
tv_spt_map, apod, /norms, scale=0.5, /forcesize, winnum=1, psfile=fname
;device, /close & set_plot, 'x'
stop
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Jackknives
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro run_lr_jack_0312
compile_opt IDL2, HIDDEN
;for ii=0, 19 do begin
for ii=1, 19 do begin
    jack_0312, ii
endfor
end

;...................................................................
; First attempt at jackknives
pro jack_0312, field_idx, $
               rand=rand, seed=seed, stopit=stopit
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = (field_arr_[field_idx])
field_name = fst.name

; make some noise
print, "JACK_0312: field ", field_name

sfreq='150'
freqi=0
freq=150

if n_elements(rand) ne 1 then rand = 0
if n_elements(seed) ne 1 then seed=7*(rand+1)

ddir='/data18/kstory/scratch/'

reso_arcmin=1.0
reso=reso_arcmin/60.*!dtor

info = get_lps12_fieldinfo(field_idx)
npix = info.npix
nbig = info.nbig

; get mask: apod+ptsrc
maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_'+field_name+'.sav'
restore, maskfile
mask = final_mask

nmask1=n_elements(mask[*,0])
nmask2=n_elements(mask[0,*])

maskb = fltarr(nbig,nbig)
maskb[0:nmask1-1,0:nmask2-1]=mask
mask=maskb


;;; apply kweight
kweight = get_lps12_kweight(field_idx)

; get all maps files
files = get_lps12_runlist(0, /xspec_maps)
mapfiles = [files]

nn = n_elements(files)
if nn lt 3 then stop

; get jackknife definitions
stub = 'lr'
jackdef=get_lps12_jack_defs(field_name,stub,files)
setdef=jackdef.setdef
dataname=jackdef.dataname
jackflag=jackdef.jackflag

tempdir=ddir+'xpec_'+field_name+'_'+stub+'_temp_'+strtrim(string(fix(freq)),2)

setdeforig=setdef
banddef = (1+findgen(6))*500.
setdeforig=setdef


;stop
unbiased_multiset_pspec, mapfiles, mask, reso_arcmin, $
  spec=spec, cov=cov, persistdir=tempdir, $
  maxell=3.2e3,mapname=dataname,$
  banddef=banddef, /resume, /cmbweighting,$
  setdef=setdef,fftrdonly=frd,allspectra=allspectra,kmask=kweight,$
  npix=npix,jackknife=jackflag,est1_cov=cov1, est2_cov=cov2

cl=spec
fstub=stub

; save output
fstub = 'lr'
sfile = '/home/kstory/lps12/jacks/0312/jack_'+sfreq+'_'+field_name+'_'+fstub+'_info.sav'

save,banddef,setdeforig,setdef,allspectra,cl,cov,cov1,cov2,filename=sfile

if keyword_set(stopit) then stop
end




;;;;;;;;;;;;;
; Make runlist of only obs that are missing
; obsolete

; pro make_missing_runlist
; idx = 3
; field_arr = lps12_fieldstruct()
; fst = field_arr[idx]
; field_name = fst.name
; field_dir_name = fst.dir_name

; maps  = get_lps12_runlist(3, /obs_maps)
; dates = get_lps12_runlist(3, /obs_dates)

; map_dir = '/data/kstory/projects/lps12/maps/20120305/ra21hdec-60/'

; runlist = fst.obs_runlist
; map_dir = fst.obs_map_dir

; ; Read the runlist of dates
; readcol,/silent,runlist,date_list,format='a'
; ndates_all = n_elements(date_list)

; files  = strarr(ndates_all)     ; map names
; has_map = intarr(ndates_all) + 1

; new_maps = ['']

; for ii=0, ndates_all-1 do begin
;     files[ii] = map_dir+'/map_'+field_dir_name+'_150_'+date_list[ii]+'.fits'
    
;     ; check for missing files
;     if ~file_test(files[ii]) then begin
;         print, "Missing file: ", files[ii]
;         new_maps = [new_maps, dates[ii]]
;     endif
; endfor

; new_maps = new_maps[1:*]
; stop

; ; print out new runlist
; new_runlist = '/data/kstory/projects/lps12/runlists/runlist_fixedIDFs_ra21hdec-60.txt'
; get_lun, lun1
; openw, lun1, new_runlist
; for ii=0, n_elements(new_maps)-1 do begin
;     printf, lun1, new_maps[ii]
; endfor
; close, lun1
; free_lun,lun1

; end
