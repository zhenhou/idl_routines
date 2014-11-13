;;;
; NAME: script_0308
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) try jackknives.  reference: ~cr/code/spt/jackknives/spt_lowell_jack.pro
; 2) make nmiss arrays
;
; MODIFICATION HISTORY:
;  03/08/2012: (KTS) Created
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Jackknives
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;...................................................................
; First attempt at jackknives
pro jack_0308;, field_idx
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_idx = 0
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name

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
restore, '/home/kstory/lps12/masks/apod_0307/apod_test_ra5h30dec-55_2008.sav'
restore, '/home/kstory/lps12/masks/ptsrc_0305/ptsrc_mask_ra5h30dec-55_2008_50mJy_proj5.sav'
;restore, maskfile
mask = apod*ptsrc_mask

nmask1=n_elements(mask[*,0])
nmask2=n_elements(mask[0,*])

maskb = fltarr(nbig,nbig)
maskb[0:nmask1-1,0:nmask2-1]=mask
mask=maskb


;;; apply kweight
;dell = 2*!pi/(nbig*reso)
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
  maxell=4.5e3,mapname=dataname,$
  banddef=banddef, /resume, /cmbweighting,$
  setdef=setdef,fftrdonly=frd,allspectra=allspectra,kmask=kweight,$
  npix=npix,jackknife=jackflag,est1_cov=cov1, est2_cov=cov2

cl=spec
fstub=stub

; save output
fstub = 'lr'
sfile = '/home/kstory/lps12/jacks/0308/jack_'+sfreq+'_'+field_name+'_'+fstub+'_info.sav'

save,banddef,setdeforig,setdef,allspectra,cl,cov,cov1,cov2,filename=sfile

stop
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; apod masks
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;...................................................................
; Check overlap in apod masks between fields
pro plot_overlap
compile_opt IDL2, HIDDEN

restore, '/home/kstory/lps12/scripts/apod_test1_0305.sav'

;------------------------
; Setup
;------------------------
field_arr = lps12_fieldstruct()
reso_arcmin = 1.0

ii = 0 ; ra5h30dec-55
fst1 = field_arr[ii] & field_name1 = fst1.name & npixels1 = get_lps12_map_npix(ii) & radec0_1  = [ fst1.ra0, fst1.dec0 ]
;apod1 = 1 & make_apod_mask_0306, ii, ret_apod=apod1
;restore, '/home/kstory/lps12/masks/apod_0305/apod_'+field_name1+'_60_0.0500_40.sav' & apod1 = apod

ii = 11 ; ra6h30dec-55
fst2 = field_arr[ii] & field_name2 = fst2.name & npixels2 = get_lps12_map_npix(ii) & radec0_2  = [ fst2.ra0, fst2.dec0 ]
;apod2 = 1 & make_apod_mask_0306, ii, ret_apod=apod2
;restore, '/home/kstory/lps12/masks/apod_0305/apod_ra6h30dec-55_60_0.0500_40.sav' & apod2 = apod

ii = 10 ; ra5h30dec-45
fst3 = field_arr[ii] & field_name3 = fst3.name & npixels3 = get_lps12_map_npix(ii) & radec0_3  = [ fst3.ra0, fst3.dec0 ]
;apod3 = 1 & make_apod_mask_0306, ii, overlap=[-1, 0., -1, 0.], ret_apod=apod3
apod3 = 1 & make_apod_mask_0306, ii, ret_apod=apod3
;restore, '/home/kstory/lps12/masks/apod_0305/apod_ra6hdec-62.5_60_0.0500_40.sav' & apod3 = apod

ii = 19 ; ra6h30dec-45
fst4 = field_arr[ii] & field_name4 = fst4.name & npixels4 = get_lps12_map_npix(ii) & radec0_4  = [ fst4.ra0, fst4.dec0 ]
;apod4 = 1 & make_apod_mask_0306, ii, overlap=[-1, 0., -1, 0.], ret_apod=apod4
apod4 = 1 & make_apod_mask_0306, ii, ret_apod=apod4
;restore, '/home/kstory/lps12/masks/apod_0305/apod_ra6hdec-62.5_60_0.0500_40.sav' & apod4 = apod

; convert to ra,dec
pix2ang_proj5, npixels1, radec0_1, reso_arcmin, ra1, dec1
pix2ang_proj5, npixels2, radec0_2, reso_arcmin, ra2, dec2
pix2ang_proj5, npixels3, radec0_3, reso_arcmin, ra3, dec3
pix2ang_proj5, npixels4, radec0_4, reso_arcmin, ra4, dec4


;------------------------
; Inset into big map
;------------------------

radec0 = [90, -50]
npixels = [floor(960*1.7), floor(960*1.7)]
bigmap = fltarr(npixels[0], npixels[1])

print, "loop 1"
npix_sm = long(npixels1[0]) * long(npixels1[1])
for ii=0, npix_sm-1 do begin
    if (apod1[ii] ne 0) then begin
        ang2pix_proj5, ra1[ii], dec1[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod1[ii]
    endif
endfor

print, "loop 2"
npix_sm = long(npixels2[0]) * long(npixels2[1])
for ii=0, npix_sm-1 do begin
    if (apod2[ii] ne 0) then begin
        ang2pix_proj5, ra2[ii], dec2[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod2[ii]
    endif
endfor

print, "loop 3"
npix_sm = long(npixels3[0]) * long(npixels3[1])
for ii=0, npix_sm-1 do begin
    if (apod3[ii] ne 0) then begin
        ang2pix_proj5, ra3[ii], dec3[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod3[ii]
    endif
endfor

print, "loop 4"
npix_sm = long(npixels4[0]) * long(npixels4[1])
for ii=0, npix_sm-1 do begin
    if (apod4[ii] ne 0) then begin
        ang2pix_proj5, ra4[ii], dec4[ii], npixels, radec0, 1.0, ipix
        bigmap[ipix] += apod4[ii]
    endif
endfor


;------------------------
; Plots
;------------------------

; Plot 1: show overlap greater than 1
ss = bigmap*0
wh = where(bigmap gt 1)
ss[wh] = 1
tv_spt_map, ss[*,150:*], /norms, scale=0.6, /forcesize, winnum=2, title='show where above 1'

; Plot 2: show where equal to 0
; ss0 = bigmap*0
; wh = where(bigmap le 0)
; ss0[wh] = 1
; tv_spt_map, ss0[*,150:*], /norms, scale=0.6, /forcesize, winnum=2, title='show where equal to 0'

; Plot 3: apod masks
tv_spt_map, bigmap[*,150:*], /norms, scale=0.6, /forcesize, winnum=1, title='5h-7h, dec-65 to -50'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_0.5overlap_0306'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_Nooverlap_0306'
figname = '/home/kstory/lps12/masks/figs_0305/apod_1deg_0306'
;err= tvread(/png, filename=figname, /nodialog)

; save output for easy recycling
save, apod1, apod2, apod3, bigmap, ss, filename='apod_test1_0305.sav'
stop
end

