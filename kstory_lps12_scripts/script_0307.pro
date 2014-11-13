;;;
; NAME: script_0302
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) try jackknives.  reference: ~cr/code/spt/jackknives/spt_lowell_jack.pro
; 2) make nmiss arrays
;
; MODIFICATION HISTORY:
;  03/07/2012: (KTS) Created
;;;

;...................................................................
; First attempt at jackknives
pro jack_0307;, field_idx
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


stop
unbiased_multiset_pspec, mapfiles, mask, reso_arcmin, $
  spec=spec, cov=cov, persistdir=tempdir, $
  maxell=4.5e3,mapname=dataname,$
  banddef=banddef, /resume, /cmbweighting,$
  setdef=setdef,fftrdonly=frd,allspectra=allspectra,kmask=kweight,$
  npix=npix,jackknife=jackflag,est1_cov=cov1, est2_cov=cov2




end







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Make nmiss arrays
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;...................................................................
; Make nmiss arrays for all fields
pro make_nmiss_0307
compile_opt IDL2, HIDDEN

for ii=0, 19 do begin

    ; skip ra21hdec-60
    if ii eq 3 then ii=4

    make_nmiss_array, ii

endfor
end


