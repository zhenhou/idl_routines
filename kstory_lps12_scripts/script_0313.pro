;;;
; NAME: script_0313
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Check LT lists for ra21hdec-60
; 2) 
;
; MODIFICATION HISTORY:
;  03/13/2012: (KTS) Created
;;;


;...................................................................
; Check overlap in apod masks between fields
pro Check_lt
compile_opt IDL2, HIDDEN

restore, '/data/kstory/projects/lps12/runlists/lt/list_ra21hdec-60_150.sav'

nobs = n_elements(lead_list)
lead_dates  = extract_date_from_filename(lead_list)
trail_dates = extract_date_from_filename(trail_list)

js=read_jitterlist_config(fieldname='ra21hdec-60')
js_dates = js.date

L_sum = intarr(2)
T_sum = intarr(2)

for ii=0, nobs-1 do begin
    whL = where(js_dates eq lead_dates[ii])
    if js[whL].islead then begin 
        L_sum[0]++
    endif else L_sum[1]++

    whT = where(js_dates eq trail_dates[ii])
    if ~js[whT].islead then begin
        T_sum[0]++
    endif else T_sum[1]++
endfor

stop    


end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Jackknives
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; temp
;...................................................................
; First attempt at jackknives
pro play_jack, field_idx, $
               dataname=dataname, $
               rand=rand, seed=seed, stopit=stopit
compile_opt IDL2, HIDDEN

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
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
print, 'get maskfile, ', maskfile
restore, maskfile

nmask1=n_elements(mask[*,0])
nmask2=n_elements(mask[0,*])

maskb = fltarr(nbig,nbig)
maskb[0:nmask1-1,0:nmask2-1]=mask
mask=maskb


;;; apply kweight
kweight = get_lps12_kweight(field_idx)
;kweight = fltarr(nbig,nbig)+1

; get all maps files
files = get_lps12_runlist(field_idx, /xspec_maps)
; --------reduce the number of files, for testing purposes
; files = files[0:10]
;-----------
mapfiles = [[files]]

nn = n_elements(files)
if nn lt 3 then stop

; get jackknife definitions
stub = 'lr_kweightBox'
;jackdef=get_lps12_jack_defs(field_name,stub,files)

;---------------------
;;; DEBUGGING, make LR definition from map.map, which should fail the jackknife
jackflag=0
i1=lindgen(nn)
i2=0
if ~keyword_set(dataname) then begin
    ;dataname='D.DMAP.MAP'
    ;dataname='MAP.MAP'
    dataname='DMAP.MAP'
endif
setdef=[[i1]]
jackdef = {dataname:dataname,setdef:setdef,jackflag:jackflag}

;---------------------

setdef=jackdef.setdef
dataname=jackdef.dataname
jackflag=jackdef.jackflag

tempdir=ddir+'xpec_'+field_name+'_'+stub+'_'+dataname+'_temp_'+strtrim(string(fix(freq)),2)
;tempdir=ddir+'xpec_'+field_name+'_'+stub+'_temp_'+strtrim(string(fix(freq)),2)

setdeforig=setdef
banddef = (1+findgen(6))*500.
setdeforig=setdef
stop
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Look at take_and_reformat_ffts, in unbiased_multiset_pspec.pro
;...................................................................
; 
pro sss
;;; My setup ------------------------
restore, '/home/kstory/lps12/masks/masks_50mJy/mask_ra5h30dec-55_2008.sav'
window   = mask
mapfiles = get_lps12_runlist(0, /xspec_maps)
mapfiles = mapfiles[0:10]
mapname  = 'DMAP.MAP'
npix = get_lps12_map_npix(0)
kmask = dblarr(npix[0], npix[1])+1.
nmapsperfile = 1

winsize=size(window)
winsize=winsize[1]
reso_arcmin = 1.0
reso=double(reso_arcmin*!dtor/60)
deltau=1./reso/(winsize)
winfactor=1.0d

ell=shift(dindgen(winsize)-winsize/2, winsize/2)*deltau*2*!PI
ellgridx=(ell#replicate(1.0, winsize))
ellgrid=ellgridx^2+transpose(ellgridx)^2
ellgrid=sqrt(ellgrid)
ellsidx=sort(ellgrid, /L64)
;------------------------

; default ram limit: 16GB
if n_elements(ramlimit) eq 0 then ramlimit=16*(ulong64(2)^30)

winsize=(size(window))[1]

nfiles=n_elements(mapfiles)

;; number of bytes in a Dcomplex: 16
;; number of arrays we need to make to do this efficiently: 6 or less
;; number of pixels in an fft: winsize^2
comp_unit=16*6*ulong64(winsize)^2

parallelism=(fix(ramlimit/comp_unit) > 1)

bufferin=dblarr(winsize, winsize, parallelism)
smallbuffer=fltarr(npix[0], npix[1])
ellfact=ellgrid*(ellgrid+1)/2/!PI

scalingfact=kmask/winfactor

cmbweighting = 1
if keyword_set(cmbweighting) then begin
    scalingfact*=ellfact
endif

scalingfact=sqrt(abs(scalingfact))

fileidx=0
mapidx=0

print, 'DB1'
scaling3D=rebin(scalingfact, winsize, winsize, parallelism) 
print, 'DB2'
window3D=rebin(window, winsize, winsize, parallelism)
print, 'DB3'
; make a sorting-index array to sort multiple ffts by ell at once
ellsidx3D=rebin(ellsidx, ulong64(winsize)^2, parallelism)
; as it stands ellsidx maps all of the sorted arrays to the
; first array.  We need to add an offset to each row
; (ie add dindgen(parallelism)[i] to the ith row) of this array
print, 'DB4'
ellsidx3D+=transpose(rebin(ul64indgen(parallelism), parallelism, ulong64(winsize)^2))$
  *ulong64(winsize)^2

print, 'DB5'
code=reverse_linefeed_code()
print,''

totalmaps=nfiles*nmapsperfile
idx=0

stop
  
;while ((fileidx lt nfiles) and (mapidx lt nmapsperfile)) do begin
;end
end


;...................................................................
; run lr jackknives for remaining fields
pro run_lr_jack_9to13
for ii=9, 13 do begin
    test_jack, ii
endfor
end

pro run_lr_jack_14to19
for ii=14, 19 do begin
    test_jack, ii
endfor
end
