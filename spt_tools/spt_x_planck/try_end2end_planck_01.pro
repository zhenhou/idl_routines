;;;
; NAME: 

pro try_end2end_planck_01, field_idx, banddef=banddef, spectrum=spectrum, jackhalf=jackhalf, resume=resume

freqs = 143
type  = 'nominal_ringfull'

f = lps12_fieldstruct()
fst = f[field_idx]
field_name = fst.name

reso_arcmin = 1.0
info = get_lps12_fieldinfo(field_idx)
npix = info.npix

; directories
mask_dir    = '/home/kstory/lps12/masks/masks_50mJy/'
kern_dir    = mask_dir
workdir     = '/data/hou/projects/spt_x_planck/powspec/planck/'+field_name+'/'
spawn, ['mkdir', workdir], /noshell
    
;;-----------------------
;; get Planck map files
;;-----------------------
mapfiles = get_planck_runlist(field_idx, freqs=freqs, type=type)
if (keyword_set(jackhalf)) then $
mapname = 'DMAP.MAP' else $
mapname = 'MAP.MAP'

;-------------------------
; get MASK and KERNEL
;   restores 'kernel' and 'mask_padded'
;-------------------------
kernfile = kern_dir+'kern_'+field_name+'.sav'
restore, kernfile, /verbose
mask_padded = get_lps12_mask(field_idx, /padded)

;-------------------------
; get BANDDEF
;-------------------------
delta_l = 50.
min_l = 250.
max_l = 3150.
nl = floor((max_l-min_l)/delta_l)
lhi = dindgen(nl)*delta_l + min_l
banddef = lhi
; this can possibly be much lower:
maxell = max(banddef)*1.5

;-------------------------
; get BEAMFILES
;-------------------------
year = fst.year
beamfiles = '/data23/hou/planck_data/2013/beams/hfi_beam_143x143_wpix_R1.10.txt'
print, 'beamfiles = ', beamfiles
print, 'beamfiles should be changed to planck ones'

;-------------------------
; if desired, get a KWEIGHT
;-------------------------
if keyword_set(use_kweight) then begin
    kmask = get_lps12_kweight(field_idx)
endif

print, "mapfiles - "
print, mapfiles

;-------------------------
; Run end_to_end
;-------------------------
end_to_end_powspec_planck, mapfiles, mask_padded, reso_arcmin, $
    banddef=banddef, mapname=mapname, $
    kernel=kernel, ellkern=ellkern, $
    maxell=maxell, $
    beamfiles=beamfiles, $
    pixfunc=pixfunc, $
    spectrum=spectrum, $
    npix=npix, workdir=workdir, $
    resume=resume, $ 
    curlyw=curlyw, $
    raw_spectrum_data=spectrum_data_raw, $
    cov_data_raw=cov_data_raw, $
    calib=calib
    
end
