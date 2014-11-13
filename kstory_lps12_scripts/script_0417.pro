;;;
; NAME: script_0417
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) coadd jacklr_*.dat files
; 2) coadd run2 maps for ra1hdec-60
;
; MODIFICATION HISTORY:
;  04/17/2012: (KTS) Created
;;;

;...................................................................
; plot kmaps
PRO sss, $
         stopit=stopit
compile_opt IDL2, HIDDEN

;--------------------
; setup
;--------------------
f = lps12_fieldstruct()
idx = 7
field_name = f[idx].name

; get simulated map size
info = get_lps12_fieldinfo(idx)
npix = info.npix
nbig = info.nbig
reso_arcmin = 1.0

; location of simulated data.
sim_dir = '/data19/hou/lps12/sim_mapmaking/output/'+field_name+'/'
spawn, 'ls '+sim_dir + 'jacklr_*.dat', files
nfiles = n_elements(files)

; output sav file dir
outputdir = '/data/kstory/projects/lps12/sims/jacklr/'

; coadded map
cmap = fltarr(npix[0],npix[1])

;--------------------
; loop over all files
;--------------------

for ii = 0, nfiles-1 do begin

    ; make some noise
    ;print, "reading map ", files[ii]

    get_lun,u
    openr,u,files[ii]
    sim_map = fltarr(npix[0],npix[1])
    readu,u,sim_map
    free_lun,u

    ; make sav files to pass into unbiased_multiset_pspec
    file_date = extract_date_from_filename(files[ii])
    savname = outputdir+field_name+'/jacklr_150_'+file_date+'.sav'
    save, sim_map, filename=savname

    cmap += sim_map
endfor

cmap /= nfiles

; make a sav file to pass into unbiased_multiset_pspec
savname = '/home/kstory/lps12/scripts/sav_files/cmap_'+field_name+'.sav'
save, cmap, filename=savname

;------------------------
; get mask: apod+ptsrc
maskfile = '/home/kstory/lps12/masks/masks_50mJy/mask_'+field_name+'.sav'
print, 'get maskfile, ', maskfile
restore, maskfile

nmask1=n_elements(mask[*,0])
nmask2=n_elements(mask[0,*])

maskb = fltarr(nbig,nbig)
maskb[0:nmask1-1,0:nmask2-1]=mask
mask=maskb

; calculate Dl's
mapfiles = [[files]]
stop
unbiased_multiset_pspec, [[savname]], mask, reso_arcmin, $
  spec=spec, cov=cov, persistdir=tempdir, $
  maxell=3.0e3,mapname='cmap',$
  banddef=banddef, /cmbweighting,$
  kmask=kweight,$
  npix=npix,jackknife=0,est1_cov=cov1, est2_cov=cov2, $
  resume=0

;  setdef=setdef,fftrdonly=frd,allspectra=allspectra, 
stop
END


;...................................................................
; coadd full 0.25 arcmin maps for ra1hdec-60
PRO coadd_fullres, remake_runlist=remake_runlist

f = lps12_fieldstruct()
fst = f[9]
auto_dir = fst.autoprocessed_map_dirs

fname ='tmp_runlist.txt' 

if keyword_set(remake_runlist) then begin
    dates = get_lps12_runlist(9, /xspec_dates)
    ndates = n_elements(dates)

    get_lun, lun1
    openw, lun1, fname
    for ii=0, ndates-1 do begin

        spawn, 'ls '+auto_dir + 'map_'+fst.name+'_150_* | grep ' + dates[ii], map
        printf, lun1, map
    endfor
    close, lun1
    free_lun, lun1
endif

coadd_name = '/data/kstory/projects/lps12/maps/tmp/coadd_'+fst.name+'_fullres.fits'
print, 'COADD_MAPS_0420: output = ', coadd_name
coadd_fits_maps, fits_list=fname, fileout=coadd_name

stop
END
