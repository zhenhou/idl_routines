;;;
; NAME: lps12_jack
; PURPOSE:
;   Run lps12 jackknives
;
; INPUTS:
;   field_idx,                index in lps12_fieldstruct
;   stub,                     jack name, i.e. 'lr', '12'
;   savdir,                   directory to save jackknifes in
;   sim,                      Set this to calculate the expected Dls from sims.
;
; NOTES:
;
; MODIFICATION HISTORY:
;  03/14/2012: (KTS) Created
;  05/11/2012: (KTS) Add ability to deal with sims
;  06/05/2012: (KTS) Change sims to use lmax8000
;  06/13/2012: (KTS) Add run keyword, now writes to jacks/run_XX/
;  06/16/2012: (KTS) Add ability to analyze randcut jacks
;;;


;...................................................................
; 
pro lps12_jack, field_idx, $
                stub, $          ; specify the jack
                run=run, $       ; Mandatory, specify the run
                sim=sim, $       ; use simulated lr jacks
                savdir=savdir, $
                mapdir=mapdir, $ ; read all maps from this directory
                dir_sim=dir_sim, $
                allowwrite=allowwrite, $
                save_temp=save_temp, $
                rand=rand, seed=seed, stopit=stopit
compile_opt IDL2, HIDDEN

if n_elements(run) eq 0 then begin
    print, 'LPS12_JACK: run keyword must be set.  Returning...'
    RETURN
endif

;------------------------
; Setup
;------------------------
if ~keyword_set(savdir) then begin
    savdir = '/home/kstory/lps12/jacks/run_'+run+'/'
    if keyword_set(sim) then savdir = '/home/kstory/lps12/jacks/sims/run_'+run+'/'
endif
print, 'savdir = ', savdir

f = lps12_fieldstruct()
field_name = f[field_idx].name
info = get_lps12_fieldinfo(field_idx)
npix = info.npix
nbig = info.nbig

reso_arcmin=1.0
reso=reso_arcmin/60.*!dtor

; make some noise
print, "LPS12_JACK: field ", field_name

sfreq='150'
freqi=0
freq=150

; read command line args
if n_elements(rand) ne 1 then rand = 0
if n_elements(seed) ne 1 then seed=7*(rand+1)
frd=1
if keyword_set(allowwrite) then frd=0

ddir='/home/kstory/lps12/scratch/' ; soft-linked to /data23/
if n_elements(dir_sim) eq 0 then dir_sim = '/data/kstory/projects/lps12/sims/run_'+run+'/'

;------------------------
; get mask: apod+ptsrc
mask_padded = get_lps12_mask(field_idx, /padded)


;------------------------
; apply kweight
kweight = get_lps12_kweight(field_idx)


;------------------------
; get all maps files
;------------------------

; Simulation mode
if keyword_set(sim) then begin
    dates = get_lps12_runlist(field_idx, /xspec_dates)
    files = file_search(dir_sim+'coaddsim_lmax8000_'+stub+'_'+field_name+'.dat')
    ;files = file_search(dir_sim+'coaddsim_'+stub+'_'+field_name+'.dat') ;; lmax4500
endif else begin

; Data mode
    if keyword_set(mapdir) then begin
        if (fst.lead_trail) then begin
            spawn, 'ls ' + mapdir+field_name+'_lt/*.fits', files
        endif else begin
            spawn, 'ls ' + mapdir+field_name+'/*.fits', files
        endelse
    endif else begin
        files = get_lps12_runlist(field_idx, /xspec_maps)

        ; Special case for randcut
        if (stub eq 'randcut_95_seed1_azrms') then files = get_lps12_runlist(field_idx, /xspec_maps, cut_name='randcut_95_seed1')
        if (stub eq 'randcut_95_seed2_azrms') then files = get_lps12_runlist(field_idx, /xspec_maps, cut_name='randcut_95_seed2')
        if (stub eq 'randcut_95_seed3_azrms') then files = get_lps12_runlist(field_idx, /xspec_maps, cut_name='randcut_95_seed3')

    endelse    
endelse
print, 'first map file: ' + files[0]

; --------reduce the number of files, for testing purposes
;files = files[0:5]
;-----------
mapfiles = [[files]]

nn = n_elements(files)
if (nn lt 3 and ~keyword_set(sim)) then stop

;------------------------
; get jackknife definitions
if ~keyword_set(sim) then begin
    jackdef=get_lps12_jack_defs(field_idx,stub,files)
    setdef=jackdef.setdef
    dataname=jackdef.dataname
    jackflag=jackdef.jackflag
    setdeforig=setdef
endif else begin
    dataname = (stub eq 'lr') ? 'DMAP.MAP' : 'MAP.MAP'
endelse


tempdir=ddir+'xpec_'+field_name+'_'+stub+'_'+dataname+'_temp_'+strtrim(string(fix(freq)),2)
if keyword_set(sim) then tempdir=ddir+'xpec_sim_'+field_name+'_'+stub+'_'+dataname+'_temp_'+strtrim(string(fix(freq)),2)
print, 'temp dir: ', tempdir

banddef = (1+findgen(6))*500.

;------------------------
; Calculate the jackknife spectra
;------------------------
if keyword_set(stopit) then stop

; simulated lr jacks
if keyword_set(sim) then begin
    nmapsperfile=100 ; sims have 100 maps.

    unbiased_multiset_pspec, mapfiles, mask_padded, reso_arcmin, $
      spec=spec, cov=cov, persistdir=tempdir, $
      maxell=3.0e3,$;mapname=dataname, $
      banddef=banddef, /cmbweighting,$
      fftrdonly=frd,allspectra=allspectra, kmask=kweight,$
      nmapsperfile=nmapsperfile, $
      npix=npix,est1_cov=cov1, est2_cov=cov2, $
      /auto, $
      resume=0

; real data
endif else begin
    unbiased_multiset_pspec, mapfiles, mask_padded, reso_arcmin, $
      spec=spec, cov=cov, persistdir=tempdir, $
      maxell=3.0e3,mapname=dataname,$
      banddef=banddef, /cmbweighting,$
      setdef=setdef,fftrdonly=frd,allspectra=allspectra, kmask=kweight,$
      npix=npix,jackknife=jackflag,est1_cov=cov1, est2_cov=cov2, $
      resume=0
endelse

dl=spec


;------------------------
; Save output
;------------------------

fstub=stub
if keyword_set(sim) then fstub +='_sim'

; save output
sfile = savdir+'jack_'+sfreq+'_'+field_name+'_'+fstub+'_info.sav'
print, 'Save jackknife in: ', sfile
save,banddef,setdeforig,setdef,allspectra,dl,cov,cov1,cov2,filename=sfile

; remove the temporary directory
if ~keyword_set(save_temp) then begin
    print, 'rm -rf ' + tempdir
    spawn, 'rm -rf ' + tempdir
endif

if keyword_set(stopit) then stop
end
