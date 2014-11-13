;;;
; NAME: calculate_sample_variance_error
; PURPOSE: 
;   calculate the sample variance error bars directly from sims.
;
; NOTES:
;
; MODIFICATION HISTORY:
;  06/01/2012: (KTS) Created
;  06/08/2012: (KTS) Add sv_err_y for by-year sample variance
;  06/11/2012: (KTS) Calculate the full cov_sv
;  06/13/2012: (KTS) Change variable names to this_dl_all, this_banddef
;  06/29/2012: (KTS) Fix cov = X/(N) bug ; [ not X/(N-1) ]
;;;

; compare errors between run_03 and run_04
PRO calculate_sample_variance_error, run=run
;if n_elements(run) eq 0 then run='04'

;-------------------------
; setup
;-------------------------
f = lps12_fieldstruct()
end_dir = '/home/kstory/lps12/end2end/'
nbin = 58
nsim = 100
reso = 1.0

;-------------------------
; Get list of fields
;-------------------------
idx_list = [0,1,indgen(17)+3]
nfields = n_elements(idx_list)

; Variables for covariance matricies
nsets   = 1
nbands  = nbin
setsize = nsim
nspectra=(nsets*(nsets+1))/2

spec_sim    = dblarr(nbin, nsim, nfields) ; raw spectra
corspec_sim = dblarr(nbin, nsim, nfields) ; corrected spectra
invkern_sim = dblarr(nbin, nbin, nfields) ; inverse coupling kernel
area        = dblarr(nfields)
years = fltarr(nfields)


;-------------------------
; retrieve the information        
;-------------------------
for ilist=0, nfields-1 do begin
    idx = idx_list[ilist]
    print, 'process field ',idx, ', '+f[idx].name
    restore, end_dir+'end_'+f[idx].name+'_'+run+'_kweight.sav'

    ; get the spectra
    spec_sim[*,*, ilist] = all_mc_spectra
    invkern_sim[*,*, ilist] = reform(sim_invkern)
    area[ilist] = total(mask_padded)*(reso/60.)^2.
    years[ilist] = f[idx].year

    for isim=0, nsim-1 do begin
        corspec_sim[*,isim,ilist] = reform(invkern_sim[*,*,ilist]) ## transpose(spec_sim[*,isim, ilist])
    endfor
endfor

;-------------------------
; get field weights
;-------------------------
w0 = dblarr(nfields) ; weight array
for i=0,nfields-1 do w0[i] = area[i]
w0 /= total(w0)

; Get the field weights by year
yy=[2008, 2009, 2010, 2011]
nuniq_year = n_elements(yy)
wyear = dblarr(nuniq_year,nbin) ; weights by year
wh2008 = where(years eq 2008, n08)
wh2009 = where(years eq 2009, n09)
wh2010 = where(years eq 2010, n10)
wh2011 = where(years eq 2011, n11)
wyear[0] = (n08 ne 0) ? total(w0[wh2008]) : 0
wyear[1] = (n09 ne 0) ? total(w0[wh2009]) : 0
wyear[2] = (n10 ne 0) ? total(w0[wh2010]) : 0
wyear[3] = (n11 ne 0) ? total(w0[wh2011]) : 0
for j=0,nbin-1 do wyear[*,j] /= total(wyear[*,j])


;-------------------------
; Make arrays to be returned
;-------------------------
dl_y = dblarr(nuniq_year, nbin, nsim)
cov_sv=dblarr(nbands*nspectra, nbands*nspectra)
diag_sv = dblarr(nbin)
cov_sv_y=dblarr(nuniq_year, nbands*nspectra, nbands*nspectra)
diag_sv_y = dblarr(nuniq_year, nbin)

;-------------------------
; Loop over all years
;-------------------------
for iy=0, nuniq_year-1 do begin
    this_y = yy[iy]
    whyear = where(years eq this_y, nyear)
    print, yy[iy], ', nyear = ', nyear

    for jbin=0, nl-1 do begin
        these_wj = reform(w0[whyear])
        these_wj /= total(these_wj)

        for isim=0, nsim-1 do begin
            dl_y[iy, jbin, isim] = total( corspec_sim[jbin,isim,whyear]*these_wj )
        endfor
    endfor

    ; calculate the variance by year
    spectrumreformed=dblarr(nbands,nspectra)
    nrealizations=(setsize)
    
    allspectra=reform(dl_y[iy,*,*], nbands*nspectra, nrealizations)
    
    spectrum=total(/double,allspectra, 2)/nrealizations
    spectrum_2d=rebin(spectrum, nbands*nspectra, nrealizations)
    cov_sv1=(transpose(allspectra-spectrum_2d)##(allspectra-spectrum_2d))/nrealizations/(nrealizations-1)
    
    cov_sv_y[iy,*,*]=cov_sv1*(nrealizations)

    ; Get diagonals
    for i=0, nbin-1 do diag_sv_y[iy,i] = sqrt(cov_sv_y[iy,i,i])

endfor


;-------------------------
; get combined spectrum
;-------------------------
this_dl_all = dblarr(nbin, nsim)
this_banddef = banddef
for ilist=0, nfields-1 do begin
    this_dl_all += corspec_sim[*,*,ilist] * w0[ilist]
endfor    

;-------------------------
; calculate variance
;-------------------------
nspectra=(nsets*(nsets+1))/2
spectrumreformed=dblarr(nbands,nspectra)
nrealizations=(setsize)

allspectra=reform(this_dl_all, nbands*nspectra, nrealizations)

spectrum=total(/double,allspectra, 2)/nrealizations
spectrum_2d=rebin(spectrum, nbands*nspectra, nrealizations)
cov_sv1=(transpose(allspectra-spectrum_2d)##(allspectra-spectrum_2d))/nrealizations/(nrealizations-1)

cov_sv=cov_sv1*(nrealizations)

; make return array of diagonal errors
for i=0, nbin-1 do diag_sv[i] = sqrt(cov_sv[i,i])

;-------------------------
; save
;-------------------------
savfile = end_dir+'cov_sv_'+run+'_kweight.sav'
print, 'Saving errors in file: ' + savfile
save, diag_sv, cov_sv, diag_sv_y, cov_sv_y, filename=savfile

savfile1 = end_dir+'cov_sv_extra_'+run+'_kweight.sav'
save, diag_sv, cov_sv, diag_sv_y, cov_sv_y, this_dl_all, this_banddef, filename=savfile1

stop
END
