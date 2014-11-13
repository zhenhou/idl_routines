;;;
; this program returns COV, which is an [NL,NL] matrix of 
; delta_CL/CL (NB: not delta_BL/BL!) taking into account all of the
; different beam errors.
;
; Used in combine_fields_lps12.pro
;
; MODIFICATION HISTORY:
; Copied from /home/rkeisler/ps09/get_beam_cov_matrix_lps12.pro
;  06/11/2012: (KTS) Change to use rev3.0 beams.
;  06/11/2012: (KTS) Change to use rev3.1 beams.
;;;

pro get_beam_cov_matrix_lps12, l=l, year=year, cov=cov, scalefactor=scalefactor


if n_elements(scalefactor) eq 0 then scalefactor=1.0

syear = year
if datatype(syear) ne 'STR' then syear = strtrim(floor(syear),2)

;beamdir = '/data/sptdat/beams/rev2.2/'
;beamdir = '/home/rkeisler/b11/rev3.0/'
beamdir = '/home/rkeisler/b11/rev3.1/'
types = ['inner','venus','alpha','dc','wobble','xtalk','outer']
names = beamdir+'errgrid_'+types+'_'+syear+'_150.txt'
ntypes = n_elements(types)

for i=0,ntypes-1 do begin
    readcol,/silent,names[i],l_grid,err_grid,format='d,d'
    this_err = interpol(err_grid,l_grid,l);delta_BL/BL
    this_err *= scalefactor
;    this_err = (1.+this_err)^2.-1. ; OLD/SLIGHTLY WRONG,convert to delta_CL/CL

    this_err = 1.-(1.+this_err)^(-2.) ; convert to delta_CL/CL
    this_cov = this_err # this_err
    if i eq 0 then begin
        nl = n_elements(l)
        cov = this_cov*0.
        errs = fltarr(ntypes,nl)
    endif
    errs[i,*] = this_err
    cov += this_cov
endfor


; now let's generate a bunch of fake experiments with the original 7
; beam types and calculate the actual covariance matrix and make sure
; it agrees with the one we just calculated.
if 0 then begin
nexp=1e4
seed=17
bigerrs = fltarr(nexp,nl)
for i=0,nexp-1 do begin
    amps = randomn(seed,ntypes)
    these_errs = errs
    for j=0,ntypes-1 do these_errs[j,*] *= (amps[j])
    this_err = total(these_errs,1)
    bigerrs[i,*] = this_err
endfor

cov2 = fltarr(nl,nl)
for i=0,nl-1 do begin
    for j=0,nl-1 do begin
        cov2[i,j] = mean((bigerrs[*,i]-0.)*(bigerrs[*,j]-0.))
    endfor
endfor
endif


;stop
end

