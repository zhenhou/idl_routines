;;;
; NAME: combine_jacks.pro
; PURPOSE:
;   function to combine jackknives
;
; INPUTS
;   jack,                 jackknife name, i.e. 'lr', '12' 'azrms'
;   fields,               optional list of fields
;   jack_dir,             directory where jackknife .sav files are
; NOTES:
;
; MODIFICATION HISTORY:
;  03/20/2012: (KTS) Created from /home/cr/code/spt/jackknives/spt_get_jack_k11.pro
;;;

FUNCTION combine_jacks,jack, $
                       fields=fields, $
                       jack_dir=jack_dir, $
                       pr=pr, $ ; print stuff
                       use_end2end=use_end2end, $ ; use output from end_to_end_powspec?
                       eqwt=eqwt,diagkernel=diagkernel,loff=loff,l250=l250,l1000=l1000

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()

;calib = 0.75^2
mask_dir = '/data/kstory/projects/lps12/masks/masks_50mJy/'
nbin=6 ; same as lps12_jack.pro

if ~keyword_set(jack_dir) then jack_dir = '/data/kstory/projects/lps12/jacks/'

; get a list of fields
if ~keyword_set(fields) then begin
    fields = strarr(20)
    for ii=0, 19 do begin
        fields[ii] = field_arr_[ii].name
    endfor
endif
nfield=n_elements(fields)

;;; Get kernel, etc from end-to-end powspec - > this is still to come
if keyword_set(use_end2end) then begin
    files = '/data/rkeisler/ps09/1.24/end_'+fields+'_1.24_kweight.sav'
    invkerns=dblarr(nbin,nbin,nfield)
endif

; arrays for all fields
diagcovs=fltarr(nbin,nfield)
cls=fltarr(nbin,nfield)
chisq = fltarr(nfield)
pte   = fltarr(nfield)

; combined field values
cov_tot = dblarr(nbin,nbin)
cl_tot  = dblarr(nbin)

; get field weights
restore, mask_dir+'field_weights.sav'
w = field_weights
; the moon jack requires a different weighting, since sets are
; significantly unbalanced
if jack eq 'moon' then begin
    restore, mask_dir+'moon_weights.sav'
    w = moon_weights
endif


;------------------------
;now for the jackknives
;------------------------

for ii=0,nfield-1 do begin
;for ii=0,0 do begin

    ; get the jackknife file (restores cl, cov)
    jack_file = jack_dir+'jack_150_'+fields[ii]+'_'+jack+'_info.sav'
    restore, jack_file
    cls[*,ii] = cl
    cov2use = double(cov)

    ; condition the covariance matrix
    offdiag = -1
    cov2use = condition_cov(cov2use, offdiag, diag=diag)
    for jj=0, nbin-1 do diagcovs[jj, ii] = cov2use[jj,jj]

    neff = n_elements(setdef)
    npair = n_elements(allspectra[0,*])
    
    ; calculate the chisq for each field
    icov2use = invert(cov2use)
    chisq[ii] = cls[*,ii] ## icov2use ## transpose(cls[*,ii])

    ; add these to cl_tot and cov_tot
    cl_tot  += cls[*,ii] * w[ii]
    cov_tot += cov2use * (w[ii])^2.

    if keyword_set(pr) then begin
        print, 'COMBINE_JACKS: restore ' + jack_file
        print,jack+' '+fields[ii],chisq
        print, 'pte - > ', 1-chisqr_pdf(chisq[ii],nbin)
        print, 'cl/err - > ',cls[*,ii]/sqrt(diagcovs[*,ii])
        print, 'cl - > ', cls[*,ii]
        print, 'err - ', sqrt(diagcovs[*,ii])
        print, ' '
    endif

endfor

plot,banddef-250.,cls[*,0]*1d12
oplot,banddef-250.,total(allspectra,2)*1./n_elements(allspectra[0,*])*1d12,color=!red

stop

; calculate the pte
icov_tot = invert(cov_tot)
chisq_tot = cl_tot ## icov_tot ## transpose(cl_tot)

pte = 1-chisqr_pdf(chisq_tot,nbin)

; calculate error on cl
cl_err = fltarr(nbin)
for jj=0, nbin-1 do cl_err[jj] = sqrt(cov_tot[jj, jj])

print,'comb fields chisq = ',chisq_tot, ' for bins = ',nbin
print,'pte = ',pte

return,{cl:cl_tot*1e12,err:1e12*cl_err,chisq:chisq_tot,pte:pte}
end
