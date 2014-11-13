;;;
; NAME: Plot_lps12_jacks_wcorrmat.pro
; PURPOSE: Plot the combined jack chisq for lps12 jacks
;   
; NOTES:
; 1) Make 2d kweight for reso=0.5, npix=4320
;
; MODIFICATION HISTORY:
;  06/21/2012: (KTS) Created from /home/cr/code/spt/jacks/plot_s10_jacks_wcorrmat.pro
;;;

PRO plot_lps12_jacks_wcorrmat,offdiag=offdiag,skipcor=skipcor,minset=minset,usevar=usevar

sfreq='150'
run='07'
;corrres = calc_s10_jacks_corrmat(minset=minset)
corrres = jack_corrmat_lps12(minset=minset)
jacks = corrres.list
njack = n_elements(jacks)
icor = invert(corrres.corrmat)
if keyword_set(skipcor) then begin
icor *=0.
icor(indgen(njack),indgen(njack))=1.
endif



;neff *=1.2
keep=fltarr(njack)+1
;keep[0]=0
;keep[5]=0
;keep[6]=0
;keep[10]=0
;keep[11]=0
;keep[13]=0
normfac = fltarr(njack)+1.
normfac[1]=1
calibfac = .825^2
normfac *= calibfac
;ll = 250+500*findgen(30)
;jj=indgen(30)

covsig=0
;clsig=0
dlsig=0


;1000-10,000
i0=1
i1=5;was 19 for 10000
nbin=i1-i0+1

globalsn = fltarr(nbin,njack)

cumchisq=0
for i=0,njack-1 do begin

    ;aaa=spt_get_jack_s10(jacks[i],sfreq)
    aaa=read_jk_ret(jacks[i],run,/expected)
    chisq = aaa.chisq
    ;globalsn[*,i]=aaa.cl/aaa.err
    globalsn[*,i]=aaa.dl/aaa.err


    pte = 1-chisqr_pdf(chisq,nbin)
    print,jacks[i],pte,chisq

endfor

binchisq=fltarr(nbin)
for i=0,nbin-1 do begin
    binchisq[i]=globalsn[i,*] # icor # transpose(globalsn[i,*])
endfor

cumchisq=total(binchisq)
pte = 1-chisqr_pdf(cumchisq,n_elements(globalsn))
print,'cum PTE: ',pte
print,'cum chisq',cumchisq,n_elements(globalsn)

stop
END
