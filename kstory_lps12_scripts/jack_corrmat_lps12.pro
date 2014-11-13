;;;
; NAME: jack_corrmat_lps12
; PURPOSE: Calculate the combined chisq for all jacks in lps12
;   
; NOTES:
; 1) Make 2d kweight for reso=0.5, npix=4320
;
; MODIFICATION HISTORY:
;  06/21/2012: (KTS) Created from /home/cr/code/spt/jacks/spt_run3_jack_cormat.pro, calc_s10_jacks_cormat.pro
;;;

FUNCTION jack_corrmat_2jacks,stub1,stub2,idx=idx,run=run

if n_elements(rand) ne 1 then rand = 0
if n_elements(seed) ne 1 then seed=7*(rand+1)

f = lps12_fieldstruct()
field = f[idx].name

sfreq='150'
info = get_lps12_fieldinfo(idx);spt_get_run3_fieldinfo(field,estub='_timeordered')
npix = info.npix
nbig = info.nbig

files = get_lps12_runlist(idx, /xspec_maps)

jackdef=get_lps12_jack_defs(idx,stub1,files)
setdef1=jackdef.setdef
jackdef=get_lps12_jack_defs(idx,stub2,files)
setdef2=jackdef.setdef

i1a=setdef1[*,0]
i2a=setdef2[*,0]
i1b=setdef1[*,1]
i2b=setdef2[*,1]
n1a=n_elements(i1a)
n2a=n_elements(i2a)
n1b=n_elements(i1b)
n2b=n_elements(i2b)
;stop
;stop
;define with respect to set 1:

q = [i1a,i1b]
q=q(sort(q))
z1=(n1b - 2*(n_elements(uniq(q))-n1a))/float(n1a)
q = [i1a,i2b]
q=q(sort(q))
z2=(n2b - 2*(n_elements(uniq(q))-n1a))/float(n1a)
q = [i2a,i1b]
q=q(sort(q))
z3=(n1b - 2*(n_elements(uniq(q))-n2a))/float(n2a)
q = [i2a,i2b]
q=q(sort(q))
z4=(n2b - 2*(n_elements(uniq(q))-n2a))/float(n2a)
RETURN,[max([z1,z2]),max([z3,z4])]

END


;-------------------------------------------------------
; Main function
FUNCTION jack_corrmat_lps12,minset=minset


jacks = ['12','azrms','tweight','moon']
njack = n_elements(jacks)
nfields = 19
idx_list = [0,1,indgen(17)+3]

corrmats=fltarr(nfields,njack+1,njack+1)
for i=0,njack-1 do for j=0,i-1 do begin
    for ifield=0,nfields-1 do begin
        idx = idx_list[ifield]
        print, 'Get corrmat for field ', idx
        xx=jack_corrmat_2jacks(jacks[i],jacks[j],idx=idx,run='07')
        corrmats[ifield,i,j]=xx[0]
        corrmats[ifield,j,i]=xx[0]
    endfor
endfor

for ifield=0, nfields-1 do begin
    for j=0, njack-1 do corrmats(ifield,j,j)=1.
endfor

; assume lr jack is uncorrelated
jacks=[jacks,'lr']
njacks = n_elements(jacks)
corrmats[*,njack,njack]=1

corrmat_tot = fltarr(njacks, njacks)
for ifield=0, nfields-1 do corrmat_tot += corrmats[ifield,*,*]
corrmat_tot = corrmat_tot / nfields

return,{corrmat:corrmat_tot,corrmats:corrmats,list:jacks}
END
