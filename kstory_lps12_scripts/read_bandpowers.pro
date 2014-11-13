;;;
; NAME: read_bandpowers.pro
; PURPOSE:
;   Easily plot bandpowers and errors for individual fields, from
;   output of end2end
;
; EXAMPLE CALL:
;   IDL> read_bandpowers, 0, run='04'
;
; INPUTS:
;   
; NOTES:
;  1) output from try_end2end.pro
;
; TODO:
;
; MODIFICATION HISTORY:
;  05/24/2012: (KTS) Created
;;;

PRO read_bandpowers, idx, run=run, use_kweight=use_kweight

;----------------------
; Setup
;----------------------
f = lps12_fieldstruct()
reso = 1.0

; directories
end_dir = '/home/kstory/lps12/end2end/'
rundir  = '/home/kstory/lps12/end2end/run_'+run+'/'
if (run eq '03') then begin  ; hack to still look at run_03
    end_dir = '/home/kstory/lps12/end2end/run_03/'
    rundir = '/home/kstory/lps12/end2end/run_03/'
endif

;----------------------
; get dls, errors
;----------------------

field = f[idx].name
fname = end_dir+'end_'+field+'_'+run
if keyword_set(use_kweight) then fname = fname + '_kweight'
fname = fname + '.sav'
print,idx, ',  '+fname
restore,fname

; variables
l = banddef-25.
nl = n_elements(l)
dls = dblarr(nl)
covs = dblarr(nl,nl)
corrs = dblarr(nl,nl)
uncorrs = dblarr(nl,nl)
diags = dblarr(nl)
diags2 = dblarr(nl)

nwf = n_elements(windowfunc[*,0])
wfs = windowfunc

area = total(mask_padded)*(reso/60.)^2.
cov = reform(cov)


; Condition the covariance matrix
nbins_off_diag = 5
cov_cond = condition_cov_matrix(cov, nbins_off_diag, $
                                corr=this_corr, uncondcorr=this_uncorr, $
                                banddef=banddef, $
                                ellmin=ellmin, ellmax=ellmax, $
                                noaverage=noaverage,$
                                knowncorr=knowncorr)

corrs = this_corr
uncorrs = this_uncorr

; Get the spectrum
dls=spectrum
covs=cov
these_diags2 = dblarr(nl)
for j=0,nl-1 do these_diags2[j] = sqrt(cov_cond[j,j])
diags2 = these_diags2

; cal hack
hack_factor=1.
dls *= (hack_factor)
;covs *= (hack_factor^2.)
;diags *= (hack_factor)



;----------------------
; plotting
;----------------------
xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8

min_l_plot = 700.
max_l_plot = 3000.
whplot = where(l ge min_l_plot and l le max_l_plot, n_plot)

; get dl_th
readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_th_uK2

plot,l[whplot],dls[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars
oplot,l_vec,dl_th_uK2,color=!red
errplot,l[whplot],(dls-diags2)[whplot]*1d12,(dls+diags2)[whplot]*1d12

stop
END

