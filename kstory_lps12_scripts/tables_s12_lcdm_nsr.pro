;;;
; NAME: tables_s12_lcdm_nsr
; PURPOSE:
;   Print numbers for tables in L12
;
; NOTES:
; 1) lcdm_table3,            Use for LCDM parameter constraints
; 2) ml_lcdm_table3,         ML values for LCDM parameters
; 3) lcdm_inflation,         Use for Inflation constraints
; 4) lcdm_r,                 Inflation
; 5) lcdm_nrun,              Inflation
; 6) lcdm_nrun_r,            Inflation
; 7) fields_table,           Table with fields information

; MODIFICATION HISTORY:
;  08/06/2012: (KTS) Created
;  09/10/2012: (KTS) Use 0828 chains
;  09/12/2012: (KTS) Add procedures: print_lcdm_table3, alens_constraint
;;;


;;;;;;;;;;;;;;;;;;;;;;
; LCDM parameters
;   LCDM table
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_table3, param=param, scale=scale
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

;name = 'ns'
subsamp=100000
nskip = 1000
if n_elements(scale) eq 0 then scale=1.

; S12-only
dir_spt   = cdir+'c27_lcdm_camb_s12tau/chains/'
files_spt = file_search(dir_spt+'c27_lcdm_camb_s12tau*.txt')
pname_spt = dir_spt+'c27_lcdm_camb_s12tau.paramnames'

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7.paramnames'

; cmb
dir_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c2_lcdm_pico_w7s12*.txt')
pname_cmb = dir_cmb+'c2_lcdm_pico_w7s12.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c5_lcdm_pico_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c5_lcdm_pico_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c5_lcdm_pico_w7s12_BAO_H0.paramnames'

print, '' & print, '************ WMAP *****************' & print, pname_w7
plot_like1dname,files_w7,pname_w7,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_w7
print, '' & print, '************ SPT *****************' & print, pname_spt
plot_like1dname,files_spt,pname_spt,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_spt
print, '' & print, '************ CMB *****************' & print, pname_cmb
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
print, '' & print, '************ CMB+H0+BAO *****************' & print, pname_cmb_h0_bao
plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'param = ', param
print, 'WMAP: ', results_w7.median, results_w7.err_plus, results_w7.err_minus, results_w7.avg_err
print, 'SPT: ', results_spt.median, results_spt.err_plus, results_spt.err_minus, results_spt.avg_err
print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END



;;;;;;;;;;;;;;;;;;;;;;
; LCDM parameters
;   Print the full LCDM table
;;;;;;;;;;;;;;;;;;;;;;
PRO print_lcdm_table3, use_oldAs=use_oldAs

subsamp=100000
nskip = 1000
extra = ''

; Set the parameters in the table
pset=['omegabh2','omegadmh2','1e-9As','ns','theta_s','tau',  'omegal*','H0*','sigma8*','z_eq','Dvp35ors','Dvp57ors'] ; names in chain.paramname
pformat = ['F5.3', 'F6.4', 'F5.3', 'F6.4', 'F7.5', 'F5.3',   'F6.4', 'F5.2', 'F5.3',   'F5.0','F5.2',    'F6.3'] ; print format
pprint_name = ['$10^2\Omega_bh^2$','$\Omega_ch^2$    ','$10^{9} \deltaR$ ','$n_s$            ','$10^2\theta_s$   ','$\tau$           ', $
               '$\Omega_\Lambda$ ','\ho              ','$\sigma_8$       ','$z_{\rm EQ}$     ','$r_s/D_V(z=0.35)$', '$r_s/D_V(z=0.57)$']
dv_fac = (150.82/154.66)/100.
pscale  = [100., 1, 1, 1, 100., 1, 1, 1, 1, 1, dv_fac, dv_fac] ; scale
ppower  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1., -1.] ; scale
nparams = n_elements(pset)


; Get the chains
cdir = '/data23/kstory/lps12/chains/'
chains  = ['c1_lcdm_pico_w7_newKp','c27_lcdm_camb_s12tau_newKp','c2_lcdm_pico_w7s12_newKp','c4_lcdm_pico_w7s12_H0_newKp','c3_lcdm_pico_w7s12_BAO_newKp','c5_lcdm_pico_w7s12_BAO_H0_newKp']
;chains  = ['c5_lcdm_pico_w7s12_BAO_H0_newKp']

;;; Obsolete, For k0=0.002 Mpc^-1
if keyword_set(use_oldAS) then begin
    cdir = '/data23/hou/lps12/paramfits/chains_0828/'
    chains  = ['c1_lcdm_pico_w7', 'c27_lcdm_camb_s12tau', 'c2_lcdm_pico_w7s12', 'c4_lcdm_pico_w7s12_H0','c3_lcdm_pico_w7s12_BAO','c5_lcdm_pico_w7s12_BAO_H0']
    extra = '_oldAs'
endif

nchains = n_elements(chains)
dirs    = cdir+chains+'/chains/'
nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
ndsets = n_elements(chains)

; ouput matricies
param_med = dblarr(ndsets, nparams)
param_err = dblarr(ndsets, nparams)

; ; Get parameter indicies in the chains
; pindex = intarr(ndsets, nparams)
; for i=0, ndsets-1 do begin
;     readcol,pnames[i], nn,n2,format='a,a'
;     for j=0, nparams-1 do pindex[i,j] = where(nn eq pset[j])
; endfor

; get ML values
;ml = {ml0:get_ml_chain_all(files.f0), ml1:get_ml_chain_all(files.f1), ml2:get_ml_chain_all(files.f2), ml3:get_ml_chain_all(files.f3), ml4:get_ml_chain_all(files.f4), ml5:get_ml_chain_all(files.f5)}

; Get parameter values
for i=0, ndsets-1 do begin
    for j=0, nparams-1  do begin
        param = pset[j]
        print, pnames[i]
        plot_like1dname,files[i].ff,pnames[i],param,subsamp=subsamp,nskip=nskip,scale=pscale[j],power=ppower[j],results=results
        param_med[i,j] = results.median
        param_err[i,j] = results.avg_err
    endfor
endfor




;-------------------------
; PRINT the table
;-------------------------
; openw,lun,/get_lun,'tmp_lcdm_tab3.txt'
; for j=0, nparams-1 do $;begin
;     printf, lun, format='("'+pprint_name[j]+'       & ",'+ $
;       pformat[j]+',"  & $",'+pformat[j]+'," \pm ",'+pformat[j]+',"$ &    ",'+$
;       pformat[j]+',"  & $",'+pformat[j]+'," \pm ",'+pformat[j]+',"$ &    ",'+$
;       pformat[j]+',"  & $",'+pformat[j]+'," \pm ",'+pformat[j]+',"$ &    ",'+$
;       pformat[j]+',"  & $",'+pformat[j]+'," \pm ",'+pformat[j]+',"$ &    ",'+$
;       pformat[j]+',"  & $",'+pformat[j]+'," \pm ",'+pformat[j]+',"$ \\")',$
;       ml.ml0[pindex[0,j]+2]*pscale[j], param_med[0,j], param_err[0,j], $
;       ml.ml1[pindex[1,j]+2]*pscale[j], param_med[1,j], param_err[1,j], $
;       ml.ml2[pindex[2,j]+2]*pscale[j], param_med[2,j], param_err[2,j], $
;       ml.ml3[pindex[3,j]+2]*pscale[j], param_med[3,j], param_err[3,j], $
;       ml.ml4[pindex[4,j]+2]*pscale[j], param_med[4,j], param_err[4,j]
; close, lun                         


; No ML values

; --- DEBUGGING ----
if nchains eq 1 then begin
openw,lun,/get_lun,'tmp.txt'
for j=0, nparams-1 do $
    printf, lun, format='("'+pprint_name[j]+'       & $",'+ $
      pformat[j]+'," \pm ",'+pformat[j]+',"$ \\")',$
      param_med[0,j], param_err[0,j]
close, lun                         
;--------

endif else begin
;openw,lun,/get_lun,'tmp_lcdm_tab3'+extra+'.txt'
openw,lun,/get_lun,'tmp_a'+extra+'.txt'
for j=0, nparams-1 do $
    printf, lun, format='("'+pprint_name[j]+'     & $",'+ $
      pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
      pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
      pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
      pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
      pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
      pformat[j]+'," \pm ",'+pformat[j]+',"$ \\")',$
      param_med[0,j], param_err[0,j], $
      param_med[1,j], param_err[1,j], $
      param_med[2,j], param_err[2,j], $
      param_med[3,j], param_err[3,j], $
      param_med[4,j], param_err[4,j], $
      param_med[5,j], param_err[5,j]
close, lun                         
endelse

stop
END



;;;;;;;;;;;;;;;;;;;;;;
; LCDM parameters
;    Currently unused
;;;;;;;;;;;;;;;;;;;;;;
PRO ml_lcdm_table3
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; S12-only
dir_spt   = cdir+'c27_lcdm_camb_s12tau/chains/'
files_spt = file_search(dir_spt+'c27_lcdm_camb_s12tau*.txt')
pname_spt = dir_spt+'c27_lcdm_camb_s12tau.paramnames'

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7.paramnames'

; cmb
dir_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c2_lcdm_pico_w7s12*.txt')
pname_cmb = dir_cmb+'c2_lcdm_pico_w7s12.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c5_lcdm_pico_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c5_lcdm_pico_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c5_lcdm_pico_w7s12_BAO_H0.paramnames'

print, '' & print, '************ WMAP *****************'
print_ml_chain_all, files_w7
print, '' & print, '************ SPT *****************'
print_ml_chain_all, files_spt
print, '' & print, '************ CMB *****************'
print_ml_chain_all, files_cmb
print, '' & print, '************ CMB+H0+BAO *****************'
print_ml_chain_all, files_cmb_h0_bao
END




;;;;;;;;;;;;;;;;;
; Table 2, bandpowers
;;;;;;;;;;;;;;;;;
PRO tab_bandpower
restore, '/home/kstory/lps12/end2end/run_09/combined_spectrum_20120828_170101_kweight.sav'

; calculate l_effective
leff = fltarr(58)
l2 = fltarr(58)
;for i=0, 57 do leff[i] = total(wf_all[*,i]*l_wf) ;/ total(wf_all[*,i])
;leff = l_wf # wf_all

for i=0, 58-1 do begin
    wh = where(l_wf ge l[i]-25 and l_wf lt l[i]+25)
    leff[i] = total(wf_all[wh,i] * l_wf[wh])
endfor
stop

istart=9
istop=55
dl_all_lps12 = dl_all[istart:istop] * 1d12
diag_nobeam_lps12 = diag_nobeam[istart:istop] * 1d12
l = l[istart:istop]
leff = leff[istart:istop]

nbin = n_elements(l)
n2bin = nbin/2
for i=0, n2bin-1 do begin
    i2 = i+n2bin+1
    print, strtrim(string(round(l[i])-24),2), ' - ', strtrim(string(round(l[i])+25),2), ' & ', round(leff[i]), ' & ', dl_all_lps12[i], ' & ', diag_nobeam_lps12[i], ' & ', $
      strtrim(string(round(l[i2])-24),2), ' - ', strtrim(string(round(l[i2])+25),2), ' & ', round(leff[i2]), ' & ', dl_all_lps12[i2], ' & ', diag_nobeam_lps12[i2], ' \\', $
      format = '(A4, A3, A5, A3, i6, A3, f6.1, A3, f4.1, A3,   A6, A3, A5, A3, i6, A3, f6.1, A3, f4.1, A3)'

endfor
i=n2bin
print, strtrim(string(round(l[i])-24),2), ' - ', strtrim(string(round(l[i])+25),2), ' & ', round(leff[i]), ' & ', dl_all_lps12[i], ' & ', sigfig(diag_nobeam_lps12[i], 3), ' & & & \\', $
  format = '(A4, A3, A5, A3, i6, A3, f6.1, A3, f4.1, A9)'

stop
END




;;;;;;;;;;;;;;;;;;;;;;
; Get CDF limits on ns
; Power-law [LCDM]
;   Print constraints for all LCDM parameters using this procedure
;   type in ['s12', 'w7', 'cmb', 'extra']
;;;;;;;;;;;;;;;;;;;;;;
PRO ns_cdf, type=type
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

if n_elements(type) eq 0 then type='cmb'
temperature = 1.0 ; use for high-T chains

;;; LCDM
case type of
    's12' : begin
        ; S12-only
        dir   = cdir+'c27_lcdm_camb_s12tau/chains/'
        files = file_search(dir+'c27_lcdm_camb_s12tau*.txt')
        pname = dir+'c27_lcdm_camb_s12tau.paramnames'
    endcase

    'w7' : begin
        ; WMAP-only
        dir = cdir+'c1_lcdm_pico_w7/chains/'
        files = file_search(dir+'c1_lcdm_pico_w7*.txt')
        pname = dir+'c1_lcdm_pico_w7.paramnames'
    endcase

    'cmb' : begin
        ; cmb
        dir = cdir+'c2_lcdm_pico_w7s12/chains/'
        files = file_search(dir+'c2_lcdm_pico_w7s12*.txt')
        pname = dir+'c2_lcdm_pico_w7s12.paramnames'
    endcase

    'cmb_h0' : begin
        ; CMB+H0
        dir = cdir+'c4_lcdm_pico_w7s12_H0/chains/'
        files = file_search(dir+'c4_lcdm_pico_w7s12_H0*.txt')
        pname = dir+'c4_lcdm_pico_w7s12_H0.paramnames'
    endcase

    'cmb_h0_highT' : begin
        ; CMB+H0
        dir = cdir+'c98_lcdm_pico_w7s12_H0/chains/'
        files = file_search(dir+'c98_lcdm_pico_w7s12_H0*.txt')
        pname = dir+'c98_lcdm_pico_w7s12_H0.paramnames'
        temperature = 6.0
    endcase

    'cmb_bao' : begin
       ; CMB+BAO 
        dir = cdir+'c3_lcdm_pico_w7s12_BAO/chains/'
        files = file_search(dir+'c3_lcdm_pico_w7s12_BAO*.txt')
        pname = dir+'c3_lcdm_pico_w7s12_BAO.paramnames'
    endcase

    'cmb_bao_highT' : begin
       ; CMB+BAO 
        dir = cdir+'c109_lcdm_pico_w7s12_BAO/chains/'
        files = file_search(dir+'c109_lcdm_pico_w7s12_BAO*.txt')
        pname = dir+'c109_lcdm_pico_w7s12_BAO.paramnames'
        temperature = 6.0 
    endcase

    'cmb_h0_bao' : begin
        ; CMB+H0+BAO
        dir = cdir+'c6_lcdm_pico_w7s12_BAO_H0/chains/'
        files = file_search(dir+'c6_lcdm_pico_w7s12_BAO_H0*.txt')
        pname = dir+'c6_lcdm_pico_w7s12_BAO_H0.paramnames'
        temperature = 6.0 
    endcase

    'cmb_h0_neff': begin
        type = 'c14_lcdm_neff_pico_w7s12_H0'
        dir = cdir+type+'/chains/'
        files = file_search(dir+type+'_*.txt')
        pname = dir+type+'.paramnames'
    endcase

    'cmb_h0_bao_neff': begin
        type = 'c15_lcdm_neff_pico_w7s12_BAO_H0'
        dir = cdir+type+'/chains/'
        files = file_search(dir+type+'_*.txt')
        pname = dir+type+'.paramnames'
    endcase

    else : begin
        print, 'un-recognized case. Returning'
        RETURN
    endcase

endcase


;;; LCDM+r
; ; cmb
; dir_cmb = cdir+'c26_lcdm_r_camb_w7s12/chains/'
; files_cmb = file_search(dir_cmb+'c26_lcdm_r_camb_w7s12_*.txt')
; pname_cmb = dir_cmb+'c26_lcdm_r_camb_w7s12.paramnames'

; ; CMB+BAO
; dir_cmb_bao = cdir+'c28_lcdm_r_camb_w7s12_BAO/chains/'
; files_cmb_bao = file_search(dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO_*.txt')
; pname_cmb_bao = dir_cmb_bao+'c28_lcdm_r_camb_w7s12_BAO.paramnames'

print, 'NS constraints for dataset, ' + type
subsamp=25
;subsamp=10000
nskip=1000
plot_like1dname,files,pname,'ns',subsamp=subsamp,nskip=nskip,scale=scale,temp=temperature,/cdf,gtx=1,/expinterp,/stopit

stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Get CDF limits on alens
;   type in ['s12', 'cmb']
;;;;;;;;;;;;;;;;;;;;;;
PRO alens_cdf, type=type
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

if n_elements(type) eq 0 then type='cmb'
temperature = 6.0 ; both are high-T chains

;;; LCDM
case type of
    's12' : begin
        ; S12-only
        dir   = cdir+'c53_lcdm_alens_camb_s12tau/chains/'
        files = file_search(dir+'c53_lcdm_alens_camb_s12tau*.txt')
        pname = dir+'c53_lcdm_alens_camb_s12tau.paramnames'
    endcase

    'cmb' : begin
        ; cmb
        dir = cdir+'c55_lcdm_alens_camb_w7s12/chains/'
        files = file_search(dir+'c55_lcdm_alens_camb_w7s12*.txt')
        pname = dir+'c55_lcdm_alens_camb_w7s12.paramnames'
    endcase

    else : begin
        print, 'un-recognized case. Returning'
        RETURN
    endcase
endcase

print, 'Alens constraints for dataset, ' + type
subsamp=25
nskip=1000
nparam = 14
;plot_like1dname,files,pname,'alens',subsamp=subsamp,nskip=nskip,scale=scale,temp=temperature,/cdf,ltx=0,/stopit
cdf_highT,files,nparam,nskip=nskip,temp=temperature,ltx=0,/stopit

;; Run the following on the command line when the program stops
;wh = where(pp ge 1.0, nwh) & frac = nwh/float(n_elements(pp)) & print, frac & print, gauss_cvf(frac)
;wh = where(bins ge 1.0) & print, ff[wh[0]] & print, gauss_cvf(ff[wh[0]])

stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Get constraint on alens^0.6
;;;;;;;;;;;;;;;;;;;;;;
PRO alens_constraint
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; cmb
dir = cdir+'c54_lcdm_alens_camb_w7s12/chains/'
files = file_search(dir+'c54_lcdm_alens_camb_w7s12*.txt')
pname = dir+'c54_lcdm_alens_camb_w7s12.paramnames'

print, 'Alens constraints for dataset, ' + pname
subsamp=25
nskip=1000
;power = 0.5
;plot_like1dname,files,pname,'alens',subsamp=subsamp,nskip=nskip,scale=scale,power=power,results=results
;print, 'Alens^'+strtrim(string(power),2)+' = ', results.median, ' \pm ', results.avg_err

plot_like1dname,files,pname,'alens',subsamp=subsamp,nskip=nskip,scale=scale,results=results
print, 'Alens = ', results.median, ' + ', results.err_plus, '  - ', results.err_minus

stop
END






;;;;;;;;;;;;;;;;;;;;;;
; Get CDF limits on oml
;;;;;;;;;;;;;;;;;;;;;;
PRO oml_cdf, type=type
cdir = '/data23/kstory/lps12/chains/'

if n_elements(type) eq 0 then type='cmb'
temperature = 1.0 ; use for high-T chains
nparam = 26 ; for ns
scale=1

;;; LCDM
case type of
    'cmb' : begin
        ; cmb
        dir = cdir+'c220_t8_lcdm_omk_camb_w7s12/chains/'
        files = file_search(dir+'c220_t8_lcdm_omk_camb_w7s12*.txt')
        pname = dir+'c220_t8_lcdm_omk_camb_w7s12.paramnames'
        temperature = 8.0
    endcase

    else : begin
        print, 'un-recognized case. Returning'
        RETURN
    endcase

endcase


print, 'OML constraints for dataset, ' + type, ', temperature = ', temperature
nskip=1000
;plot_like1dname,files,pname,'omegal*',subsamp=subsamp,nskip=nskip,scale=scale,temp=temperature,/cdf,gtx=1,/expinterp,/stopit
cdf_highT,files,nparam,nskip=nskip,temp=temperature,ltx=0,/stopit
stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Get constraint on omk
;;;;;;;;;;;;;;;;;;;;;;
PRO omk_constraint, type=type, param=param

if n_elements(param) eq 0 then param='omegak'

cdir = '/data23/hou/lps12/paramfits/chains_0828/'
case type of
    'cmb' : begin
        chain = 'c46_lcdm_omk_camb_w7s12'
        dir = cdir+chain+'/chains/'
        files = file_search(dir+chain+'*.txt')
        pname = dir+chain+'.paramnames'
    endcase

    'cmb_h0' : begin
        chain = 'c48_lcdm_omk_camb_w7s12_H0'
        dir = cdir+chain+'/chains/'
        files = file_search(dir+chain+'*.txt')
        pname = dir+chain+'.paramnames'
    endcase

    'cmb_bao' : begin
        chain = 'c47_lcdm_omk_camb_w7s12_BAO'
        dir = cdir+chain+'/chains/'
        files = file_search(dir+chain+'*.txt')
        pname = dir+chain+'.paramnames'
    endcase

    'cmb_h0_bao' : begin
        chain = 'c120_lcdm_omk_camb_w7s12_BAO_H0'
        cdir = '/data/cr/chains_lps12/chains_20121005/'
        dir = cdir
        files = file_search(dir+chain+'*.txt')
        pname = dir+chain+'.paramnames'
    endcase

    else : begin
        print, 'un-recognized case. Returning'
        RETURN
    endcase

endcase


print, param+' constraints for dataset, ' + pname
subsamp=100000
nskip=1000
plot_like1dname,files,pname,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results
print, param+' = ', results.median, ' \pm ', results.avg_err
print, param+' = ', results.median, ' + ', results.err_plus, ' - ', results.err_minus

END



;;;;;;;;;;;;;;;;;;;;;;
; Get CDF limits on nrun
; Power-law [LCDM]
;   Print constraints for all LCDM parameters using this procedure
;   type in ['s12', 'w7', 'cmb', 'extra']
;;;;;;;;;;;;;;;;;;;;;;
PRO nrun_cdf, type=type, use_r=use_r
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; CMB+H0+BAO
dir = cdir+'c10_lcdm_nrun_pico_w7s12_BAO_H0/post_chains/'
files = file_search(dir+'c10_lcdm_nrun_pico_w7s12_BAO_H0*.txt')
pname = dir+'c10_lcdm_nrun_pico_w7s12_BAO_H0.paramnames'

subsamp=10000
nskip=5000
stale=1
plot_like1dname,files,pname,'nrun',subsamp=subsamp,nskip=nskip,scale=scale,/cdf,/stopit

;; Run the following on the command line when the program stops
;wh = where(bins ge 0.0) & print, ff[wh[0]] & print, gauss_cvf(ff[wh[0]])

stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Print the full inflation
;;;;;;;;;;;;;;;;;;;;;;
PRO print_inflation_table4

subsamp=100000
nskip = 1000

; ; Set the parameters in the table
pset=['ns','r']
pformat = ['F6.4', 'F5.3']
pprint_name = ['$n_s$       ','$r$           ']
pscale  = [1, 1] ; scale
nparams = n_elements(pset)


; Get the chains
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

chains  = ['c2_lcdm_pico_w7s12', 'c4_lcdm_pico_w7s12_H0','c3_lcdm_pico_w7s12_BAO','c5_lcdm_pico_w7s12_BAO_H0',$
           'c56_lcdm_r_camb_w7s12','c58_lcdm_r_camb_w7s12_H0','c57_lcdm_r_camb_w7s12_BAO','c59_lcdm_r_camb_w7s12_BAO_H0']
nchains = n_elements(chains)
dirs    = cdir+chains+'/chains/'
nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
ndsets = n_elements(chains)

; ouput matricies
param_med = dblarr(ndsets, nparams)
param_err = dblarr(ndsets, nparams)
param_uppers = dblarr(ndsets, nparams)

; Get parameter indicies in the chains
pindex = intarr(ndsets, nparams)
for i=0, ndsets-1 do begin
    readcol,pnames[i], nn,n2,format='a,a'
    for j=0, nparams-1 do pindex[i,j] = where(nn eq pset[j])
endfor

; ; get ML values
; ml = {ml0:get_ml_chain_all(files.f0), ml1:get_ml_chain_all(files.f1), ml2:get_ml_chain_all(files.f2), ml3:get_ml_chain_all(files.f3), ml4:get_ml_chain_all(files.f4), ml5:get_ml_chain_all(files.f5)}

; Get parameter values
for i=0, ndsets-1 do begin
    for j=0, nparams-1  do begin
        param = pset[j]
        print, pnames[i]

        ; skip 'r' for the LCDM set
        if (param eq 'r' and i lt 4) then continue
        plot_like1dname,files[i].ff,pnames[i],param,subsamp=subsamp,nskip=nskip,scale=pscale[j],results=results
        param_med[i,j] = results.median
        param_err[i,j] = results.avg_err
        param_uppers[i,j] = results.uppers[1]
    endfor
endfor

;-------------------------
; PRINT the table
;-------------------------
openw,lun,/get_lun,'tmp_inflation_tab4.txt'
j=0; 'ns'
printf, lun, format='("'+pprint_name[j]+'       & $",'+ $
  pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
  pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
  pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
  pformat[j]+'," \pm ",'+pformat[j]+',"$ \\")',$
  param_med[0,j], param_err[0,j], $
  param_med[1,j], param_err[1,j], $
  param_med[2,j], param_err[2,j], $
  param_med[3,j], param_err[3,j]

printf, lun, format='("'+pprint_name[j]+'       & $",'+ $
  pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
  pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
  pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
  pformat[j]+'," \pm ",'+pformat[j]+',"$ \\")',$
  param_med[4,j], param_err[4,j], $
  param_med[5,j], param_err[5,j], $
  param_med[6,j], param_err[6,j], $
  param_med[7,j], param_err[7,j]


j=1; 'r'
printf, lun, format='("'+pprint_name[j]+' (95\% CL) & $< ",'+ $
  pformat[j]+',"$ & $< ",'+$
  pformat[j]+',"$ & $< ",'+$
  pformat[j]+',"$ & $< ",'+$
  pformat[j]+',"$ \\")',$
  param_uppers[4,j],$
  param_uppers[5,j],$
  param_uppers[6,j],$
  param_uppers[7,j]

close, lun                         

stop
END




;;;;;;;;;;;;;;;;;;;;;;
; Power-law [LCDM]
;   Print constraints for all LCDM parameters using this procedure
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_inflation, param=param
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

;name = 'ns'
subsamp=100000
nskip = 2000
scale=1.

; cmb
dir_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c2_lcdm_pico_w7s12*.txt')
pname_cmb = dir_cmb+'c2_lcdm_pico_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = cdir+'c4_lcdm_pico_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c4_lcdm_pico_w7s12_H0*.txt')
pname_cmb_h0 = dir_cmb_h0+'c4_lcdm_pico_w7s12_H0.paramnames'
stop
; CMB+BAO
dir_cmb_bao = cdir+'c3_lcdm_pico_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c3_lcdm_pico_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c3_lcdm_pico_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c5_lcdm_pico_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c5_lcdm_pico_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c5_lcdm_pico_w7s12_BAO_H0.paramnames'

print, '' & print, '************ CMB *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
print, '' & print, '************ CMB+H0 *****************'
plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
print, '' & print, '************ CMB+H0+BAO *****************'
plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
print, 'CMB+H0: ', results_cmb_h0.median, results_cmb_h0.err_plus, results_cmb_h0.err_minus, results_cmb_h0.avg_err
print, 'CMB+BAO: ', results_cmb_bao.median, results_cmb_bao.err_plus, results_cmb_bao.err_minus, results_cmb_bao.avg_err
print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END


;;;;;;;;;;;;;;;;;;;;;;
; Tensor [LCDM+r]
;   param = 'ns' or 'r'
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_r, param=param
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; cmb
dir_cmb = cdir+'c56_lcdm_r_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c56_lcdm_r_camb_w7s12_*.txt')
pname_cmb = dir_cmb+'c56_lcdm_r_camb_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = cdir+'c58_lcdm_r_camb_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c58_lcdm_r_camb_w7s12_H0_*.txt')
pname_cmb_h0 = dir_cmb_h0+'c58_lcdm_r_camb_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = cdir+'c57_lcdm_r_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c57_lcdm_r_camb_w7s12_BAO_*.txt')
pname_cmb_bao = dir_cmb_bao+'c57_lcdm_r_camb_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c59_lcdm_r_camb_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c59_lcdm_r_camb_w7s12_BAO_H0_*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c59_lcdm_r_camb_w7s12_BAO_H0.paramnames'

;param = 'ns'
subsamp=100000
nskip = 1000
scale=1.

print, '' & print, '************ CMB *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
print, '' & print, '************ CMB+H0 *****************'
plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
print, '' & print, '************ CMB+H0+BAO *****************'
plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
print, 'CMB+H0: ', results_cmb_h0.median, results_cmb_h0.err_plus, results_cmb_h0.err_minus, results_cmb_h0.avg_err
print, 'CMB+BAO: ', results_cmb_bao.median, results_cmb_bao.err_plus, results_cmb_bao.err_minus, results_cmb_bao.avg_err
print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Running [LCDM+nrun]
;   param = 'nrun' or 'ns'
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_nrun, param=param
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; cmb
dir_cmb = cdir+'c7_lcdm_nrun_pico_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c7_lcdm_nrun_pico_w7s12*.txt')
pname_cmb = dir_cmb+'c7_lcdm_nrun_pico_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = cdir+'c9_lcdm_nrun_pico_w7s12_H0/post_chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c9_lcdm_nrun_pico_w7s12_H0*.txt')
pname_cmb_h0 = dir_cmb_h0+'c9_lcdm_nrun_pico_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = cdir+'c8_lcdm_nrun_pico_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c8_lcdm_nrun_pico_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c8_lcdm_nrun_pico_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c10_lcdm_nrun_pico_w7s12_BAO_H0/post_chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c10_lcdm_nrun_pico_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c10_lcdm_nrun_pico_w7s12_BAO_H0.paramnames'


subsamp=100000
nskip = 2000
scale=1.

print, '' & print, '************ CMB *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
print, '' & print, '************ CMB+H0 *****************'
plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
print, '' & print, '************ CMB+H0+BAO *****************'
plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
print, 'CMB+H0: ', results_cmb_h0.median, results_cmb_h0.err_plus, results_cmb_h0.err_minus, results_cmb_h0.avg_err
print, 'CMB+BAO: ', results_cmb_bao.median, results_cmb_bao.err_plus, results_cmb_bao.err_minus, results_cmb_bao.avg_err
print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Running + tensor [LCDM+nrun+r
;   param = 'ns', 'r', or 'nrun'
;;;;;;;;;;;;;;;;;;;;;;
PRO lcdm_nrun_r, param=param
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; cmb
dir_cmb = cdir+'c60_lcdm_nrun_r_camb_w7s12/post_chains/'
files_cmb = file_search(dir_cmb+'c60_lcdm_nrun_r_camb_w7s12_*.txt')
pname_cmb = dir_cmb+'c60_lcdm_nrun_r_camb_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = cdir+'c62_lcdm_nrun_r_camb_w7s12_H0/post_chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c62_lcdm_nrun_r_camb_w7s12_H0_*.txt')
pname_cmb_h0 = dir_cmb_h0+'c62_lcdm_nrun_r_camb_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = cdir+'c61_lcdm_nrun_r_camb_w7s12_BAO/post_chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c61_lcdm_nrun_r_camb_w7s12_BAO_*.txt')
pname_cmb_bao = dir_cmb_bao+'c61_lcdm_nrun_r_camb_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c63_lcdm_nrun_r_camb_w7s12_BAO_H0/post_chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c63_lcdm_nrun_r_camb_w7s12_BAO_H0_*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c63_lcdm_nrun_r_camb_w7s12_BAO_H0.paramnames'


subsamp=100000
nskip = 2000
scale=1.

print, '' & print, '************ CMB *****************'
plot_like1dname,files_cmb,pname_cmb,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb
print, '' & print, '************ CMB+H0 *****************'
plot_like1dname,files_cmb_h0,pname_cmb_h0,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
plot_like1dname,files_cmb_bao,pname_cmb_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_bao
print, '' & print, '************ CMB+H0+BAO *****************'
plot_like1dname,files_cmb_h0_bao,pname_cmb_h0_bao,param,subsamp=subsamp,nskip=nskip,scale=scale,results=results_cmb_h0_bao

print, 'CMB: ', results_cmb.median, results_cmb.err_plus, results_cmb.err_minus, results_cmb.avg_err
print, 'CMB+H0: ', results_cmb_h0.median, results_cmb_h0.err_plus, results_cmb_h0.err_minus, results_cmb_h0.avg_err
print, 'CMB+BAO: ', results_cmb_bao.median, results_cmb_bao.err_plus, results_cmb_bao.err_minus, results_cmb_bao.avg_err
print, 'CMB+H0+BAO: ', results_cmb_h0_bao.median, results_cmb_h0_bao.err_plus, results_cmb_h0_bao.err_minus, results_cmb_h0_bao.avg_err
stop
END








;;;;;;;;;;;;;;;;;;;;;;
; ML: LCDM
;;;;;;;;;;;;;;;;;;;;;;
PRO ml_lcdm
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; cmb
dir_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c2_lcdm_pico_w7s12*.txt')
pname_cmb = dir_cmb+'c2_lcdm_pico_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = cdir+'c4_lcdm_pico_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c4_lcdm_pico_w7s12_H0_*.txt')
pname_cmb_h0 = dir_cmb_h0+'c4_lcdm_pico_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = cdir+'c3_lcdm_pico_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c3_lcdm_pico_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c3_lcdm_pico_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c5_lcdm_pico_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c5_lcdm_pico_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c5_lcdm_pico_w7s12_BAO_H0.paramnames'

; print, '' & print, '************ CMB *****************'
; print_ml_chain_all, files_cmb
print, '' & print, '************ CMB+H0 *****************'
print_ml_chain_all, files_cmb_h0
; print, '' & print, '************ CMB+BAO *****************'
; print_ml_chain_all, files_cmb_bao
; print, '' & print, '************ CMB+H0+BAO *****************'
; print_ml_chain_all, files_cmb_h0_bao

END


;;;;;;;;;;;;;;;;;;;;;;
; ML: LCDM+r
;;;;;;;;;;;;;;;;;;;;;;
PRO ml_lcdm_r
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; cmb
dir_cmb = cdir+'c56_lcdm_r_camb_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c56_lcdm_r_camb_w7s12_*.txt')
pname_cmb = dir_cmb+'c56_lcdm_r_camb_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = cdir+'c58_lcdm_r_camb_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c58_lcdm_r_camb_w7s12_H0_*.txt')
pname_cmb_h0 = dir_cmb_h0+'c58_lcdm_r_camb_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = cdir+'c57_lcdm_r_camb_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c57_lcdm_r_camb_w7s12_BAO_*.txt')
pname_cmb_bao = dir_cmb_bao+'c57_lcdm_r_camb_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c59_lcdm_r_camb_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c59_lcdm_r_camb_w7s12_BAO_H0_*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c59_lcdm_r_camb_w7s12_BAO_H0.paramnames'

print, '' & print, '************ CMB *****************'
print_ml_chain_all, files_cmb
print, '' & print, '************ CMB+H0 *****************'
print_ml_chain_all, files_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
print_ml_chain_all, files_cmb_bao
print, '' & print, '************ CMB+H0+BAO *****************'
print_ml_chain_all, files_cmb_h0_bao

END


;;;;;;;;;;;;;;;;;;;;;;
; ML: LCDM+nrun
;;;;;;;;;;;;;;;;;;;;;;
PRO ml_lcdm_nrun
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

; cmb
dir_cmb = cdir+'c7_lcdm_nrun_pico_w7s12/chains/'
files_cmb = file_search(dir_cmb+'c7_lcdm_nrun_pico_w7s12*.txt')
pname_cmb = dir_cmb+'c7_lcdm_nrun_pico_w7s12.paramnames'

; CMB+H0
dir_cmb_h0 = cdir+'c9_lcdm_nrun_pico_w7s12_H0/chains/'
files_cmb_h0 = file_search(dir_cmb_h0+'c9_lcdm_nrun_pico_w7s12_H0*.txt')
pname_cmb_h0 = dir_cmb_h0+'c9_lcdm_nrun_pico_w7s12_H0.paramnames'

; CMB+BAO
dir_cmb_bao = cdir+'c8_lcdm_nrun_pico_w7s12_BAO/chains/'
files_cmb_bao = file_search(dir_cmb_bao+'c8_lcdm_nrun_pico_w7s12_BAO*.txt')
pname_cmb_bao = dir_cmb_bao+'c8_lcdm_nrun_pico_w7s12_BAO.paramnames'

; CMB+H0+BAO
dir_cmb_h0_bao = cdir+'c10_lcdm_nrun_pico_w7s12_BAO_H0/chains/'
files_cmb_h0_bao = file_search(dir_cmb_h0_bao+'c10_lcdm_nrun_pico_w7s12_BAO_H0*.txt')
pname_cmb_h0_bao = dir_cmb_h0_bao+'c10_lcdm_nrun_pico_w7s12_BAO_H0.paramnames'

print, '' & print, '************ CMB *****************'
print_ml_chain_all, files_cmb
;print, '' & print, '************ CMB+H0 *****************'
;print_ml_chain_all, files_cmb_h0
print, '' & print, '************ CMB+BAO *****************'
print_ml_chain_all, files_cmb_bao
;print, '' & print, '************ CMB+H0+BAO *****************'
;print_ml_chain_all, files_cmb_h0_bao
END



;;;;;;;;;;;;;;;;;
;
; Table 1, Field information
;
;;;;;;;;;;;;;;;;;
PRO fields_table
f = lps12_fieldstruct()
edir='/home/kstory/lps12/end2end/'
adir = '/home/kstory/lps12/masks/apod/'

dra    = fltarr(20)
ddec   = fltarr(20)
thresh = 0.5

total_area = 0.0
ita = 0

for i=0, 19 do begin
    if i eq 2 then i=3

    case i of
        0: ss='\zero' & 1: ss='\one' & 2: ss='\two' & 3: ss='\three' & 4: ss='\four' & 5: ss='\five' & 6: ss='\six' & 7: ss='\seven' & 8: ss='\eight' & 9: ss='\nine' & 10: ss='\ten' & $
          11: ss='\eleven' & 12: ss='\twelve' & 13: ss='\thirteen' & 14: ss='\fourteen' & 15: ss='\fifteen' & 16: ss='\sixteen' & 17: ss='\seventeen' & 18: ss='\eighteen' & 19: ss='\nineteen'
    endcase

    restore, '/home/kstory/lps12/masks/apod/apod_'+f[i].name+'_60_0.0500_30.sav' ; get apod
    mask = get_lps12_mask(i)    ; has ptsrc masks

    nx = n_elements(mask[*,0])
    ny = n_elements(mask[0,*])
    radec0 = [f[i].ra0, f[i].dec0]
    radec0[0] = 180.
    xc = round(0.5*nx)
    yc = round(0.5*ny)
    xind = findgen(nx)
    yind = findgen(ny)

; dRA
    tmp = apod[*,yc]
    wh=where(xind lt 0.8*xc)
    xmin = wh[where_closest(tmp[wh],thresh*max(tmp))]
    wh=where(xind gt 1.2*xc)
    xmax = wh[where_closest(tmp[wh],thresh*max(tmp))]

    pix2ang_proj5, [nx,ny], radec0, 1., ra, dec, $
                   xpix=[xmin,xmax], ypix=[yc,yc]

    dra[i] = abs(ra[1]-ra[0])
    dra[i] = 5* round(dra[i]/5.)

; dDEC
    tmp = apod[xc,*]
    wh=where(yind lt 0.8*yc)
    ymin = wh[where_closest(tmp[wh],thresh*max(tmp))]
    wh=where(yind gt 1.2*yc)
    ymax = wh[where_closest(tmp[wh],thresh*max(tmp))]

    pix2ang_proj5, [nx,ny], radec0, 1., ra, dec, $
                   xpix=[xc,xc], ypix=[ymin,ymax]

    ddec[i] = abs(dec[1]-dec[0])
    ddec[i] = 5* round(ddec[i]/5.)

; area
    area = total(mask) / (60^2.)
     ;print,f[i].name
    pformat = ['F5.1', 'F5.1', 'I', 'I', 'I3']

    print, format='("'+ss+' & ",'+ $
      pformat[0]+'," &   ",'+$
      pformat[1]+'," &   ",'+$
      pformat[2]+'," &   ",'+$
      pformat[3]+'," &   ",'+$
      pformat[4]+'," \\")',$
      f[i].ra0, f[i].dec0, dra[i], ddec[i], round(area)

;     print, ss+'    & ', sigfig(f[i].ra0,3),' & ',sigfig(f[i].dec0,1),' & ',$
;       sigfig(dra[i],3),' & ',sigfig(ddec[i],3), ' & ', round(area), ' \\'
    
    total_area += area
    ita += round(area)
endfor
print, "total area = ", total_area
print, "integer total area = ", ita

stop
END
