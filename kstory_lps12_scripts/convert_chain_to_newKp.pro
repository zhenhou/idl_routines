;;;
; NAME: convert_chain_to_newKp.pro
; PURPOSE:
;  Write new lcdm chains with different pivot point
;
; NOTES:
;  If the old chain has "logA" instead of 1e-9As, convert that as well
;
; MODIFICATION HISTORY:
;  04/27/2012: (KTS) Created from script13_1011.pro
;  05/15/2012: (KTS) Modify to work with Planck chains
;;;


;;;;;;;;;;;;;;;;;;;
;
; 1) change pivot point in chains from 0.002 to 0.05
;
;;;;;;;;;;;;;;;;;;;
PRO convert_chain_to_newKp, mkdir=mkdir, logA=logA, doPlanck=doPlanck

;;;SPT chains
k0_old = 0.002
k0_new = 0.05

; Set the parameters in the table
pset=['1e-9As','ns'] ; names in chain.paramname
if keyword_set(logA) then pset[0] = 'logA'
nparams = n_elements(pset)


; Get the chains
cdir = '/data23/hou/lps12/paramfits/chains_0828/'
;chains  = ['c1_lcdm_pico_w7', 'c27_lcdm_camb_s12tau', 'c2_lcdm_pico_w7s12', 'c4_lcdm_pico_w7s12_H0','c3_lcdm_pico_w7s12_BAO','c5_lcdm_pico_w7s12_BAO_H0']
chains  = ['c2_lcdm_pico_w7s12_lmax1500', 'c2_lcdm_pico_w7s12_1500ell3000']

;;;Planck chains
if keyword_set(doPlanck) then begin
    print, "*** Running on planck chains"
    k0_old = 0.05
    k0_new = 0.05
    cdir = '/data/kstory/projects/planck_chains/planck_lowl_lowLike/' 
    chains = ['base_planck_lowl_lowLike']
endif


nchains = n_elements(chains)
dirs    = cdir+chains+'/chains/'
nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
;ndsets = n_elements(chains)

if keyword_set(doPlanck) then begin
    nfiles = n_elements(file_search(cdir+chains[0]+'_*.txt'))
    st = {ff:strarr(nfiles)}
    files = replicate(st, nchains)
    for i=0, nchains -1 do files[i].ff = file_search(cdir+chains[i]+'_*.txt')
    pnames  = cdir + chains+'.paramnames'
endif

; New output chains
base_dir = '/data23/kstory/lps12/chains/'
ochains = chains+'_newKp'
odirs = base_dir+ochains+'/chains/'
ofiles = replicate(st, nchains)
opnames  = odirs + ochains+'.paramnames'

if keyword_set(doPlanck) then begin
    odirs = '/data/kstory/projects/planck_chains/planck_lowl_lowLike_newKp/' 
    ochains = chains+'_newKp'
    ofiles = replicate(st, nchains)
    opnames  = odirs + ochains+'.paramnames'
endif

; setup output filenames 
for i=0, nchains-1 do begin
    if keyword_set(mkdir) then begin 
        spawn, 'mkdir '+base_dir+ochains[i]
        spawn, 'mkdir '+odirs[i]
    endif

    for j=1, nfiles do begin
        ssj = strtrim(string(j),2)
        ofiles[i].ff[j-1] = odirs[i]+ochains[i]+'_'+ssj+'.txt'
    endfor
endfor

; Get parameter indicies in the chains
pindex = intarr(nchains, nparams)
for i=0, nchains-1 do begin
    readcol,pnames[i], nn,n2,format='a,a'
    for j=0, nparams-1 do pindex[i,j] = where(nn eq pset[j])
endfor
pindex += 2 ; skip Likelihood and ML column

;--------------------------
; Calculate new chains
for i=0, nchains-1 do begin
    print, ''
    print, 'process chain: '+chains[i]
    for j=0, nfiles-1 do begin ;for j=0, 0 do begin

        print, format='($,I,:,",")', j
        b=read_ascii(files[i].ff[j])

        ;Read out As and ns
        As = double(reform(b.field01[pindex[i,0],*]))
        ns = double(reform(b.field01[pindex[i,1],*]))
        
        ;Convert to 1e-9As
        if keyword_set(logA) then As = 0.1*exp(As)

        As_new = double(As)*0
        for k=0L, n_elements(As)-1 do As_new[k] = As[k] * (k0_new / k0_old) ^ (ns[k]-1)

        b.field01[pindex[i,0],*] = float(As_new)

        ; Write new file
        ncol = n_elements(b.field01[*,0])
        openw,lun1,/get_lun,ofiles[i].ff[j]
        for k=0L, n_elements(As)-1 do begin
            printf, lun1, format='(52(E15, :, " "))', b.field01[*,k]
        endfor
        close, lun1

    endfor

    spawn, 'cp '+pnames[i]+' '+opnames[i]

    ; I'm being lazy...
    if keyword_set(logA) then print, "*********** YOU NEED TO MANUALLY CHANGE THE .pnames FILE TO USE 1e-9As **************"

endfor

stop
END






;;;;;;;;;;;;;;;;;;;
;
; 2) Check that the new chains are right
;
;;;;;;;;;;;;;;;;;;;
PRO check_chain, mkdir=mkdir
subsamp=100000
nskip = 1000

; Set the parameters in the table
pset=['omegabh2','omegadmh2','1e-9As','ns','theta_s','tau',  'omegal*','H0*','sigma8*','z_eq','Dvp35ors','Dvp57ors'] ; names in chain.paramname
pformat = ['F5.3', 'F6.4', 'F5.3', 'F6.4', 'F6.4', 'F5.3',   'F5.3', 'F4.1', 'F5.3',   'F5.0','F6.4',    'F7.4'] ; print format
pprint_name = ['$10^2\Omega_bh^2$ ','$\Omega_ch^2$    ','$10^{9} \deltaR$  ','$n_s$             ','$10^2\theta_s$  ','$\tau$   ', $
               '$\Omega_\Lambda$   ','\ho     ','$\sigma_8$      ','$z_{\rm EQ}$     ','$r_s/D_V(z=0.35)$', '$r_s/D_V(z=0.57)$']
dv_fac = (150.82/154.66)/100.
pscale  = [100., 1, 1, 1, 100., 1, 1, 1, 1, 1, dv_fac, dv_fac] ; scale
ppower  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1., -1.] ; scale
nparams = n_elements(pset)


; Get the chains
cdir = ['/data23/hou/lps12/paramfits/chains_0828/', '/data23/kstory/lps12/chains/']
;chains  = ['c1_lcdm_pico_w7', 'c27_lcdm_camb_s12tau', 'c2_lcdm_pico_w7s12', 'c4_lcdm_pico_w7s12_H0','c3_lcdm_pico_w7s12_BAO','c5_lcdm_pico_w7s12_BAO_H0']
;chains  = ['c1_lcdm_pico_w7','c1_lcdm_pico_w7_newKp']
chains  = ['c27_lcdm_camb_s12tau', 'c27_lcdm_camb_s12tau_newKp']
nchains = n_elements(chains)
dirs    = cdir[0]+chains+'/chains/'
dirs[1] = cdir[1]+chains[1]+'/chains/'
nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
ndsets = n_elements(chains)

; ouput matricies
param_med = dblarr(ndsets, nparams)
param_err = dblarr(ndsets, nparams)

; Get parameter indicies in the chains
pindex = intarr(ndsets, nparams)
for i=0, ndsets-1 do begin
    readcol,pnames[i], nn,n2,format='a,a'
    for j=0, nparams-1 do pindex[i,j] = where(nn eq pset[j])
endfor

; get ML values
;ml = {ml0:get_ml_chain_all(files.f0), ml1:get_ml_chain_all(files.f1), ml2:get_ml_chain_all(files.f2), ml3:get_ml_chain_all(files.f3), ml4:get_ml_chain_all(files.f4), ml5:get_ml_chain_all(files.f5)}

subset=[2,3] ; indicies of pset

; Get parameter values
for i=0, ndsets-1 do begin
    for k=0, n_elements(subset)-1 do begin
    ;for j=0, nparams-1  do begin
        j = subset[k]
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

for k=0, n_elements(subset)-1 do $ 
  print, format='("'+pprint_name[subset[k]]+'       & $",'+ $
  pformat[subset[k]]+'," \pm ",'+pformat[subset[k]]+',"$ &    $",'+$
  pformat[subset[k]]+'," \pm ",'+pformat[subset[k]]+',"$ \\")',$
  param_med[0,subset[k]], param_err[0,subset[k]], $
  param_med[1,subset[k]], param_err[1,subset[k]]

openw,lun,/get_lun,'tmp.txt'
for k=0, n_elements(subset)-1 do $ 
  printf, lun, format='("'+pprint_name[subset[k]]+'       & $",'+ $
  pformat[subset[k]]+'," \pm ",'+pformat[subset[k]]+',"$ &    $",'+$
  pformat[subset[k]]+'," \pm ",'+pformat[subset[k]]+',"$ \\")',$
  param_med[0,subset[k]], param_err[0,subset[k]], $
  param_med[1,subset[k]], param_err[1,subset[k]]
close, lun

; ; No ML valuse
; openw,lun,/get_lun,'tmp.txt'
; for j=0, nparams-1 do $
;     printf, lun, format='("'+pprint_name[j]+'       & $",'+ $
;       pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
;       pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
;       pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
;       pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
;       pformat[j]+'," \pm ",'+pformat[j]+',"$ &    $",'+$
;       pformat[j]+'," \pm ",'+pformat[j]+',"$ \\")',$
;       param_med[0,j], param_err[0,j], $
;       param_med[1,j], param_err[1,j], $
;       param_med[2,j], param_err[2,j], $
;       param_med[3,j], param_err[3,j], $
;       param_med[4,j], param_err[4,j], $
;       param_med[5,j], param_err[5,j]
; close, lun                         

stop



END
