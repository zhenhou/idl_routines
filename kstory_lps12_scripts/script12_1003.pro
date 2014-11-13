;;;
; NAME: script12_0927
;
; NOTES:
;  1) test constraints
; 2) make mode coupling kernels
;;;


;;;;;;;;;;;;;;;;;;;
;
; Make mode coupling kernels
;
;;;;;;;;;;;;;;;;;;;
PRO make_kerns, interp=interp, oversamp=oversamp, stopit=stopit

if n_elements(interp) eq 0 then interp=1000
if n_elements(oversamp) eq 0 then oversamp=8

idx = 4
f = lps12_fieldstruct()
mask = get_lps12_mask(idx, /padded)
reso_arcmin = 1.0

; get BANDDEF
delta_l = 50.
min_l = 250.
max_l = 3150.
nl = floor((max_l-min_l)/delta_l)
lhi = dindgen(nl)*delta_l + min_l
banddef = lhi
; this can possibly be much lower:
maxell = max(banddef)*1.5

; make kern
kernel=coupling_kernel(mask, reso_arcmin, maxell=maxell, /changevar, $
                       interp=1000, oversamp=8, /cheby, ellkern=ellkern, $
                       curlyw=curlyw)
reso=double(reso_arcmin)/60*!dtor
kernsize=(size(kernel))[1]
u=(dindgen(kernsize)+0.5)#replicate(1./(reso*winsize)^4, kernsize)
kernel*=u    

; Save output
dir = '/data/kstory/projects/lps12/masks/test_kerns_1004/'
sinterp = strtrim(string(interp),2)
soversamp = strtrim(string(oversamp),2)
save, ellkern, kernel, filename=dir+'kern_'+f[idx].name+'_interp'+sinterp+'_oversamp'+soversamp+'.sav'

if keyword_set(stopit) then stop
END




PRO a

; ------------------------
; CHAINS
; ------------------------
; Get the new chains
cdir = '/data23/hou/lps12/paramfits/chains_0828_pipetest/'

spawn, 'ls '+cdir+' | grep _lcdm', chains
nchains = n_elements(chains)
dirs = cdir+'/'+chains+'/chains/'

nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
;ndsets = n_elements(chains)

; Specify which chains to run
ch_beamerr = [3,5,9,10,13]
ch_nofg = [0,4,6,11,12]

; Get the baseline chains
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

bchains  = ['c2_lcdm_pico_w7s12', $
            'c54_lcdm_alens_camb_w7s12', $
            'c12_lcdm_neff_pico_w7s12', $
            'c7_lcdm_nrun_pico_w7s12' ]
           

bnchains = n_elements(chains)
bdirs = cdir+'/'+bchains+'/chains/'

nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
bfiles = replicate(st, nchains)
for i=0, bnchains -1 do bfiles[i].ff = file_search(bdirs[i]+bchains[i]+'*.txt')
bpnames  = dirs + chains+'.paramnames'

; Specify which chains to run
;ch_beamerr = [3,5,9,10,13]
;ch_nofg = [0,4,6,11,12]



; ------------------------
; PARAMETERS
; ------------------------
pset=['omegabh2','omegadmh2','1e-9As','ns','theta_s','tau',  'omegal*','H0*','sigma8*','z_eq','Dvp35ors','Dvp57ors'] ; names in chain.paramname
dv_fac = (150.82/154.66)/100.
pscale  = [100., 1, 1, 1, 100., 1, 1, 1, 1, 1, dv_fac, dv_fac] ; scale
ppower  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1., -1.] ; scale
nparams = n_elements(pset)

; ouput matricies
param_med = dblarr(ndsets, nparams)
param_err = dblarr(ndsets, nparams)

; Get parameter indicies in the chains
pindex = intarr(ndsets, nparams)
for i=0, ndsets-1 do begin
    readcol,pnames[i], nn,n2,format='a,a'
    for j=0, nparams-1 do pindex[i,j] = where(nn eq pset[j])
endfor


; Specify which chain and params to look at
;ch_set = ch_beamerr ; specify which chains to use
ch_set = [3]
nch = n_elements(ch_set)
puse_idx=[3] ; parameters to use
nparam = n_elements(puse_idx)

subsamp=100000
nskip = 1000

for chidx=0, nch-1 do begin
    i = ch_set[chidx]
    for pidx=0, nparam-1 do begin
        j = puse_idx[pidx]
        print, chains[i], pset[j]
        plot_like1dname,files[i].ff,pnames[i],pset[j],subsamp=subsamp,nskip=nskip,scale=pscale[j],power=ppower[pidx],results=results

        param_med[i,j] = results.median
        param_err[i,j] = results.avg_err
    endfor
endfor



stop
END







;;;;;;;;;;;;;;
;
; Compare two chains
;
;;;;;;;;;;;;;;
PRO comp

; ------------------------
; CHAINS
; ------------------------
; beamerr
;chains = ['c2_lcdm_pico_w7s12','ptc1_lcdm_pico_w7s12_beamerr'] ; lcdm
chains = ['c27_lcdm_camb_s12tau','ptc9_lcdm_camb_s12tau_beamerr'] ;lcdm_s12
; chains = ['c54_lcdm_alens_camb_w7s12', 'ptc13_lcdm_alens_camb_w7s12_beamerr'] ; lcdm_alens
; chains = ['c12_lcdm_neff_pico_w7s12', 'ptc5_lcdm_neff_pico_w7s12_beamerr'] ; lcdm_neff
; chains = ['c7_lcdm_nrun_pico_w7s12', 'ptc6_lcdm_nrun_pico_w7s12_beamerr'] ; lcdm_nrun

; ; nofgprior
; chains = ['c2_lcdm_pico_w7s12','ptc2_lcdm_pico_w7s12_nofgprior'] ; lcdm
; chains = ['c27_lcdm_camb_s12tau','ptc10_lcdm_camb_s12tau_nofgprior'] ;lcdm_s12
; chains = ['c54_lcdm_alens_camb_w7s12', 'ptc14_lcdm_alens_camb_w7s12_nofgprior'] ; lcdm_alens
; chains = ['c12_lcdm_neff_pico_w7s12', 'ptc7_lcdm_neff_pico_w7s12_nofgprior'] ; lcdm_neff
; chains = ['c7_lcdm_nrun_pico_w7s12', 'ptc8_lcdm_nrun_pico_w7s12_nofgprior'] ; lcdm_nrun

cdir = ['/data23/hou/lps12/paramfits/chains_0828/','/data23/hou/lps12/paramfits/chains_0828_pipetest/']

nchains = 2
dirs = [cdir[0]+'/'+chains[0]+'/chains/', cdir[1]+'/'+chains[1]+'/chains/']

nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
ndsets = 2

; ------------------------
; PARAMETERS
; ------------------------
pset=['omegabh2','omegadmh2','1e-9As','ns','theta_s','tau',  'omegal*','H0*','sigma8*','z_eq','Dvp35ors','Dvp57ors'] ; names in chain.paramname
dv_fac = (150.82/154.66)/100.
pscale  = [100., 1, 1, 1, 100., 1, 1, 1, 1, 1, dv_fac, dv_fac] ; scale
ppower  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1., -1.] ; scale
nparams = n_elements(pset)

; Get parameter indicies in the chains
pindex = intarr(ndsets, nparams)
for i=0, ndsets-1 do begin
    readcol,pnames[i], nn,n2,format='a,a'
    for j=0, nparams-1 do pindex[i,j] = where(nn eq pset[j])
endfor


; Specify which chain and params to look at
puse_idx=[3] ; parameters to use
nparam = n_elements(puse_idx)

; ouput matricies
param_med = dblarr(ndsets, nparam)
param_err = dblarr(ndsets, nparam)

subsamp=100000
nskip = 1000

for i=0, nchains-1 do begin
    for pidx=0, nparam-1 do begin
        j = puse_idx[pidx]
        print, chains[i], ', ', pset[j]
        plot_like1dname,files[i].ff,pnames[i],pset[j],subsamp=subsamp,nskip=nskip,scale=pscale[j],power=ppower[pidx],results=results

        param_med[i,pidx] = results.median
        param_err[i,pidx] = results.avg_err
    endfor
endfor

print, '---- ----'
print, chains, ', ', pset[puse_idx]
print, param_med
print, param_err

stop
END






;;;;;;;;;;;;;;
;
; Plot two chains
;
;;;;;;;;;;;;;;
PRO plot_ch, ii, stopit=stopit

; ------------------------
; CHAINS
; ------------------------
; beamerr
case ii of
    0: chains = ['c2_lcdm_pico_w7s12','ptc1_lcdm_pico_w7s12_beamerr'] ; lcdm
    1: chains = ['c27_lcdm_camb_s12tau','ptc9_lcdm_camb_s12tau_beamerr'] ;lcdm_s12
    2: chains = ['c54_lcdm_alens_camb_w7s12', 'ptc13_lcdm_alens_camb_w7s12_beamerr'] ; lcdm_alens
    3: chains = ['c12_lcdm_neff_pico_w7s12', 'ptc5_lcdm_neff_pico_w7s12_beamerr'] ; lcdm_neff
    4: chains = ['c7_lcdm_nrun_pico_w7s12', 'ptc6_lcdm_nrun_pico_w7s12_beamerr'] ; lcdm_nrun

; ; nofgprior
    5: chains = ['c2_lcdm_pico_w7s12','ptc2_lcdm_pico_w7s12_nofgprior'] ; lcdm
    6: chains = ['c27_lcdm_camb_s12tau','ptc10_lcdm_camb_s12tau_nofgprior'] ;lcdm_s12
    7: chains = ['c54_lcdm_alens_camb_w7s12', 'ptc14_lcdm_alens_camb_w7s12_nofgprior'] ; lcdm_alens
    8: chains = ['c12_lcdm_neff_pico_w7s12', 'ptc7_lcdm_neff_pico_w7s12_nofgprior'] ; lcdm_neff
    9: chains = ['c7_lcdm_nrun_pico_w7s12', 'ptc8_lcdm_nrun_pico_w7s12_nofgprior'] ; lcdm_nrun
endcase

cdir = ['/data23/hou/lps12/paramfits/chains_0828/','/data23/hou/lps12/paramfits/chains_0828_pipetest/']

nchains = 2
dirs = [cdir[0]+'/'+chains[0]+'/chains/', cdir[1]+'/'+chains[1]+'/chains/']

nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
ndsets = 2

; ------------------------
; PARAMETERS
; ------------------------
pset=['omegabh2','omegadmh2','1e-9As','ns','theta_s','tau',  'omegal*','H0*','sigma8*','z_eq','Dvp35ors','Dvp57ors'] ; names in chain.paramname
dv_fac = (150.82/154.66)/100.
pscale  = [100., 1, 1, 1, 100., 1, 1, 1, 1, 1, dv_fac, dv_fac] ; scale
ppower  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1., -1.] ; scale
nparams = n_elements(pset)

; Get parameter indicies in the chains
pindex = intarr(ndsets, nparams)
for i=0, ndsets-1 do begin
    readcol,pnames[i], nn,n2,format='a,a'
    for j=0, nparams-1 do pindex[i,j] = where(nn eq pset[j])
endfor


; Specify which chain and params to look at
puse_idx=[4] ; parameters to use
nparam = n_elements(puse_idx)

; ouput matricies
param_med = dblarr(ndsets, nparam)
param_err = dblarr(ndsets, nparam)

subsamp=5
nskip = 1000

for i=0, nchains-1 do begin
    for pidx=0, nparam-1 do begin
        j = puse_idx[pidx]
        print, chains[i], ', ', pset[j]
        if (i+pidx eq 0) then begin
            plot_like1dname,files[i].ff,pnames[i],pset[j],title=chains[0],subsamp=subsamp,nskip=nskip,scale=pscale[j],power=ppower[pidx],results=results
        endif else $
          plot_like1dname,files[i].ff,pnames[i],pset[j],/oplot,subsamp=subsamp,nskip=nskip,scale=pscale[j],power=ppower[pidx],color=!red,results=results
        
        param_med[i,pidx] = results.median
        param_err[i,pidx] = results.avg_err
    endfor
endfor

print, '---- ----'
print, chains, ', ', pset[puse_idx]
print, param_med
print, param_err

if keyword_set(stopit) then stop
END







;;;;;;;;;;;;;;
;
; Plot two chains
;
;;;;;;;;;;;;;;
PRO plot_sne, ii, pidx=pidx, stopit=stopit

; ------------------------
; CHAINS
; ------------------------

cdir = ['/data23/hou/lps12/paramfits/chains_0828/','/data23/hou/lps12/paramfits/chains_0828/']
case ii of
    0: chains = ['c2_lcdm_pico_w7s12','c95_lcdm_pico_w7s12_SNe'] ; lcdm
    1: chains = ['c27_lcdm_camb_s12tau','c94_lcdm_camb_s12tau_SNe'] ;lcdm_s12
    2: begin
        chains = ['c5_lcdm_pico_w7s12_BAO_H0', 'c205_lcdm_pico_w7s12_BAO_H0_SNe'] ; lcdm_cmb_h0_bao
        cdir = ['/data23/hou/lps12/paramfits/chains_0828/','/data23/kstory/lps12/chains/']
    end
endcase


nchains = 2
dirs = [cdir[0]+'/'+chains[0]+'/chains/', $
        cdir[1]+'/'+chains[1]+'/chains/']

nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))
st = {ff:strarr(nfiles)}
files = replicate(st, nchains)
for i=0, nchains -1 do files[i].ff = file_search(dirs[i]+chains[i]+'*.txt')
pnames  = dirs + chains+'.paramnames'
ndsets = 2

; ------------------------
; PARAMETERS
; ------------------------
pset=['omegabh2','omegadmh2','1e-9As','ns','theta_s','tau',  'omegal*','H0*','sigma8*','z_eq','Dvp35ors','Dvp57ors'] ; names in chain.paramname
dv_fac = (150.82/154.66)/100.
pscale  = [100., 1, 1, 1, 100., 1, 1, 1, 1, 1, dv_fac, dv_fac] ; scale
ppower  = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1., -1.] ; scale
nparams = n_elements(pset)

; Get parameter indicies in the chains
pindex = intarr(ndsets, nparams)
for i=0, ndsets-1 do begin
    readcol,pnames[i], nn,n2,format='a,a'
    for j=0, nparams-1 do pindex[i,j] = where(nn eq pset[j])
endfor


; Specify which chain and params to look at
if n_elements(pidx) eq 0 then pidx = 0
puse_idx=[pidx] ; parameters to use
nparam = n_elements(puse_idx)

; ouput matricies
param_med = dblarr(ndsets, nparam)
param_err = dblarr(ndsets, nparam)

subsamp=5
nskip = 1000

for i=0, nchains-1 do begin
    for pidx=0, nparam-1 do begin
        j = puse_idx[pidx]
        print, chains[i], ', ', pset[j]
        if (i+pidx eq 0) then begin
            plot_like1dname,files[i].ff,pnames[i],pset[j],title=chains[1],subsamp=subsamp,nskip=nskip,scale=pscale[j],power=ppower[pidx],results=results
        endif else $
          plot_like1dname,files[i].ff,pnames[i],pset[j],/oplot,subsamp=subsamp,nskip=nskip,scale=pscale[j],power=ppower[pidx],color=!red,results=results
        
        param_med[i,pidx] = results.median
        param_err[i,pidx] = results.avg_err
    endfor
endfor

print, '---- ----'
print, chains, ', ', pset[puse_idx]
print, param_med
print, param_err

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
if keyword_set(stopit) then stop
END

