;;;
; NAME: script12_1214
;
; NOTES:
;  1) Check nrun constraints from beamerr, nofgprior
;;;


;;;;;;;;;;;;;;
;
; Plot two chains
;
;;;;;;;;;;;;;;
PRO check_nrun, ii, stopit=stopit

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

subsamp=1000
nskip = 1000

i=0
print, chains[i]
plot_like1dname,files[i].ff,pnames[i],'nrun',title=chains[i],subsamp=subsamp,nskip=nskip,results=results

i=1
print, chains[i]
plot_like1dname,files[i].ff,pnames[i],'nrun',title=chains[i],subsamp=subsamp,nskip=nskip,results=results

stop
print, '---- ----'
print, chains, ', ', pset[puse_idx]
print, param_med
print, param_err

if keyword_set(stopit) then stop
END

