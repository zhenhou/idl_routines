;;;
; NAME: script13_0512
;
; NOTES:
;  1) Check how much params move with no foreground priors
;;;

PRO check_nofg,pnum
chains = ['c2_lcdm_pico_w7s12','ptc2_lcdm_pico_w7s12_nofgprior'] ; lcdm
cdir = ['/data23/hou/lps12/paramfits/chains_0828/','/data23/hou/lps12/paramfits/chains_0828_pipetest/']
dirs = [cdir[0]+'/'+chains[0]+'/chains/', cdir[1]+'/'+chains[1]+'/chains/']
nfiles = n_elements(file_search(dirs[0]+chains[0]+'*.txt'))

spawn, 'ls '+dirs[0]+chains[0]+'_*.txt', files0
params0 = dirs[0]+chains[0]+'.paramnames'

spawn, 'ls '+dirs[1]+chains[1]+'_*.txt', files1
params1 = dirs[1]+chains[1]+'.paramnames'

pset=['omegabh2','omegadmh2','1e-9As','ns','theta_s','tau',  'omegal*','H0*','sigma8*','z_eq'] ; names in chain.paramname
pscale  = [100., 1, 1, 1, 100., 1, 1, 1, 1, 1] ; scale

stop
plot_like1dname, files0,params0,pset[pnum],subsamp=1000,nskip=1000,scale=pscale[pnum], results=results0
plot_like1dname, files1,params1,pset[pnum],subsamp=1000,nskip=1000,scale=pscale[pnum], results=results1

x = (results1.median - results0.median) / (sqrt(results1.avg_err^2. + results0.avg_err^2.))
print, "*** "+pset[pnum]
print, results0.median, results0.avg_err
print, results1.median, results1.avg_err
print, "change: ", x

stop
END
