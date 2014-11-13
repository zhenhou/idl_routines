;;;
; NAME: script12_0927
;
; NOTES:
;  1) test constraints
;;;


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
PRO plot_sne, ii, stopit=stopit

; ------------------------
; CHAINS
; ------------------------

case ii of
    0: chains = ['c2_lcdm_pico_w7s12','c95_lcdm_pico_w7s12_SNe'] ; lcdm
    1: chains = ['c27_lcdm_camb_s12tau','c94_lcdm_camb_s12tau_SNe'] ;lcdm_s12
endcase

cdir = ['/data23/hou/lps12/paramfits/chains_0828/','/data23/hou/lps12/paramfits/chains_0828/']

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
puse_idx=[6] ; parameters to use
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







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 8-panel LCDM plot
PRO plot_lcdm_8panel
cdir = ['/data23/hou/lps12/paramfits/chains_0828/','/data23/hou/lps12/paramfits/chains_0828_pipetest/']

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data
ii = 9
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


; baseline
i = 0
bdir   = cdir[i]+'/'+chains[i]+'/chains/'
bfiles = file_search(bdir+chains[i]+'*.txt')
bpname = bdir + chains[i]+'.paramnames'

; extension
i = 1
dir   = cdir[i]+'/'+chains[i]+'/chains/'
files = file_search(dir+chains[i]+'*.txt')
pname = dir + chains[i]+'.paramnames'

stop

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 8
csize = 0.8
yb = xb * 0.5

filename_stub = '~kstory/lps12/scripts/figs/'+chains[1]+'_0927'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.095
ddx=.22
y0=.12
ddy=.35

; Define the output plot and font sizes
axlCthick = 3
axlCsize  = 0.8
ytxt = [-0.23, .45] ; y position of [x,y]-axis label
;bigspt = 1.4
;thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_w7s12 = pscolors[0] ; black
color_s12   = pscolors[8] ; blue
color_w7    = pscolors[2] ; !red

ls      = [3,2,0] ; [SPT-only, WMAP-only, SPT+WMAP]
mycolor = [color_s12, color_w7, color_w7s12] ; [SPT-only, WMAP-only, SPT+WMAP]


!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

; axis properties
!y.crange=[0,1]
ytitle='!3Likelihood!X!N'

; TMP
!x.crange=[0,9.5]
xtmp=indgen(10)
ytmp=0.1*xtmp

;;; Top Row

for i=0,3 do begin
    !x.window=[x0+i*ddx,x0+(i+1)*ddx]
    !y.window=[(y0+ddy)+y0,(y0+ddy)+y0+ddy]

    case i of
        0: begin ; omegahb plot
            !x.crange=[1.8, 2.92]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xyouts, 1.61, ytxt[1], alignment=0.5, orientation=90, charsize=0.9, font=-1, charthick=axlCthick, ytitle
            xtitle='!4X!3!Db!Nh!U2!X!N'
            xyouts, 2.3, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,bfiles,bpname,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files,pname,'omegabh2',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
;             plot_like1dname,files_w7s12,pname_w7s12,'omegabh2',subsamp=3,nskip=1000,scale=100,$
;               thick=4,linestyle=ls[2], color=mycolor[2], /oplot

            oplot, [2.6,2.9], [0.9+0.02,0.9+0.02], thick=4,linestyle=ls[2],color=mycolor[2]
            oplot, [2.6,2.9], [0.7+0.02,0.7+0.02], thick=4,linestyle=ls[1],color=mycolor[1]
;             oplot, [2.6,2.9], [0.6+0.02,0.6+0.02], thick=4,linestyle=ls[0],color=mycolor[0]
            
        endcase
        1: begin ; omegadmh2
            !x.crange=[0.075, 0.139]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.02
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!4X!3!Dc!Nh!U2!X!N'
            xyouts, 0.105, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,bfiles,bpname,'omegadmh2',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files,pname,'omegadmh2',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
;             plot_like1dname,files_w7s12,pname_w7s12,'omegadmh2',subsamp=3,nskip=1000,$
;               thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
            xyouts, 0.076, 0.9, 'baseline', charsize=0.7, charthick=2, font=0, color=mycolor[2]
            xyouts, 0.076, 0.8, 'extension', charsize=0.7, charthick=2, font=0, color=mycolor[1]

        endcase
        2: begin ; theta_s
            !x.crange=[1.028, 1.053]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.01
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!3 100 !4h!D!3s!X!N'
            xyouts, 1.039, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,bfiles,bpname,'theta_s',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files,pname,'theta_s',subsamp=3,nskip=1000,scale=100,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
;             plot_like1dname,files_w7s12,pname_w7s12,'theta_s',subsamp=3,nskip=1000,scale=100,$
;               thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        3: begin ; omegal
            !x.crange=[0.56, 0.89]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=0.1 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!U*!N!4X!DK!X!N'
            xyouts, 0.73, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,bfiles,bpname,'omegal*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files,pname,'omegal*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
;             plot_like1dname,files_w7s12,pname_w7s12,'omegal*',subsamp=3,nskip=1000,$
;               thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        else: begin
            print, 'i = ', i
        endcase
    endcase
endfor

;;; Bottom Row
for i=0,3 do begin
    !x.window=[x0+i*ddx,x0+(i+1)*ddx]
    !y.window=[y0,y0+ddy]

    case i of
        0: begin ; tau
            !x.crange=[0.041, 0.149]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xyouts, 0.022, ytxt[1], alignment=0.5, orientation=90, charsize=0.9, font=-1, charthick=axlCthick, ytitle
            xtitle='!4s!X!N'
            xyouts, 0.095, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            ;;   Do not plot s12, since we use the WMAP prior on tau
            plot_like1dname,bfiles,bpname,'tau',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files,pname,'tau',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
;             plot_like1dname,files_w7s12,pname_w7s12,'tau',subsamp=3,nskip=1000,$
;               thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        1: begin ; ns
            !x.crange=[0.81, 1.05]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!3n!Ds!X!N'
            xyouts, 0.93, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,bfiles,bpname,'ns',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files,pname,'ns',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
;             plot_like1dname,files_w7s12,pname_w7s12,'ns',subsamp=3,nskip=1000,$
;               thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        2: begin ; As
            ;!x.crange=[2.81, 3.9]
            !x.crange=[1.5, 4.2]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            ;xtitle='!3log[10!U10!N!4D!3!U2!DR!N]!X'
            xtitle='!310!U9!N!4D!3!U2!DR!N!X'
            xyouts, 2.5, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            scale=1.

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,bfiles,bpname,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files,pname,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
;             plot_like1dname,files_w7s12,pname_w7s12,'1e-9As',subsamp=3,nskip=1000,scale=scale,$
;               thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        3: begin ; H0
            !x.crange=[61, 89.9]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            xtitle='!U*!N!3H!D0!X!N'
            xyouts, 75, ytxt[0], charsize=axlCsize, charthick=axlCthick, font=-1, xtitle

            oplot,[0,0],[0,0],color=pscolors[0]
            plot_like1dname,bfiles,bpname,'H0*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
            plot_like1dname,files,pname,'H0*',subsamp=3,nskip=1000,$
              thick=4,linestyle=ls[1], color=mycolor[1], /oplot
;             plot_like1dname,files_w7s12,pname_w7s12,'H0*',subsamp=3,nskip=1000,$
;               thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            
        endcase
        else: begin
            !x.crange=[0,9.5]
            axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
            axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
            axis,yaxis=0,ytickname=empty,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 
            axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

            oplot,[0,0],[0,0],color=pscolors[0]
            oplot, xtmp, ytmp
        endcase
    endcase
            
endfor
;*******************************
ps_close

print, 'output file: figs/'+filename_stub+'.pdf'
spawn,'epstopdf '+filename_stub+'.ps'

stop
END

