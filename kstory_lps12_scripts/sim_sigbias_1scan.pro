;;;
; NAME: sim_sigbias_map.pro
; PURPOSE:
;   Simulate the signal-bias question from the S12 referee report.
;   Calculate the bias in map-space.
;
; INPUTS: None
;
; OUTPUTS: Sav file
;
; UTIL FUNCTIONS:
;   1) get_weights,     return the relative weights of detectors
;   2) get_map_name,    return the name of a sim map to use
;   3) plot_sim,        plotting routine
;   4) plot_bias,       second plotting routine
;
; NOTES:
;
; MODIFICATION HISTORY:
;  03/04/2013: (KTS) Created
;  03/09/2012: (KTS) Created from sim_sigbias_pow.pro
;  03/17/2012: (KTS) Get rid of nscans from 03/16
;;;


;;;;;;;;;;;;;;
;
; Get the simulated map in units of [K]
;   You must set one (and only one) of the three flags
;
;;;;;;;;;;;;;;
FUNCTION get_sim_map, npix=npix, $
                      cmb_flatsky=cmb_flatsky, $
                      coaddsim=coaddsim, $
                      single_sim=single_sim
                      ;reso

if keyword_set(cmb_flatsky) then begin
    print, ' *** Using simulated map from cmb_flatsky *** '
    restore, '/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_0.sav'
    ss = size(sim_map)
    ret = dblarr(ss[1],ss[2]*7)
    i=0 & ret[*,i*ss[2]:i*ss[2]+ss[2]-1] = sim_map

    restore, '/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_1.sav'
    i=1 & ret[*,i*ss[2]:i*ss[2]+ss[2]-1] = sim_map
    restore, '/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_2.sav'
    i=2 & ret[*,i*ss[2]:i*ss[2]+ss[2]-1] = sim_map
    restore, '/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_3.sav'
    i=3 & ret[*,i*ss[2]:i*ss[2]+ss[2]-1] = sim_map
    restore, '/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_4.sav'
    i=4 & ret[*,i*ss[2]:i*ss[2]+ss[2]-1] = sim_map
    restore, '/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_5.sav'
    i=5 & ret[*,i*ss[2]:i*ss[2]+ss[2]-1] = sim_map
    restore, '/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_6.sav'
    i=6 & ret[*,i*ss[2]:i*ss[2]+ss[2]-1] = sim_map

    sim_map = ret
;     camb2cl_k11,cl_all,ell=ellcmb
;     clcmb=cl_all[*,0]
;     sim_map=cmb_flatsky(ellcmb,clcmb,npix[0],reso*!dtor/60.)
;     save, sim_map, filename='/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_x.sav'
endif

if keyword_set(coaddsim) then begin
    print, ' *** Reading simulated map from coadd *** '
    sim_map = read_dat_map('/home/kstory/lps12/sims/coaddsim_lmax8000_ra4h10dec-50.dat',npix[0],npix[1],1)
endif

if keyword_set(single_sim) then begin
    ; sim_dir     = '/home/kstory/lps12/sims/' ; soft-linked to /data23/
    ; bdir = '/data25/hou/lps12/sim_mapmaking/output/'
    ; spawn, 'ls ' + bdir+fname+'/map_*lmax8000*.bin', list
    ; map_name = list[10] ; arbitrarily take 10'th map
    map_name = '/data25/hou/lps12/sim_mapmaking/output/ra4h10dec-50/map_150_lmax8000_20100212_112516.bin'
    maps = read_bin_map(map_name, /nowt)
    sim_map = maps[*,*,0] * 1d-6
endif

RETURN, sim_map
END



;;;;;;;;;;;;;;
;
; Util function: Get the weights
;
;;;;;;;;;;;;;;
FUNCTION get_weights, wtype, nbins, ncoadd, cor, hz, n_tstream, sn_tstream,$
                      stopit=stopit, narrow_band=narrow_band

wts = fltarr(ncoadd)
psd_1to3 = fltarr(ncoadd)
wh_1to3 = where(hz ge 1 and hz lt 3)

; use narrow-band noise, as a test
if keyword_set(narrow_band) then begin
    wh_1to3 = wh_1to3[0:5] 
endif

case wtype of 
    0: begin                    ; no weights
        wts += 1.
    end
    
    1: begin                    ; noise
;        wts = findgen(ncoadd)+1.
        for jj=0L, ncoadd-1 do begin
            pp = abs(fft(n_tstream[*,jj]))^(2.) * cor
            psd_1to3[jj] = mean(pp[wh_1to3])
            wts[jj] = 1/psd_1to3[jj]
        endfor
    end

    2: begin                    ; signal + noise
        for jj=0L, ncoadd-1 do begin
            pp = abs(fft(sn_tstream[*,jj]))^(2.) * cor
            psd_1to3[jj] = mean(pp[wh_1to3])
            wts[jj] = 1/psd_1to3[jj]
        endfor
    end
end
wts /= total(wts)

if keyword_set(stopit) then stop
RETURN, wts
END



;;;;;;;;;;;;;;
;
; Get the beams
;
;;;;;;;;;;;;;;
FUNCTION get_bl, year, l
simbeamfiles = get_lps12_beams(year, l_beam, bl_beam, /s09)
simbeams=dblarr(n_elements(l))
simbeam_interp=dblarr(n_elements(l))

simbeam=read_ascii(simbeamfiles)
simbeam_interp=interpol(simbeam.field1[1, *], simbeam.field1[0, *], l)
idx=where(simbeam_interp ge 0)
simbeams[idx]= abs(simbeam_interp[idx])
bl = simbeams

RETURN, bl
END



;;;;;;;;;;;;;;
;
; Make high-pass filter
;
;;;;;;;;;;;;;;
FUNCTION get_hpfilt, npix, sample_rate, tau_hz
hpfilt = complexarr(npix) + 1.
hpfilt = hpfilt*single_pole_filter(npix,sample_rate,1./2./!pi/tau_hz,1)
RETURN, hpfilt
END




;;;;;;;;;;;;;;
;
; Calculate the noise rms
;   Returns [ncoadd, nsims]
;   s,   should be sim_st
;
;;;;;;;;;;;;;;
FUNCTION get_sigma_noise, s, use_1to3=use_1to3

dc_level = 500 * sqrt(s.sample_rate) ; units uK

; draw sigma from a distribution
undefine, seed
sigma = randomn(seed, s.ncoadd, s.nsims) * s.bolo_dist + dc_level
sigma = abs(sigma) ; this will mess up the distribution a bit, but that's ok.  We don't want negative sigmas.

if keyword_set(use_1to3) then begin
    bdir = '/home/kstory/lps12/lowellfits/'+fst.name+'/'
    fname = bdir+'field_scan_stub_150_20100212_112516.fits'
    data = krf(fname)
    whgood = where(data.observation.bolo_flags eq 0)
    bolo_id = 100               ; random choice for now
    wn = data.observation.bolo_psd_1to3[bolo_id] ; noise psd between 1 and 3 Hz?
    sigma = fltarr(ncoadd, nsims) + 1*(wn * 1d6 * 8d5) ; Should only be one factor of 1d6
endif

RETURN, sigma
END



;;;;;;;;;;;;;;
;
; Make plots
;
;;;;;;;;;;;;;;
PRO plot_sim, s, sim, tag=tag, pn=pn, wnum=wnum
if n_elements(wnum) eq 0 then wnum=0
hzbins = s.lbins * s.hz[1]/s.l[1]

wh = where(s.lbins lt 4000, nwh)

wset, 0+wnum
!p.multi=[0,2,2]

plot, hzbins[wh], sim.pow_mb[wh,0], /ylog, title=tag+', "measured" power, 1sim',xtitle='Hz',ytitle='Power, [uK^2]'
errplot, hzbins[wh], sim.pow_mb[wh,0]-sim.err_mb, sim.pow_mb[*,0]+sim.err_mb

plot, hzbins[wh], sim.pow_rb[wh,0], title='mean residual power, 1sim',xtitle='Hz',ytitle='Power, [uK^2]'
errplot, hzbins[wh], sim.pow_rb[wh,0]-sim.err_rb, sim.pow_rb[*,0]+sim.err_rb

plot, s.hz[s.vec],pn[s.vec],/ylog,xr=[0,5],title='Power of single timestream',xtitle='Hz',ytitle='Power, [uK^2]'
oplot, s.hz[s.vec],sim.pow_s[s.vec,0],color=!green

plot, hzbins, sim.pow_mb,/ylog,xr=[0,hzbins[wh[nwh-1]]],title='Signal and Measured power',xtitle='Hz',ytitle='Power, [uK^2]'
oplot, hzbins, sim.pow_sb, color=!green
legend,['Signal','Noise'],colors=[!green,!black],linestyle=[0,0]

!p.multi=0
END


;;;;;;;;;;;;;;
;
; Make plots
; INPUTS:
;   sim,        sim structure
;
;;;;;;;;;;;;;;
PRO plot_bias, s, sim0, sim1, sim2, s0,s1,s2,whmask, tag=tag, yr=yr, wnum=wnum
if n_elements(wnum) eq 0 then wnum=0
hzbins = s.lbins * s.hz[1]/s.l[1]

npts = 1e5 ; number of points to plot
if (npts gt n_elements(whmask)) then npts = n_elements(whmask)

!p.multi=[0,2,2]
wset, 0+wnum

; PLOT 1
plot, s0.ms,xtitle='sim number',ytitle='<T_resid>_pix [uK]',title=tag+', Bias by sim',xr=[0,200]
oplot, s1.ms,color=!red
oplot, s2.ms,color=!blue

; PTE plot
; h0 = histogram(sim0.pte_rb, locations=xb0,nbins=11)
; xb0 = xb0+(xb0[1]-xb0[0])/2.
; h1 = histogram(sim1.pte_rb, locations=xb1,nbins=11)
; h2 = histogram(sim2.pte_rb, locations=xb2,nbins=11)
; mx = max([h0,h1,h2])
; plot, xb0, h0, psym=10,xr=[0,1],/xst,yr=[0,mx],/yst,title='PTE distribution'
; oplot, xb0,h1,psym=10,color=!red
; oplot, xb0,h2,psym=10,color=!blue

; PLOT 2
wh = where(s.l ge 4000)
plot, hzbins,sim0.pow_rb[*,0],title='Residual power for single sim',xtitle='Hz',ytitle='pow_resid [uK^2]',xr=[0,s.hz[wh[0]]]
oplot, hzbins,sim1.pow_rb[*,0],color=!red
oplot, hzbins,sim2.pow_rb[*,0],color=!blue
legend,['sim0','sim1','sim2'],colors=[!black,!red,!blue],linestyle=[0,0,0],pos=[2500,50]

; PLOT 3
if n_elements(yr) eq 0 then yr=[min(sim2.coadd_sn[*,0]/sim2.signal[*,0]), max(sim2.coadd_sn[*,0]/sim2.signal[*,0])]
plot, sim2.signal[*,0], sim2.coadd_sn[*,0]/sim2.signal[*,0],psym=3,yr=yr,$
  title='Weight 2',xtitle='T_true',ytitle='T_meas/T_true'
oplot, sim2.signal[whmask[0:npts-1]], sim2.coadd_sn[whmask[0:npts-1]]/sim2.signal[whmask[0:npts-1]],psym=3,color=!blue

; PLOT 4
plot, sim2.signal[*,0], sim2.resid[*,0], psym=3, $
  title='Weight 2',xtitle='T_true',ytitle='T_meas - T_true'
oplot, sim2.signal[whmask[0:npts-1]], sim2.resid[whmask[0:npts-1]],psym=3,color=!blue
; plot, sim2.signal[*,0], sim2.coadd_sn[*,0]-sim2.signal[*,0], psym=3, $
;   xtitle='T_true',ytitle='T_meas - T_true'
; oplot, sim2.signal, sim2.coadd_sn-sim2.signal,psym=3,color=!blue

!p.multi=0
END



;;;;;;;;;;;;;;
;
; Make power-bias plots
;
;;;;;;;;;;;;;;
PRO plot_pbias, s,s0,s1,s2,tag=tag
if n_elements(tag) eq 0 then tag=''
hzbins = s.lbins * s.hz[1]/s.l[1]
wh = where(s.lbins lt 4000, nwh)
xr=[0, hzbins[max(wh)]]

wset, 5
!p.multi=[0,1,3]
chsize = 3

plot, hzbins, s0.rpb,xr=xr,charsize=chsize,title=tag+', Weight 0',xtitle='Hz',ytitle='Power_resid [uK^2]'
errplot, hzbins, s0.rpb-s0.rpb_err, s0.rpb+s0.rpb_err

plot, hzbins, s1.rpb,xr=xr,charsize=chsize,title='Weight 1',xtitle='Hz',ytitle='Power_resid [uK^2]'
errplot, hzbins, s1.rpb-s1.rpb_err, s1.rpb+s1.rpb_err

plot, hzbins, s2.rpb,xr=xr,charsize=chsize,title='Weight 2',xtitle='Hz',ytitle='Power_resid [uK^2]'
errplot, hzbins, s2.rpb-s2.rpb_err, s2.rpb+s2.rpb_err

!p.multi=0

wset, 0
!p.multi=[0,3,2]

plot, s0.r_psd_inband,xr=[0,1000],charsize=chsize,title=tag+', Power in-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, s1.r_psd_inband, color=!red
oplot, s2.r_psd_inband, color=!blue

plot, s0.r_psd_low,xr=[0,1000],charsize=chsize,title='Power low-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, s1.r_psd_low, color=!red
oplot, s2.r_psd_low, color=!blue

plot, s0.r_psd_high,xr=[0,1000],charsize=chsize,title='Power high-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, s1.r_psd_high, color=!red
oplot, s2.r_psd_high, color=!blue

h0 = histogram(s0.r_psd_inband, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(s1.r_psd_inband, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(s2.r_psd_inband, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue


h0 = histogram(s0.r_psd_low, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in low-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(s1.r_psd_low, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(s2.r_psd_low, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue

h0 = histogram(s0.r_psd_high, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in high-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(s1.r_psd_high, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(s2.r_psd_high, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue

!p.multi=0

; Print stats
print, ' --- Sim 0 ---'
print, 'Detect bias at xx sigma: ', (mean(s0.r_psd_inband) / mean(s0.s_psd_inband)) / (stddev(s0.r_psd_inband) / sqrt(s.nsims) / mean(s0.s_psd_inband))
print, 'mean, stddev of r_psd_inband = ', mean(s0.r_psd_inband), stddev(s0.r_psd_inband)
print, 'sigma_r_inband / sqrt(nsims) = ', stddev(s0.r_psd_inband) / sqrt(s.nsims)
print, 'mean_r / mean_s, inband      = ', mean(s0.r_psd_inband) / mean(s0.s_psd_inband)
print, 'sigma(mean_r / mean_s)       = ', stddev(s0.r_psd_inband) / sqrt(s.nsims) / mean(s0.s_psd_inband)

print, ' --- Sim 1 ---'
print, 'Detect bias at xx sigma: ', (mean(s1.r_psd_inband) / mean(s1.s_psd_inband)) / (stddev(s1.r_psd_inband) / sqrt(s.nsims) / mean(s1.s_psd_inband))
print, 'mean, stddev of r_psd_inband = ', mean(s1.r_psd_inband), stddev(s1.r_psd_inband)
print, 'sigma_r_inband / sqrt(nsims) = ', stddev(s1.r_psd_inband) / sqrt(s.nsims)
print, 'mean_r / mean_s, inband      = ', mean(s1.r_psd_inband) / mean(s1.s_psd_inband)
print, 'sigma(mean_r / mean_s)       = ', stddev(s1.r_psd_inband) / sqrt(s.nsims) / mean(s1.s_psd_inband)

print, ' --- Sim 2 ---'
print, 'Detect bias at xx sigma: ', (mean(s2.r_psd_inband) / mean(s2.s_psd_inband)) / (stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband))
print, 'mean, stddev of r_psd_inband = ', mean(s2.r_psd_inband), stddev(s2.r_psd_inband)
print, 'sigma_r_inband / sqrt(nsims) = ', stddev(s2.r_psd_inband) / sqrt(s.nsims)
print, 'mean_r / mean_s, inband      = ', mean(s2.r_psd_inband) / mean(s2.s_psd_inband)
print, 'sigma(mean_r / mean_s)       = ', stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband)
END


;;;;;;;;;;;;;;
;
; Print a summary of bias statistics
;
;;;;;;;;;;;;;;
PRO print_summary, tag=tag
print, tag
restore, 'sav_files/sim_sigbias_map_'+tag+'_stub.sav'
s = sim_st


; print, ' --- Sim 0 ---'
; print, 'Detect Power bias at xx sigma: ', (mean(s0.r_psd_inband) / mean(s0.s_psd_inband)) / (stddev(s0.r_psd_inband) / sqrt(s.nsims) / mean(s0.s_psd_inband))
; print, 'Inband: mean, err             = ', mean(s0.r_psd_inband), stddev(s0.r_psd_inband)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s0.r_psd_inband) / mean(s0.s_psd_inband), stddev(s0.r_psd_inband) / sqrt(s.nsims) / mean(s0.s_psd_inband)
; print, 'Low:    mean, err             = ', mean(s0.r_psd_low), stddev(s0.r_psd_low)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s0.r_psd_low) / mean(s0.s_psd_low), stddev(s0.r_psd_low) / sqrt(s.nsims) / mean(s0.s_psd_low)
; print, 'High:   mean, err             = ', mean(s0.r_psd_high), stddev(s0.r_psd_high)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s0.r_psd_high) / mean(s0.s_psd_high), stddev(s0.r_psd_high) / sqrt(s.nsims) / mean(s0.s_psd_high)
; print, 'Pixel Bias [uK]               = ', s0.m, s0.sd/sqrt(sim_st.nsims)
; print, ''
; print, ' --- Sim 1 ---'
; print, 'Detect Power bias at xx sigma: ', (mean(s1.r_psd_inband) / mean(s1.s_psd_inband)) / (stddev(s1.r_psd_inband) / sqrt(s.nsims) / mean(s1.s_psd_inband))
; print, 'Inband: mean, err             = ', mean(s1.r_psd_inband), stddev(s1.r_psd_inband)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s1.r_psd_inband) / mean(s1.s_psd_inband), stddev(s1.r_psd_inband) / sqrt(s.nsims) / mean(s1.s_psd_inband)
; print, 'Low:    mean, err             = ', mean(s1.r_psd_low), stddev(s1.r_psd_low)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s1.r_psd_low) / mean(s1.s_psd_low), stddev(s1.r_psd_low) / sqrt(s.nsims) / mean(s1.s_psd_low)
; print, 'High:   mean, err             = ', mean(s1.r_psd_high), stddev(s1.r_psd_high)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s1.r_psd_high) / mean(s1.s_psd_high), stddev(s1.r_psd_high) / sqrt(s.nsims) / mean(s1.s_psd_high)
; print, 'Pixel Bias [uK]               = ', s1.m, s1.sd/sqrt(sim_st.nsims)
; print, ''
; print, ' --- Sim 2 ---'
; print, 'Detect Power bias at xx sigma: ', (mean(s2.r_psd_inband) / mean(s2.s_psd_inband)) / (stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband))
; print, 'Inband: mean, err             = ', mean(s2.r_psd_inband), stddev(s2.r_psd_inband)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s2.r_psd_inband) / mean(s2.s_psd_inband), stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband)
; print, 'Low:    mean, err             = ', mean(s2.r_psd_low), stddev(s2.r_psd_low)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s2.r_psd_low) / mean(s2.s_psd_low), stddev(s2.r_psd_low) / sqrt(s.nsims) / mean(s2.s_psd_low)
; print, 'High:   mean, err             = ', mean(s2.r_psd_high), stddev(s2.r_psd_high)/sqrt(s.nsims)
; print, '        mean_r / mean_s, err  = ', mean(s2.r_psd_high) / mean(s2.s_psd_high), stddev(s2.r_psd_high) / sqrt(s.nsims) / mean(s2.s_psd_high)
; print, 'Pixel Bias [uK]               = ', s2.m, s2.sd/sqrt(sim_st.nsims)
; print, ''

print, ' --- Sim 0 ---'
print, 'Detect Power bias at xx sigma: ', (mean(s0.r_psd_inband) / mean(s0.s_psd_inband)) / (stddev(s0.r_psd_inband) / sqrt(s.nsims) / mean(s0.s_psd_inband))
print, 'Pixel Bias [uK]              : ', s0.m, s0.sd/sqrt(sim_st.nsims)
print, '                           Inband                     Low-band                  High-band'
print, 'mean, err             = ', mean(s0.r_psd_inband), stddev(s0.r_psd_inband)/sqrt(s.nsims), mean(s0.r_psd_low), stddev(s0.r_psd_low)/sqrt(s.nsims), mean(s0.r_psd_high), stddev(s0.r_psd_high)/sqrt(s.nsims)
print, 'mean_r / mean_s, err  = ', mean(s0.r_psd_inband) / mean(s0.s_psd_inband), stddev(s0.r_psd_inband) / sqrt(s.nsims) / mean(s0.s_psd_inband), mean(s0.r_psd_low) / mean(s0.s_psd_low), stddev(s0.r_psd_low) / sqrt(s.nsims) / mean(s0.s_psd_low), mean(s0.r_psd_high) / mean(s0.s_psd_high), stddev(s0.r_psd_high) / sqrt(s.nsims) / mean(s0.s_psd_high)
print, ''

print, ' --- Sim 1 ---'
print, 'Detect Power bias at xx sigma: ', (mean(s1.r_psd_inband) / mean(s1.s_psd_inband)) / (stddev(s1.r_psd_inband) / sqrt(s.nsims) / mean(s1.s_psd_inband))
print, 'Pixel Bias [uK]              : ', s1.m, s1.sd/sqrt(sim_st.nsims)
print, '                           Inband                     Low-band                  High-band'
print, 'mean, err             = ', mean(s1.r_psd_inband), stddev(s1.r_psd_inband)/sqrt(s.nsims), mean(s1.r_psd_low), stddev(s1.r_psd_low)/sqrt(s.nsims), mean(s1.r_psd_high), stddev(s1.r_psd_high)/sqrt(s.nsims)
print, 'mean_r / mean_s, err  = ', mean(s1.r_psd_inband) / mean(s1.s_psd_inband), stddev(s1.r_psd_inband) / sqrt(s.nsims) / mean(s1.s_psd_inband), mean(s1.r_psd_low) / mean(s1.s_psd_low), stddev(s1.r_psd_low) / sqrt(s.nsims) / mean(s1.s_psd_low), mean(s1.r_psd_high) / mean(s1.s_psd_high), stddev(s1.r_psd_high) / sqrt(s.nsims) / mean(s1.s_psd_high)
print, ''

print, ' --- Sim 2 ---'
print, 'Detect Power bias at xx sigma: ', (mean(s2.r_psd_inband) / mean(s2.s_psd_inband)) / (stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband))
print, 'Pixel Bias [uK]              : ', s2.m, s2.sd/sqrt(sim_st.nsims)
print, '                           Inband                     Low-band                  High-band'
print, 'mean, err             = ', mean(s2.r_psd_inband), stddev(s2.r_psd_inband)/sqrt(s.nsims), mean(s2.r_psd_low), stddev(s2.r_psd_low)/sqrt(s.nsims), mean(s2.r_psd_high), stddev(s2.r_psd_high)/sqrt(s.nsims)
print, 'mean_r / mean_s, err  = ', mean(s2.r_psd_inband) / mean(s2.s_psd_inband), stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband), mean(s2.r_psd_low) / mean(s2.s_psd_low), stddev(s2.r_psd_low) / sqrt(s.nsims) / mean(s2.s_psd_low), mean(s2.r_psd_high) / mean(s2.s_psd_high), stddev(s2.r_psd_high) / sqrt(s.nsims) / mean(s2.s_psd_high)
print, ''

END





;;;;;;;;;;;;;;
;
; Signal bias sims, map level
;
;;;;;;;;;;;;;;
PRO sim_sigbias_map, ncoadd=ncoadd, nsims=nsims, stopit=stopit, $
                     bolo_dist=bolo_dist, $ ; use a distribution of bolo noises
                     sigfac=sigfac, sav_name=sav_name, verbose=verbose

if n_elements(ncoadd) eq 0 then ncoadd = 50 ; number of coadds per sim
if n_elements(nsims) eq 0  then nsims  = 1000   ; number of sims

if n_elements(bolo_dist) eq 0 then bolo_dist = 2e2
if n_elements(sav_name) eq 0 then print, 'Not saving any output files!'
if n_elements(sigfac) eq 0 then sigfac = 1.
print, ' --- Scaling the signal by a factor of '+strtrim(string(sigfac),2)+' ---' 

; ------------------------
; Timing
; ------------------------
t_start = systime(/seconds)
t_stop = double(0.)

; ------------------------
; Setup
; ------------------------
field_idx = 6  ; ra4h10dec-50
scan_speed = 0.417812 ; dps on-sky
reso = 1.0
reso_rad = reso/60.*!dtor
sample_rate = scan_speed*60./reso

f = lps12_fieldstruct()
fst = f[field_idx]
info = get_lps12_fieldinfo(field_idx)
npix = info.npix
fname = fst.name

l = make_fft_grid(reso/60.*!dtor/2./!pi,npix[0],fx=lx)
hz = make_fft_grid(reso/60./scan_speed,npix[0])
vec = indgen(npix[0]/2 - 1)

; For binning the residuals
dl = 50
max_l = max(l[vec])
nbins = floor(max_l / dl)
lbins = findgen(nbins) * dl + dl/2.
hzbins = lbins * hz[1]/l[1]

; Convenient structure for passing to plotting routines
sim_st = {field_idx:field_idx, $  ; ra4h10dec-50
          scan_speed: scan_speed, $ ; dps on-sky
          reso: reso, $
          reso_rad: reso_rad, $
          sample_rate: sample_rate, $
          npix:info.npix, $
          ncoadd : ncoadd, $
          nsims : nsims, $
          l: l, $
          hz: hz, $
          vec: vec, $
          bolo_dist: bolo_dist, $
          sigfac: sigfac, $
          dl : dl, $
          max_l : max_l, $
          nbins : nbins, $
          lbins : lbins, $
          hzbins : hzbins}

; ------------------------
; Signal Timestream setup
; ------------------------
; if (nsims gt npix[0]) then begin
;     print, 'nsims gt npix.  Fix this.'
;     stop
; endif
sim_map = get_sim_map(/cmb_flatsky)
ss = size(sim_map)
sig_idx = scramble(indgen(ss[2]))
num_sig_idx = n_elements(sig_idx)

; get the mask
restore, '/home/kstory/lps12/masks/apod/apod_'+fname+'_60_0.0500_30.sav'
mask = apod
mask_row = mask[*,npix[1]/2]
mask_fac = mean(mask_row^2.)
cor = (1/mask_fac) ; correction factor for PS

; get high-pass filter
wh = where(l gt 200)
tau_hz = hz[ min(wh)]
hpfilt = get_hpfilt(npix[0],sample_rate,tau_hz)

; Get the beams
year = fst.year
bl = get_bl(year, l)
bl *= 0
bl += 1. ; No beams for cmb_flatsky

; ------------------------
; Noise timetsreams
; ------------------------
sigma = get_sigma_noise(sim_st) ; 500 uK \pm bolo_dist

; sim structure
tmp = {coadd_sn:fltarr(npix[0], nsims), $ ; coadded signal+noise timestream
       coadd_n:fltarr(npix[0], nsims), $ ; coadded noise timestream
       signal:fltarr(npix[0], nsims), $ ; signal timestream
       resid:fltarr(npix[0], nsims), $ ; T_meas - T_true
       wts:fltarr(ncoadd, nsims), $ ; weights
       pow_m :fltarr(npix[0], nsims), $ ; measured power
       pow_n :fltarr(npix[0], nsims), $ ; noise power
       pow_r :fltarr(npix[0], nsims), $ ; residual power
       pow_s :fltarr(npix[0], nsims), $ ; signal power
       pow_mb:fltarr(nbins, nsims), $ ; pow_m, bin-averaged
       pow_nb:fltarr(nbins, nsims), $ ; pow_n, bin-averaged
       pow_rb:fltarr(nbins, nsims), $ ; pow_r, bin-averaged
       pow_sb:fltarr(nbins, nsims), $ ; pow_s, bin-averaged
       err_mb:fltarr(nbins), $  ; measured power errors
       err_nb:fltarr(nbins), $  ; measured power errors
       err_rb:fltarr(nbins), $  ; residual power errors
       chisq_rb:fltarr(nsims), $
       pte_rb:fltarr(nsims) $
      }



n_tstream = fltarr(npix[0], ncoadd) ; noise timestreams
sn_tstream = fltarr(npix[0], ncoadd) ; signal+noise timestreams

sim0 = replicate(tmp,1)
sim1 = replicate(tmp,1)
sim2 = replicate(tmp,1)


; ------------------------
; Coadd noise and calculate power
; ------------------------

print, "isim = ", format="(A,$)"
for isim=0L, nsims-1 do begin

    if keyword_set(verbose) then begin 
        print, "isim = ", isim
    endif else begin 
        if (isim mod 10 eq 0) then print, isim, format='((x, I0),$)' 
    endelse
    
    ; ------------------------
    ; Signal timestream
    idx = isim mod num_sig_idx ; repeat use of signal timestreams
    signal0 = reform( sim_map[*, idx] ) * 1d6 ; units of uK
    signal = signal0 * mask_row
    signal *= sigfac            ; scale the signal, if desired
    
    signal = float(fft(fft(signal)*hpfilt,1)) ; filter the timestream
    
    
    ; Calculate signal power
    sim0.pow_s[*,isim]   = abs(fft(signal))^(2.) * cor / (bl^2.)
    sim0.pow_sb[*,isim]  = fltarr(nbins)
    for ii=0, nbins-1 do begin 
        wh = where(l ge dl*ii and l lt dl*(ii+1))
        sim0.pow_sb[ii,isim] = mean(sim0.pow_s[wh,isim])
    endfor
    
    ; copy power to all sims
    sim1.pow_s[*,isim]  = sim0.pow_s[*,isim]
    sim1.pow_sb[*,isim] = sim0.pow_sb[*,isim]
    sim2.pow_s[*,isim]  = sim0.pow_s[*,isim]
    sim2.pow_sb[*,isim] = sim0.pow_sb[*,isim]
    
    sim0.signal[*,isim] = signal
    sim1.signal[*,isim] = signal
    sim2.signal[*,isim] = signal
    
    
    ; ------------------------
    ; Noise timestreams; coadd many "bolometers," realizations of noise
    
    n_tstream   *= 0
    sn_tstream  *= 0
    
    for jj=0L, ncoadd-1 do begin
        undefine, seed
        noise_tmp = randomn(seed, npix[0]) * sigma[jj,isim]
        noise_tmp *= mask_row   ; apply mask
        noise_tmp = float(fft(fft(noise_tmp)*hpfilt,1)) ; high-pass filter
        
        n_tstream[*,jj]  = noise_tmp
        sn_tstream[*,jj] = signal + noise_tmp
    endfor
    
    ; ------------------------
    ; Get the weights
    ; ------------------------
    sim0.wts[*,isim] = get_weights(0, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)
    sim1.wts[*,isim] = get_weights(1, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)
    sim2.wts[*,isim] = get_weights(2, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)


    ; coadd the timestreams
    for k=0, npix[0]-1 do begin
        sim0.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim0.wts[*,isim])
        sim1.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim1.wts[*,isim])
        sim2.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim2.wts[*,isim])
        
        sim0.coadd_n[k,isim] += total(n_tstream[k,*] * sim0.wts[*,isim])
        sim1.coadd_n[k,isim] += total(n_tstream[k,*] * sim1.wts[*,isim])
        sim2.coadd_n[k,isim] += total(n_tstream[k,*] * sim2.wts[*,isim])
    endfor
        
    ; Calculate map residuals
    sim0.resid[*,isim] = sim0.coadd_sn[*,isim] - sim0.signal[*,isim]
    sim1.resid[*,isim] = sim1.coadd_sn[*,isim] - sim1.signal[*,isim]
    sim2.resid[*,isim] = sim2.coadd_sn[*,isim] - sim2.signal[*,isim]
    
    
    
    ; Calculate the measured power, and residuals from expected power
    sim0.pow_m[*,isim] = abs(fft(sim0.coadd_sn[*,isim]))^(2.) * cor
    sim0.pow_n[*,isim] = abs(fft(sim0.coadd_n[*,isim]))^(2.) * cor
    sim0.pow_r[*,isim] = sim0.pow_m[*,isim] - (sim0.pow_s[*,isim] + sim0.pow_n[*,isim])
    
    sim1.pow_m[*,isim] = abs(fft(sim1.coadd_sn[*,isim]))^(2.) * cor
    sim1.pow_n[*,isim] = abs(fft(sim1.coadd_n[*,isim]))^(2.) * cor
    sim1.pow_r[*,isim] = sim1.pow_m[*,isim] - (sim1.pow_s[*,isim] + sim1.pow_n[*,isim])
    
    sim2.pow_m[*,isim] = abs(fft(sim2.coadd_sn[*,isim]))^(2.) * cor
    sim2.pow_n[*,isim] = abs(fft(sim2.coadd_n[*,isim]))^(2.) * cor
    sim2.pow_r[*,isim] = sim2.pow_m[*,isim] - (sim2.pow_s[*,isim] + sim2.pow_n[*,isim])
    
    ; Bin the power residuals
    for ii=0, nbins-1 do begin
        wh = where(l ge dl*ii and l lt dl*(ii+1))
        sim0.pow_mb[ii,isim] = mean(sim0.pow_m[wh,isim])
        sim0.pow_nb[ii,isim] = mean(sim0.pow_n[wh,isim])
        sim0.pow_rb[ii,isim] = mean(sim0.pow_r[wh,isim])
        sim1.pow_mb[ii,isim] = mean(sim1.pow_m[wh,isim])
        sim1.pow_nb[ii,isim] = mean(sim1.pow_n[wh,isim])
        sim1.pow_rb[ii,isim] = mean(sim1.pow_r[wh,isim])
        sim2.pow_mb[ii,isim] = mean(sim2.pow_m[wh,isim])
        sim2.pow_nb[ii,isim] = mean(sim2.pow_n[wh,isim])
        sim2.pow_rb[ii,isim] = mean(sim2.pow_r[wh,isim])
    end
    
endfor
print, ''

; ------------------------
; Calculate the power error in each bin
; ------------------------
print, '--- calculate error ---'
for ii=0, nbins-1 do begin
    sim0.err_mb[ii] = stddev(sim0.pow_mb[ii,*])
    sim0.err_nb[ii] = stddev(sim0.pow_nb[ii,*])
    sim0.err_rb[ii] = stddev(sim0.pow_rb[ii,*])
    sim1.err_mb[ii] = stddev(sim1.pow_mb[ii,*])
    sim1.err_nb[ii] = stddev(sim1.pow_nb[ii,*])
    sim1.err_rb[ii] = stddev(sim1.pow_rb[ii,*])
    sim2.err_mb[ii] = stddev(sim2.pow_mb[ii,*])
    sim2.err_nb[ii] = stddev(sim2.pow_nb[ii,*])
    sim2.err_rb[ii] = stddev(sim2.pow_rb[ii,*])
endfor

print, '--- calculate chisq ---'
wh = where(lbins ge 500 and lbins le 3000, nwh)

for jj=0L, nsims-1 do begin
    sim0.chisq_rb[jj] = total( (sim0.pow_rb[wh,jj]/sim0.err_rb[wh])^2.)
    sim0.pte_rb[jj]   = mpchitest(sim0.chisq_rb[jj], nwh)
    sim1.chisq_rb[jj] = total( (sim1.pow_rb[wh,jj]/sim1.err_rb[wh])^2.)
    sim1.pte_rb[jj]   = mpchitest(sim1.chisq_rb[jj], nwh)
    sim2.chisq_rb[jj] = total( (sim2.pow_rb[wh,jj]/sim2.err_rb[wh])^2.)
    sim2.pte_rb[jj]   = mpchitest(sim2.chisq_rb[jj], nwh)
endfor

; ------------------------
; Calculate bias
; ------------------------

; m "mean <T_resid>"
; ms "mean <T_resid> in each sim"
; sd "stddev of ms over sims"
; rpb "residual power binned, averaged over all sims"
; rpb_err "error on the above"
; psd_inband "power per Hz in the 1-3 Hz band [uK^2/Hz]"
; psd_low "power per Hz in the 0-1 Hz band [uK^2/Hz]"
; psd_high "power per Hz in the 3Hz-10,000(l) band [uK^2/Hz]"

tmp1 = {m:0., ms:fltarr(nsims), sd:0., $
        rpb:fltarr(nbins), rpb_err:fltarr(nbins), $
        r_psd_inband:fltarr(nsims), r_psd_low:fltarr(nsims), r_psd_high:fltarr(nsims), $
        n_psd_inband:fltarr(nsims), n_psd_low:fltarr(nsims), n_psd_high:fltarr(nsims), $
        s_psd_inband:fltarr(nsims), s_psd_low:fltarr(nsims), s_psd_high:fltarr(nsims)}

s0 = replicate(tmp1, 1)
s1 = replicate(tmp1, 1)
s2 = replicate(tmp1, 1)

; ------------------------
; Bias in pixel values
mask_expanded = [0]
for i=0L,nsims-1 do mask_expanded = [mask_expanded, mask_row]
mask_expanded = mask_expanded[1:*]
whmask = where(mask_expanded gt 0.9)
whp = where((sim0.signal gt 0) and (mask_expanded gt 0.9) )
whn = where((sim0.signal lt 0) and (mask_expanded gt 0.9) )

for isim=0L, nsims-1 do begin
    whp = where((sim0.signal[*,isim] gt 0) and (mask_row gt 0.9) )
    whn = where((sim0.signal[*,isim] lt 0) and (mask_row gt 0.9) )

    s0.ms[isim]   = mean( [sim0.resid[whp,isim], -1*sim0.resid[whn,isim]] )
    s1.ms[isim]   = mean( [sim1.resid[whp,isim], -1*sim1.resid[whn,isim]] )
    s2.ms[isim]   = mean( [sim2.resid[whp,isim], -1*sim2.resid[whn,isim]] )
endfor

s0.m   = mean(   s0.ms )
s0.sd  = stddev( s0.ms )
s1.m   = mean(   s1.ms )
s1.sd  = stddev( s1.ms )
s2.m   = mean(   s2.ms )
s2.sd  = stddev( s2.ms )


; ------------------------
; Bias in power

; get 1 to 3 Hz bins
wh = where(hz ge 1 and hz lt 3)
lmin = l[wh[0]]
lmax = l[wh[n_elements(wh)-1]]
wh_1to3 = where((lbins ge lmin) and (lbins le lmax), complement=whnot)

; get bandwidth outside of 1-3Hz, with lmax=10,000
l1 = lmax
l2 = 10000
wh_low  = where(lbins lt lmin) ; 0-1Hz
wh_high = where((lbins gt lmax) and (lbins lt 10000)) ; above 3Hz

wh = where((l gt l1) and (l lt l2)) ; indexes l
hz_out = [hz[wh[0]], hz[wh[n_elements(wh)-1]]]

nHz_inband = 3-1.
nHz_low    = 1.
nHz_high   = hz_out[1] - hz_out[0]
dHz = dl * (hz[1]/l[1]) ; size of a bin in Hz

;integrate in-band and out-of-band power in the residuals
for j=0L,nsims-1 do begin
    s0.r_psd_inband[j] = total(sim0.pow_rb[wh_1to3,j]) * dHz / nHz_inband
    s0.r_psd_low[j]    = total(sim0.pow_rb[wh_low,j])  * dHz / nHz_low
    s0.r_psd_high[j]   = total(sim0.pow_rb[wh_high,j]) * dHz / nHz_high

    s0.s_psd_inband[j] = total(sim0.pow_sb[wh_1to3,j]) * dHz / nHz_inband
    s0.s_psd_low[j]    = total(sim0.pow_sb[wh_low,j])  * dHz / nHz_low
    s0.s_psd_high[j]   = total(sim0.pow_sb[wh_high,j]) * dHz / nHz_high

    s0.n_psd_inband[j] = total(sim0.pow_nb[wh_1to3,j]) * dHz / nHz_inband
    s0.n_psd_low[j]    = total(sim0.pow_nb[wh_low,j])  * dHz / nHz_low
    s0.n_psd_high[j]   = total(sim0.pow_nb[wh_high,j]) * dHz / nHz_high


    s1.r_psd_inband[j] = total(sim1.pow_rb[wh_1to3,j]) * dHz / nHz_inband
    s1.r_psd_low[j]    = total(sim1.pow_rb[wh_low,j])  * dHz / nHz_low
    s1.r_psd_high[j]   = total(sim1.pow_rb[wh_high,j]) * dHz / nHz_high

    s1.s_psd_inband[j] = total(sim1.pow_sb[wh_1to3,j]) * dHz / nHz_inband
    s1.s_psd_low[j]    = total(sim1.pow_sb[wh_low,j])  * dHz / nHz_low
    s1.s_psd_high[j]   = total(sim1.pow_sb[wh_high,j]) * dHz / nHz_high

    s1.n_psd_inband[j] = total(sim1.pow_nb[wh_1to3,j]) * dHz / nHz_inband
    s1.n_psd_low[j]    = total(sim1.pow_nb[wh_low,j])  * dHz / nHz_low
    s1.n_psd_high[j]   = total(sim1.pow_nb[wh_high,j]) * dHz / nHz_high


    s2.r_psd_inband[j] = total(sim2.pow_rb[wh_1to3,j]) * dHz / nHz_inband
    s2.r_psd_low[j]    = total(sim2.pow_rb[wh_low,j])  * dHz / nHz_low
    s2.r_psd_high[j]   = total(sim2.pow_rb[wh_high,j]) * dHz / nHz_high

    s2.s_psd_inband[j] = total(sim2.pow_sb[wh_1to3,j]) * dHz / nHz_inband
    s2.s_psd_low[j]    = total(sim2.pow_sb[wh_low,j])  * dHz / nHz_low
    s2.s_psd_high[j]   = total(sim2.pow_sb[wh_high,j]) * dHz / nHz_high

    s2.n_psd_inband[j] = total(sim2.pow_nb[wh_1to3,j]) * dHz / nHz_inband
    s2.n_psd_low[j]    = total(sim2.pow_nb[wh_low,j])  * dHz / nHz_low
    s2.n_psd_high[j]   = total(sim2.pow_nb[wh_high,j]) * dHz / nHz_high
endfor

for i=0, nbins-1 do begin
    s0.rpb[i] = mean(sim0.pow_rb[i,0:nsims-1])
    s1.rpb[i] = mean(sim1.pow_rb[i,0:nsims-1])
    s2.rpb[i] = mean(sim2.pow_rb[i,0:nsims-1])
endfor

s0.rpb_err = sim0.err_rb / sqrt(nsims)
s1.rpb_err = sim1.err_rb / sqrt(nsims)
s2.rpb_err = sim2.err_rb / sqrt(nsims)


; Plot stuff
pn = abs(fft(noise_tmp))^(2.) * cor ; for plotting
if n_elements(sav_name) eq 0 then begin
    plot_sim, sim_st, sim2, tag='sim_sigbias_map', pn=pn
    plot_bias, sim_st, sim0, sim1, sim2, s0,s1,s2,whmask,tag='sim_sigbias_map'
    ;plot_pbias,sim_st,s0,s1,s2
endif


; Timing
t_stop = systime(/seconds)
t_tot = t_stop - t_start
min = floor(t_tot / 60.)
print,  '*** Simulation Time ***'
print, min, ' minutes, ', t_tot - 60*min, ' sec.'

if n_elements(sav_name) ne 0 then save, sim_st,wh_1to3,hz,mask_row,whmask,mask_fac,cor,hpfilt,sigma,n_tstream,sn_tstream,pn,$
  sim0,sim1,sim2,s0,s1,s2,$
  filename=sav_name
if keyword_set(stopit) then stop
END



