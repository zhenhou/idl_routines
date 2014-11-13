;;;
; NAME: script13_0318
;
; NOTES:
;  1) re-process 1e5 sim
;;;



;;;;;;;;;;;;;;
;
; Fix 0317_1e5_realistic
;
;;;;;;;;;;;;;;
PRO fix_sim
stub = '0317_1e5_realistic'

sav_name = 'sav_files/sim_sigbias_map_'+stub+'.sav'
restore, sav_name

nsims = sim_st.nsims
nbins = sim_st.nbins
l = sim_st.l
hz = sim_st.hz
dl = sim_st.dl
lbins = sim_st.lbins

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
    print, "isim = ", j
    
    
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
;pn = abs(fft(noise_tmp))^(2.) * cor ; for plotting
if n_elements(sav_name) eq 0 then begin
    plot_sim, sim_st, sim2, stub='sim_sigbias_map', pn=pn
    plot_bias, sim_st, sim0, sim1, sim2, s0,s1,s2,whmask,stub='sim_sigbias_map'
    ;plot_pbias,sim_st,s0,s1,s2
endif


; Timing
t_stop = systime(/seconds)
t_tot = t_stop - t_start
min = floor(t_tot / 60.)
print,  '*** Simulation Time ***'
print, min, ' minutes, ', t_tot - 60*min, ' sec.'

save, sim_st,wh_1to3,hz,mask_row,whmask,mask_fac,cor,hpfilt,sigma,n_tstream,sn_tstream,pn,$
  sim0,sim1,sim2,s0,s1,s2,$
  filename=sav_name

stop
END





;;;;;;;;;;;;;;
;
; Make sim stub
;
;;;;;;;;;;;;;;
PRO make_sim_stub, tag=tag
;tag = '0317_1e5_realistic'
;sav_name = 
print, 'restore, sav_files/sim_sigbias_map_'+tag+'.sav'
restore, 'sav_files/sim_sigbias_map_'+tag+'.sav'

sav_name_stub = 'sav_files/sim_sigbias_map_'+tag+'_stub.sav'


; Make sim_st if it doesn't exist.
if n_elements(sim_st) eq 0 then begin
    print, '*** make sim_st ***'

    if n_elements(bolo_dist) eq 0 then bolo_dist = 2e2
    if n_elements(sigfac) eq 0 then sigfac = 1.
    nsims = n_elements(sim0.chisq_rb)
    ncoadd = n_elements(sim0.wts[*,0])

    field_idx = 6               ; ra4h10dec-50
    scan_speed = 0.417812       ; dps on-sky
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
    
    ; Convenient structure for passing to plotting routines
    sim_st = {field_idx:field_idx, $ ; ra4h10dec-50
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
              lbins : lbins}
endif

print, 'save file to '+sav_name_stub
save, sim_st,wh_1to3,hz,mask_row,whmask,mask_fac,cor,hpfilt,sigma,n_tstream,sn_tstream,pn,$
  s0,s1,s2,$
  filename=sav_name_stub
END




;;;;;;;;;;;;;;
;
; Plot stuff
;
;;;;;;;;;;;;;;
PRO ps1, tag=tag, stopit=stopit
restore, 'sav_files/sim_sigbias_map_'+tag+'_stub.sav'
plot_pbias,sim_st,s0,s1,s2,tag=tag

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
wset, 0
err = tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+tag+'_presid')

wset, 5
err = tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+tag+'_pbias')
if keyword_set(stopit) then stop
END



;;;;;;;;;;;;;;
;
; Look at how the bias averages down with nscans
;
;;;;;;;;;;;;;;
PRO nscans, npix=npix

if (n_elements(npix) eq 0) then npix = 1.e3
sig_n = 5.
sig_s = 5.

undefine, seed
n = randomn(seed, npix) * sig_n
undefine, seed
s = randomn(seed, npix) * sig_s

pow_n  = abs(fft(n))^2.
pow_s  = abs(fft(s))^2.
pow_sn = abs(fft(s+n))^2.
pow_r = pow_sn - pow_n - pow_s

; print, 'Noise : '
; print, rms(n) & print, rms(pow_n) & print, ''
; print, 'Signal: '
; print, rms(s) & print, rms(pow_s) & print, ''
; print, 'S+N   : '
; print, rms(s+n) & print, rms(pow_sn) & print, ''

print, 'Resid : '
print, rms(pow_r) & print, ''

stop



END





PRO scratch
tag = '0318_1e4_realistic'
;tag = '0318_1e4_sn1'
;tag = '0317_1e5_realistic'

print, tag
restore, 'sav_files/sim_sigbias_map_'+tag+'_stub.sav'
;restore, 'sav_files/sim_sigbias_map_'+tag+'.sav'
s = sim_st

if 0 then begin
print, ' --- Sim 0 ---'
print, 'Detect bias at xx sigma: ', (mean(s0.r_psd_inband) / mean(s0.s_psd_inband)) / (stddev(s0.r_psd_inband) / sqrt(s.nsims) / mean(s0.s_psd_inband))
print, 'mean, err r_psd_inband = ', mean(s0.r_psd_inband), stddev(s0.r_psd_inband)/sqrt(s.nsims)
print, 'mean_r / mean_s, inband      = ', mean(s0.r_psd_inband) / mean(s0.s_psd_inband)
print, 'sigma(mean_r / mean_s)       = ', stddev(s0.r_psd_inband) / sqrt(s.nsims) / mean(s0.s_psd_inband)

print, ' --- Sim 1 ---'
print, 'Detect bias at xx sigma: ', (mean(s1.r_psd_inband) / mean(s1.s_psd_inband)) / (stddev(s1.r_psd_inband) / sqrt(s.nsims) / mean(s1.s_psd_inband))
print, 'mean, err r_psd_inband = ', mean(s1.r_psd_inband), stddev(s1.r_psd_inband)/sqrt(s.nsims)
print, 'mean_r / mean_s, inband      = ', mean(s1.r_psd_inband) / mean(s1.s_psd_inband)
print, 'sigma(mean_r / mean_s)       = ', stddev(s1.r_psd_inband) / sqrt(s.nsims) / mean(s1.s_psd_inband)

print, ' --- Sim 2 ---'
print, 'Detect bias at xx sigma: ', (mean(s2.r_psd_inband) / mean(s2.s_psd_inband)) / (stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband))
print, 'mean, err r_psd_inband = ', mean(s2.r_psd_inband), stddev(s2.r_psd_inband)/sqrt(s.nsims)
print, 'mean_r / mean_s, inband      = ', mean(s2.r_psd_inband) / mean(s2.s_psd_inband)
print, 'sigma(mean_r / mean_s)       = ', stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband)
endif

if 1 then begin
;     nsamp = 0L
;     for isim=0L, s.nsims-1 do begin
;         whp = where((sim0.signal[*,isim] gt 0) and (mask_row gt 0.9) )
;         whn = where((sim0.signal[*,isim] lt 0) and (mask_row gt 0.9) )
;         nsamp += n_elements(whp) + n_elements(whn) 
;     endfor

    print, ' --- Sim 0 ---'
    print, s0.m, s0.sd/sqrt(sim_st.nsims)
    print, ' --- Sim 1 ---'
    print, s1.m, s1.sd/sqrt(sim_st.nsims)
    print, ' --- Sim 2 ---'
    print, s2.m, s2.sd/sqrt(sim_st.nsims)

    stop
endif

if 0 then begin
print, ' --- Sim 2 ---'
print, ' inband '
print, 'Detect bias at xx sigma: ', (mean(s2.r_psd_inband) / mean(s2.s_psd_inband)) / (stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband))
print, 'mean, err r_psd_inband = ', mean(s2.r_psd_inband), stddev(s2.r_psd_inband)/sqrt(s.nsims)
print, 'mean_r / mean_s, inband      = ', mean(s2.r_psd_inband) / mean(s2.s_psd_inband)
print, 'sigma(mean_r / mean_s)       = ', stddev(s2.r_psd_inband) / sqrt(s.nsims) / mean(s2.s_psd_inband)
print, ''

print, ' low '
print, 'Detect bias at xx sigma: ', (mean(s2.r_psd_low) / mean(s2.s_psd_low)) / (stddev(s2.r_psd_low) / sqrt(s.nsims) / mean(s2.s_psd_low))
print, 'mean, err r_psd_inband = ', mean(s2.r_psd_low), stddev(s2.r_psd_low)/sqrt(s.nsims)
print, 'mean_r / mean_s, low      = ', mean(s2.r_psd_low) / mean(s2.s_psd_low)
print, 'sigma(mean_r / mean_s)       = ', stddev(s2.r_psd_low) / sqrt(s.nsims) / mean(s2.s_psd_low)
print, ''

print, ' high '
print, 'Detect bias at xx sigma: ', (mean(s2.r_psd_high) / mean(s2.s_psd_high)) / (stddev(s2.r_psd_high) / sqrt(s.nsims) / mean(s2.s_psd_high))
print, 'mean, err r_psd_inband = ', mean(s2.r_psd_high), stddev(s2.r_psd_high)/sqrt(s.nsims)
print, 'mean_r / mean_s, high      = ', mean(s2.r_psd_high) / mean(s2.s_psd_high)
print, 'sigma(mean_r / mean_s)       = ', stddev(s2.r_psd_high) / sqrt(s.nsims) / mean(s2.s_psd_high)
endif

END





PRO plot_pix
tag = '0318_1e4_realistic'
;tag = '0318_1e4_sn1'
;tag = '0317_1e5_realistic'

print, tag
restore, 'sav_files/sim_sigbias_map_'+tag+'_stub.sav'
;restore, 'sav_files/sim_sigbias_map_'+tag+'.sav'
s = sim_st

print, ' --- Sim 0 ---'
print, s0.m, s0.sd/sqrt(sim_st.nsims)
print, ' --- Sim 1 ---'
print, s1.m, s1.sd/sqrt(sim_st.nsims)
print, ' --- Sim 2 ---'
print, s2.m, s2.sd/sqrt(sim_st.nsims)


h0 = histogram(s0.ms, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title=tag+', Temperature Residuals',xtitle='T_resid [uK]'
h1 = histogram(s1.ms, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(s2.ms, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue

stop
END



