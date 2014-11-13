;;;
; NAME: sim_sigbias.pro
; PURPOSE:
;   Simulate the signal-bias question from the S12 referee report
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
;  03/04/2012: (KTS) Created
;  03/08/2012: (KTS) Change nbig from 2160 to 4320
;;;


;;;;;;;;;;;;;;
;
; Get the simulated map in units of [K]
;   You must set one (and only one) of the three flags
;
;;;;;;;;;;;;;;
FUNCTION get_sim_map, npix, $
                      cmb_flatsky=cmb_flatsky, $
                      coaddsim=coaddsim, $
                      single_sim=single_sim
                      ;reso

if keyword_set(cmb_flatsky) then begin
    print, ' *** Using simulated map from cmb_flatsky *** '
    restore, '/home/kstory/lps12/scripts/sav_files/cmb_flatsky_map_20130308.sav'
;     camb2cl_k11,cl_all,ell=ellcmb
;     clcmb=cl_all[*,0]
;     sim_map=cmb_flatsky(ellcmb,clcmb,npix[0],reso*!dtor/60.)
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
;
;;;;;;;;;;;;;;
FUNCTION get_sigma_noise, ncoadd, sample_rate, use_1to3=use_1to3
sigma = fltarr(ncoadd) + 500 * sqrt(sample_rate) ; units uK

if keyword_set(use_1to3) then begin
    bdir = '/home/kstory/lps12/lowellfits/'+fst.name+'/'
    fname = bdir+'field_scan_stub_150_20100212_112516.fits'
    data = krf(fname)
    whgood = where(data.observation.bolo_flags eq 0)
    bolo_id = 100               ; random choice for now
    wn = data.observation.bolo_psd_1to3[bolo_id] ; noise psd between 1 and 3 Hz?
    sigma = fltarr(ncoadd) + 1*(wn * 1d6 * 8d5) ; Should only be one factor of 1d6
endif

RETURN, sigma
END



;;;;;;;;;;;;;;
;
; Make plots
; INPUTS:
;   sim,        sim structure
;
;;;;;;;;;;;;;;
PRO plot_sim, sim, lbins, l, vec, pn, wnum=wnum
if n_elements(wnum) eq 0 then wnum=0

wh = where(lbins lt 4000)

wset, 0+wnum
!p.multi=[0,2,2]

plot, lbins[wh], sim.pow_mb[wh,0], /ylog, title='mean pow_meas'
errplot, lbins[wh], sim.pow_mb[wh,0]-sim.err_mb, sim.pow_mb[*,0]+sim.err_mb

;wset, 1+wnum
plot, lbins[wh], sim.pow_rb[wh,0], title='mean pow_resid'
errplot, lbins[wh], sim.pow_rb[wh,0]-sim.err_rb, sim.pow_rb[*,0]+sim.err_rb

;wset, 2+wnum
plot, l[vec],pn[vec],/ylog,xr=[0,4000],title='FFT of single timestream'
oplot, l[vec],sim.pow_s[vec,0],color=!green

; h = histogram(sim.pte_rb, locations=xb,nbins=11)
; xb = xb+(xb[1]-xb[0])/2.
; plot, xb, h, psym=10,xr=[0,1],/xst

;wset, 3+wnum
plot, lbins, sim.pow_mb, title='Signal and Measured power',/ylog,xr=[0,4000]
oplot, lbins, sim.pow_sb, color=!green

!p.multi=0
END


;;;;;;;;;;;;;;
;
; Make plots
; INPUTS:
;   sim,        sim structure
;
;;;;;;;;;;;;;;
PRO plot_bias, sim0, sim1, sim2, signal, lbins, wnum=wnum
if n_elements(wnum) eq 0 then wnum=5

!p.multi=[0,2,2]
wset, 0+wnum

h0 = histogram(sim0.pte_rb, locations=xb0,nbins=11)
xb0 = xb0+(xb0[1]-xb0[0])/2.
h1 = histogram(sim1.pte_rb, locations=xb1,nbins=11)
h2 = histogram(sim2.pte_rb, locations=xb2,nbins=11)
mx = max([h0,h1,h2])
plot, xb0, h0, psym=10,xr=[0,1],/xst,yr=[0,mx],/yst,title='PTE distribution'
oplot, xb0,h1,psym=10,color=!red
oplot, xb0,h2,psym=10,color=!blue

plot, lbins,sim0.pow_rb[*,0],xr=[0,4000]
oplot, lbins,sim1.pow_rb[*,0],color=!red
oplot, lbins,sim2.pow_rb[*,0],color=!blue

plot, signal, sim2.coadd_sn[*,0]/signal, psym=3,yr=[0.5,1.5]
nsims = n_elements(sim0.pte_rb)
for i=0,nsims-1 do oplot,signal,sim2.coadd_sn[*,i]/signal,psym=3,color=!blue

plot, signal, sim2.coadd_sn[*,0]- signal, psym=3
for i=0,nsims-1 do oplot, signal,sim2.coadd_sn[*,i]-signal, psym=3,color=!blue

!p.multi=0
END



;;;;;;;;;;;;;;
;
; Make plots
; INPUTS:
;   sim,        sim structure
;
;;;;;;;;;;;;;;
PRO calc_bias
restore, 'sav_files/sig_bias_0309a.sav'

whg = where(mask_row gt 0.9) ; only use values where the mask is above 90%
whp = where(signal gt 0)
whn = where(signal lt 0)

s0 = {resid_p:fltarr(nsims*whp), resid_n:fltarr(nsims*whn),ms_p:fltarr(nsims), ms_n:fltarr(nsims), m_n:0., m_p:0., m:0.}
s1 = {ms_p:fltarr(nsims), ms_n:fltarr(nsims), m_n:0., m_p:0., m:0.}
s2 = {ms_p:fltarr(nsims), ms_n:fltarr(nsims), m_n:0., m_p:0., m:0.}

for i=0,nsims-1 do begin
    s0.ms_p[i] = mean(sim0.coadd_sn[whp,i]-signal[whp])
    s0.ms_n[i] = mean(sim0.coadd_sn[whn,i]-signal[whn])
    s1.ms_p[i] = mean(sim1.coadd_sn[whp,i]-signal[whp])
    s1.ms_n[i] = mean(sim1.coadd_sn[whn,i]-signal[whn])
    s2.ms_p[i] = mean(sim2.coadd_sn[whp,i]-signal[whp])
    s2.ms_n[i] = mean(sim2.coadd_sn[whn,i]-signal[whn])
endfor

s0.m_n = mean(s0.ms_n) & s0.m_p = mean(s0.ms_p)
s1.m_n = mean(s1.ms_n) & s1.m_p = mean(s1.ms_p)
s2.m_n = mean(s2.ms_n) & s2.m_p = mean(s2.ms_p)

s0.m = (s0.m_p - s0.m_n)/2.
s1.m = (s1.m_p - s1.m_n)/2.
s2.m = (s2.m_p - s2.m_n)/2.

stop

END




;;;;;;;;;;;;;;
;
; Signal bias sims
;
;;;;;;;;;;;;;;
PRO sim_sigbias, ncoadd=ncoadd, nsims=nsims, stopit=stopit, $
                 sigfac=sigfac, mk_sav=mk_sav, vebose=verbose

if n_elements(ncoadd) eq 0 then ncoadd = 1000;0 ; number of coadds per sim
if n_elements(nsims) eq 0  then nsims  = 100   ; number of sims

if not keyword_set(mk_sav) then print, 'Not saving any output files!'

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
    

; ------------------------
; Signal Timestream
; ------------------------

sim_map = get_sim_map(npix,/cmb_flatsky)
sig_idx = scramble(indgen(nsims))
signal0 = reform( sim_map[*, sig_idx ) * 1d6; units of uK

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
bl += 1.

; Mask the signal
signal = signal0 * mask_row
if n_elements(sigfac) ne 0 then begin
    print, ' --- Scaling the signal by a factor of '+strtrim(string(sigfac),2)+' ---'
    signal *= sigfac            ; no signal, for now
endif

; filter the timestream
signal = float(fft(fft(signal)*hpfilt,1))

; Calculate the signal power
pow_s    = abs(fft(signal))^(2.) * cor / (bl^2.)
pow_sb = fltarr(nbins)
for ii=0, nbins-1 do begin 
    wh = where(l ge dl*ii and l lt dl*(ii+1))
    pow_sb[ii] = mean(pow_s[wh])
endfor

for isim=0L, nsims-1 do begin
    sim0.pow_s[*,isim] = pow_s   & sim1.pow_s[*,isim] = pow_s   & sim2.pow_s[*,isim] = pow_s
    sim0.pow_sb[*,isim] = pow_sb & sim1.pow_sb[*,isim] = pow_sb & sim2.pow_sb[*,isim] = pow_sb
endfor

; ------------------------
; Noise timetsreams
; ------------------------
sigma = get_sigma_noise(ncoadd,sample_rate) ; currently just 500 uK

; sim structure
tmp = {coadd_sn:fltarr(npix[0], nsims), $ ; coadded signal+noise timestream
       coadd_n:fltarr(npix[0], nsims), $  ; coadded noise timestream
       pow_m :fltarr(npix[0], nsims), $ ; measured power
       pow_n :fltarr(npix[0], nsims), $ ; noise power
       pow_r :fltarr(npix[0], nsims), $ ; residual power
       pow_s :fltarr(npix[0], nsims), $ ; signal power
       pow_mb:fltarr(nbins, nsims), $ ; pow_m, bin-averaged
       pow_nb:fltarr(nbins, nsims), $ ; pow_n, bin-averaged
       pow_rb:fltarr(nbins, nsims), $ ; pow_r, bin-averaged
       pow_sb:fltarr(nbins, nsims), $ ; pow_s, bin-averaged
       pow_Nexp:fltarr(nsims), $ ; expected noise power
       err_mb:fltarr(nbins), $  ; measured power errors
       err_nb:fltarr(nbins), $  ; measured power errors
       err_rb:fltarr(nbins), $  ; residual power errors
       wts:fltarr(ncoadd, nsims), $ ; weights
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

    n_tstream   *= 0
    sn_tstream  *= 0

    ; Loop over "bolometer timestreams"
    for jj=0L, ncoadd-1 do begin
        ;seed = isim*ncoadd + jj
        undefine, seed
        noise_tmp = randomn(seed, npix[0]) * sigma[jj]
        noise_tmp *= mask_row   ; apply mask
        noise_tmp = float(fft(fft(noise_tmp)*hpfilt,1)) ; high-pass filter
        
        n_tstream[*,jj]  = noise_tmp
        sn_tstream[*,jj] = signal + noise_tmp
    endfor

    ; Get the weights
    sim0.wts[*,isim] = get_weights(0, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)
    sim1.wts[*,isim] = get_weights(1, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)
    sim2.wts[*,isim] = get_weights(2, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)

    ; calculate expected noise power
    for jj=0L, ncoadd-1 do begin
        sim0.pow_Nexp[isim] += (sigma[jj])^2. * sim0.wts[jj,isim]^2. * (1./npix[0]) * reso_rad^2.
        sim1.pow_Nexp[isim] += (sigma[jj])^2. * sim1.wts[jj,isim]^2. * (1./npix[0]) * reso_rad^2.
        sim2.pow_Nexp[isim] += (sigma[jj])^2. * sim2.wts[jj,isim]^2. * (1./npix[0]) * reso_rad^2.
    endfor

    ; coadd the timestreams
    for k=0, npix[0]-1 do begin
        sim0.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim0.wts[*,isim])
        sim1.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim1.wts[*,isim])
        sim2.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim2.wts[*,isim])

        sim0.coadd_n[k,isim] += total(n_tstream[k,*] * sim0.wts[*,isim])
        sim1.coadd_n[k,isim] += total(n_tstream[k,*] * sim1.wts[*,isim])
        sim2.coadd_n[k,isim] += total(n_tstream[k,*] * sim2.wts[*,isim])
    endfor

    ; Calculate the measured power, and residuals from expected power
    sim0.pow_m[*,isim] = abs(fft(sim0.coadd_sn[*,isim]))^(2.) * cor
    sim0.pow_n[*,isim] = abs(fft(sim0.coadd_n[*,isim]))^(2.) * cor
    ;sim0.pow_r[*,isim] = sim0.pow_m[*,isim] - (sim0.pow_s + sim0.pow_Nexp[isim])
    sim0.pow_r[*,isim] = sim0.pow_m[*,isim] - (sim0.pow_s[*,isim] + sim0.pow_n[*,isim])

    sim1.pow_m[*,isim] = abs(fft(sim1.coadd_sn[*,isim]))^(2.) * cor
    sim1.pow_n[*,isim] = abs(fft(sim1.coadd_n[*,isim]))^(2.) * cor
    sim1.pow_r[*,isim] = sim1.pow_m[*,isim] - (sim1.pow_s[*,isim] + sim1.pow_n[*,isim])

    sim2.pow_m[*,isim] = abs(fft(sim2.coadd_sn[*,isim]))^(2.) * cor
    sim2.pow_n[*,isim] = abs(fft(sim2.coadd_n[*,isim]))^(2.) * cor
    sim2.pow_r[*,isim] = sim2.pow_m[*,isim] - (sim2.pow_s[*,isim] + sim2.pow_n[*,isim])


    ; Bin the residuals
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
; Calculate the error in each bin
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

;print, "isim = ", format="(A,$)"
for jj=0L, nsims-1 do begin
    ;if (isim mod 10 eq 0) then print, isim, format='((x, I0),$)' 
    sim0.chisq_rb[jj] = total( (sim0.pow_rb[wh,jj]/sim0.err_rb[wh])^2.)
    sim0.pte_rb[jj]   = mpchitest(sim0.chisq_rb[jj], nwh)
    sim1.chisq_rb[jj] = total( (sim1.pow_rb[wh,jj]/sim1.err_rb[wh])^2.)
    sim1.pte_rb[jj]   = mpchitest(sim1.chisq_rb[jj], nwh)
    sim2.chisq_rb[jj] = total( (sim2.pow_rb[wh,jj]/sim2.err_rb[wh])^2.)
    sim2.pte_rb[jj]   = mpchitest(sim2.chisq_rb[jj], nwh)
endfor


; DEBUGGING
pn = abs(fft(noise_tmp))^(2.) * cor

; Plot stuff
plot_sim, sim2, lbins, l, vec, pn
plot_bias, sim0, sim1, sim2, signal, lbins

if keyword_set(mk_sav) then save, map_name, ncoadd, nsims, signal,dl,lbins,mask_row,mask_fac,cor,hpfilt,sigma,n_tstream,sn_tstream,sim0,sim1,sim2,$
  filename='sav_files/sig_bias_0309a.sav'
if keyword_set(stopit) then stop
END



