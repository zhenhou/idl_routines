;;;
; NAME: script13_0306
;
; NOTES:
;  1) first attempt at signal-bias sims
;;;


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
if keyword_set(narrow_band) then wh_1to3 = wh_1to3[0:1] 

case wtype of 
    0: begin                    ; no weights
        wts += 1.
    end
    
    1: begin                    ; noise
;        wts = findgen(ncoadd)+1.
        for jj=0, ncoadd-1 do begin
            pp = abs(fft(n_tstream[*,jj]))^(2.) * cor
            psd_1to3[jj] = mean(pp[wh_1to3])
            wts[jj] = 1/psd_1to3[jj]
        endfor
    end

    2: begin                    ; signal + noise
        for jj=0, ncoadd-1 do begin
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
; Util function: Get the name of the simulated map we want to use
;
;;;;;;;;;;;;;;
FUNCTION get_map_name
; sim_dir     = '/home/kstory/lps12/sims/' ; soft-linked to /data23/
; bdir = '/data25/hou/lps12/sim_mapmaking/output/'
; spawn, 'ls ' + bdir+fname+'/map_*lmax8000*.bin', list
; map_name = list[10] ; arbitrarily take 10'th map
map_name = '/data25/hou/lps12/sim_mapmaking/output/ra4h10dec-50/map_150_lmax8000_20100212_112516.bin'
RETURN, map_name
END



;;;;;;;;;;;;;;
;
; Make plots
; INPUTS:
;   sim,        sim structure
;
;;;;;;;;;;;;;;
PRO plot_sim, sim, lbins, wnum=wnum
if n_elements(wnum) eq 0 then wnum=0

wh = where(lbins lt 4000)

wset, 0+wnum
plot, lbins[wh], sim.pow_mb[wh,0], /ylog, title='mean pow_meas'
errplot, lbins[wh], sim.pow_mb[wh,0]-sim.err_mb, sim.pow_mb[*,0]+sim.err_mb

wset, 1+wnum
plot, lbins[wh], sim.pow_rb[wh,0], title='mean pow_resid'
errplot, lbins[wh], sim.pow_rb[wh,0]-sim.err_rb, sim.pow_rb[*,0]+sim.err_rb

wset, 2+wnum
h = histogram(sim.pte_rb, locations=xb,nbins=11)
xb = xb+(xb[1]-xb[0])/2.
plot, xb, h, psym=10,xr=[0,1],/xst

wset, 3+wnum
plot, lbins, sim.pow_mb, title='Signal and Measured power'
oplot, lbins, sim.pow_sb, color=!green

END

;;;;;;;;;;;;;;
;
; Make plots
; INPUTS:
;   sim,        sim structure
;
;;;;;;;;;;;;;;
PRO plot_bias, sim0, sim1, sim2, lbins, wnum=wnum
if n_elements(wnum) eq 0 then wnum=5

wset, 0+wnum
h0 = histogram(sim0.pte_rb, locations=xb0,nbins=11)
xb0 = xb0+(xb0[1]-xb0[0])/2.
h1 = histogram(sim1.pte_rb, locations=xb1,nbins=11)
h2 = histogram(sim2.pte_rb, locations=xb2,nbins=11)
mx = max([h0,h1,h2])
plot, xb0, h0, psym=10,xr=[0,1],/xst,yr=[0,mx],/yst,title='PTE distribution'
oplot, xb0,h1,psym=10,color=!red
oplot, xb0,h2,psym=10,color=!blue

wset, 1+wnum
plot, lbins,sim0.pow_rb[*,0],xr=[0,4000]
oplot, lbins,sim1.pow_rb[*,0],color=!red
oplot, lbins,sim2.pow_rb[*,0],color=!blue

END



;;;;;;;;;;;;;;
;
; Signal bias sims
;   wtype: 0 = no weight
;          1 = psd in noise timestreams
;          2 = psd in signal+noise timestreams
;
;;;;;;;;;;;;;;
PRO sig_bias, stopit=stopit, nosig=nosig

field_idx = 6  ; ra4h10dec-50
scan_speed = 0.417812 ; dps on-sky

; ------------------------
; Setup
; ------------------------
reso = 1.0
reso_rad = reso/60.*!dtor
f = lps12_fieldstruct()
fst = f[field_idx]
fname = fst.name

map_name = get_map_name()

ncoadd = 10000 ; number of coadds per sim
nsims  = 100                    ; number of sims


; ------------------------
; Signal Timestream
; ------------------------

; Read in Simulated Map
sim_map_arr = read_bin_map(map_name)
sim_map = reform(sim_map_arr[*,*,0])

ss = size(sim_map)
npix = [ss[1], ss[2]]
l = make_fft_grid(reso/60.*!dtor/2./!pi,npix[0],fx=lx)
hz = make_fft_grid(reso/60./scan_speed,npix[0])
vec = indgen(npix[0]/2 - 1)

; For binning the residuals
dl = 50
max_l = max(l[vec])
nbins = floor(max_l / dl)
lbins = findgen(nbins) * dl + dl/2.
    
; get "Signal" timestream
signal = reform( sim_map[*, npix[1]/2] )

mask = get_lps12_mask(field_idx)
mask_row = mask[*,npix[1]/2]
mask_fac = mean(mask_row^2.)
cor = (1/mask_fac) * reso_rad^2. ; correction factor for PS
signal *= mask_row
if keyword_set(nosig) then begin
    print, ' --- Zeroing the signal ---'
    signal *= 0.                ; no signal, for now
endif

; Calibrate simulated signal timestream ??

; ------------------------
; Noise timetsreams
; ------------------------

bdir = '/home/kstory/lps12/lowellfits/'+fst.name+'/'
fname = bdir+'field_scan_stub_150_20100212_112516.fits'
data = krf(fname)
whgood = where(data.observation.bolo_flags eq 0)
bolo_id = 100 ; random choice for now
wn = data.observation.bolo_wnoise[bolo_id] ; white-noise level?
sigma = fltarr(ncoadd) + 1*(wn * 1d6 * 8d5) ; Should only be one factor of 1d6

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


; Calculate the signal power
pow_s    = abs(fft(signal))^(2.) * cor
pow_sb = fltarr(nbins)
for ii=0, nbins-1 do begin 
    wh = where(l ge dl*ii and l lt dl*(ii+1))
    pow_sb[ii] = mean(pow_s[wh])
endfor

sim0.pow_s = pow_s & sim1.pow_s = pow_s & sim2.pow_s = pow_s
sim0.pow_sb = pow_sb & sim1.pow_sb = pow_sb & sim2.pow_sb = pow_sb

; ------------------------
; Coadd noise and calculate power
; ------------------------

print, "isim = ", format="(A,$)"
for isim=0, nsims-1 do begin
    if (isim mod 10 eq 0) then print, isim, format='((x, I0),$)' 
    n_tstream   *= 0
    sn_tstream  *= 0

    ; Loop over "bolometer timestreams"
    for jj=0, ncoadd-1 do begin
        ;seed = isim*ncoadd + jj
        undefine, seed
        noise_tmp = randomn(seed, npix[0]) * sigma[jj]
        noise_tmp *= mask_row   ; apply mask
        
        n_tstream[*,jj]  = noise_tmp
        sn_tstream[*,jj] = signal + noise_tmp
    endfor

    ; Get the weights
    sim0.wts[*,isim] = get_weights(0, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)
    sim1.wts[*,isim] = get_weights(1, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)
    sim2.wts[*,isim] = get_weights(2, nbins, ncoadd, cor, hz, n_tstream, sn_tstream)

    ; calculate expected noise power
    for jj=0, ncoadd-1 do begin
        ;sim0.pow_Nexp[isim] += (sigma[jj])^2./float(ncoadd) * wts0[jj] * (1./npix[0]) * reso_rad^2.
        sim0.pow_Nexp[isim] += (sigma[jj])^2. * sim0.wts[jj,isim]^2. * (1./npix[0]) * reso_rad^2.
        sim1.pow_Nexp[isim] += (sigma[jj])^2. * sim1.wts[jj,isim]^2. * (1./npix[0]) * reso_rad^2.
        sim2.pow_Nexp[isim] += (sigma[jj])^2. * sim2.wts[jj,isim]^2. * (1./npix[0]) * reso_rad^2.
    endfor

    ; coadd the timestreams
    for k=0, npix[0]-1 do begin
        sim0.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim0.wts[*,isim])
        sim1.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim1.wts[*,isim])
        sim2.coadd_sn[k,isim] += total(sn_tstream[k,*] * sim2.wts[*,isim])

        sim0.coadd_n[k,isim] += total(n_tstream[k,*] * sim0.wts.[*,isim])
        sim1.coadd_n[k,isim] += total(n_tstream[k,*] * sim1.wts.[*,isim])
        sim2.coadd_n[k,isim] += total(n_tstream[k,*] * sim2.wts.[*,isim])
    endfor

    ; Calculate the measured power, and residuals from expected power
    sim0.pow_m[*,isim] = abs(fft(sim0.coadd_sn[*,isim]))^(2.) * cor
    sim0.pow_n[*,isim] = abs(fft(sim0.coadd_n[*,isim]))^(2.) * cor
    ;sim0.pow_r[*,isim] = sim0.pow_m[*,isim] - (sim0.pow_s + sim0.pow_Nexp[isim])
    sim0.pow_r[*,isim] = sim0.pow_m[*,isim] - (sim0.pow_s + sim0.pow_n[*,isim])

    sim1.pow_m[*,isim] = abs(fft(sim1.coadd_sn[*,isim]))^(2.) * cor
    sim1.pow_n[*,isim] = abs(fft(sim1.coadd_n[*,isim]))^(2.) * cor
    sim1.pow_r[*,isim] = sim1.pow_m[*,isim] - (sim1.pow_s + sim1.pow_n[*,isim])

    sim2.pow_m[*,isim] = abs(fft(sim2.coadd_sn[*,isim]))^(2.) * cor
    sim2.pow_n[*,isim] = abs(fft(sim2.coadd_n[*,isim]))^(2.) * cor
    sim2.pow_r[*,isim] = sim2.pow_m[*,isim] - (sim2.pow_s + sim2.pow_n[*,isim])


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
wh = where(lbins ge 500and lbins le 3000, nwh)

;print, "isim = ", format="(A,$)"
for jj=0, nsims-1 do begin
    ;if (isim mod 10 eq 0) then print, isim, format='((x, I0),$)' 
    sim0.chisq_rb[jj] = total( (sim0.pow_rb[wh,jj]/sim0.err_rb[wh])^2.)
    sim0.pte_rb[jj]   = mpchitest(sim0.chisq_rb[jj], nwh)
    sim1.chisq_rb[jj] = total( (sim1.pow_rb[wh,jj]/sim1.err_rb[wh])^2.)
    sim1.pte_rb[jj]   = mpchitest(sim1.chisq_rb[jj], nwh)
    sim2.chisq_rb[jj] = total( (sim2.pow_rb[wh,jj]/sim2.err_rb[wh])^2.)
    sim2.pte_rb[jj]   = mpchitest(sim2.chisq_rb[jj], nwh)
endfor


; Plot stuff
plot_sim, sim2, lbins
plot_bias, sim0, sim1, sim2, lbins

if keyword_set(stopit) then stop
END



;;;;;;;;;;;;;;
;
; Toy test program
;
;;;;;;;;;;;;;;
PRO sc1
sigma = .18209569*1d6
nsims = 1000
ncoadd = 1;00
r = fltarr(nsims)
npts = 2^14.
reso = 1.
l = make_fft_grid(reso/60.*!dtor/2./!pi,npts,fx=lx)
vec=indgen(npts/2-1)
coadd = fltarr(npts)

; make mask
mask = fltarr(npts) + 1
cut = 4000
mask[0:2*cut-1] = 0 & mask = shift(mask, -cut)
mask_fac = mean(mask^2.)
cor = (1/mask_fac) ; *reso_rad^2

; Bin the residuals
dl = 500
max_l = max(l[vec])
nbins = floor(max_l / 500)
lbins = indgen(nbins) * 500 + 250

pow_mb = fltarr(nbins, nsims) ; pow_meas, bin-averaged
pow_rb = fltarr(nbins, nsims) ; pow_resid, bin-averaged

pow_Nexp = (sigma/sqrt(ncoadd))^2. * (1/npts)

for ii=0, nsims-1 do begin
    if (ii mod 20 eq 0) then print, ii, format='((x, I0),$)' 

    ; make coadded vector
    coadd *= 0.
    for jj=0, ncoadd-1 do begin
        x = randomn(seed, npts) * sigma * mask
        coadd += x/ncoadd
    endfor

    px = abs(fft(coadd))^2. * cor
    r[ii] = mean(px) / pow_Nexp

    ; Bin the residuals
    for ii=0, nbins-1 do begin
        wh = where(l ge dl*ii and l lt dl*(ii+1))
        pow_mb[ii,isim] = mean(pow_meas[wh])
        pow_rb[ii,isim] = mean(pow_resid[wh])
    end
endfor

h = histogram(r, locations=xvec, nbins=40)
xmid = xvec + (xvec[1] - xvec[0])/2.

wset, 0
plot, xmid, h, psym=10

wset, 1
plot, (abs(fft(x))^2. * cor)[vec],/ylog
oplot, px[vec],color=!red

stop
END
