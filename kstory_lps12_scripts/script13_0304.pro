;;;
; NAME: script13_0304
;
; NOTES:
;  1) first attempt at signal-bias sims
;;;


;;;;;;;;;;;;;;
;
; Signal bias sims
;   wtype: 0 = no weight
;          1 = psd in noise timestreams
;          2 = psd in signal+noise timestreams
;
;;;;;;;;;;;;;;
PRO sig_bias, wtype, stopit=stopit, nosig=nosig

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

; sim_dir     = '/home/kstory/lps12/sims/' ; soft-linked to /data23/
; bdir = '/data25/hou/lps12/sim_mapmaking/output/'
; spawn, 'ls ' + bdir+fname+'/map_*lmax8000*.bin', list
; map_name = list[10] ; arbitrarily take 10'th map
map_name = '/data25/hou/lps12/sim_mapmaking/output/ra4h10dec-50/map_150_lmax8000_20100212_112516.bin'


ncoadd = 20;0 ; number of coadds per sim
nsims  = 100 ; number of sims


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
dl = 200
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
if keyword_set(nosig) then signal *= 0. ; no signal, for now

; Calibrate simulated signal timestream ??

; ------------------------
; Noise timetsreams
; ------------------------

bdir = '/home/kstory/lps12/lowellfits/'+fst.name+'/'
fname = bdir+'field_scan_stub_150_20100212_112516.fits'
data = krf(fname)
bolo_id = 100 ; random choice for now
wn = data.observation.bolo_wnoise[bolo_id] ; white-noise level?
sigma = wn * 1d6 * 5d5 ; Should only be one factor of 1d6

n_tstream = fltarr(npix[0], ncoadd) ; tmp array for coadding timestreams
sn_tstream = fltarr(npix[0], ncoadd) ; tmp array for coadding timestreams
noise_coadd = fltarr(npix[0])
pow_m  = fltarr(npix[0], nsims)
pow_r  = fltarr(npix[0], nsims)
pow_mb = fltarr(nbins, nsims) ; pow_meas, bin-averaged
pow_rb = fltarr(nbins, nsims) ; pow_resid, bin-averaged


; ------------------------
; Coadd noise and calculate power
; ------------------------
pow_Nexp = (sigma/sqrt(ncoadd))^2. * (1./npix[0]) * reso_rad^2. 
pow_s    = abs(fft(signal))^(2.) * cor
pow_sb = fltarr(nbins)
for ii=0, nbins-1 do begin 
    wh = where(l ge dl*ii and l lt dl*(ii+1))
    pow_sb[ii] = mean(pow_s[wh])
endfor

print, "isim = ", format="(A,$)"
for isim=0, nsims-1 do begin
    if (isim mod 10 eq 0) then print, isim, format='((x, I0),$)' 
    n_tstream   *= 0
    sn_tstream  *= 0
    noise_coadd *= 0

    ; Loop over "bolometer timestreams"
    for ii=0, ncoadd-1 do begin
        ;seed = isim*ncoadd + ii
        undefine, seed
        noise_tmp = randomn(seed, npix[0]) * sigma
        noise_tmp *= mask_row   ; apply mask
        
        n_tstream[*,ii]  = noise_tmp
        sn_tstream[*,ii] = signal + noise_tmp
    endfor

    ; coadd timestreams with weight scheme
    case wtype of 
        0 begin: ; no weights
            noise_coadd = total(sn_tstream) / ncoadd
        end
        1 begin:
            wts = fltarr(nbins)
            psd_1to3 = fltarr(nbins)
            wh_1to3 = where(hz ge 1 and hz lt 3)
            for ii=0, ncoadd-1 do begin
                pp = abs(fft(n_tstream))^(2.) * cor
                psd_1to3[ii] = mean(pp[wh_1to3])
                wts[ii] = 1/psd_1to3[ii]
            endfor

    pow_m[*,isim] = abs(fft(noise_coadd))^(2.) * cor
    pow_r[*,isim] = pow_m[*,isim] - (pow_s + pow_Nexp)

    ; Bin the residuals
    for ii=0, nbins-1 do begin
        wh = where(l ge dl*ii and l lt dl*(ii+1))
        pow_mb[ii,isim] = mean(pow_m[wh,isim])
        pow_rb[ii,isim] = mean(pow_r[wh,isim])
    end
endfor

; ------------------------
; Calculate the error in each bin
; ------------------------
print, '--- calculate error ---'
err_mb = fltarr(nbins) 
err_rb = fltarr(nbins) 
for ii=0, nbins-1 do begin
    err_mb[ii] = stddev(pow_mb[ii,*])
    err_rb[ii] = stddev(pow_rb[ii,*])
endfor

print, '--- calculate chisq ---'
chisq_rb = fltarr(nsims)
pte_rb   = fltarr(nsims)

print, "isim = ", format="(A,$)"
for jj=0, nsims-1 do begin
    if (isim mod 10 eq 0) then print, isim, format='((x, I0),$)' 
    chisq_rb[jj] = total( (pow_rb[*,jj]/err_rb)^2.)
    pte_rb[jj]   = mpchitest(chisq_rb[jj], nbins)
endfor



; print, " --- Chisq --- "
; print, mpchitest(chisq_rb, n_elements(pow_rb)), ", ", mpchitest(chisq_rb, n_elements(pow_rb),/sigma), " sigma"

wset, 0
plot, lbins, pow_mb[*,0], /ylog, title='mean pow_meas'
errplot, lbins, pow_mb[*,0]-err_mb, pow_mb[*,0]+err_mb

wset, 1
plot, lbins, pow_rb[*,0], title='mean pow_resid'
errplot, lbins, pow_rb[*,0]-err_rb, pow_rb[*,0]+err_rb

wset, 2
h = histogram(pte_rb, locations=xb,nbins=11)
xb = xb+(xb[1]-xb[0])/2.
plot, xb, h, psym=10,xr=[0,1],/xst

wset, 3
plot, lbins, pow_mb, title='Signal and Measured power'
oplot, lbins, pow_sb, color=!green

; h_r = histogram(pow_resid, locations=xbins_r,nbins=50)
; xbins_r = xbins_r + (xbins_r[1]-xbins_r[0])/2.
; h_m = histogram(pow_meas, locations=xbins_m,nbins=50)
; xbins_m = xbins_m + (xbins_m[1]-xbins_m[0])/2.


; wset, 0
; plot, noise_coadd, title='Coadded Noise timestream'

; wset, 1
; plot, noise_tmp, title='Single Noise Timestream'
; oplot, noise_coadd, color=!green
; oplot, signal, color=!red

; wset, 3
; plot, l[vec], pow_meas[vec], /ylog, title='Measured Power'

; wset, 4
; plot, xbins_r, h_r, psym=10, title='Power Residuals'

; wset, 5
; plot, xbins_m, h_m, psym=10, title='Power Measured'



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
