;;;
; NAME: script13_0314
;
; NOTES:
;  1) calculate nhits per pixel in individual maps
;  1) re-process sims
;;;



;;;;;;;;;;;;;;
;
; Plot stuff
;
;;;;;;;;;;;;;;
PRO s1_plot, save_figs=save_figs

;stub = '0309a_realistic' & yr=[-50,50]
;stub = '0309a_maxBias' & yr=[0.8,1.2]
;stub = '0309a_manyCoadd' & yr=[-9,11]

;stub = '0313_realistic' & yr=[-50,50]
;stub = '0313_maxBias' & yr=[0.8,1.2]
stub = '0313_manyCoadd' & yr=[-9,11]

sav_file = 'sav_files/sim_sigbias_map_'+stub+'.sav'
restore, sav_file

vec = indgen(npix[0]/2 - 1)

; In-band v.s. out-of-band power:
print, '*** '+stub+' ***'
print, "--- sim 0 ---"
print, 'inband : ', mean(s0.psd_inband), stddev(s0.psd_inband)
print, 'low    : ', mean(s0.psd_low), stddev(s0.psd_low)
print, 'high    : ', mean(s0.psd_high), stddev(s0.psd_high)
print, "--- sim 1 ---"
print, 'inband : ', mean(s1.psd_inband), stddev(s1.psd_inband)
print, 'low    : ', mean(s1.psd_low), stddev(s1.psd_low)
print, 'high    : ', mean(s1.psd_high), stddev(s1.psd_high)
print, "--- sim 2 ---"
print, 'inband : ', mean(s2.psd_inband), stddev(s2.psd_inband)
print, 'low    : ', mean(s2.psd_low), stddev(s2.psd_low)
print, 'high    : ', mean(s2.psd_high), stddev(s2.psd_high)


; PLOTS
wset, 0
!p.multi=[0,1,3]
chsize = 3

plot, lbins, s0.rpb,xr=[0,4000],charsize=chsize,title='Sim 0',xtitle='ell',ytitle='Power_resid [uK^2]'
errplot, lbins, s0.rpb-s0.rpb_err, s0.rpb+s0.rpb_err

plot, lbins, s1.rpb,xr=[0,4000],charsize=chsize,title='Sim 1',xtitle='ell',ytitle='Power_resid [uK^2]'
errplot, lbins, s1.rpb-s1.rpb_err, s1.rpb+s1.rpb_err

plot, lbins, s2.rpb,xr=[0,4000],charsize=chsize,title='Sim 2',xtitle='ell',ytitle='Power_resid [uK^2]'
errplot, lbins, s2.rpb-s2.rpb_err, s2.rpb+s2.rpb_err

!p.multi=0

wset, 5
!p.multi=[0,3,2]

plot, s0.psd_inband,charsize=chsize,title='Power in-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, s1.psd_inband, color=!red
oplot, s2.psd_inband, color=!blue

plot, s0.psd_low,charsize=chsize,title='Power low-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, s1.psd_low, color=!red
oplot, s2.psd_low, color=!blue

plot, s0.psd_high,charsize=chsize,title='Power high-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, s1.psd_high, color=!red
oplot, s2.psd_high, color=!blue

h0 = histogram(s0.psd_inband, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(s1.psd_inband, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(s2.psd_inband, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue


h0 = histogram(s0.psd_low, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in low-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(s1.psd_low, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(s2.psd_low, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue

h0 = histogram(s0.psd_high, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in high-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(s1.psd_high, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(s2.psd_high, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue

!p.multi=0

if keyword_set(save_figs) then begin
    print, 'saving figures...'
    fdir = '/home/kstory/public_html/notebook/spt_lps12/'

    wset, 0
    err=tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+stub+'_3tsplot')

    wset, 5
    err=tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+stub+'_hzbandBias')
endif

stop
END




;;;;;;;;;;;;;;
;
; Make plots for wiki post
;
;;;;;;;;;;;;;;
PRO make_plots,save_figs=save_figs

;stub = '0309a_realistic' & yr=[-50,50]
;stub = '0309a_maxBias' & yr=[0.8,1.2]
;stub = '0309a_manyCoadd' & yr=[-9,11]

stub = '0313_realistic' & yr=[-50,50]
;stub = '0313_maxBias' & yr=[0.8,1.2]
;stub = '0313_manyCoadd' & yr=[-9,11]

sav_file = 'sav_files/sim_sigbias_map_'+stub+'.sav'
restore, sav_file

vec = indgen(npix[0]/2 - 1)

print, '--- make plot ---'
if n_elements(pn) eq 0 then begin
    noise_tmp = randomn(seed, npix[0]) * sigma[0]
    noise_tmp *= mask_row       ; apply mask
    noise_tmp = float(fft(fft(noise_tmp)*hpfilt,1)) ; high-pass filter
    pn = abs(fft(noise_tmp))^(2.) * cor
endif

wset, 0
plot_sim, sim2, lbins, l, vec, stub, pn=pn

wset, 5
plot_bias, sim0, sim1, sim2, s0,s1,s2,lbins, whmask,stub, yr=yr


help, s0,s1,s2,/st

if keyword_set(save_figs) then begin
    print, 'saving figures...'
    fdir = '/home/kstory/public_html/notebook/spt_lps12/'

    wset, 0
    err=tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+stub+'_tsplot')

    wset, 5
    err=tvread(/png,/nodialog,filename=fdir+'sim_sigbias_'+stub+'_biasplot')
endif

stop

END




;;;;;;;;;;;;;;
;
; Restuff power bias information
;
;;;;;;;;;;;;;;
PRO restuff_sims_2, input, output
restore, input
nbins = n_elements(sim0.err_rb)

st0 = {rpb:fltarr(nbins), rpb_err:fltarr(nbins), psd_inband:fltarr(nsims), psd_low:fltarr(nsims), psd_high:fltarr(nsims)}
st1 = {rpb:fltarr(nbins), rpb_err:fltarr(nbins), psd_inband:fltarr(nsims), psd_low:fltarr(nsims), psd_high:fltarr(nsims)}
st2 = {rpb:fltarr(nbins), rpb_err:fltarr(nbins), psd_inband:fltarr(nsims), psd_low:fltarr(nsims), psd_high:fltarr(nsims)}

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

;integrate in-band and out-of-band power in the residuals
for j=0,nsims-1 do begin
    st0.psd_inband[j] = total(sim0.pow_rb[wh_1to3,j])*dl / nsims / nHz_inband
    st0.psd_low[j]    = total(sim0.pow_rb[wh_low,j])*dl  / nsims / nHz_low
    st0.psd_high[j]   = total(sim0.pow_rb[wh_high,j])*dl / nsims / nHz_high

    st1.psd_inband[j] = total(sim1.pow_rb[wh_1to3,j])*dl / nsims / nHz_inband
    st1.psd_low[j]    = total(sim1.pow_rb[wh_low,j])*dl  / nsims / nHz_low
    st1.psd_high[j]   = total(sim1.pow_rb[wh_high,j])*dl / nsims / nHz_high

    st2.psd_inband[j] = total(sim2.pow_rb[wh_1to3,j])*dl / nsims / nHz_inband
    st2.psd_low[j]    = total(sim2.pow_rb[wh_low,j])*dl  / nsims / nHz_low
    st2.psd_high[j]   = total(sim2.pow_rb[wh_high,j])*dl / nsims / nHz_high

endfor

for i=0, nbins-1 do begin
    st0.rpb[i] = mean(sim0.pow_rb[i,0:nsims-1])
    st1.rpb[i] = mean(sim1.pow_rb[i,0:nsims-1])
    st2.rpb[i] = mean(sim2.pow_rb[i,0:nsims-1])
endfor

st0.rpb_err = sim0.err_rb / sqrt(nsims)
st1.rpb_err = sim1.err_rb / sqrt(nsims)
st2.rpb_err = sim2.err_rb / sqrt(nsims)

s0 = st0
s1 = st1
s2 = st2

save, ncoadd,nsims,npix,nbins,dl,l,wh_1to3,hz,lbins,mask_row,whmask,mask_fac,cor,hpfilt,sigfac,sigma,n_tstream,sn_tstream,pn,$
  sim0,sim1,sim2,s0,s1,s2,$
  filename=output
END


;;;;;;;;;;;;;;
;
; play around
;
;;;;;;;;;;;;;;
PRO s1

;stub = '0309a_realistic' & yr=[-50,50]
;stub = '0309a_maxBias' & yr=[0.8,1.2]
;stub = '0309a_manyCoadd' & yr=[-9,11]

;stub = '0313_realistic' & yr=[-50,50]
;stub = '0313_maxBias' & yr=[0.8,1.2]
stub = '0313_manyCoadd' & yr=[-9,11]

sav_file = 'sav_files/sim_sigbias_map_'+stub+'.sav'
restore, sav_file

vec = indgen(npix[0]/2 - 1)

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

nbins = n_elements(sim0.err_rb)
ns = nsims ;100
st0 = {rpb:fltarr(nbins), rpb_err:fltarr(nbins), psd_inband:fltarr(ns), psd_low:fltarr(ns), psd_high:fltarr(ns)}
st1 = {rpb:fltarr(nbins), rpb_err:fltarr(nbins), psd_inband:fltarr(ns), psd_low:fltarr(ns), psd_high:fltarr(ns)}
st2 = {rpb:fltarr(nbins), rpb_err:fltarr(nbins), psd_inband:fltarr(ns), psd_low:fltarr(ns), psd_high:fltarr(ns)}

;integrate in-band and out-of-band power in the residuals
for j=0,ns-1 do begin
    st0.psd_inband[j] = total(sim0.pow_rb[wh_1to3,j])*dl / ns / nHz_inband
    st0.psd_low[j]    = total(sim0.pow_rb[wh_low,j])*dl  / ns / nHz_low
    st0.psd_high[j]   = total(sim0.pow_rb[wh_high,j])*dl / ns / nHz_high

    st1.psd_inband[j] = total(sim1.pow_rb[wh_1to3,j])*dl / ns / nHz_inband
    st1.psd_low[j]    = total(sim1.pow_rb[wh_low,j])*dl  / ns / nHz_low
    st1.psd_high[j]   = total(sim1.pow_rb[wh_high,j])*dl / ns / nHz_high

    st2.psd_inband[j] = total(sim2.pow_rb[wh_1to3,j])*dl / ns / nHz_inband
    st2.psd_low[j]    = total(sim2.pow_rb[wh_low,j])*dl  / ns / nHz_low
    st2.psd_high[j]   = total(sim2.pow_rb[wh_high,j])*dl / ns / nHz_high

endfor

for i=0, nbins-1 do begin
    st0.rpb[i] = mean(sim0.pow_rb[i,0:ns-1])
    st1.rpb[i] = mean(sim1.pow_rb[i,0:ns-1])
    st2.rpb[i] = mean(sim2.pow_rb[i,0:ns-1])
endfor

st0.rpb_err = sim0.err_rb / sqrt(ns)
st1.rpb_err = sim1.err_rb / sqrt(ns)
st2.rpb_err = sim2.err_rb / sqrt(ns)


; In-band v.s. out-of-band power:
print, "--- sim 0 ---"
print, 'inband : ', mean(st0.psd_inband), stddev(st0.psd_inband)
print, 'low    : ', mean(st0.psd_low), stddev(st0.psd_low)
print, 'high    : ', mean(st0.psd_high), stddev(st0.psd_high)
print, "--- sim 1 ---"
print, 'inband : ', mean(st1.psd_inband), stddev(st1.psd_inband)
print, 'low    : ', mean(st1.psd_low), stddev(st1.psd_low)
print, 'high    : ', mean(st1.psd_high), stddev(st1.psd_high)
print, "--- sim 2 ---"
print, 'inband : ', mean(st2.psd_inband), stddev(st2.psd_inband)
print, 'low    : ', mean(st2.psd_low), stddev(st2.psd_low)
print, 'high    : ', mean(st2.psd_high), stddev(st2.psd_high)


; PLOTS
wset, 0
!p.multi=[0,1,3]
chsize = 3

plot, lbins, st0.rpb,xr=[0,4000],charsize=chsize,title='Sim 0',xtitle='ell',ytitle='Power_resid [uK^2]'
errplot, lbins, st0.rpb-st0.rpb_err, st0.rpb+st0.rpb_err

plot, lbins, st1.rpb,xr=[0,4000],charsize=chsize,title='Sim 1',xtitle='ell',ytitle='Power_resid [uK^2]'
errplot, lbins, st1.rpb-st1.rpb_err, st1.rpb+st1.rpb_err

plot, lbins, st2.rpb,xr=[0,4000],charsize=chsize,title='Sim 2',xtitle='ell',ytitle='Power_resid [uK^2]'
errplot, lbins, st2.rpb-st2.rpb_err, st2.rpb+st2.rpb_err

!p.multi=0

wset, 5
!p.multi=[0,3,2]

plot, st0.psd_inband,charsize=chsize,title='Power in-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, st1.psd_inband, color=!red
oplot, st2.psd_inband, color=!blue

plot, st0.psd_low,charsize=chsize,title='Power low-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, st1.psd_low, color=!red
oplot, st2.psd_low, color=!blue

plot, st0.psd_high,charsize=chsize,title='Power high-band',xtitle='sim number',ytitle='PSD [uK^2/Hz]'
oplot, st1.psd_high, color=!red
oplot, st2.psd_high, color=!blue

h0 = histogram(st0.psd_inband, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(st1.psd_inband, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(st2.psd_inband, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue


h0 = histogram(st0.psd_low, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in low-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(st1.psd_low, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(st2.psd_low, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue

h0 = histogram(st0.psd_high, locations=xb0,nbins=50)
xb0 = xb0+(xb0[1]-xb0[0])/2.
plot, xb0, h0, psym=10,charsize=chsize,title='Power in high-band',xtitle='PSD [uK^2/Hz]'
h1 = histogram(st1.psd_high, locations=xb1,nbins=50)
xb1 = xb1+(xb1[1]-xb1[0])/2.
oplot, xb1, h1, psym=10,color=!red
h2 = histogram(st2.psd_high, locations=xb2,nbins=50)
xb2 = xb2+(xb2[1]-xb2[0])/2.
oplot, xb2, h2, psym=10,color=!blue

!p.multi=0



stop
END







PRO read_my_data
;restore, '/path_to_file/myfile.sav'
x = findgen(100)
y = findgen(100) * 5

npts = n_elements(x)
get_lun, lun1
openw,lun1,'output.txt'
for i=0, npts-1 do begin
   printf,lun1,x[i],y[i]
endfor
close, lun1
free_lun,lun1
END
