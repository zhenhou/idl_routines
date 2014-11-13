;;;
; NAME: script13_0315
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
PRO s1

;stub = '0309a_realistic' & yr=[-50,50]
;stub = '0309a_maxBias' & yr=[0.8,1.2]
;stub = '0309a_manyCoadd' & yr=[-9,11]

stub = '0313_realistic' & yr=[-50,50]
;stub = '0313_maxBias' & yr=[0.8,1.2]
;stub = '0313_manyCoadd' & yr=[-9,11]

sav_file = 'sav_files/sim_sigbias_map_'+stub+'.sav'
restore, sav_file

vec = indgen(npix[0]/2 - 1)

pow_s = fltarr(npix[0])
; average power over all sims
for i=0, npix[0]-1 do pow_s[i] = mean(sim0.pow_s[i,*])

pow_sb = fltarr(nbins)
for i=0, nbins-1 do pow_sb[i] = mean(sim0.pow_sb[i,*])

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


tmp = {m:0., ms:fltarr(nsims), $
       rpb:fltarr(nbins), rpb_err:fltarr(nbins), $
       r_psd_inband:fltarr(nsims), r_psd_low:fltarr(nsims), r_psd_high:fltarr(nsims), $
       n_psd_inband:fltarr(nsims), n_psd_low:fltarr(nsims), n_psd_high:fltarr(nsims), $
       s_psd_inband:fltarr(nsims), s_psd_low:fltarr(nsims), s_psd_high:fltarr(nsims)}

st0 = replicate(tmp,1)
st1 = replicate(tmp,1)
st2 = replicate(tmp,1)

;integrate in-band and out-of-band power in the residuals
for j=0,nsims-1 do begin
    st0.r_psd_inband[j] = total(sim0.pow_rb[wh_1to3,j]) / nsims / nHz_inband
    st0.r_psd_low[j]    = total(sim0.pow_rb[wh_low,j])  / nsims / nHz_low
    st0.r_psd_high[j]   = total(sim0.pow_rb[wh_high,j]) / nsims / nHz_high

    st0.s_psd_inband[j] = total(sim0.pow_sb[wh_1to3,j]) / nsims / nHz_inband
    st0.s_psd_low[j]    = total(sim0.pow_sb[wh_low,j])  / nsims / nHz_low
    st0.s_psd_high[j]   = total(sim0.pow_sb[wh_high,j]) / nsims / nHz_high

    st0.n_psd_inband[j] = total(sim0.pow_nb[wh_1to3,j]) / nsims / nHz_inband
    st0.n_psd_low[j]    = total(sim0.pow_nb[wh_low,j])  / nsims / nHz_low
    st0.n_psd_high[j]   = total(sim0.pow_nb[wh_high,j]) / nsims / nHz_high

    st1.r_psd_inband[j] = total(sim1.pow_rb[wh_1to3,j]) / nsims / nHz_inband
    st1.r_psd_low[j]    = total(sim1.pow_rb[wh_low,j])  / nsims / nHz_low
    st1.r_psd_high[j]   = total(sim1.pow_rb[wh_high,j]) / nsims / nHz_high

    st1.s_psd_inband[j] = total(sim1.pow_sb[wh_1to3,j]) / nsims / nHz_inband
    st1.s_psd_low[j]    = total(sim1.pow_sb[wh_low,j])  / nsims / nHz_low
    st1.s_psd_high[j]   = total(sim1.pow_sb[wh_high,j]) / nsims / nHz_high

    st1.n_psd_inband[j] = total(sim1.pow_nb[wh_1to3,j]) / nsims / nHz_inband
    st1.n_psd_low[j]    = total(sim1.pow_nb[wh_low,j])  / nsims / nHz_low
    st1.n_psd_high[j]   = total(sim1.pow_nb[wh_high,j]) / nsims / nHz_high

    st2.r_psd_inband[j] = total(sim2.pow_rb[wh_1to3,j]) / nsims / nHz_inband
    st2.r_psd_low[j]    = total(sim2.pow_rb[wh_low,j])  / nsims / nHz_low
    st2.r_psd_high[j]   = total(sim2.pow_rb[wh_high,j]) / nsims / nHz_high

    st2.s_psd_inband[j] = total(sim2.pow_sb[wh_1to3,j]) / nsims / nHz_inband
    st2.s_psd_low[j]    = total(sim2.pow_sb[wh_low,j])  / nsims / nHz_low
    st2.s_psd_high[j]   = total(sim2.pow_sb[wh_high,j]) / nsims / nHz_high

    st2.n_psd_inband[j] = total(sim2.pow_nb[wh_1to3,j]) / nsims / nHz_inband
    st2.n_psd_low[j]    = total(sim2.pow_nb[wh_low,j])  / nsims / nHz_low
    st2.n_psd_high[j]   = total(sim2.pow_nb[wh_high,j]) / nsims / nHz_high
endfor

stop


END
