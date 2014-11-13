;;;
; NAME: script13_0321
;
; NOTES:
;  1) check power levels
;;;


PRO re_process, tag=tag,stopit=stopit,save_file=save_file
;tag = '0318_1e4_realistic' & yr=[-20,20]
;tag = '0318_1e4_sn1'
;tag = '0317_1e5_realistic' & yr=[-20,20]

print, tag
;restore, 'sav_files/sim_sigbias_map_'+tag+'_stub.sav'
dir = '/data23/kstory/lps12/sim_sigbias/'
restore, dir+'sim_sigbias_map_'+tag+'.sav'
s = sim_st

nbins = s.nbins
nsims = s.nsims
l = s.l
hz = s.hz
dl = s.dl
dHz = dl * (hz[1]/l[1]) ; size of a bin in Hz
lbins = s.lbins
hzbins = lbins * hz[1]/l[1]

tmp = {psd_mb:fltarr(nbins, nsims), $ ; psd [uK^2/Hz], binned
       psd_nb:fltarr(nbins, nsims), $ ; psd [uK^2/Hz], binned
       psd_rb:fltarr(nbins, nsims), $ ; psd [uK^2/Hz], binned
       psd_sb:fltarr(nbins, nsims), $
       rpb:fltarr(nbins), rpb_err:fltarr(nbins), $
       r_psd_inband:fltarr(nsims), r_psd_low:fltarr(nsims), r_psd_high:fltarr(nsims), $
       n_psd_inband:fltarr(nsims), n_psd_low:fltarr(nsims), n_psd_high:fltarr(nsims), $
       s_psd_inband:fltarr(nsims), s_psd_low:fltarr(nsims), s_psd_high:fltarr(nsims)}

st0 = replicate(tmp,1)
st1 = replicate(tmp,1)
st2 = replicate(tmp,1)

;for isim=0L, 9 do begin
for isim=0, nsims-1 do begin
    if (isim mod 10 eq 0) then print, isim, format='((x, I0),$)' 
    for ii=0, nbins-1 do begin
        wh = where(l ge dl*ii and l lt dl*(ii+1))
        st0.psd_sb[ii,isim] = total(sim0.pow_s[wh,isim])/dHz

        st0.psd_mb[ii,isim] = total(sim0.pow_m[wh,isim])/dHz
        st0.psd_nb[ii,isim] = total(sim0.pow_n[wh,isim])/dHz
        st0.psd_rb[ii,isim] = total(sim0.pow_r[wh,isim])/dHz
        st1.psd_mb[ii,isim] = total(sim1.pow_m[wh,isim])/dHz
        st1.psd_nb[ii,isim] = total(sim1.pow_n[wh,isim])/dHz
        st1.psd_rb[ii,isim] = total(sim1.pow_r[wh,isim])/dHz
        st2.psd_mb[ii,isim] = total(sim2.pow_m[wh,isim])/dHz
        st2.psd_nb[ii,isim] = total(sim2.pow_n[wh,isim])/dHz
        st2.psd_rb[ii,isim] = total(sim2.pow_r[wh,isim])/dHz
    endfor
endfor

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

print, 'processing step 2'
;for j=0L, 9 do begin
for j=0L,nsims-1 do begin
    st0.r_psd_inband[j] = total(st0.psd_rb[wh_1to3,j]) * dHz / nHz_inband
    st0.r_psd_low[j]    = total(st0.psd_rb[wh_low,j])  * dHz / nHz_low
    st0.r_psd_high[j]   = total(st0.psd_rb[wh_high,j]) * dHz / nHz_high

    st0.s_psd_inband[j] = total(st0.psd_sb[wh_1to3,j]) * dHz / nHz_inband
    st0.s_psd_low[j]    = total(st0.psd_sb[wh_low,j])  * dHz / nHz_low
    st0.s_psd_high[j]   = total(st0.psd_sb[wh_high,j]) * dHz / nHz_high

    st0.n_psd_inband[j] = total(st0.psd_nb[wh_1to3,j]) * dHz / nHz_inband
    st0.n_psd_low[j]    = total(st0.psd_nb[wh_low,j])  * dHz / nHz_low
    st0.n_psd_high[j]   = total(st0.psd_nb[wh_high,j]) * dHz / nHz_high


    st1.r_psd_inband[j] = total(st1.psd_rb[wh_1to3,j]) * dHz / nHz_inband
    st1.r_psd_low[j]    = total(st1.psd_rb[wh_low,j])  * dHz / nHz_low
    st1.r_psd_high[j]   = total(st1.psd_rb[wh_high,j]) * dHz / nHz_high

    st1.s_psd_inband[j] = total(st1.psd_sb[wh_1to3,j]) * dHz / nHz_inband
    st1.s_psd_low[j]    = total(st1.psd_sb[wh_low,j])  * dHz / nHz_low
    st1.s_psd_high[j]   = total(st1.psd_sb[wh_high,j]) * dHz / nHz_high

    st1.n_psd_inband[j] = total(st1.psd_nb[wh_1to3,j]) * dHz / nHz_inband
    st1.n_psd_low[j]    = total(st1.psd_nb[wh_low,j])  * dHz / nHz_low
    st1.n_psd_high[j]   = total(st1.psd_nb[wh_high,j]) * dHz / nHz_high


    st2.r_psd_inband[j] = total(st2.psd_rb[wh_1to3,j]) * dHz / nHz_inband
    st2.r_psd_low[j]    = total(st2.psd_rb[wh_low,j])  * dHz / nHz_low
    st2.r_psd_high[j]   = total(st2.psd_rb[wh_high,j]) * dHz / nHz_high

    st2.s_psd_inband[j] = total(st2.psd_sb[wh_1to3,j]) * dHz / nHz_inband
    st2.s_psd_low[j]    = total(st2.psd_sb[wh_low,j])  * dHz / nHz_low
    st2.s_psd_high[j]   = total(st2.psd_sb[wh_high,j]) * dHz / nHz_high

    st2.n_psd_inband[j] = total(st2.psd_nb[wh_1to3,j]) * dHz / nHz_inband
    st2.n_psd_low[j]    = total(st2.psd_nb[wh_low,j])  * dHz / nHz_low
    st2.n_psd_high[j]   = total(st2.psd_nb[wh_high,j]) * dHz / nHz_high
endfor


if keyword_set(save_file) then begin
    save, st0, st1,st2,filename=dir+'sim_sigbias_map_'+tag+'_psd.sav'
endif

if keyword_set(stopit) then stop
END


