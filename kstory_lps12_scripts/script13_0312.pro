;;;
; NAME: script13_0312
;
; NOTES:
;  1) calculate nhits per pixel in individual maps
;  1) re-process sims
;;;


;;;;;;;;;;;;;;
;
; Util function: Get the weights
;
;;;;;;;;;;;;;;
PRO nhits
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

; get the map files
maps = get_lps12_runlist(field_idx, /xspec_maps)
;nmaps = n_elements(maps)
nmaps = 10
info = expand_fits_struct(read_spt_fits(maps[0])) ; get map shape

ss = size(info.map.map)

nhit = fltarr(nmaps, ss[1],ss[2])

tmp = 0*info.map.map
for ii=0, nmaps-1 do begin
    tmp *= 0.

    ; make some noise
    if (ii mod 20 eq 2) then print, "processing " + strtrim(string(ii),2) + "/" + strtrim(string(nmaps),2)

    d = expand_fits_struct(read_spt_fits(maps[ii]))
    weight = d.weight.map
    weight /=max(weight)

    nhit[ii,*,*] = weight
endfor


whg = where(nhit ne 0)
; plot stuff
h = histogram(nhit[whg], locations=xb,nbins=1000)
xb = xb+(xb[1]-xb[0])/2.
plot, xb, h, psym=10,xr=[0,1],/xst

stop



END


;;;;;;;;;;;;;;
;
; Util function: Get the weights
;
;;;;;;;;;;;;;;
PRO restuff_sims, input, output
restore, input

npix = [1536,960]
vec = indgen(npix[0]/2 - 1)


mask_expanded = [0]
for i=0,nsims-1 do mask_expanded = [mask_expanded, mask_row]
mask_expanded = mask_expanded[1:*]
whmask = where(mask_expanded gt 0.9)

; ------------------------
; Calculate bias
; ------------------------

s0 = {m:0., ms:fltarr(nsims), sd:0.}
s1 = {m:0., ms:fltarr(nsims), sd:0.}
s2 = {m:0., ms:fltarr(nsims), sd:0.}

for isim=0, nsims-1 do begin
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

save, ncoadd,nsims,npix,dl,l,hz,lbins,mask_row,whmask,mask_fac,cor,hpfilt,sigfac,sigma,n_tstream,sn_tstream,pn,$
  sim0,sim1,sim2,s0,s1,s2,$
  filename=output
END




;;;;;;;;;;;;;;
;
; Make plots for wiki post
;
;;;;;;;;;;;;;;
PRO make_plots
;restore, 'sav_files/sim_sigbias_map_0312_realistic.sav'
;yr=[-50,50]

;restore, 'sav_files/sim_sigbias_map_0312_maxBias.sav'
;yr=[0.8,1.2]

restore, 'sav_files/sim_sigbias_map_0312_manyCoadd.sav'
yr=[-9,11]

vec = indgen(npix[0]/2 - 1)

print, '--- make plot ---'
noise_tmp = randomn(seed, npix[0]) * sigma[0]
noise_tmp *= mask_row           ; apply mask
noise_tmp = float(fft(fft(noise_tmp)*hpfilt,1)) ; high-pass filter
pn = abs(fft(noise_tmp))^(2.) * cor

wset, 0
plot_sim, sim2, lbins, l, vec, pn=pn

wset, 5
plot_bias, sim0, sim1, sim2, s0,s1,s2,lbins, whmask,yr=yr

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
;err=tvread(/png,/nodialog,filename=fdir+'sim2_tsplots_0312')
stop

END
