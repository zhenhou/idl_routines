;;;
; NAME: plots_s12_lrange_MLspectra.pro
; PURPOSE:
;   Plot the differences between the best-fit models of different ell-ranges
; 
; Notes:
;  1) 
;
; MODIFICATION HISTORY:
;  08/27/2012: (KTS) Created
;;;

PRO plot_s12_lrange_MLspectra
readcol, '/data/kstory/projects/lps12/best_fit/lcdm_w7s12_lmax1500_ML_dl.txt',format='d,d',l_x15, dl_x15
readcol, '/data/kstory/projects/lps12/best_fit/lcdm_w7s12_1500ell3000_ML_dl.txt',format='d,d',l_x30, dl_x30

; subtract the foregrounds
; Get ML values, from params_c2_lcdm_pico_w7s12_lmax1500_ML.ini
poisson3000_x15 = 17.6350
;poisson3000_x15 = 21
clust_3000_x15  = 4.21565
sz3000_x15      = 5.17675

; Get ML values, from params_c2_lcdm_pico_w7s12_1500ell3000_ML.ini
poisson3000_x30 = 19.2641
clust_3000_x30  = 4.59950
sz3000_x30      = 4.73804

;---------------------
; get the residual Poisson power
help,poisson3000_x15
dl_poisson_x15 = poisson3000_x15*((l_x15/3000.)^2.)

help,poisson3000_x30
dl_poisson_x30 = poisson3000_x30*((l_x30/3000.)^2.)

;---------------------
; get the power from SZ
b=read_ascii('/home/cr/paramfits/cosmomc.r11/ptsrc/dl_shaw_tsz_s10_153ghz.txt')
yy = b.field1
; match the ell range with that of the cmb
istart = (where_closest(reform(yy[0,*]),2))[0]
istop = istart + n_elements(l_x15) -1
yy = yy[*,istart:istop]
l_sz = reform(yy[0,*])
dl_sz = reform(yy[1,*])
wh_3000 = (where_closest(l_sz,3000))[0]

dl_sz_x15 = dl_sz/dl_sz[wh_3000]*sz3000_x15
dl_sz_x30 = dl_sz/dl_sz[wh_3000]*sz3000_x30

;---------------------
; get the power from spatially correlated galaxies and tSZ and kSZ
help, clust_3000_x15
dl_flat_x15 = clust_3000_x15*((l_x15/3000.)^(0.8))
dl_flat_x30 = clust_3000_x30*((l_x30/3000.)^(0.8))


;---------------------
; subtract foregrounds
;---------------------
dl_cmb15 = dl_x15 - dl_poisson_x15 - dl_sz_x15 - dl_flat_x15
dl_cmb30 = dl_x30 - dl_poisson_x30 - dl_sz_x30 - dl_flat_x30

;dl_poisson_x15b = 21*((l_x15/3000.)^2.)
;dl_cmb15b = dl_x15 - dl_poisson_x15b - dl_sz_x15 - dl_flat_x15
dl_cmb15b = dl_x15 - dl_poisson_x30 - dl_sz_x30 - dl_flat_x30

; wset, 1
; xtxt='!12l!X!N'
; ;plot, l_x30, (dl_x15-dl_x30)/dl_x15,yr=[0.9, 1.1],/yst,title='lrange',xtitle=xtxt,ytitle='dl_30/dl_15'
; plot, l_x30, (dl_x30 - dl_x15)/dl_x30,yr=[-0.1,0.1],/yst,title='lrange comparison',xtitle=xtxt,ytitle='(dl_30 - dl_15)/dl_30'

; wset, 2
plot, l_x30, (dl_cmb30 - dl_cmb15)/dl_cmb30,yr=[-0.2,0.2],/yst,xr=[650,3000],/xst,$
  title='lrange comparison',xtitle=xtxt,ytitle='(dl_30 - dl_15)/dl_30'
oplot, l_x30, (dl_cmb30 - dl_cmb15b)/dl_cmb30,color=!purple
;err=tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/lrange_mlRatio')


stop
END

