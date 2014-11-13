;;;
; NOTES:
;  1) plot best-fit spectra from lrange test
;;;

PRO plot_lrange
;--------------------
; get spectra
;--------------------
readcol, '/data/kstory/projects/lps12/best_fit/lcdm_w7s12_lmax1500_ML_dl.txt',format='d,d',l_x15, dl_x15
readcol, '/data/kstory/projects/lps12/best_fit/lcdm_w7s12_1500ell3000_ML_dl.txt',format='d,d',l_x30, dl_x30

window, 1
xtxt='!12l!X!N'
plot, l_x30, dl_x30/dl_x15,yr=[0.9, 1.1],/yst,title='lrange',xtitle=xtxt,ytitle='dl_30/dl_15'
legend,['dl_(1500<ell<3000) / dl_(650<ell<1500)'],linestyle=0
;err=tvread(/png,/nodialog,filename='/home/kstory/public_html/notebook/spt_lps12/lrange_mlRatio_0826')

stop
END



;;;;;;;;;;;;;;;;
; 2) Make ini files lrange, to make ML spectra
;;;;;;;;;;;;;;;;
PRO make_ML_lrange
chain = 'c2_lcdm_pico_w7s12_lmax1500'
chain_dir = '/data23/hou/lps12/paramfits/chains_0828/c2_lcdm_pico_w7s12_lmax1500/chains/'
make_ML_ini, chain, chain_dir=chain_dir

chain = 'c2_lcdm_pico_w7s12_1500ell3000'
chain_dir = '/data23/hou/lps12/paramfits/chains_0828/c2_lcdm_pico_w7s12_1500ell3000/chains/'
make_ML_ini, chain, chain_dir=chain_dir

END
