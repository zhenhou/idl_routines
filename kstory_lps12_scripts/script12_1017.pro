;;;
; Daily script
;
; Notes:
;  1) check CalErr chains
;
;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; How strongly do we rule out oml<0 ?
PRO oml

;----------------------
; Get the data
cdir = '/data23/kstory/lps12/chains/'
chains = 'c220_t8_lcdm_omk_camb_w7s12'

dir   = cdir+chains+'/chains/'
;files = file_search(dir+chains+'*.txt')
;pname = dir + chains+'.paramnames'
files = file_search(cdir+chains+'/chains_1.txt')
pname = cdir + chains+'/chains.paramnames'
temperature = 8

subsamp=25
;subsamp=10000
nskip=1000

stop
plot_like1dname,files,pname,'omegal*',subsamp=subsamp,nskip=nskip,scale=scale,temp=temperature,/cdf,ltx=0,/expinterp,/stopit


stop
END

