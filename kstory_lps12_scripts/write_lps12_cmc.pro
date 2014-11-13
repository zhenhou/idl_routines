;;;
; NAME: write_lps12_cmc.pro
; PURPOSE:
;   Write 'newdat' files for mcmc parameter estimations.
;
; CALLING SEQUENCE: 
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   fake,                    Replace the real spectum with a fake input spectrum
;
; KEYWORD PARAMETERS:
;   k11,                     Use this to use only 2008, 2009 data.
;
; OUTPUTS:
;
; NOTES:
;   1) 
;
; MODIFICATION HISTORY:
;   05/11/2012 (KTS) : Copied from Ryan's directory, when it was named write_lowell_cmc.pro
;   06/11/2012 (KTS) : Add run_05, re-name window functions
;   06/13/2012 (KTS) : Apply correction to window function to comply with newdat convention.
;   06/14/2012 (KTS) : Add run keyword argument
;   07/17/2012 (KTS) : Add run_08 sav file
;   10/30/2012 (KTS) : Commit final version, with run_09 files
;   11/15/2012 (KTS) : Add "nocal" keyword
;;;

PRO write_lps12_cmc, run, fake=fake, k11=k11, only1011=only1011, beamErr=beamErr, halfcal=halfcal, doublecal=doublecal, nocal=nocal
;, halfbeam=halfbeam, doublebeam=doublebeam, area2010=area2010, area2011=area2011, zerobeam=zerobeam, zerocal=zerocal

if keyword_set(halfcal) or $
  keyword_set(doublecal) or $
  keyword_set(halfbeam) or $
  keyword_set(doublebeam) or $
  keyword_set(fake) or $
  keyword_set(area2010) or $
  keyword_set(area2011) or $
  keyword_set(zerobeam) or $
  keyword_set(nocal) $
  then dowf=0 else dowf=1

case run of

    ; FINAL bandpowers
    '09' : begin
        file='/home/kstory/lps12/end2end/run_09/combined_spectrum_20120828_170101_kweight.sav'
        odir='/home/kstory/lps12/end2end/run_09/spt_lps12_20120828/'
        wdir = 'window_lps12'
        if keyword_set(k11) then begin
             file='/home/kstory/lps12/end2end/run_09/combined_spectrum_20120830_154612_kweight_0809.sav'
             odir='/home/kstory/lps12/end2end/run_09/spt_lps12_0809_20120829/'
             wdir = 'window_0809'
        endif
        if keyword_set(only1011) then begin
             file='/home/kstory/lps12/end2end/run_09/combined_spectrum_20120830_210134_kweight_1011.sav'
             odir='/home/kstory/lps12/end2end/run_09/spt_lps12_1011_20120831/'
             wdir = 'window_1011'
        endif
        if keyword_set(beamErr) then begin
             file='/home/kstory/lps12/end2end/run_09/combined_spectrum_20120829_225610_kweight_beamErr.sav'
             odir='/home/kstory/lps12/end2end/run_09/spt_lps12_beamErr_20120831/'
             wdir = 'window_beamErr'
        endif
        if keyword_set(halfcal) then begin
             file='/home/kstory/lps12/end2end/run_09/combined_spectrum_20121015_122000_kweight_calErr_small.sav'
             odir='/home/kstory/lps12/end2end/run_09/spt_lps12_calErr_halfcal_20121015/'
             wdir = 'window_beamErr'
        endif
        if keyword_set(doublecal) then begin
             file='/home/kstory/lps12/end2end/run_09/combined_spectrum_20121015_122000_kweight_calErr_big.sav'
             odir='/home/kstory/lps12/end2end/run_09/spt_lps12_calErr_doublecal_20121015/'
             wdir = 'window_beamErr'
        endif
        if keyword_set(nocal) then begin
             file='/home/kstory/lps12/end2end/run_09/combined_spectrum_20121115_101000_kweight_NOcalErr.sav'
             odir='/home/kstory/lps12/end2end/run_09/spt_lps12_NOcalErr_20121115/'
             wdir = 'window_beamErr'
        endif
    end 

    '08' : begin
        file='/home/kstory/lps12/end2end/run_08/combined_spectrum_20120717_174249_kweight.sav'
        odir='/home/kstory/lps12/end2end/run_08/spt_lps12_20120717/'
        wdir = 'window_lps12'
        if keyword_set(k11) then begin
             file='/home/kstory/lps12/end2end/run_08/combined_spectrum_20120718_132910_kweight_0809.sav'
             odir='/home/kstory/lps12/end2end/run_08/spt_lps12_20120718/'
             wdir = 'window_0809'
        endif
    end 

    '07' : begin
        file='/home/kstory/lps12/end2end/run_07/combined_spectrum_20120617_144542_kweight.sav'
        odir='/home/kstory/lps12/end2end/run_07/spt_lps12_20120617/'
        wdir = 'window_lps12'
        if keyword_set(k11) then begin
;             file='/home/kstory/lps12/end2end/run_07/'
;             odir='/home/kstory/lps12/end2end/run_07/spt_lps12_20120614/'
;             wdir = 'window_0809'
        endif
    end 

    '05' : begin
        file='/home/kstory/lps12/end2end/run_05/combined_spectrum_20120614_192137_kweight.sav'
        odir='/home/kstory/lps12/end2end/run_05/spt_lps12_20120614/'
        wdir = 'window_lps12'
        if keyword_set(k11) then begin
            file='/home/kstory/lps12/end2end/run_05/combined_spectrum_20120614_193416_kweight_0809.sav'
            odir='/home/kstory/lps12/end2end/run_05/spt_lps12_20120614/'
            wdir = 'window_0809'
        endif
    end 
    else : begin
        print, "Not set up for anything but run='05'.  Returning..."
        RETURN
    end
endcase

print, 'use input file: '+file
restore,file

; if desired, replace the real spectrum with a fake spectrum (a
; perfect theory spectrum).
if keyword_set(fake) then begin
    ntmp = n_elements(dl_all)
    dl_all_fake = fltarr(ntmp)
    readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',lth,dl_uK2
    dlth_wf = interpol(dlth,lth,l_wf)/1d12
    for i=0,ntmp-1 do dl_all_fake[i] = total(dlth_wf*wf_all[*,i])
    dl_all = dl_all_fake
endif

if dowf then begin
    spawn,'mkdir '+odir,soutput
    spawn,'mkdir '+odir+'windows/',soutput
    spawn,'mkdir '+wdir
endif

window = wdir+'/window_'
comments='lps12 run_05'
nn=n_elements(diag)



;dump=9
;nn=nn-dump
istart=9
istop=55
nn = istop-istart+1
nbands = [nn,0,0,0,0,0]
bandsel=intarr(2,6)
bandsel[0,0]=1
bandsel[1,0]=nn
ls=lonarr(2,nn)
l=l[istart:istop]
ls[0,*]=long(l+.0001)-25
ls[1,*]=long(l+.0001)+24
beamerrtype=0
fwhm=0.
fwhmerr=0.
cal=1.0
calerr=0.00

cov = cov_all[istart:istop,istart:istop]*1e24
bands = fltarr(4,nn)+1e6
bands[0,*]=dl_all[istart:istop]*1e12
err = diag_nobeam[istart:istop]*1e12
;err = sqrt(diag[istart:istop])*1e12
bands[1,*]=err
bands[2,*]=err

correl=fltarr(nn,nn)
for i=0,nn-1 do for j=0,nn-1 do $
  correl[j,i]=cov[j,i]/sqrt(cov[j,j]*cov[i,i])
cmc={window:window,$
     comments:comments,$
     bands:bands,$
     ls:ls,$
     cal:cal,$
     calerr:calerr,$
     fwhm:fwhm,$
     fwhmerr:fwhmerr,$
     nbands:nbands,$
     bandsel:bandsel,$
     correl:correl,$
     cov:cov,$
     beamerrtype:beamerrtype}

suffix = ''
if keyword_set(halfcal) then suffix += '_halfcal'
if keyword_set(doublecal) then suffix += '_doublecal'
if keyword_set(nocal) then suffix += '_nocal'
if keyword_set(halfbeam) then suffix += '_halfbeam'
if keyword_set(doublebeam) then suffix += '_doublebeam'
if keyword_set(fake) then suffix += '_fake'
if keyword_set(area2010) then suffix += '_area2010'
if keyword_set(area2011) then suffix += '_area2011'
if keyword_set(zerobeam) then suffix += '_zerobeam'
if keyword_set(zerocal) then suffix += '_zerocal'
if keyword_set(k11) then suffix += '_0809'
if keyword_set(only1011) then suffix += '_1011'
if keyword_set(only1011) then suffix += '_beamErr'


print, 'write output file to: '+odir+'Spectrum_spt2500deg2_lps12'+suffix+'.newdat'

write_cosmomc,odir+'Spectrum_spt2500deg2_lps12'+suffix+'.newdat',cmc, $
  comment='SPT bandpowers - Kyle Story'

if dowf then begin
    winfuncs = wf_all[0:3250,istart:istop]

    ; Apply the correction for the newdat convention
    ells = dindgen(4000-50-1)+50
    ells = ells[0:3250]
    nl = n_elements(ells)
    nkept = istop - istart + 1
    correction = rebin((ells+1.0)/(ells+0.5),nl,nkept)

    winfuncs *= correction

    lwin=l_wf[0:3250]
    banddef = reform(ls[1,*])

    window_functions_to_oz_cosmomc, banddef, winfuncs, window, $
      ells=lwin, ellmin=0, ellmax=3300
    spawn,'mv '+wdir+' '+odir+'windows/',soutput
endif






;stop


end
