;;;
; Daily script
;
; Notes:
;  1) check high-temperature chain treatement
;;;




;;;;;;;;;;
;
; Simplified treatment
;
;;;;;;;;;;
PRO limit,files,nparam,temp=temp,gtx=gtx,ltx=ltx,nskip=nskip,stopit=stopit

if n_elements(nskip) ne 1 then nskip = 1000
if n_elements(temp) eq 0 then temp=1.
if (n_elements(gtx) eq 0 and n_elements(ltx) eq 0) then begin
    print, 'must set either gtx or ltx.  Returning.'
    RETURN
endif

a = read_ascii(files[0])
pp = a.field01[nparam,nskip:*]
wt = long(a.field01[0,nskip:*])
like = a.field01[1,nskip:*]

; read all files
nf = n_elements(files)
print, format='($,%"Read file: ")'
for i=1, nf-1 do begin
    print, FORMAT='($,%"%d, ")', i
    b = read_ascii(files[i])

    pp2 = b.field01[nparam,nskip:*]
    wt2 = long(b.field01[0,nskip:*])
    like2 = b.field01[1,nskip:*]

    pp = [reform(pp),reform(pp2)]
    wt = [reform(wt),reform(wt2)]
    like = [reform(like),reform(like2)]

endfor
print, ''

; Find ML value
ind = (where(like eq min(like)))[0]

; calculate new weights
dLnl = like - like[ind]
ww = wt * double(EXP( -(temp-1)*dLnl ))

if n_elements(ltx) ne 0 then begin
    wh = where(pp lt ltx, nwh)

    px = total(ww[wh]) / total(ww)
    spx = string(px,format='(e)')
    print,'P(<x) ',ltx,' is: ',spx,', ',gauss_cvf(double(px)),' sigma'
endif
if n_elements(gtx) ne 0 then begin
    wh = where(pp gt gtx, nwh)

    px = total(ww[wh]) / total(ww)
    spx = string(px,format='(e)')
    print,'P(>x) ',gtx,' is: ',spx,', ',gauss_cvf(double(px)),' sigma'
endif

if keyword_set(stopit) then stop
END




;;;;;;;;;;
;
; attempt 1
;
;;;;;;;;;;
PRO read_t_chain,files,nparam,nskip=nskip,stopit=stopit,subsamp=subsamp,offset=offset,results=results,_extra=extras,lincomb=lcomb, power=power,cdf=cdf

if not tag_exist(/top_level,extras,'zeropad') then $
  zeropad=0 $
else $
  zeropad=extras.zeropad

if not tag_exist(/top_level,extras,'temp') then begin
    temp=1
endif else begin
    temp=extras.temp
endelse

if not tag_exist(/top_level,extras,'xtitle') then $
  if n_elements(extras) ne 0 then $
  extras = create_struct('xtitle',strtrim(string(fix(nparam)),2),extras) $
  else $
  extras = {xtitle:strtrim(string(fix(nparam)),2)}
if not tag_exist(/top_level,extras,'ytitle') then begin
    if keyword_set(cdf) then begin
        extras = create_struct('ytitle','CDF',extras) 
    endif else begin
        extras = create_struct('ytitle','Likelihood',extras) 
    endelse
endif


;if n_elements(xtitle) ne 1 then xtitle='param_'+strtrim(string(fix(nparam)),2)
;if n_elements(ytitle) ne 1 then ytitle='Likelihood'

if n_elements(subsamp) ne 1 then subsamp=5.
if not tag_exist(/top_level,extras,'scale') then $
  scale=1.0 $
else $
  scale=extras.scale
;print,scale

if n_elements(power) eq 0 then power = 1.0
if n_elements(nskip) ne 1 then nskip=4000
if n_elements(offset) ne 1 then offset =0.

a=read_ascii(files[0])
pp = a.field01[nparam,nskip:*]*scale
wt = long(a.field01[0,nskip:*])
like = double(a.field01[1,nskip:*])

if  tag_exist(/top_level,extras,'priorrng') then begin
    pwt = total(wt)
    prng = extras.priorrng
    prip = reform(prng[0,*])
    primn = reform(prng[1,*])
    primx = reform(prng[2,*])
    nprip = n_elements(prip)
    for k=0,nprip-1 do begin
        prp = a.field01[prip[k],nskip:*]
        ind = where(prp lt primn[k] or prp gt primx[k],nbad)

        if nbad gt 0 then wt(ind) = 0
    endfor
    print,'prior kept: ',float(total(wt))/total(pwt)
endif

nf = n_elements(files)

for i=1,nf-1 do begin
    b=read_ascii(files[i])
    pp2 = b.field01[nparam,nskip:*]*scale
    wt2 = long(b.field01[0,nskip:*])
    like2 = b.field01[1,nskip:*]

    if  tag_exist(/top_level,extras,'priorrng') then begin
        for k=0,nprip-1 do begin
            prp = b.field01[prip[k],nskip:*]
            ind = where(prp lt primn[k] or prp gt primx[k],nbad)
            if nbad gt 0 then wt2(ind) = 0
        endfor
    endif
    
    pp = [reform(pp),reform(pp2)]
    wt = [reform(wt),reform(wt2)]
    like = [reform(like),reform(like2)]
    stop
endfor

ind = where(wt gt 0,ngood)
if ngood lt 2000 then begin
    print,'only X samples, stopping',ngood
    stop
endif
pp = pp(ind)
wt = wt(ind)

;first expand pp into fully sampled array
np = total(wt)

pp2 = fltarr(np)
npp = n_elements(pp)
tnp = [0,reform(total(wt,/cum))-1]

for i=0l,npp-1 do begin
    pp2[tnp[i]:tnp[i+1]]=pp[i]
endfor
oldpp=pp
pp=pp2+offset

pp = pp^(power)

if temp ne 1 then begin
    print, 'high temp'
endif

sig = stddev(pp)
if (sig le 0) or (finite(sig) eq 0) then begin
    print,'data constant, no histogram possible'
    results={fail:1}
    return
endif
print,sig
binsize=sig/subsamp
; if keyword_set(constraint) then begin
;     binsize = (max(pp) - min(pp)) / n_elements(pp)
; endif
h=double(histogram(pp,binsize=binsize,omin=omn,omax=omx))
if n_elements(h) eq 1 then begin
    print,'data constant2, no histogram possible'
    results={fail:1}
    return
endif

h/=max(h)
;h = h^temp

bins = (findgen(n_elements(h))+.5)*binsize +omn

if zeropad eq 1 then begin
    h=[0,h,0]
    bins=[min(bins)-binsize,bins,max(bins)+binsize]
endif
boplot=0
if tag_exist(extras,'oplot',/top_level) then boplot =  extras.oplot 
bnoplot=0
if tag_exist(extras,'noplot',/top_level) then bnoplot =  extras.noplot 

; plot the CDF if requested
if keyword_set(cdf) then begin
    nbins = n_elements(h)
    ff = dblarr(nbins)
    for i=0L, nbins-1 do ff[i] = double(total(h[0:i]))
    ff = ff/max(ff)
    ff = 1-ff

    if bnoplot eq 0 then begin
        if boplot eq 0 then begin
            plot,bins,ff,_extra=extras
;,xtitle=xtitle,ytitle=ytitle,title=title,color=color,linestyle=linestyle,psym=psym
        endif else begin
            inds = where(bins ge !x.crange[0] and bins le !x.crange[1])
            if (extras.yl ne 0) then inds = where(bins ge !x.crange[0] and bins le !x.crange[1] and ff ge 10^(!y.crange[0]) and ff le 10^(!y.crange[1]) )
            oplot,bins(inds),ff(inds),_extra=extras
        endelse
    endif
    bnoplot = 1

; ; --- For ns>1 ---------------
;     if keyword_set(ns_gt_1) then begin
;         if temp eq 1 then begin
;             wh = where(pp ge 1.0, nwh) & frac = nwh/double(n_elements(pp)) 
;             print, 'Frac above 1 = ', frac, ', sigma = ', gauss_cvf(frac)
;         endif else begin
;             plot, bins, h
;             wh = where(bins ge 1.0)
;             d1 = 1./abs(1.-bins[wh[0]-1])
;             d2 = 1./abs(1.-bins[wh[0]])
;             frac = double( (d1*ff[wh[0]-1] + d2*ff[wh[0]]) / (d1+d2) )
;             print, 'Frac above 1 = ', frac, ', sigma = ', gauss_cvf(frac)
;             ;print, val, gauss_cvf(val)
;             stop
;         endelse
;     endif
; ;---- end ---------------------

; ; --- For alens<0 ---------------
;     if keyword_set(alens_lt_0) then begin
;         if temp eq 1 then begin
;             wh = where(pp le 0.0, nwh) & frac = nwh/double(n_elements(pp)) 
;             print, 'Frac above 1 = ', frac, ', sigma = ', gauss_cvf(frac)
;         endif else begin
;             plot, bins, h
;             wh = where(bins gt 0.0)
;             jj0 = wh[0]-1
;             jj1 = jj0+1
;             bb0 = bins[jj0]
;             bb1 = bins[jj1]

;             d1 = 1./abs(bins[wh+1])
;             d2 = 1./abs(bins[wh])
;             frac = 1. - double( (d1*ff[wh+1] + d2*ff[wh]) / (d1+d2) )
;             print, 'Frac below 0 = ', frac, ', sigma = ', gauss_cvf(frac)
;             ;print, val, gauss_cvf(val)
;             stop
;         endelse
;     endif
; ;---- end ---------------------

endif

if bnoplot eq 0 then begin
    if boplot eq 0 then begin
        plot,bins,h,_extra=extras
;,xtitle=xtitle,ytitle=ytitle,title=title,color=color,linestyle=linestyle,psym=psym
    endif else begin
        inds = where(bins ge !x.crange[0] and bins le !x.crange[1])
        oplot,bins(inds),h(inds),_extra=extras
    endelse
endif
;,color=color,linestyle=linestyle,psym=psym
res = gaussfit(bins,h,a,nterms=3)
print,a
z = total(h,/cum)/total(h)

if tag_exist(/top_level,extras,'ltx') then begin
    junkind = where(bins gt extras.ltx,njunk)
    if njunk le 0 then print,'no elements above limit:',extras.ltx
    jj0=junkind[0]-1
    jj1=jj0+1
    bb0=bins[jj0]
    bb1=bins[jj1]
    p0 = z[jj0]
    p1 = z[jj1]
    if tag_exist(/top_level,extras,'expinterp') then $
      px = p0 + exp( alog(p1-p0)/alog(bb1-bb0) * alog(extras.ltx-bb0)) $
    else $
      px = p0 + (p1-p0)/(bb1-bb0) * (extras.ltx - bb0) 

    spx = string(px,format='(e)')
    print,'P(<x) ',extras.ltx,' is: ',spx,', ',gauss_cvf(double(px)),' sigma'
    stop
endif
if tag_exist(/top_level,extras,'gtx') then begin
    junkind = where(bins gt extras.gtx,njunk)
    if njunk le 0 then print,'no elements below limit:',extras.gtx
;    help,pp,where(pp gt extras.gtx)
    jj0=junkind[0]-1
    jj1=jj0+1
    bb0=bins[jj0]
    bb1=bins[jj1]
    p0 = z[jj0]
    p1 = z[jj1]
    if tag_exist(/top_level,extras,'expinterp') then $
      px = p0 + exp( alog(p1-p0)/alog(bb1-bb0) * alog(extras.gtx-bb0)) $
    else $
      px = p0 + (p1-p0)/(bb1-bb0) * (extras.gtx - bb0) 

    spx = string(1-px,format='(e)')
    print,'P(>x) ',extras.gtx,' is: ',spx,', ',gauss_cvf(double(1.-px)),' sigma'
    stop
endif



ind = where(z gt 0.95)
ind2 = where(z gt 0.05,n)
ind3 = where(z gt 0.67)
ind4 = where(z gt 0.023,n)
ind5 = where(z gt 0.977)
ind6 = where(z gt 0.159)
ind7 = where(z gt 0.841)

print,'95% limit: ',bins[ind[0]]
print,'67% limit: ',bins[ind3[0]]
uppers=[bins[ind3[0]],bins[ind[0]]]
print,'5% limit: ',bins[ind2[0]]
print,'2sigma: 2.3 - 97.7% limit: ',bins[ind4[0]],bins[ind5[0]]
print,'1sigma: 15.9 - 84.1% limit: ',bins[ind6[0]],bins[ind7[0]]

print,'2.3% limit: ',bins[ind5[0]]

mean = total(bins*h)/total(h)
print,'mean ',mean
ind = where(z gt .5-.683/2)
ind2 = where(z gt .5+.683/2)
ind3 = where(z gt .5)
print,bins[ind[0]],bins[ind2[0]]
print,'median ',bins(ind3[0])
print,'err+, err-, avg_err ',bins[ind2[0]]-mean, bins[ind[0]]-mean, (bins[ind2[0]]-bins[ind[0]])/2.

results={bins:bins,prob:h,rng:[bins[ind[0]],bins[ind2[0]]],$
         avg_err:(bins[ind2[0]]-bins[ind[0]])/2.,median:bins(ind3[0]),$
         err_plus:bins[ind2[0]]-mean, err_minus:bins[ind[0]]-mean,$
         uppers:uppers,$
         fail:0}

; Use this for getting the most Gaussain power on alens
if tag_exist(extras,'get_chain',/top_level) then $
  results={bins:bins,prob:h,rng:[bins[ind[0]],bins[ind2[0]]],$
           avg_err:(bins[ind2[0]]-bins[ind[0]])/2.,median:bins(ind3[0]),$
           err_plus:bins[ind2[0]]-mean, err_minus:bins[ind[0]]-mean,$
           uppers:uppers,$
           fail:0, $
           pp:pp}

if keyword_set(stopit) then stop


END



;;;;;;;;;;;;;;;;;
;
; Test program
;
;;;;;;;;;;;;;;;;;
PRO al, type=type
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

if n_elements(type) eq 0 then type='cmb'
temperature = 6.0 ; both are high-T chains

;;; LCDM
case type of
    's12' : begin
        ; S12-only
        dir   = cdir+'c53_lcdm_alens_camb_s12tau/chains/'
        files = file_search(dir+'c53_lcdm_alens_camb_s12tau*.txt')
        pname = dir+'c53_lcdm_alens_camb_s12tau.paramnames'
    endcase

    'cmb' : begin
        ; cmb
        dir = cdir+'c55_lcdm_alens_camb_w7s12/chains/'
        files = file_search(dir+'c55_lcdm_alens_camb_w7s12*.txt')
        pname = dir+'c55_lcdm_alens_camb_w7s12.paramnames'
    endcase

    else : begin
        print, 'un-recognized case. Returning'
        RETURN
    endcase
endcase

print, 'Alens constraints for dataset, ' + type
subsamp=25
nskip=1000

limit,files,14,temp=temperature,ltx=0,nskip=nskip,/stopit
;read_t_chain,files,14,subsamp=subsamp,nskip=nskip,scale=scale,temp=temperature,/cdf,ltx=0,/stopit
;plot_like1dname,files,pname,'alens',subsamp=subsamp,nskip=nskip,scale=scale,temp=temperature,/cdf,ltx=0,/stopit

stop
END


;;;;;;;;;;;;;;;;;;;;;;
; Get CDF limits on ns
; Power-law [LCDM]
;   Print constraints for all LCDM parameters using this procedure
;   type in ['s12', 'w7', 'cmb', 'extra']
;;;;;;;;;;;;;;;;;;;;;;
PRO ns_cdf, type=type
cdir = '/data23/hou/lps12/paramfits/chains_0828/'

if n_elements(type) eq 0 then type='cmb'
temperature = 1.0 ; use for high-T chains
nparam = 15 ; for ns

;;; LCDM
case type of
    's12' : begin
        ; S12-only
        dir   = cdir+'c27_lcdm_camb_s12tau/chains/'
        files = file_search(dir+'c27_lcdm_camb_s12tau*.txt')
        pname = dir+'c27_lcdm_camb_s12tau.paramnames'
    endcase

    'w7' : begin
        ; WMAP-only
        dir = cdir+'c1_lcdm_pico_w7/chains/'
        files = file_search(dir+'c1_lcdm_pico_w7*.txt')
        pname = dir+'c1_lcdm_pico_w7.paramnames'
    endcase

    'cmb' : begin
        ; cmb
        dir = cdir+'c2_lcdm_pico_w7s12/chains/'
        files = file_search(dir+'c2_lcdm_pico_w7s12*.txt')
        pname = dir+'c2_lcdm_pico_w7s12.paramnames'
    endcase

    'cmb_h0' : begin
        ; CMB+H0
        dir = cdir+'c4_lcdm_pico_w7s12_H0/chains/'
        files = file_search(dir+'c4_lcdm_pico_w7s12_H0*.txt')
        pname = dir+'c4_lcdm_pico_w7s12_H0.paramnames'
    endcase

    'cmb_h0_highT' : begin
        ; CMB+H0
        dir = cdir+'c98_lcdm_pico_w7s12_H0/chains/'
        files = file_search(dir+'c98_lcdm_pico_w7s12_H0*.txt')
        pname = dir+'c98_lcdm_pico_w7s12_H0.paramnames'
        temperature = 6.0
    endcase

    'cmb_bao' : begin
       ; CMB+BAO 
        dir = cdir+'c3_lcdm_pico_w7s12_BAO/chains/'
        files = file_search(dir+'c3_lcdm_pico_w7s12_BAO*.txt')
        pname = dir+'c3_lcdm_pico_w7s12_BAO.paramnames'
    endcase

    'cmb_bao_highT' : begin
       ; CMB+BAO 
        dir = cdir+'c109_lcdm_pico_w7s12_BAO/chains/'
        files = file_search(dir+'c109_lcdm_pico_w7s12_BAO*.txt')
        pname = dir+'c109_lcdm_pico_w7s12_BAO.paramnames'
        temperature = 6.0 
    endcase

    'cmb_h0_bao' : begin
        ; CMB+H0+BAO
        dir = cdir+'c6_lcdm_pico_w7s12_BAO_H0/chains/'
        files = file_search(dir+'c6_lcdm_pico_w7s12_BAO_H0*.txt')
        pname = dir+'c6_lcdm_pico_w7s12_BAO_H0.paramnames'
        temperature = 6.0 
    endcase

    'cmb_h0_neff': begin
        type = 'c14_lcdm_neff_pico_w7s12_H0'
        dir = cdir+type+'/chains/'
        files = file_search(dir+type+'_*.txt')
        pname = dir+type+'.paramnames'
    endcase

    else : begin
        print, 'un-recognized case. Returning'
        RETURN
    endcase

endcase

;stop

print, 'NS constraints for dataset, ' + type
nskip=1000
;plot_like1dname,files,pname,'ns',subsamp=subsamp,nskip=nskip,scale=scale,temp=temperature,/cdf,gtx=1,/expinterp,/stopit
limit,files,nparam,nskip=nskip,temp=temperature,gtx=1;,/stopit


stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Get CDF limits on oml
;;;;;;;;;;;;;;;;;;;;;;
PRO oml_cdf, type=type
cdir = '/data23/kstory/lps12/chains/'

if n_elements(type) eq 0 then type='cmb'
temperature = 1.0 ; use for high-T chains
nparam = 26 ; for ns
scale=1

;;; LCDM
case type of
    'cmb' : begin
        ; cmb
        dir = cdir+'c220_t8_lcdm_omk_camb_w7s12/chains/'
        files = file_search(dir+'c220_t8_lcdm_omk_camb_w7s12*.txt')
        pname = dir+'c220_t8_lcdm_omk_camb_w7s12.paramnames'
        temperature = 8.0
    endcase

    else : begin
        print, 'un-recognized case. Returning'
        RETURN
    endcase

endcase


print, 'OML constraints for dataset, ' + type, ', temperature = ', temperature
nskip=1000
;plot_like1dname,files,pname,'omegal*',subsamp=subsamp,nskip=nskip,scale=scale,temp=temperature,/cdf,gtx=1,/expinterp,/stopit
limit,files,nparam,nskip=nskip,temp=temperature,ltx=0,/stopit


stop
END


