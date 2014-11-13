;code to plot 2d likelihood contours for cosmomc outputs

pro plot_like2d,files,nparam1,nparam2,nskip=nskip,sigma=sigma,scale1=scale1,$
                scale2=scale2,stopit=stopit,fracwt=fracwt,$
                power1=power1, power2=power2,$
                lincombx=lcombx,lincomby=lcomby,_extra=extras
if n_elements(scale1) ne 1 then scale1=1.
if n_elements(scale2) ne 1 then scale2=1.
if n_elements(power1) eq 0 then power1 = 1.0
if n_elements(power2) eq 0 then power2 = 1.0
nf = n_elements(files)


if not tag_exist(/top_level,extras,'xtitle') then $
  if n_elements(extras) ne 0 then $
  extras = create_struct('xtitle',strtrim(string(fix(nparam1)),2),extras) $
  else $
  extras = {xtitle:strtrim(string(fix(nparam2)),2)}
if not tag_exist(/top_level,extras,'ytitle') then $
    extras = create_struct('ytitle',strtrim(string(fix(nparam2)),2),extras)
if not tag_exist(/top_level,extras,'ycharsize') then $
    extras = create_struct('ycharsize',1.5,extras)
if not tag_exist(/top_level,extras,'xcharsize') then $
    extras = create_struct('xcharsize',1.5,extras)

aypos=0
if  tag_exist(/top_level,extras,'ypos') then $
  ypos = extras.ypos


colors = [!red,!green,!blue,!yellow,!black,!white,!skyblue,!purple]
if tag_exist(/top_level,extras,'icolors') then $
  colors=extras.icolors

smooth=0.0
if tag_exist(/top_level,extras,'smooth') then $
  smooth = extras.smooth > 0

subsamp=5.0
if tag_exist(/top_level,extras,'subsamp') then $
  subsamp = extras.subsamp > 0

if n_elements(nskip) ne 1 then nskip=4000

a=read_ascii(files[0])
if n_elements(lcombx) gt 0 then begin
    pp1=0
    for j=0,n_elements(lcombx.inds)-1 do begin
        pp1 += lcombx.factors[j]*a.field01[lcombx.inds[j],nskip:*]*scale1
    endfor
endif else $
  pp1 = a.field01[nparam1,nskip:*]*scale1
if n_elements(lcomby) gt 0 then begin
    pp2=0
    for j=0,n_elements(lcomby.inds)-1 do begin
        pp2 += lcomby.factors[j]*a.field01[lcomby.inds[j],nskip:*]*scale2
    endfor
endif else $
  pp2 = a.field01[nparam2,nskip:*]*scale2
twt = a.field01[0,nskip:*]

if  tag_exist(/top_level,extras,'priorrng') then begin
    prng = extras.priorrng
    prip = reform(prng[0,*])
    primn = reform(prng[1,*])
    primx = reform(prng[2,*])
    nprip = n_elements(prip)
    for k=0,nprip-1 do begin
        prp = a.field01[prip[k],nskip:*]
        ind = where(prp lt primn[k] or prp gt primx[k],nbad)
        if nbad gt 0 then twt(ind) = 0
    endfor
endif



for i=1,nf-1 do begin
    a=read_ascii(files[0])
    if n_elements(lcombx) gt 0 then begin
        app1=0
        for j=0,n_elements(lcombx.inds)-1 do begin
            app1 += lcombx.factors[j]*a.field01[lcombx.inds[j],nskip:*]*scale1
        endfor
    endif else $
      app1 = a.field01[nparam1,nskip:*]*scale1
    if n_elements(lcomby) gt 0 then begin
        app2=0
        for j=0,n_elements(lcomby.inds)-1 do begin
            app2 += lcomby.factors[j]*a.field01[lcomby.inds[j],nskip:*]*scale2
        endfor
    endif else $
      app2 = a.field01[nparam2,nskip:*]*scale2
    atwt = a.field01[0,nskip:*]


    if  tag_exist(/top_level,extras,'priorrng') then begin
        for k=0,nprip-1 do begin
            prp = a.field01[prip[k],nskip:*]
            ind = where(prp lt primn[k] or prp gt primx[k],nbad)
            if nbad gt 0 then atwt(ind) = 0
        endfor
    endif


    pp1 = [reform(pp1),reform(app1)]
    pp2 = [reform(pp2),reform(app2)]
    twt = [reform(twt),reform(atwt)]
endfor

if n_elements(fracwt) eq 1 then begin
    mwt=max(twt)
    twt*= (float(fracwt)/mwt)
    wt=long(twt)
endif else $
  wt = long(twt)

ind = where(wt gt 0)
pp1=pp1(ind)
pp2=pp2(ind)
wt=wt(ind)

np = total(wt)
print,np
qq1 = fltarr(np)
qq2 = fltarr(np)
npp = n_elements(pp1)
tnp = [0,reform(total(wt,/cum))-1]

for i=0l,npp-1 do begin
    qq2[tnp[i]:tnp[i+1]]=pp2[i]
    qq1[tnp[i]:tnp[i+1]]=pp1[i]
endfor
oldpp1=pp1
oldpp2=pp2
pp2=qq2
pp1=qq1

pp1 = pp1^(power1)
pp2 = pp2^(power2)

nn=n_elements(pp1)
sig1 = stddev(pp1)
sig2 = stddev(pp2)
if min([sig1,sig2]) le 0. then begin
print,'zero stddev',sig1,sig2
;stop
return
endif
omn1 = min(pp1,max=omx1)
omn2 = min(pp2,max=omx2)
if smooth gt 0 then begin
    omn1-=2*smooth
    omn2-=2*smooth
    omx1+=2*smooth
    omx2+=2*smooth
endif

bs1=sig1/subsamp
bs2=sig2/subsamp
nb1 = fix((omx1-omn1)/bs1)+1
nb2 = fix((omx2-omn2)/bs2)+1
;h1=(histogram(pp1,binsize=sig1/5.,omin=omn1,omax=omx1,REVERSE_INDICES=ri1))
;h2=(histogram(pp2,binsize=sig2/5.,omin=omn1,omax=omx2,REVERSE_INDICES=ri2))
;nb1=n_elements(h1)
;nb2=n_elements(h2)
bins1 = (findgen(nb1)+.5)*sig1/subsamp +omn1
bins2 = (findgen(nb2)+.5)*sig2/subsamp +omn2
nhits = dblarr(nb1,nb2)
for i=0l,nn-1 do begin
   nhits[fix((pp1[i]-omn1)/bs1),fix((pp2[i]-omn2)/bs2)]++
endfor

; now order to get significance contours
if n_elements(sigma) lt 1 then sigma = [1,2,3]

;colors = [!red,!green,!blue,!yellow,!black,!white,!skyblue,!purple]
if n_elements(sigma) le n_elements(colors) then begin
    ccolors=colors[0:n_elements(sigma)-1]

endif
ntot = double(total(nhits))
prob = double(nhits)/ntot
;possibly smooth prob
if smooth gt 0 then begin
    sl1 = smooth/(sig1/subsamp)
    sl2 = smooth/(sig2/subsamp)
    nxx = fix(6*sl1)+1
    nyy = fix(6*sl2)+1
    xx = (dindgen(nxx)-nxx/2) # (dblarr(nyy)+1.)
    yy = (dblarr(nxx)+1.) # (dindgen(nyy)-nyy/2)
    
    kernel = exp(-0.5 * ((xx/sl1)^2 + (yy/sl2)^2))
    kernel /= total(kernel)
;    prob = convol(prob,kernel,/edge_truncate)
    prob = convol(prob,kernel,/edge_zero)
    prob /= total(prob)
endif


ptes = erf(sigma/sqrt(2.))

ind = reverse(sort(prob)) ; in descending order
psort = prob(ind)
pcum=total(psort,/cum)
nline = n_elements(sigma)
levels = fltarr(nline)
for i=0,nline-1 do begin
    ii = min(where(pcum gt ptes[i],nn))
    if nn gt 0 then $
      levels[i]=psort(ii) $
    else stop
endfor
clevels = levels(sort(levels))

if keyword_set(ypos) then begin
ind = where(bins2 ge 0)
bins2=bins2(ind)
prob=prob(*,ind)
endif
;stop
;xcharsize=1.5,ycharsize=1.5,_extra=extras $
if tag_exist(extras,'useline',/top_level) then $
  contour,prob,bins1,bins2,levels=clevels,c_colors=ccolors,_extra=extras $
else $
  contour,prob,bins1,bins2,levels=clevels,/fill,c_colors=ccolors,_extra=extras 
;xcharsize=1.5,ycharsize=1.5,_extra=extras
ii=(where(nhits eq max(nhits)))[0]
ii1= ii mod nb1
ii2= ii / nb1
;print,nhits[ii1,ii2],max(nhits)
print,bins1[ii1],bins2[ii2]
if keyword_set(stopit) then stop
end 
