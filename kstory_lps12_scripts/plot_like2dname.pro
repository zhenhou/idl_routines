pro plot_like2dname,chain,params,paramname1,paramname2,nskip=nskip,$
                    stopit=stopit,sigma=sigma,scale1=scale1,scale2=scale2,$
                    fracwt=fracwt,lincombx=lincombx,lincomby=lincomby,$
                    _extra=extras

if not tag_exist(extras,/top_level,'xtitle') then $
  if n_elements(extras) ne 0 then $
  extras = create_struct('xtitle',paramname1,extras) $
else $
  extras = {xtitle:paramname1}

if not tag_exist(/top_level,extras,'ytitle') then $
  if n_elements(extras) ne 0 then $
  extras = create_struct('ytitle',paramname2,extras) $
else $
  extras = {ytitle:paramname2}


readcol,params,names,format='a'

ii1=-1 
ii2=-2

if n_elements(lincombx) gt 0 then begin
    nc=n_elements(lincombx.names)
    linind=intarr(nc)
    for i=0,nc-1 do begin
        ind = where(lincombx.names[i] eq names,ni)
        if ni ne 1 then begin
            print,'found ',ni,' matches when looking for ',lincombx.names[i], ' in file ',params
            print,'list was: ',names
;    stop                                                                       
            return
        endif
        linind[i] = ind[0]+2
    endfor

    lcombx={names:lincombx.names,inds:linind,factors:lincombx.factors}
endif    else begin
    ind = where(paramname1 eq names,ni)
    if ni ne 1 then begin
        print,'found ',ni,' matches when looking for ',paramname1, ' in file ',params
        print,'list was: ',names
;    stop
        return
    endif
    ii1 = ind[0]+2
print,paramname1,' found at chain column ',ii1
endelse
if n_elements(lincomby) gt 0 then begin
    nc=n_elements(lincomby.names)
    linind=intarr(nc)
    for i=0,nc-1 do begin
        ind = where(lincomby.names[i] eq names,ni)
        if ni ne 1 then begin
            print,'found ',ni,' matches when looking for ',lincomby.names[i], ' in file ',params
            print,'list was: ',names
;    stop                                                                       
            return
        endif
        linind[i] = ind[0]+2
    endfor

    lcomby={names:lincomby.names,inds:linind,factors:lincomby.factors}

endif    else begin
    ind = where(paramname2 eq names,ni)
    if ni ne 1 then begin
        print,'found ',ni,' matches when looking for ',paramname2, ' in file ',params
        print,'list was: ',names
;    stop
        return
    endif

ii2 = ind[0]+2
print,paramname2,' found at chain column ',ii2
endelse

if tag_exist(/top_level,extras,'priorname') then begin
    pname=extras.priorname
    n = n_elements(pname)
    pall = fltarr(3,n)
    pall[1,*]=extras.priormin
    pall[2,*]=extras.priormax
    for i=0,n-1 do begin
        ind = where(pname[i] eq names,ni)
        if ni lt 1 then begin
            print,'missing prior var ',pname
            print,names
            return
        endif
        jj=ind[0]+2
        pall[0,i]=jj
    endfor
    extras = create_struct('priorrng',pall,extras)
endif



plot_like2d,chain,ii1,ii2,nskip=nskip,stopit=stopit,$
  sigma=sigma,scale1=scale1,scale2=scale2,fracwt=fracwt,$
  lincombx=lcombx,lincomby=lcomby,$
  _extra=extras


end
