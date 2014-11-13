pro plot_like1dname,chain,params,paramname,nskip=nskip,stopit=stopit,subsamp=subsamp,offset=offset,results=results,lincomb=lincomb,_extra=extras,power=power
if not tag_exist(/top_level,extras,'xtitle') then $
  if n_elements(extras) ne 0 then $
  extras = create_struct('xtitle',paramname,extras) $
  else $
  extras = {xtitle:paramname}

if n_elements(power) eq 0 then power=1.0
readcol,params,names,format='a'
if n_elements(lincomb) gt 0 then begin
    nc=n_elements(lincomb.names)
    linind=intarr(nc)
    for i=0,nc-1 do begin
        ind = where(lincomb.names[i] eq names,ni)
        if ni ne 1 then begin
            print,'found ',ni,' matches when looking for ',lincomb.names[i], ' in file ',params
            print,'list was: ',names
;    stop                                                                                   
            return
        endif
        linind[i] = ind[0]+2
    endfor

    lcomb={names:lincomb.names,inds:linind,factors:lincomb.factors}
    ii=-1
endif    else begin

    ind = where(paramname eq names,ni)
    if ni ne 1 then begin
        print,'found ',ni,' matches when looking for ',paramname, ' in file ',params
        print,'list was: ',names
;    stop
        return
    endif
    ii = ind[0]+2
    print,paramname,' found at chain column ',ii
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

plot_like1d,chain,ii,nskip=nskip,stopit=stopit,subsamp=subsamp, $
  offset=offset,results=results,lincomb=lcomb,_extra=extras,power=power

end
