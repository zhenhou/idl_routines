;;;
; NAME: cdf_highT
; PURPOSE:
;   Calculate the probability to exceed some limit from cosmomc chains
;
; KEYWORDS:
;   files,            cosmomc files
;   nparam,           the column number of the parameter
;   temp,             the temperature
;   nskip,            skip burn-in period
; OPTIONAL: (pick one or the other)
;   gtx,              Calculate P(variable > gtx)
;   ltx,              Calculate P(variable < ltx)
;
; MODIFICATION HISTORY:
;  10/19/2012: (KTS) Created
;;;

PRO cdf_highT,files,nparam,temp=temp,gtx=gtx,ltx=ltx,nskip=nskip,stopit=stopit

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
