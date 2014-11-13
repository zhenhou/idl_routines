;;;
; NAME: study_1hz_bolos
; PURPOSE:
;   Plotting program to study bad 1hz bolos
;
; CALLING SEQUENCE: study_1hz_bolos, 17
;
; INPUTS: 
;   field index,        index in lps12_fieldstruct()
;   yr,                 yrange
;   xr,                 xrange
;
; OUTPUTS:
;   saves a plot if keyword 'save_plot' is set
;
; NOTES: (from RK)
;   Expects a save file at: '/data/kstory/projects/lps12/1hzBolo/lines_'+field_name+'.sav'
;
; MODIFICATION HISTORY:
;   02/08/2012: (KTS) Created from /data/rkeisler/lines/study3.pro
;   04/16/2012: (KTS) get bad1hzbolos list from function
;;;

pro study_1hz_bolos, field_idx, stopit=stopit, xr=xr, yr=yr, chars=chars, lchars=lchars_in, save_plot=save_plot
compile_opt IDL2, HIDDEN

field_arr = lps12_fieldstruct()
fitsdir = field_arr[field_idx].idf_dirs[0]
field_name = field_arr[field_idx].name

restore, '/data/kstory/projects/lps12/1hzBolo/lines_'+field_name+'.sav'
;restore, '/data/kstory/projects/lps12/1hzBolo/lines_ra5h30dec-55.sav'
nfreq = n_elements(width[0,*,0])
nobs = n_elements(date)

; explode bolo flag array and bolowts array
bigflag = power*0.
bigbolowts = power*0.
wh_inf = where(~finite(bolowts),ninf)
if ninf gt 0 then bolowts[wh_inf]=0.
for i=0,nfreq-1 do begin
    bigflag[*,i,*] = flag
    bigbolowts[*,i,*] = bolowts
endfor
; multiply by bolowts
power *= bigbolowts

; take into account the fact that we'll be lpf in the low ell
; analysis.
tf = exp(-(center/5.)^6. -(center/7.5)^6.)
power *= tf

; explode ptc array
uniqptc = uniqval(ptc)
nuniqptc = n_elements(uniqptc)
bigptc = power*0.
for i=0,nobs-1 do begin
    bigptc[i,*,*] = (ptc[i])[0]
endfor

;bad1hzbolos_rk = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768]

;; 2008
if( intersect(field_idx, [0,1]) ne -1 ) then begin
    bad1hzbolos = bad_1hz_bolo_array(2008)
;; 2009, 2010, or 2011
endif else begin
    bad1hzbolos = bad_1hz_bolo_array(2009)
endelse

nbad1hz = n_elements(bad1hzbolos)
bigflag2   = bigflag ; bad1hzbolos
for i=0,nbad1hz-1 do begin
    wh=where(BOLOID eq bad1hzbolos[i],nwh)
    if nwh ne 1 then begin
        print,'uh oh'
        stop
    endif
    bigflag2[*,*,wh[0]] += (2.^14.)
endfor

;;; RK bolo set
; bigflag2_rk = bigflag ; bad1hzbolos_rk
; for i=0, n_elements(bad1hzbolos_rk)-1 do begin
;     wh=where(BOLOID eq bad1hzbolos_rk[i],nwh)
;     bigflag2_rk[*,*,wh[0]] += (2.^14.)
; endfor

;;; Find set of good bolometers
bigok = 1.*(bigflag eq 0)   ; 1 for no bolo flags
bigok2 = 1.*(bigflag2 eq 0) ; 1 for no bolo flags and not in bad1hzbolos
;bigok2_rk = 1.*(bigflag2_rk eq 0) ; 1 for no bolo flags and not in bad1hzbolos

; get per-bolo stats
hump = fltarr(nbolos)  ; power between freq = [0.5, 0.9]
inthz = fltarr(nbolos) ; power at integer frequencies??
broad = fltarr(nbolos) ; power in [f_ptc+0.2, 7]
for i=0,nbolos-1 do begin
    this_power = reform((power*bigok)[*,*,i])
    this_freq = reform((center)[*,*,i])
    
    wh_hump = where(this_freq ge 0.5 and this_freq le 0.9,nhump)
    if nhump gt 0 then hump[i] = total(this_power[wh_hump])

    wh_inthz = where(abs(this_freq mod 1.0) lt 0.05,ninthz)
    if ninthz gt 0 then inthz[i] = total(this_power[wh_inthz])

    wh_broad = where(this_freq gt (max(uniqptc)+0.2) and this_freq le 7.0,nbroad)
    if nbroad gt 0 then broad[i] = total(this_power[wh_broad])
endfor

whuse = where(broad le median(broad),compl=whnotuse)
bad_broad_bolos = boloid[whnotuse]
save,bad_broad_bolos,filename='/data/kstory/projects/lps12/1hzBolo/bad_broad_bolos_'+field_name+'.sav'

;;; Broad noise
bigflag3 = bigflag
for i=0,n_elements(whnotuse)-1 do begin
    j=whnotuse[i]
    bigflag3[*,*,j] = 2.^15.
endfor
bigok3 = 1.*(bigflag3 eq 0)
;bigok2 = bigok3;tempp


if n_elements(xr) eq 0 then xr=[0,7]
if n_elements(yr) eq 0 then yr=[1,1e6]
;if n_elements(yr) eq 0 then yr=[1e-6,1]
plot,center,power*bigok,ps=3,yr=yr,/yst,xr=xr,/xst, $
  xtitle='Frequency (Hz)',ytitle='Line Power (width*height)', $
  title=field_name,/yl,chars=chars,/nodata

;;; Plot 1hz lines
y=findgen(10)
inthz_color = !red
ptc_color = !darkgreen
for i=0,n_elements(y)-1 do oplot,[1,1]*y[i],[1e-10,1e10],lines=2,color=inthz_color

;;; Plot the PTC line
for i=0,nuniqptc-1 do begin
    oplot,[1,1]*uniqptc[i],[1e-10,1e10],lines=2,color=ptc_color
    oplot,[1,1]*uniqptc[i]*2.,[1e-10,1e10],lines=2,color=ptc_color
    oplot,[1,1]*uniqptc[i]*3.,[1e-10,1e10],lines=2,color=ptc_color
    oplot,[1,1]*uniqptc[i]*4.,[1e-10,1e10],lines=2,color=ptc_color
    oplot,[1,1]*uniqptc[i]*5.,[1e-10,1e10],lines=2,color=ptc_color
endfor

; Make the plot legend
if n_elements(lchars_in) eq 0 and n_elements(chars) gt 0 then lchars=0.5*chars else begin
    if n_elements(lchars_in) gt 0 then lchars=lchars_in
endelse
legend,/top,/left,['Integers','PTC hrmncs','have broad noise','bad1hzbolos'],textc=[inthz_color,ptc_color,!skyblue,!red],chars=lchars,box=0

oplot,center,power*(bigok),ps=3,color=!white
oplot,center,power*(bigok3),ps=3,color=!skyblue
;oplot,center,power*(bigok-bigok2_rk),ps=3,color=!orange
oplot,center,power*(bigok-bigok2),ps=3,color=!red


if 0 then begin

; let's bin in frequency and create power maps.
dhz = 0.005
maxhz = 6.
nhz = maxhz/dhz
hz = findgen(nhz)*dhz
binpower = fltarr(nhz)
binpower2 = fltarr(nhz)
wh=where(center le maxhz)
pf_center = center[wh]
pf_power_bigok = (power*bigok)[wh]
pf_power_bigok2 = (power*bigok2)[wh]
for i=0,nhz-1 do begin
;    print,strtrim(i,2),'/',strtrim(nhz-1,2)
    wh=where(pf_center ge (hz[i]-dhz) and pf_center lt (hz[i]+dhz),nwh)
    if nwh eq 0 then continue
    binpower[i] += total((pf_power_bigok)[wh])
    binpower2[i] += total((pf_power_bigok2)[wh])
endfor
plot,hz,binpower,/yl,yr=[1e-6,10],xr=[0,6],ps=3
oplot,hz,binpower
endif

; Save the plot?
if keyword_set(save_plot) then begin
    figdir = '/home/kstory/lps12/figs/'
    filename = figdir+'boloPower_0416_'+field_name
    err = tvread(/png, filename=filename, /nodialog)
endif

if keyword_set(stopit) then stop
end

