;;;
; NAME: script_0802
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Alens plot
;;;

PRO p1_alens

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

; S12
dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/'
files_s12 = [file_search(dir_s12+'c24_lcdm_alens_camb_w7s12*.txt')]
pname_s12 = dir_s12+'c24_lcdm_alens_camb_w7s12.paramnames'

;n=n_elements(files)

; WMAP
dir_w = '/data23/hou/lps12/paramfits/chains_final/lcdm_alens_camb_w7/chains/'
files_w = file_search(dir_w+'lcdm_alens_camb_w7*.txt')
pname_w = dir_w+'lcdm_alens_camb_w7.paramnames'

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait
loadct,39,ncolors=12

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1700 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.
thickspt = 1
lfont = 0
ls = [0,2]

;*******************************
; colors
color_this_work = pscolors[8] ; blue
color_wmap = pscolors[2] ; red

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[-1, 5]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, -2.1, 0.45, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 1.8, -0.16, charsize=1.5, charthick=3, font=0, '!3A!DL!X!N'

oplot,[0,0],[0,0],color=3

; S12
plot_like1dname,files_s12,pname_s12,'alens',subsamp=3,thick=7,linestyle=ls[0], color=color_this_work, /oplot
xyouts,2,0.8,'SPT+WMAP',charsize=csize*1.3, color=color_this_work,font=lfont

; WMAP
i=1
plot_like1dname,files_w,pname_w,'alens',subsamp=3,thick=4,linestyle=ls[i], color=color_wmap, /oplot
xyouts,2,0.67,'WMAP',charsize=csize*1.3, color=color_wmap,font=lfont

; vertical line
oplot, [1.,1.], [0,1], linestyle=2, thick=5

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/alens_0802.pdf'
stop
END




;;;;;;;;;;;;;;;;;;;
; WMAP only, for testing
PRO p1_alens_w

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

; S12
dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/'
files_s12 = file_search(dir_s12+'c24_lcdm_alens_camb_w7s12*.txt')
pname_s12 = dir_s12+'c24_lcdm_alens_camb_w7s12.paramnames'

; WMAP
dir_w = '/data23/hou/lps12/paramfits/chains_final/lcdm_alens_camb_w7/chains/'
files_w = file_search(dir_w+'lcdm_alens_camb_w7*.txt')
pname_w = dir_w+'lcdm_alens_camb_w7.paramnames'

; K11
dir_k11   = '/data17/rkeisler/ps09/chains_20110202/'
files_k11 = file_search(dir_k11+'chain_k10_alens_1.txt')
pname_k11 = dir_k11+'chain_k10_alens.paramnames'


;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 5
csize = 1.3
;yb = xb *.88/.9*(.44/.9)
yb = xb * 0.9

filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.20
ddx=.75
y0=.16
ddy=.77

; Define the output plot and font sizes
thick=3 ; points
xtxt=1700 ; xyouts label
ytxt=2300
dytxt = 0.63
bigspt = 1.
thickspt = 1
lfont = 0
ls = [0,2]

;*******************************
; colors
color_this_work = pscolors[8] ; blue
color_wmap = pscolors[2] ; red
color_k11 = pscolors[5] ; red

!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

;-----------------
; The plot
!x.window=[x0,x0+ddx]
!x.crange=[-1, 3]
!y.window=[y0,y0+ddy]
!y.crange=[0, 1]

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

; axes labels
xyouts, -2.1, 0.45, alignment=0.5, orientation=90, charsize=1.5, font=0, charthick=3, 'Likelihood'
xyouts, 1.8, -0.16, charsize=1.5, charthick=3, font=0, '!3A!DL!X!N'

oplot,[0,0],[0,0],color=3

; S12
plot_like1dname,files_s12,pname_s12,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=ls[0], color=color_this_work, /oplot
xyouts,2,0.9,'SPT+WMAP',charsize=csize, color=color_this_work,font=lfont

; WMAP
i=1
plot_like1dname,files_w,pname_w,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=ls[i], color=color_wmap, /oplot
xyouts,2,0.8,'WMAP',charsize=csize, color=color_wmap,font=lfont

; K11
plot_like1dname,files_k11,pname_k11,'alens',subsamp=3,nskip=1000,$
  thick=7,linestyle=3, color=color_k11, /oplot
xyouts,2,0.7,'K11',charsize=csize, color=color_K11,font=lfont

; vertical line
oplot, [1.,1.], [0,1], linestyle=1, thick=3

;*******************************

; Save the plot
ps_close

spawn,'epstopdf '+filename_stub+'.ps'
spawn, 'cp figs/tmp.pdf figs/alens_w_0802.pdf'
stop
END




;----------------------------
; ns constraints
PRO ns_constraint
nskip = 5000

dataset = 3
case dataset of 

    ; CMB only
    0: begin
        print, "---------- CMB only ------------------"
        dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
        files_s12 = file_search(dir_s12+'c4_lcdm_camb_w7s12_*.txt')
        pname_s12 = dir_s12+'c4_lcdm_camb_w7s12.paramnames'
        plot_like1dname,files_s12,pname_s12,'ns',subsamp=3,nskip=nskip
    endcase 
    
    ; CMB+H0
    1: begin
        print, "---------- CMB+H0 ------------------"
        ;dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
        ;files_s12 = file_search(dir_s12+'c4_lcdm_camb_w7s12_*.txt')
        ;pname_s12 = dir_s12+'c4_lcdm_camb_w7s12.paramnames'
        ;plot_like1dname,files_s12,pname_s12,'ns',subsamp=3,nskip=1000
    endcase

    ; CMB+BAO
    2: begin
        print, "---------- CMB+BAO ------------------"
        dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c11_lcdm_camb_w7s12_BAO/chains/'
        files_s12 = file_search(dir_s12+'c11_lcdm_camb_w7s12_BAO_*.txt')
        pname_s12 = dir_s12+'c11_lcdm_camb_w7s12_BAO.paramnames'
        plot_like1dname,files_s12,pname_s12,'ns',subsamp=3,nskip=1000
    endcase

    ; CMB+H0+BAO
    3: begin
        print, "---------- CMB+H0+BAO ------------------"
        dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c15_lcdm_camb_w7s12_BAO_H0/chains/'
        files_s12 = file_search(dir_s12+'c15_lcdm_camb_w7s12_BAO_H0_*.txt')
        pname_s12 = dir_s12+'c15_lcdm_camb_w7s12_BAO_H0.paramnames'
        plot_like1dname,files_s12,pname_s12,'ns',subsamp=3,nskip=1000
    endcase
end

stop
END


;----------------------------
; ns constraints
PRO r_constraint, dataset
nskip = 1000

case dataset of 

    ; CMB only
    0: begin
        print, "---------- CMB only ------------------"
        dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c26_lcdm_r_camb_w7s12/chains/'
        files_s12 = file_search(dir_s12+'c26_lcdm_r_camb_w7s12_*.txt')
        pname_s12 = dir_s12+'c26_lcdm_r_camb_w7s12.paramnames'
        plot_like1dname,files_s12,pname_s12,'r',subsamp=3,nskip=1000
        print, "****************** ns ********************"
        plot_like1dname,files_s12,pname_s12,'ns',subsamp=3,nskip=1000
    endcase 
    
    ; CMB+H0
    1: begin
        print, "---------- CMB+H0 ------------------"
;         dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c54_lcdm_r_camb_w7s12_H0/chains/'
;         files_s12 = file_search(dir_s12+'c54_lcdm_r_camb_w7s12_H0_*.txt')
;         pname_s12 = dir_s12+'c54_lcdm_r_camb_w7s12_H0.paramnames'
;         plot_like1dname,files_s12,pname_s12,'r',subsamp=3,nskip=1000
;         print, "****************** ns ********************"
;         plot_like1dname,files_s12,pname_s12,'ns',subsamp=3,nskip=1000
    endcase

    ; CMB+BAO
    2: begin
        print, "---------- CMB+BAO ------------------"
        dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c28_lcdm_r_camb_w7s12_BAO/chains/'
        files_s12 = file_search(dir_s12+'c28_lcdm_r_camb_w7s12_BAO_*.txt')
        pname_s12 = dir_s12+'c28_lcdm_r_camb_w7s12_BAO.paramnames'
        plot_like1dname,files_s12,pname_s12,'r',subsamp=3,nskip=1000
        print, "****************** ns ********************"
        plot_like1dname,files_s12,pname_s12,'ns',subsamp=3,nskip=1000
    endcase

    ; CMB+H0+BAO
    3: begin
        print, "---------- CMB+H0+BAO ------------------"
        dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c35_lcdm_r_camb_w7s12_BAO_SNe/chains/'
        files_s12 = file_search(dir_s12+'c35_lcdm_r_camb_w7s12_BAO_SNe_*.txt')
        pname_s12 = dir_s12+'c35_lcdm_r_camb_w7s12_BAO_SNe.paramnames'
        plot_like1dname,files_s12,pname_s12,'r',subsamp=3,nskip=1000
        print, "****************** ns ********************"
        plot_like1dname,files_s12,pname_s12,'ns',subsamp=3,nskip=1000
    endcase
end

stop
END



;----------------------------
; LCDM constraints
PRO lcdm_constraint
nskip = 5000
subsamp=1000

dataset = 0
case dataset of 

    ; CMB only
    0: begin
        print, "---------- CMB only ------------------"
        dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c4_lcdm_camb_w7s12/chains/'
        files_s12 = file_search(dir_s12+'c4_lcdm_camb_w7s12_*.txt')
        pname_s12 = dir_s12+'c4_lcdm_camb_w7s12.paramnames'

;         plot_like1dname,files_s12,pname_s12,'omegabh2',subsamp=subsamp,nskip=nskip
;         plot_like1dname,files_s12,pname_s12,'omegadmh2',subsamp=subsamp,nskip=nskip
;         plot_like1dname,files_s12,pname_s12,'omegal*',subsamp=subsamp,nskip=nskip
;         plot_like1dname,files_s12,pname_s12,'logA',subsamp=subsamp,nskip=nskip
;         plot_like1dname,files_s12,pname_s12,'ns',subsamp=subsamp,nskip=nskip
;         plot_like1dname,files_s12,pname_s12,'tau',subsamp=subsamp,nskip=nskip

        print, '- - - - Derived Parameters - - - -'
        ;Derived parameters
        plot_like1dname,files_s12,pname_s12,'H0*',subsamp=subsamp,nskip=nskip
        plot_like1dname,files_s12,pname_s12,'sigma8*',subsamp=subsamp,nskip=nskip
        plot_like1dname,files_s12,pname_s12,'z_eq',subsamp=subsamp,nskip=nskip
        plot_like1dname,files_s12,pname_s12,'theta_s',subsamp=subsamp,nskip=nskip
        plot_like1dname,files_s12,pname_s12,'zrei*',subsamp=subsamp,nskip=nskip
        ;plot_like1dname,files_s12,pname_s12,'Dv/rs?',subsamp=subsamp,nskip=nskip

    endcase 
end

stop    
END



;;;;
; Get Alens^0.6
;;;;
pro plot_like1d_0p6,files,nparam,nskip=nskip,stopit=stopit,subsamp=subsamp,offset=offset,results=results,_extra=extras,lincomb=lcomb
;,xtitle=xtitle,scale=scale,title=title,ytitle=ytitle,cz=cz,oplot=oplot,color=color,linestyle=linestyle,psym=psym,stopit=stopit,subsamp=subsamp,offset=offset,results=results,_extra=extras


if not tag_exist(/top_level,extras,'zeropad') then $
  zeropad=0 $
else $
  zeropad=extras.zeropad

if not tag_exist(/top_level,extras,'temp') then $
  temp=1 $
else $
  temp=extras.temp

if not tag_exist(/top_level,extras,'xtitle') then $
  if n_elements(extras) ne 0 then $
  extras = create_struct('xtitle',strtrim(string(fix(nparam)),2),extras) $
  else $
  extras = {xtitle:strtrim(string(fix(nparam)),2)}
if not tag_exist(/top_level,extras,'ytitle') then $
    extras = create_struct('ytitle','Likelihood',extras) 

;if n_elements(xtitle) ne 1 then xtitle='param_'+strtrim(string(fix(nparam)),2)
;if n_elements(ytitle) ne 1 then ytitle='Likelihood'

;if n_elements(subsamp) ne 1 then subsamp=5.
if n_elements(subsamp) ne 1 then subsamp=10000.
if not tag_exist(/top_level,extras,'scale') then $
  scale=1.0 $
else $
  scale=extras.scale
;print,scale

if n_elements(nskip) ne 1 then nskip=4000
if n_elements(offset) ne 1 then offset =0.

a=read_ascii(files[0])
if n_elements(lcomb) gt 0 then begin
    pp=0
    for j=0,n_elements(lcomb.inds)-1 do begin
        pp += lcomb.factors[j]*a.field01[lcomb.inds[j],nskip:*]*scale
    endfor
endif else $
  pp = a.field01[nparam,nskip:*]*scale
wt = long(a.field01[0,nskip:*])

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
    if n_elements(lcomb) gt 0 then begin
        pp2=0
        for j=0,n_elements(lcomb.inds)-1 do begin
            pp2 += lcomb.factors[j]*b.field01[lcomb.inds[j],nskip:*]*scale
        endfor
    endif else $
      pp2 = b.field01[nparam,nskip:*]*scale
    wt2 = long(b.field01[0,nskip:*])

    if  tag_exist(/top_level,extras,'priorrng') then begin
        for k=0,nprip-1 do begin
            prp = b.field01[prip[k],nskip:*]
            ind = where(prp lt primn[k] or prp gt primx[k],nbad)
            if nbad gt 0 then wt2(ind) = 0
        endfor
    endif
    
    pp = [reform(pp),reform(pp2)]
    wt = [reform(wt),reform(wt2)]
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

; set pp = pp^0.6
pp = pp^0.6

sig = stddev(pp)
if (sig le 0) or (finite(sig) eq 0) then begin
    print,'data constant, no histogram possible'
    results={fail:1}
    return
endif
print,sig
h=double(histogram(pp,binsize=sig/subsamp,omin=omn,omax=omx))
if n_elements(h) eq 1 then begin
    print,'data constant2, no histogram possible'
    results={fail:1}
    return
endif

h/=max(h)
h = h^temp

bins = (findgen(n_elements(h))+.5)*sig/subsamp +omn

if zeropad eq 1 then begin
    h=[0,h,0]
    bins=[min(bins)-sig/subsamp,bins,max(bins)+sig/subsamp]
endif
boplot=0
if tag_exist(extras,'oplot',/top_level) then boplot =  extras.oplot 
bnoplot=0
if tag_exist(extras,'noplot',/top_level) then bnoplot =  extras.noplot 

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
ind = where(z gt 0.95)
ind5 = where(z gt 0.975)
ind3 = where(z gt 0.67)
ind2 = where(z gt 0.05,n)
ind4 = where(z gt 0.025,n)
ind5 = where(z gt 0.975)
ind6 = where(z gt 0.159)
ind7 = where(z gt 0.841)

print,'95% limit: ',bins[ind[0]]
print,'67% limit: ',bins[ind3[0]]
uppers=[bins[ind3[0]],bins[ind[0]]]
print,'5% limit: ',bins[ind2[0]]
print,'2sigma: 2.5 - 97.5% limit: ',bins[ind4[0]],bins[ind5[0]]
print,'1sigma: 15.9 - 84.1% limit: ',bins[ind6[0]],bins[ind7[0]]

print,'2.5% limit: ',bins[ind5[0]]

print,'mean ',total(bins*h)/total(h)
ind = where(z gt .5-.683/2)
ind2 = where(z gt .5+.683/2)
ind3 = where(z gt .5)
print,bins[ind[0]],bins[ind2[0]]
print,'median ',bins(ind3[0])
print,'errs ',(bins[ind2[0]]-bins[ind[0]])/2.
results={bins:bins,prob:h,rng:[bins[ind[0]],bins[ind2[0]]],$
        err:(bins[ind2[0]]-bins[ind[0]])/2.,median:bins(ind3[0]),$
        uppers:uppers,$
        fail:0}
if keyword_set(stopit) then stop
end 


;;;;;-------------------------
; Run alens^0.6
; S12
PRO run_al6

dir_s12 = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/'
files_s12 = [file_search(dir_s12+'c24_lcdm_alens_camb_w7s12*.txt')]
pname_s12 = dir_s12+'c24_lcdm_alens_camb_w7s12.paramnames'

;plot_like1dname,files_s12,pname_s12,'alens',subsamp=10000,/stopit
END

