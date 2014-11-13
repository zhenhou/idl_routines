;;;
; NAME: script_0724
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) test plots
;
; MODIFICATION HISTORY:
;  07/23/2012: (KTS) Created
;;;

PRO p1_alens
dir = '/data23/hou/lps12/paramfits/chains_final/c24_lcdm_alens_camb_w7s12/chains/'
files = file_search(dir+'c24_lcdm_alens_camb_w7s12*.txt')
pname=dir+'c24_lcdm_alens_camb_w7s12.paramnames'
;n=n_elements(files)

;;; plotting options
cc=[!blue, !red, !black]
ls=[0,2]
xb = 7
csize = xb/7.
yb = 2.5
csize=2
lfont=1
loadct,39,ncolors=12
!p.multi=0

;;;;;;;;;;;;
; The plot
;;;;;;;;;;;;
window,0
oplot,[0,0],[0,0],color=3
i=0
plot_like1dname,files,pname,'alens',nskip=1000,subsamp=3,$
  thick=4,linestyle=ls[0], xr=[-1,5],yr=[0,1.],/xstyle,xtitle=textoidl('A_L');, color=cc[i],/overplot
plot_like1dname,files,pname,'alens',subsamp=3,thick=4,linestyle=ls[i], color=cc[i], /oplot
xyouts,2.5,0.85,'SPT+WMAP',charsize=csize*1.3, color=cc[0],font=lfont

i=1 ;& plot_like1dname,files[i],pname[i],'alens',subsamp=3,thick=4,linestyle=ls[i], color=cc[i], /oplot
xyouts,2.5,0.75,'WMAP',charsize=csize*1.3, color=cc[i],font=lfont

; vertical line
oplot, [1.,1.], [0,1.2], linestyle=2

;!p.multi=0
fdir='/home/kstory/lps12/scripts/figs/'
ff=fdir+'alense_0724' & fname=ff
err = tvread(/png,/nodialog,filename=fname)
spawn, 'convert '+fname+' '+ff+'.pdf '
stop
END


;--------------------------------
; Plot r, r v.s. ns
;--------------------------------
PRO p2_r
dir = '/data23/hou/lps12/paramfits/chains_final/c35_lcdm_r_camb_w7s12_BAO_SNe/chains/'
files = file_search(dir+'c35_lcdm_r_camb_w7s12_BAO_SNe*.txt')
pname=dir+'c35_lcdm_r_camb_w7s12_BAO_SNe.paramnames'


cc=[!blue, !red]
ls=[0,2]
xb = 7
csize = xb/7.
yb = 2.5
csize=2
lfont=1
!p.multi=[0,2,1]

window, 0, xsize=900,ysize=400

plot_like1dname,files,pname,'r',subsamp=3,thick=4,linestyle=ls[0], xr=[0,0.5],yr=[0,1.],/xstyle,xtitle='r'
i=0 & plot_like1dname,files,pname,'r',subsamp=3,thick=4,linestyle=ls[i], color=cc[i], /oplot
xyouts,0.3,0.85,'SPT+WMAP',charsize=csize*1.3, color=cc[0],font=lfont
xyouts,0.3,0.75,'WMAP',charsize=csize*1.3, color=cc[1],font=lfont

plot_like2dname,files,pname,'ns','r',xtitle='n_s',ytitle=textoidl('r');,icolors=[3,8,10] ;,smooth=0.15,priorname=['czero_ksz'],priormin=[0],priormax=[100]
xyouts,0.94,0.35,'SPT+WMAP',charsize=csize*1.3, color=cc[0],font=lfont


fdir='/home/kstory/lps12/scripts/figs/'
ff=fdir+'r_ns_0724' & fname=ff
err = tvread(/png,/nodialog,filename=fname)
spawn, 'convert '+fname+' '+ff+'.pdf '
stop
END


;--------------------------------
; Plot r, r v.s. ns
;   - Try using RK's stuff, definitely not working yet
;--------------------------------
PRO p2p1_r
dir = '/data23/hou/lps12/paramfits/chains_final/c35_lcdm_r_camb_w7s12_BAO_SNe/chains/'
files = file_search(dir+'c35_lcdm_r_camb_w7s12_BAO_SNe*.txt')
pname=dir+'c35_lcdm_r_camb_w7s12_BAO_SNe.paramnames'

;;; plotting options
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 6.5
csize = xb/7.
yb = xb *.88/.9*(.44/.9)

ps_open,'~kstory/lps12/scripts/figs/tmp',/color,xsize=xb,ysize=yb,/inches,/portrait
;window, 1, xsize=900, ysize=300
loadct,39,ncolors=12
ox = !x.window
cx = !x.crange
oy = !y.window
cy = !y.crange
;.11 - .99 = .44 each 
;.07 - .97  = .9 each
x0=.11
ddx=.44
y0=.07
ddy=.9

xt2='!17!12l!17 '
yt2='!13D!I!12l!N!17  (!7l!17K!E2!N)!17'

label='Modified BB'
label2='tSZ-CIB correl.'
k=0


cc=[!blue, !red]
ls=[0,2]
xb = 7
csize = xb/7.
yb = 2.5
csize=2
lfont=1

;;;;;;;;;;;;
; The plot
;;;;;;;;;;;;
for i=0,1 do begin
;for i=0,2 do for j=0,1 do begin

    !x.window=[x0+i*ddx,x0+(i+1)*ddx]        ;*xb
    !x.crange=[0,7.5]
    
    !y.window=[y0,y0+ddy]        ;*yb
    !y.crange=[0,10.5]
    


    axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2 
    axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2
    
    if i eq 0 then $
      axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 $
    else $
      axis,yaxis=0,ytickname=empty,ystyle=1,/save,ythick=2

    axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2


;    if (j eq 0 and i eq 1) then    xyouts,5000,.02,xt2,alignment=0.5
    if ( i eq 0) then  begin
        xyouts,max(!x.crange)*0.9,-1.7,alignment=0.5,'D!S!D3000!R!N!UtSZ  !N',font=0,charsize=1.2
        xyouts,max(!x.crange)*1.1,-1.7,alignment=0.5,'(!7l!XK!U2!N)',charsize=1.2
        xyouts,-.9,4.4,alignment=0.5,orientation=90,'D!S!D3000!R!N!UkSZ  !N',font=0,charsize=1.2
       xyouts,-.9,6.5,alignment=0.5,orientation=90,'(!7l!XK!U2!N)',charsize=1.2
    endif
    oplot,[0,0],[0,0],color=3

    if i eq 0 then begin
        plot_like1dname,files,pname,'r',subsamp=3,thick=4,linestyle=ls[0], xr=[0,0.5],yr=[0,1.],/xstyle,xtitle='r'
        i=0 & plot_like1dname,files[i],pname[i],'r',subsamp=3,thick=4,linestyle=ls[i], color=cc[i], /oplot
        xyouts,0.3,0.85,'SPT+WMAP',charsize=csize*1.3, color=cc[0],font=lfont
        xyouts,0.3,0.75,'WMAP',charsize=csize*1.3, color=cc[1],font=lfont
    endif else begin

   ; plot 2
        plot_like2dname,files,pname,'r','ns',/overplot,xtitle='',ytitle='',icolors=[3,8,10] ;,smooth=0.15,priorname=['czero_ksz'],priormin=[0],priormax=[100]
        xyouts,7,8.5,label,font=0,alignment=1
    endelse
endfor
xyouts,7.3,7.7,label2,font=0,alignment=1

ps_close

!p.multi=0
stop

END
