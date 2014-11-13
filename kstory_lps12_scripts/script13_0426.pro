;;;
; NAME: script13_0426.pro
; PURPOSE:
;   Make a fields plot
;
; NOTES:
;;;

;;;;;;;;;;;;;;;;;;;
; Make a triangle plot
;;;;;;;;;;;;;;;;;;;
PRO triangle_plot

; Set up the environment for ps plotting, get colors
setup_ps_plotting, pscolors=pscolors

;----------------------
; Get the data

cdir = '/data23/kstory/lps12/chains/'
chains  = ['c1_lcdm_pico_w7_newKp','c27_lcdm_camb_s12tau_newKp','c2_lcdm_pico_w7s12_newKp']



; S12-only
dir_s12   = cdir+'c27_lcdm_camb_s12tau_newKp/chains/'
files_s12 = file_search(dir_s12+'c27_lcdm_camb_s12tau_newKp*.txt')
pname_s12 = dir_s12+'c27_lcdm_camb_s12tau_newKp.paramnames'

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7_newKp/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7_newKp*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7_newKp.paramnames'

; S12+wmap
dir_w7s12 = cdir+'c2_lcdm_pico_w7s12_newKp/chains/'
files_w7s12 = file_search(dir_w7s12+'c2_lcdm_pico_w7s12_newKp*.txt')
pname_w7s12 = dir_w7s12+'c2_lcdm_pico_w7s12_newKp.paramnames'


; Params
p_names = ['omegabh2', 'omegadmh2', 'theta_s', 'tau', 'ns', '1e-9As', 'A_Pois', 'A_Clust', 'A_tSZ']
np = n_elements(p_names)

pst = {n:'',range:[0.,1.],scale:1.,tiv:1.}
p = replicate(pst,np)
p[0].n=p_names[0] & p[0].range=[1.8, 2.92]    & p[0].scale=100 & p[0].tiv=1.
p[1].n=p_names[1] & p[1].range=[0.075, 0.139] & p[1].scale=1.  & p[1].tiv=0.02
p[2].n=p_names[2] & p[2].range=[1.028, 1.053] & p[2].scale=100
p[3].n=p_names[3] & p[3].range=[0.041, 0.149]
p[4].n=p_names[4] & p[4].range=[0.81, 1.05]
p[5].n=p_names[5] & p[5].range=[1.81, 2.49]
p[6].n=p_names[6] & p[6].range=[5,35]
p[7].n=p_names[7] & p[7].range=[0,20]
p[8].n=p_names[8] & p[8].range=[0,20]

; ;TESTING
; i=1
; p_names = p_names[0:i]
; p = p[0:i]
; np = n_elements(p_names)

;Define the window size
empty = [' ']
for i=1,59 do empty = [empty,' ']
xb = 7;8
csize = 0.8
yb = xb

filename_stub = '~kstory/lps12/scripts/figs/triangle_lcdm'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

;Define the axis size relative to the window size
x0=.05
ddx=(1-2*x0)/np
y0=.05
ddy=(1-2*y0)/np
yt = y0+ddy*np

; Define the output plot and font sizes
axlCthick = 3
axlCsize  = 0.9;0.8
ytxt = [-0.23, .5] ; y position of [x,y]-axis label
;bigspt = 1.4
;thickspt = 1
lfont = 0

;*******************************
; colors; use setup_ps_plotting
color_w7s12 = pscolors[0] ; black
color_s12   = pscolors[5] ; forest
color_w7    = pscolors[10] ; purple

ls      = [3,2,0] ; [SPT-only, WMAP-only, SPT+WMAP]
mycolor = [color_s12, color_w7, color_w7s12] ; [SPT-only, WMAP-only, SPT+WMAP]


!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

; axis properties
!y.crange=[0,1]
ytitle='!X!3Likelihood!X!N'

print, '************ start loop'
for i=0, np-1 do begin ; row, top to bottom
    !y.window=[yt-(i+1)*ddy,yt-i*ddy]

    for k=0,i do begin ; column
        !x.window=[x0+k*ddx,x0+(k+1)*ddx]
        
        !x.crange=p[k].range
        !y.crange=p[i].range
        if k eq i then !y.crange=[0,1]

        print, "*** ", i,k, !x.crange, !y.crange
        ; X-axis
        case i of
            np-1: axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2;,XTICKINTERVAL=p[i].tiv ;bottom-row
            else: axis,xaxis=0,xtickname=empty,/save,xstyle=1,xthick=2
        endcase
        axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2

        ; Y-axis
        case k of
            0: axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 ;,YTICKINTERVAL=p[i].tiv  ;left-column
            else: axis,yaxis=0,ytickname=empty,/save,ystyle=1,ythick=2
        endcase                
        axis,yaxis=1,ytickname=empty,ystyle=1,ythick=2

        ;-------------------
        ; Plots
        ;-------------------
        oplot,[0,0],[0,0],color=pscolors[0]

        ; testing
        xvec = !x.crange[0] + findgen(10)*(!x.crange[1]-!x.crange[0])/10.
        yvec = !y.crange[0] + findgen(10)*(!y.crange[1]-!y.crange[0])/10.
        
        if k eq i then begin    ; diagonal
            oplot, xvec, yvec
        endif else begin
            plot_like2dname,files_s12,pname_s12,p[k].n,p[i].n,scale1=p[k].scale,scale2=p[i].scale,$
              /overplot,sigma=[1,2],$
              icolors=[pscolors[8],pscolors[3]],subsamp=3,yl=0
              ;smooth=0.01,
        endelse
        


;         if k eq i then begin ; on diagonal
;             plot_like1dname,files_s12,pname_s12,pst[i].n,subsamp=3,nskip=1000,scale=pst[i].scale,$
;               thick=4,linestyle=ls[0], color=mycolor[0]
            
    endfor
endfor
       
ps_close

print, 'output file: '+filename_stub+'.pdf'
spawn,'epstopdf '+filename_stub+'.ps'

stop
END

