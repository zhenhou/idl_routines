;;;
; NAME: plots_s12_triangle.pro
; PURPOSE:
;   Make a triangle plot
;
; NOTES:
;;;

;;;;;;;;;;;;;;;;;;;
; Make a triangle plot
;;;;;;;;;;;;;;;;;;;
PRO triangle_plot
t0=systime(/seconds)

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

pst = {n:'',range:[0.,1.],label:'',scale:1.,tiv:1.}
p = replicate(pst,np)

; SPT-only ranges
if 1 then begin 
i=0 & p[i].n=p_names[i] & p[i].range=[1.8, 2.92]    & p[i].label='100!4X!3!Db!Nh!U2!X!N' & p[i].scale=100 & p[i].tiv=0.5
i=1 & p[i].n=p_names[i] & p[i].range=[0.075, 0.139] & p[i].label='!4X!3!Dc!Nh!U2!X!N'    & p[i].scale=1.  & p[i].tiv=0.02
i=2 & p[i].n=p_names[i] & p[i].range=[1.035, 1.053] & p[i].label='!3 100 !4h!D!3s!X!N'   & p[i].scale=100 & p[i].tiv=0.01
i=3 & p[i].n=p_names[i] & p[i].range=[0.041, 0.149] & p[i].label='!4s!X!N'               & p[i].scale=1.  & p[i].tiv=0.02
i=4 & p[i].n=p_names[i] & p[i].range=[0.81, 1.05]   & p[i].label='!3n!Ds!X!N'            & p[i].scale=1.  & p[i].tiv=0.1
i=5 & p[i].n=p_names[i] & p[i].range=[1.81, 2.49]   & p[i].label='!310!U9!N!4D!3!U2!DR!N!X' & p[i].scale=1.  & p[i].tiv=0.5
i=6 & p[i].n=p_names[i] & p[i].range=[5,35]         & p[i].label='A_Pois'                & p[i].scale=1.  & p[i].tiv=10
i=7 & p[i].n=p_names[i] & p[i].range=[0,20]         & p[i].label='A_Clust'               & p[i].scale=1.  & p[i].tiv=10
i=8 & p[i].n=p_names[i] & p[i].range=[0,20]         & p[i].label='A_tSZ'                 & p[i].scale=1.  & p[i].tiv=10
endif

; w7s12 ranges
if 0 then begin
i=0 & p[i].n=p_names[i] & p[i].range=[2.0, 2.44]    & p[i].label='100!4X!3!Db!Nh!U2!X!N' & p[i].scale=100 & p[i].tiv=0.5
i=1 & p[i].n=p_names[i] & p[i].range=[0.09, 0.135]  & p[i].label='!4X!3!Dc!Nh!U2!X!N'    & p[i].scale=1.  & p[i].tiv=0.02
i=2 & p[i].n=p_names[i] & p[i].range=[1.036, 1.048] & p[i].label='!3 100 !4h!D!3s!X!N'   & p[i].scale=100 & p[i].tiv=0.01
i=3 & p[i].n=p_names[i] & p[i].range=[0.041, 0.149] & p[i].label='!4s!X!N'               & p[i].scale=1.  & p[i].tiv=0.02
i=4 & p[i].n=p_names[i] & p[i].range=[0.91,1.02]    & p[i].label='!3n!Ds!X!N'            & p[i].scale=1.  & p[i].tiv=0.1
i=5 & p[i].n=p_names[i] & p[i].range=[1.81, 2.49]   & p[i].label='!310!U9!N!4D!3!U2!DR!N!X' & p[i].scale=1.  & p[i].tiv=0.5
i=6 & p[i].n=p_names[i] & p[i].range=[11,29]        & p[i].label='A_Pois'                & p[i].scale=1.  & p[i].tiv=5
i=7 & p[i].n=p_names[i] & p[i].range=[0,15]         & p[i].label='A_Clust'               & p[i].scale=1.  & p[i].tiv=5
i=8 & p[i].n=p_names[i] & p[i].range=[0,15]          & p[i].label='A_tSZ'                 & p[i].scale=1.  & p[i].tiv=5
endif

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
x0=.09
ddx=(1-2*x0)/np
y0=.06
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
color_s12   = pscolors[8] ; blue
color_w7    = pscolors[2] ; red

ls      = [3,2,0] ; [SPT-only, WMAP-only, SPT+WMAP]
mycolor = [color_s12, color_w7, color_w7s12] ; [SPT-only, WMAP-only, SPT+WMAP]


!p.charsize=1.
!p.color = pscolors[0]
!p.background = pscolors[1]
!p.multi = [0,1,1]
xmargin=[8.5,2]
ymargin=[4,2]

print, '************ start loop'
for i=0, np-1 do begin ; row, top to bottom
    !y.window=[yt-(i+1)*ddy,yt-i*ddy]

    for k=0,i do begin ; column
        !x.window=[x0+k*ddx,x0+(k+1)*ddx]
        
        !x.crange=p[k].range
        !y.crange=p[i].range
        if k eq i then !y.crange=[0,1]

        print, "*** ", i,k, !x.crange, !y.crange, p[i].tiv
        ; X-axis
        case i of
            np-1: begin ; bottom row
                axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=2,XTICKINTERVAL=p[k].tiv ;bottom-row
                xpos = !x.crange[0] + 0.2*(!x.crange[1]-!x.crange[0])
                ypos = !y.crange[0] - 0.5*(!y.crange[1]-!y.crange[0])
                if k eq i then ypos = (p[i].range)[0] - 0.5*((p[i].range)[1]-(p[i].range)[0])
                if k eq 0 then ypos = -1.5
                ;xyouts, xpos, ypos, alignment=0.5, charsize=0.9, font=-1, charthick=axlCthick, p[k].n
                xyouts, xpos, ypos, charsize=0.9, font=-1, charthick=axlCthick, p[k].label
                ;stop
            end
            else: axis,xaxis=0,xtickname=empty,/save,xstyle=1,xthick=2
        endcase
        axis,xaxis=1,xtickname=empty,xstyle=1,xthick=2

        ; Y-axis
        case k of
            0: begin
                axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=2 ;,YTICKINTERVAL=p[i].tiv  ;left-column
                xpos = !x.crange[0] - 0.7*(!x.crange[1]-!x.crange[0])
                ypos = !y.crange[0] + 0.4*(!y.crange[1]-!y.crange[0])
                xyouts, xpos, ypos, alignment=0.5, orientation=90, charsize=0.9, font=-1, charthick=axlCthick, p[i].label
            end
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
            plot_like1dname,files_s12,pname_s12,p[k].n,subsamp=3,nskip=1000,scale=p[k].scale,$
              thick=4,linestyle=ls[0], color=mycolor[0], /oplot
;             plot_like1dname,files_w7,pname_w7,p[k].n,subsamp=3,nskip=1000,scale=p[k].scale,$
;               thick=4,linestyle=ls[1], color=mycolor[1], /oplot
            plot_like1dname,files_w7s12,pname_w7s12,p[k].n,subsamp=3,nskip=1000,scale=p[k].scale,$
              thick=4,linestyle=ls[2], color=mycolor[2], /oplot
            ;oplot, xvec, yvec
        
        endif else begin
            plot_like2dname,files_s12,pname_s12,p[k].n,p[i].n,scale1=p[k].scale,scale2=p[i].scale,$
              /overplot,sigma=[1,2],$
              icolors=[pscolors[5],pscolors[10]],subsamp=3,yl=0 ; forest, purple
;             plot_like2dname,files_w7,pname_w7,p[k].n,p[i].n,scale1=p[k].scale,scale2=p[i].scale,$
;               /overplot,sigma=[1,2],$
;               icolors=[pscolors[8],pscolors[3]],subsamp=3,yl=0 ; blue, orange
            plot_like2dname,files_w7s12,pname_w7s12,p[k].n,p[i].n,scale1=p[k].scale,scale2=p[i].scale,$
              /overplot,sigma=[1,2],$
              icolors=[pscolors[2],pscolors[6]],subsamp=3,yl=0 ; red, yellow
            ;oplot, xvec, yvec
        endelse

    endfor
endfor
       
ps_close

print, 'output file: '+filename_stub+'.pdf'
spawn,'epstopdf '+filename_stub+'.ps'

t1=systime(/seconds)
print, 'This procedure took ', (t1-t0)/60.,' seconds'
stop
END

