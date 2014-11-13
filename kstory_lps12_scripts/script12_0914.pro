;;;
; NAME: script12_0910
; PURPOSE:
;   Make WMAP-only comparison plot
;
; NOTES:
;  1) plot new theory spectra
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 8-pannel LCDM plot
PRO plot_dl_th_fromCR
readcol, '/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_0.txt',l0,dl0
readcol, '/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_1.txt',l1,dl1
readcol, '/data23/kstory/lps12/sims/pipetest3/dl_th_fromCR_2.txt',l2,dl2

!p.multi=[0,2,1]
plot, l0, dl0, /yl, xr=[650,3000],/xst & oplot, l1,dl1,color=!red & oplot, l2,dl2,color=!blue
plot, dl2/dl0 & oplot, dl2/dl0,color=!blue & oplot, dl1/dl0, color=!red
!p.multi=0
ff = '/home/kstory/public_html/notebook/spt_lps12/dl_th_fromCR_0914'
err=tvread(/png,/nodialog,filename=ff) 

END



;; Transparent plots
PRO trans
;setup_ps_plotting
cdir = '/data23/hou/lps12/paramfits/chains_0828/' & 

dir_w7 = cdir+'c1_lcdm_pico_w7/chains/' &  files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7*.txt') & pname_w7 = dir_w7+'c1_lcdm_pico_w7.paramnames'
dir_cmb = cdir+'c2_lcdm_pico_w7s12/chains/'  & files_cmb = file_search(dir_cmb+'c2_lcdm_pico_w7s12*.txt') & pname_cmb = dir_cmb+'c2_lcdm_pico_w7s12.paramnames'

x0=.20
ddx=.75
y0=.16
ddy=.77

csize = 1.0

!x.window=[x0,x0+ddx]
!x.crange=[60,80]
!y.window=[y0,y0+ddy]
!y.crange=[12, 15]

wset, 1
plot_like2dname, files_w7,pname_w7, 'H0*','Dvp57ors',sigma=[1,2,3],colors=[!black,!green,!blue], yr=[12,16],xr=[60,85],/oplot
fig1 = tvread(true=3)

wset, 2
plot_like2dname, files_cmb,pname_cmb, 'H0*','Dvp57ors',sigma=[1,2,3],colors=[!red,!orange,!yellow],yr=[12,16],xr=[60,85],/oplot
fig2 = tvread(true=3)


setup_ps_plotting
filename_stub = '~kstory/lps12/scripts/figs/tmp'
ps_open,filename_stub,/color,xsize=xb,ysize=yb,/inches,/portrait

axis,xaxis=0,font=0,xstyle=1,/save,xcharsize=csize,xthick=3
axis,xaxis=1,xtickname=empty,xstyle=1,xthick=3
axis,yaxis=0,font=0,ystyle=1,/save,ycharsize=csize,ythick=3
axis,yaxis=1,ytickname=empty,ystyle=1,ythick=3

alpha = 0.5
tv, fig1*alpha + fig2*alpha, true=3

ps_close

END






PRO make_1col_runlists
f = lps12_fieldstruct()
nfields = 21
dir = '/home/kstory/lps12/runlists/'

for ii=0, 20 do begin
;for ii=1, 2 do begin
    rin  = dir+'runlist_lps12_'+f[ii].name+'.txt'
    rout = dir+'for_oliver/runlist_lps12_1col_'+f[ii].name+'.txt'


    obslist = (read_ascii(rin)).field1
    nrow=n_elements(obslist[0,*])
    ncol=n_elements(obslist[*,0])
    if nrow eq 1 and ncol gt 1 then begin
        tmp=nrow
        nrow=ncol
        ncol=tmp
    endif

    obslist=strarr(ncol,nrow)
    case ncol of
        1: begin
            readcol,rin,a,format='a'
            obslist[0,*]=a
        end
        2: begin
            readcol,rin,a,b,format='a,a'
            obslist[0,*]=a
            obslist[1,*]=b
        end
        
        4: begin
            readcol,rin,a,b,c,d,format='a,a,a,a'
            obslist[0,*]=a
            obslist[1,*]=b
            obslist[2,*]=c
            obslist[3,*]=d
        end
        else: stop
    endcase

    obs_out = reform(obslist, ncol*nrow)

    print, 'writing out new runlist: ', rout
    openw, lun1, /get_lun,rout
    for j=0, n_elements(obs_out) -1 do printf,lun1,obs_out[j]
    close, lun1
    free_lun,lun1

endfor    

END
