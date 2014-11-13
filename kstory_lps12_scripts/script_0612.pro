;;;
; NAME: script_0612
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) Sort out runlists -> moved to make_runlist_0612
; 2) 
;
; MODIFICATION HISTORY:
;  06/12/2012: (KTS) Created
;;;

;-----------------------------
; apod maps
PRO apod_masks
f = lps12_fieldstruct()

END





;;;;;;;;;;;;;;;;;;;;;;;;;
; Runlists
;;;;;;;;;;;;;;;;;;;;;;;;;

;-----------------------------
; Make runlist_rawmap_field.txt
PRO make_runlist_rawmap
rdir_old = '/home/kstory/lps12/runlists_0612/'
rdir_new = '/home/kstory/lps12/runlists/'
f = lps12_fieldstruct()

for i=0, 20 do begin
    print, 'cp ' + rdir_old+'runlist_lps12_'+f[i].name+'.txt  ' + $
      rdir_new+'runlist_rawmap_'+f[i].name+'.txt'
endfor
END


;-----------------------------
; Make runlist_raw_field.txt
PRO make_runlist_raw
rdir_old = '/home/kstory/lps12/runlists_0612/'
rdir_new = '/home/kstory/lps12/runlists/'
f = lps12_fieldstruct()

for i=0, 20 do begin
    print, 'cp ' + rdir_old+'runlist_dateOnly_lps12_'+f[i].name+'.txt  ' + $
      rdir_new+'runlist_raw_'+f[i].name+'.txt'
endfor
END



;-----------------------------
; Make runlist_preaz_field.txt
PRO make_runlist_preaz
rdir_old = '/home/kstory/lps12/runlists_0612/'
rdir_new = '/home/kstory/lps12/runlists/'
f = lps12_fieldstruct()

for i=0, 20 do begin
    if f[i].lead_trail then begin

        ; get ncol, nrow
        fobslist = rdir_old+'lt/runlist_lt_'+f[i].name+'.txt'
        obslist = (read_ascii(fobslist)).field1
        nrow=n_elements(obslist[0,*])
        ncol=n_elements(obslist[*,0])
        if nrow eq 1 and ncol gt 1 then begin
            tmp=nrow
            nrow=ncol
            ncol=tmp
        endif
        obslist=strarr(ncol,nrow)

        ; read obs list
        case ncol of
            4: begin
                readcol,fobslist,a,b,c,d,format='a,a,a,a'
                obslist[0,*]=a
                obslist[1,*]=b
                obslist[2,*]=c
                obslist[3,*]=d
            end
             2: begin
                readcol,fobslist,a,b,format='a,a'
                obslist[0,*]=a
                obslist[1,*]=b
            end
        endcase

        ; sort
        sind = sort(obslist[0,*])
        sobslist = strarr(ncol,nrow)
        for jcol=0, ncol-1 do begin
            sobslist[jcol,*] = obslist[jcol,sind]
        endfor


        ; write new runlist
        new_runlist = rdir_new+'runlist_preaz_'+f[i].name+'.txt'
        get_lun, lun1
        openw, lun1, new_runlist
        for irow=0, nrow-1 do begin
            case ncol of 
                4: printf, lun1, sobslist[0,irow] +' '+ sobslist[1,irow] +' '+ sobslist[2,irow] +' '+ sobslist[3,irow]
                2: printf, lun1, sobslist[0,irow] +' '+ sobslist[1,irow]
            endcase
        endfor
        close, lun1
        free_lun,lun1

    ; not lead trail
    endif else begin
        new_runlist = rdir_new+'runlist_preaz_'+f[i].name+'.txt'
        print, 'new runlist: ' + new_runlist
        spawn, 'cp ' + rdir_old+'runlist_xspec_lps12_'+f[i].name+'.txt  ' + new_runlist
          
    endelse
endfor
stop
END




;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check preaz runlists
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO check_preaz_field1
ncol = 4
nrow = 272

fold = '/home/kstory/lps12/runlists_0612/lt/runlist_lt_ra23h30dec-55_2008.txt'
readcol,fold,a,b,c,d,format='a,a,a,a'

fnew = '/home/kstory/lps12/runlists/runlist_preaz_ra23h30dec-55_2008.txt'
readcol,fnew,a1,b1,c1,d1,format='a,a,a,a'

whidx_a = intarr(nrow)
whfail_a = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_a[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_a[ir] = 1
endfor

whidx_b = intarr(nrow)
whfail_b = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_b[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_b[ir] = 1
endfor

whidx_c = intarr(nrow)
whfail_c = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_c[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_c[ir] = 1
endfor

whidx_d = intarr(nrow)
whfail_d = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_d[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_d[ir] = 1
endfor


stop
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO check_preaz_field3
ncol = 2
nrow = 270

fold = '/home/kstory/lps12/runlists_0612/lt/runlist_lt_ra23h30dec-55_2008.txt'
readcol,fold,a,b,format='a,a'

fnew = '/home/kstory/lps12/runlists/runlist_preaz_ra23h30dec-55_2008.txt'
readcol,fnew,a1,b1,format='a,a'

whidx_a = intarr(nrow)
whfail_a = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_a[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_a[ir] = 1
endfor

whidx_b = intarr(nrow)
whfail_b = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_b[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_b[ir] = 1
endfor

stop
END





PRO plot_ps_1field
restore, '/home/kstory/lps12/end2end/end_ra22h30dec-55_05_kweight.sav'
whplot = indgen(45) + 10
yr = [40,8e3]
xr = [0,3.5e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
l = banddef-25
readcol,'/home/kstory/lps12/cls_theory/Dls_theory.txt',l_vec,dl_th
cov = reform(cov)
dl = spectrum*0.68
diag=dblarr(58)
for i=0, 57 do diag[i] = sqrt(cov[i,i])*0.68

; plots
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

plot, l[whplot], dl[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3,/xst,xtitle=xtitle,ytitle=yitle,chars=chars,title='PS, ra22h30dec-55'
oplot,l_vec,dl_th,color=!red
errplot,l[whplot],(dl-diag)[whplot]*1d12,(dl+diag)[whplot]*1d12


stop
END



PRO plot_err
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_022240_kweight.sav' ; to be sure.
dl_05  = dl_all
diag05 = diag_nobeam

diag05_data = dblarr(58)
diag05_sv   = dblarr(58)
;cov_e2e = cov_all_nobeam_data + cov_all_nobeam_mc
for ii=0, 57 do begin
    diag05_data[ii] = sqrt(cov_all_nobeam_data[ii, ii]); / 0.68
    diag05_sv[ii]   = sqrt(cov_sv[ii,ii])
endfor


yr = [0.5, 10e1]
xtitle= '!12l!X!N'
ytitle= textoidl('\delta')+'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10


; plots
fdir = '/home/kstory/public_html/notebook/spt_lps12/'
v1 = indgen(50)+8

plot, l[v1], diag05[v1]*1d12, /ylog,xr=xr,yr=yr,xtitle=xtitle,ytitle=ytitle,title='errors, 2500 deg^2'
oplot, l[v1], diag05_sv[v1]*1d12, color=!red
oplot, l[v1], diag05_data[v1]*1d12, color=!yellow

legend, ['dDl_tot', 'dDl_sv', 'dDl_data'], colors=[!white,!red,!yellow], linestyle=[0,0,0], $
  pos=[2000,50]

stop

END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot sv error bars
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO plot_sv_err

restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_022240_kweight.sav' ; to be sure.
restore, '/home/kstory/lps12/end2end/sv_err_05_kweight.sav' ; gets sv_err
;restore, '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav'
;restore, '/home/kstory/lps12/end2end/sv_err_04_kweight.sav' ; gets sv_err

v1 = indgen(50)+8
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

window, 1, xsize=700, ysize=500
plot, l[v1], sqrt(cov_all_nobeam[v1,v1])*1d12, yr=[0.5, 100],/ylog, xtitle='ell', ytitle='dDl [K^2]', title='Simulated SV errors'
oplot, l[v1], sqrt(cov_all_nobeam_e2e[v1,v1])*1d12, color=!red                                                    
oplot, l[v1], sqrt(cov_all_nobeam_data[v1,v1])*1d12, color=!skyblue
oplot, l[v1], sqrt(cov_all_nobeam_mc[v1,v1])*1d12, color=!lavender
oplot, l[v1], sqrt(cov_sv[v1,v1])*1d12, color=!green
oplot, l[v1], sqrt(cov_all_nobeam[v1,v1])*1d12
legend, ['Sim_total_err', 'original_err', 'data_err', 'mc_err', 'sv_err'], linestyle=[0,0,0,0,0], $
  colors=[!white, !red, !skyblue, !lavender, !green], pos=[2000, 50]

stop
END



PRO plot_beams
bf05 = get_lps12_beams(2010, lb, bl10_05)
bf04 = get_lps12_beams(2010, lb, bl10_04, /sim_lmax8000)

fdir = '/home/kstory/public_html/notebook/spt_lps12/'
window, 0
plot, lb, bl10_04, xr=[0,4000], yr=[0.6,1.1], title='Beams 2010, prelim-release v.s. rev3.0'
oplot, bl10_05, color=!red
err=tvread(/png,/nodialog,filename=fdir+'beams_0612')

window, 1
plot, lb, bl10_05 / bl10_04, xr=[0,4000], yr=[1.04,1.12], title='Beams 2010, prelim-release v.s. rev3.0'
oplot, [650,650], [0,10], linestyle=1
oplot, [3000,3000], [0,10], linestyle=1
err=tvread(/png,/nodialog,filename=fdir+'beams_frac_0612')

stop
END



;;;;;;;;;;;;;;;;;;;;;;
; Plot error bars
;;;;;;;;;;;;;;;;;;;;;;

PRO plot_ps_run_05

; plotting setup
if n_elements(xr) eq 0 then xr = [0,3.5e3]
yr = [40,8e3]
xtitle= '!12l!X!N'
ytitle= 'D!D!12l!X!N'+textoidl(' (\muK^2)')
chars=1.8
whplot = indgen(45) + 10

;run_05
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120610_190210_kweight.sav' ; calib applied to cov
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120611_175439_kweight.sav' ; no calib on cov
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120611_225834_kweight.sav' ; new cov_sv
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120611_231531_kweight.sav' ; fix by-year
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120611_235136_kweight.sav' ; fix beam errors
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_022240_kweight.sav' ; to be sure.
dl_05  = dl_all
diag05 = diag_nobeam
diag05_e2e = dblarr(58)
cov_e2e = cov_all_nobeam_data + cov_all_nobeam_mc
for ii=0, 57 do diag05_e2e[ii] = sqrt(cov_e2e[ii, ii])

; run_05; 0809 only
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_024104_kweight_0809.sav'; best guess
;restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_025500_kweight_0809.sav' ; no combined_calib applied to cov_data
restore, '/home/kstory/lps12/end2end/run_05/combined_spectrum_20120612_040328_kweight_0809.sav' ; corrected cov_sv
dl_0809     = dl_all
diag05_0809 = diag_nobeam

; run_04
restore, '/home/kstory/lps12/end2end/run_04/combined_spectrum_20120606_134756_kweight.sav'
diag04 = diag_nobeam

; K11
restore, '/data/rkeisler/ps09/1.37/combined_spectrum_20110402_164956_kweight.sav'
err_rk = diag_nobeam


; cl_theory
readcol,'/home/kstory/lps12/cls_theory/Cls_theory.txt',l_vec,cl_uK2
cl_uK2[0] = cl_uK2[1]       ; get rid of the zero at l=0
dl_th = cl_uK2 * l_vec*(l_vec+1) / (2*!pi)

;----------------
; Plots
;----------------
fdir = '/home/kstory/public_html/notebook/spt_lps12/'

; ps plot
window, 0
plot,l[whplot],dl_05[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_05'
oplot,l_vec,dl_th,color=!red
errplot,l[whplot],(dl_05-diag05)[whplot]*1d12,(dl_05+diag05)[whplot]*1d12
; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0611a')

window, 1
plot,l[whplot],dl_05[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
  xtitle=xtitle,ytitle=ytitle,chars=chars, title='Full 2500 deg^2, run_05'
errplot,l[whplot],(dl_05-diag05)[whplot]*1d12,(dl_05+diag05)[whplot]*1d12

; plot,l[whplot],dl_all[whplot]*1d12,xr=xr,yr=yr,/yst,/yl,ps=3, /xst, $
;   xtitle=xtitle,ytitle=ytitle,chars=chars, title='2008, 2009 data, run_05'
; oplot,l_vec,dl_th,color=!red
; errplot,l[whplot],(dl_all-diag_nobeam)[whplot]*1d12,(dl_all+diag_nobeam)[whplot]*1d12
; ;err=tvread(/png,/nodialog,filename=fdir+'ps_0611a')


; ;----------------
; ; make err plot comparing run 04 to 05

; window, 1
; plot, l[whplot], ((diag04 - diag05)/diag05)[whplot], xtitle=xtitle, ytitle='(diag04 - diag05)/diag05', $
;   title='Error bars, 04 vs 05'
; ;err=tvread(/png,/nodialog,filename=fdir+'err_04v05_0611a')

; ; compare against K11
window, 2
plot, l[whplot], ((err_rk )/diag05)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05', title='Error bars, K11 vs 05'
oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
;err=tvread(/png,/nodialog,filename=fdir+'err_k11v05_0611a')

; window, 3
; plot, l[whplot], ((err_rk )/diag05_e2e)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05', title='Error bars, K11 vs 05_e2e'
; oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
; ;err=tvread(/png,/nodialog,filename=fdir+'err_k11v05e2e_0611a')

; window, 4
; plot, l[whplot], (diag05/diag05_e2e)[whplot], xtitle=xtitle, ytitle='diag05/diag05_e2e', title='Error bars, 05_sim vs 05_e2e'
; oplot, [0,10000], [sqrt(2500./790), sqrt(2500./790)], linestyle=1
; ;err=tvread(/png,/nodialog,filename=fdir+'err_05simv05e2e_0611a')

window, 5
plot, l[whplot], (err_rk/diag05_0809)[whplot], xtitle=xtitle, ytitle='diag_k11/diag05_0809', title='Error bars, K11 vs 05_0809'
stop
END
