;;;
; NAME: print_lrange
; PURPOSE:
;   Print the median values for the lrange test
;
; MODIFICATION HISTORY:
;  08/12/2012: (KTS) Created
;  09/10/2012: (KTS) Use 0828 chains
;  04/19/2013: (KTS) Copied from plot_s12_lcdm
;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 8-panel LCDM plot
PRO print_lrange, zoom=zoom

;----------------------
; Get the data

cdir = '/data23/kstory/lps12/chains/'
;cdir = '/data23/hou/lps12/paramfits/chains_0828/'
;chains  = ['c1_lcdm_pico_w7_newKp','c27_lcdm_camb_s12tau_newKp','c2_lcdm_pico_w7s12_newKp']

;;; SPT+WMAP7, lmax = 1500
dir_lmx15   = cdir+'c2_lcdm_pico_w7s12_lmax1500_newKp/chains/'
files_lmx15 = file_search(dir_lmx15+'c2_lcdm_pico_w7s12_lmax1500_newKp*.txt')
pname_lmx15 = dir_lmx15+'c2_lcdm_pico_w7s12_lmax1500_newKp.paramnames'

;;; SPT+WMAP7, 1500 < l < 3000
dir_lmx30 = cdir+'c2_lcdm_pico_w7s12_1500ell3000_newKp/chains/'
files_lmx30 = file_search(dir_lmx30+'c2_lcdm_pico_w7s12_1500ell3000_newKp*.txt')
pname_lmx30 = dir_lmx30+'c2_lcdm_pico_w7s12_1500ell3000_newKp.paramnames'


; S12-only
; dir_s12   = cdir+'c27_lcdm_camb_s12tau_newKp/chains/'
; files_s12 = file_search(dir_s12+'c27_lcdm_camb_s12tau_newKp*.txt')
; pname_s12 = dir_s12+'c27_lcdm_camb_s12tau_newKp.paramnames'

; WMAP-only
dir_w7 = cdir+'c1_lcdm_pico_w7_newKp/chains/'
files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7_newKp*.txt')
pname_w7 = dir_w7+'c1_lcdm_pico_w7_newKp.paramnames'

; S12+wmap
dir_w7s12 = cdir+'c2_lcdm_pico_w7s12_newKp/chains/'
files_w7s12 = file_search(dir_w7s12+'c2_lcdm_pico_w7s12_newKp*.txt')
pname_w7s12 = dir_w7s12+'c2_lcdm_pico_w7s12_newKp.paramnames'


;;; As(k0=0.002)
;cdir = '/data23/hou/lps12/paramfits/chains_0828/'
; ; S12-only
; dir_s12   = cdir+'c27_lcdm_camb_s12tau/chains/'
; files_s12 = file_search(dir_s12+'c27_lcdm_camb_s12tau*.txt')
; pname_s12 = dir_s12+'c27_lcdm_camb_s12tau.paramnames'

; ; WMAP-only
; dir_w7 = cdir+'c1_lcdm_pico_w7/chains/'
; files_w7 = file_search(dir_w7+'c1_lcdm_pico_w7*.txt')
; pname_w7 = dir_w7+'c1_lcdm_pico_w7.paramnames'

; ; S12+wmap
; dir_w7s12 = cdir+'c2_lcdm_pico_w7s12/chains/'
; files_w7s12 = file_search(dir_w7s12+'c2_lcdm_pico_w7s12*.txt')
; pname_w7s12 = dir_w7s12+'c2_lcdm_pico_w7s12.paramnames'


; For output
params = ['omegabh2', 'omegadmh2', 'theta_s', 'tau', 'ns', 'logA']
lmx15 = {value:fltarr(6), err:fltarr(6)}
lmx30 = {value:fltarr(6), err:fltarr(6)}

;;; Top Row
for i=0,2 do begin

    case i of
        0: begin ; omegahb plot
            print, ""
            print, " ------------------- omegabh2 -------------------"
            print, "*** lmx15"
            plot_like1dname,files_lmx15,pname_lmx15,'omegabh2',subsamp=3,nskip=1000,scale=100,result=res
            lmx15.value[i]=res.median & lmx15.err[i]=res.avg_err

            print, "*** lmx30"
            plot_like1dname,files_lmx30,pname_lmx30,'omegabh2',subsamp=3,nskip=1000,scale=100,result=res
            lmx30.value[i]=res.median & lmx30.err[i]=res.avg_err

        endcase
        1: begin ; omegadmh2
            print, ""
            print, " ------------------- omegadmh2 -------------------"
            print, "*** lmx15"
            plot_like1dname,files_lmx15,pname_lmx15,'omegadmh2',subsamp=3,nskip=1000,result=res
            lmx15.value[i]=res.median & lmx15.err[i]=res.avg_err

            print, "*** lmx30"
            plot_like1dname,files_lmx30,pname_lmx30,'omegadmh2',subsamp=3,nskip=1000,result=res
            lmx30.value[i]=res.median & lmx30.err[i]=res.avg_err

        endcase
        2: begin ; theta_s
            print, ""
            print, " ------------------- theta_s -------------------"
            print, "*** lmx15"
            plot_like1dname,files_lmx15,pname_lmx15,'theta_s',subsamp=3,nskip=1000,scale=100,result=res
            lmx15.value[i]=res.median & lmx15.err[i]=res.avg_err

            print, "*** lmx30"
            plot_like1dname,files_lmx30,pname_lmx30,'theta_s',subsamp=3,nskip=1000,scale=100,result=res
            lmx30.value[i]=res.median & lmx30.err[i]=res.avg_err

        endcase
;         3: begin ; omegal
;             print, ""
;             print, " ------------------- omegal -------------------"
;             print, "*** lmx15"
;             plot_like1dname,files_lmx15,pname_lmx15,'omegal*',subsamp=3,nskip=1000,result=res
;             print, "*** lmx30"
;             plot_like1dname,files_lmx30,pname_lmx30,'omegal*',subsamp=3,nskip=1000,result=res
;             lmx30.value[i]=res.median & lmx30.err[i]=res.avg_err
;         endcase
        else: begin
            print, 'i = ', i
        endcase
    endcase
endfor

;;; Bottom Row
for i=0,2 do begin

    case i of
        0: begin ; tau
            print, ""
            print, " ------------------- tau -------------------"
            print, "*** lmx15"
            plot_like1dname,files_lmx15,pname_lmx15,'tau',subsamp=3,nskip=1000,result=res
            lmx15.value[i+3]=res.median & lmx15.err[i+3]=res.avg_err

            print, "*** lmx30"
            plot_like1dname,files_lmx30,pname_lmx30,'tau',subsamp=3,nskip=1000,result=res
            lmx30.value[i+3]=res.median & lmx30.err[i+3]=res.avg_err
            
        endcase
        1: begin ; ns
            print, ""
            print, " ------------------- ns -------------------"
            print, "*** lmx15"
            plot_like1dname,files_lmx15,pname_lmx15,'ns',subsamp=3,nskip=1000,result=res
            lmx15.value[i+3]=res.median & lmx15.err[i+3]=res.avg_err

            print, "*** lmx30"
            plot_like1dname,files_lmx30,pname_lmx30,'ns',subsamp=3,nskip=1000,result=res
            lmx30.value[i+3]=res.median & lmx30.err[i+3]=res.avg_err
            
        endcase
        2: begin ; logA
            print, ""
            print, " ------------------- 1e-9As -------------------"
            print, "*** lmx15"
            plot_like1dname,files_lmx15,pname_lmx15,'1e-9As',subsamp=3,nskip=1000,scale=scale,result=res
            lmx15.value[i+3]=res.median & lmx15.err[i+3]=res.avg_err

            print, "*** lmx30"
            plot_like1dname,files_lmx30,pname_lmx30,'1e-9As',subsamp=3,nskip=1000,scale=scale,result=res
            lmx30.value[i+3]=res.median & lmx30.err[i+3]=res.avg_err
            
        endcase
;         3: begin ; H0
;             print, ""
;             print, " ------------------- H0 -------------------"
;             print, "*** lmx15"
;             plot_like1dname,files_lmx15,pname_lmx15,'H0*',subsamp=3,nskip=1000,result=res
;             print, "*** lmx30"
;             plot_like1dname,files_lmx30,pname_lmx30,'H0*',subsamp=3,nskip=1000,result=res
;         endcase
        else: begin
            print, "else case ???????????"
        endcase
    endcase
            
endfor
;*******************************

stop
END


