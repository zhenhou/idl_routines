;;;
; NAME: read_jk_ret.pro
; PURPOSE:
;   Return chisq from jackknives, for plot_lps12_jacks_wcorrmat.pro
;
; INPUTS:
;   stub,             name of jackknife, i.e. 'lr', '12', 'azrms', 'tweight', 'moon', 'sun'
;   field_idx,        use this for checking individual fields
;   jackdir,          directory with input jack files
;   diagkernel,       diagonalize the kernel?
;   print_wiki,       print output in convenient form for spt Wiki
;
; NOTES:
; 0) This is for plot_lps12_jacks_wcorrmat.pro
; 1) expected dls saved in '/home/kstory/lps12/jacks/sims/expected_jack_dls_[stub].sav'
; 2) Returns chisq
;
; MODIFICATION HISTORY:
;  03/14/2012: (KTS) Created
;  03/26/2012: CR: adapted to include mode-mixing
;  04/01/2012: (KTS) Fix mode mixing index error
;  06/10/2012: (KTS) Add use_e2e_calib keyword
;;;


;...................................................................
; Check overlap in apod masks between fields
; from read_jk.pro, RK
FUNCTION read_jk_ret, stub, run, $
             field_idx = field_idx, $
             expected=expected, $ ; use expected Dl's from sims
             jackdir=jackdir,$
             diagkernel=diagkernel, $
             print_wiki=print_wiki, $
             k11=k11, $
             use_e2e_calib=use_e2e_calib, $
             save_plots=save_plots, $
             plot_combined=plot_combined, $ ; plot combined jacks
             stopit=stopit

compile_opt IDL2, HIDDEN

; get the calibration
mycalib = 1.0
if keyword_set(use_e2e_calib) then begin
    print, 'READ_JK: using calib from end2end'
endif else begin
    case run of
        '01' : mycalib = 0.760
        '02' : mycalib = 0.760
        '03' : mycalib = 0.760
        '04' : mycalib = 0.760
        '05' : mycalib = 0.680
        '06' : use_e2e_calib = 1 ; calib = 0.825
        '07' : use_e2e_calib = 1 ; calib = 0.825
        else : begin
            print, 'un-recognised run: '+run+', quitting.'
            RETURN, -1
        endelse
    endcase
    print, 'READ_JK: using mycalib = ', mycalib
endelse

;-------------------------
; Setup
;-------------------------
f = lps12_fieldstruct()

; directories
if ~keyword_set(jackdir) then jackdir = '/home/kstory/lps12/jacks/run_'+run+'/'
kern_dir = '/home/kstory/lps12/masks/masks_50mJy/'
end_dir  = '/home/kstory/lps12/end2end/'
mask_dir = '/data/kstory/projects/lps12/masks/masks_50mJy/'

if run eq '03' then begin
    ;if ~keyword_set(jackdir) then jackdir = '/home/kstory/lps12/jacks/run_03/'
    kern_dir = '/home/kstory/lps12/masks/0522_run03/masks_50mJy/'
    end_dir  = '/home/kstory/lps12/end2end/run_03/'
    mask_dir = '/data/kstory/projects/lps12/masks/0522_run03/masks_50mJy/'
endif

; get expected dls from sims, 'expected_dls'
if keyword_set(expected) then begin
    exp_file = '/home/kstory/lps12/jacks/sims/run_'+run+'/expected_jack_dls_'+stub+'.sav'
    print, 'Get expected Dls from file: '+exp_file
    restore, exp_file
    ;restore, '/home/kstory/lps12/jacks/sims/expected_jack_dls.sav'
endif

;-------------------------
; get weight mask dir for combining fields
;-------------------------
case stub of
    'moon': begin
        restore, mask_dir+'moon_weights.sav'
        w = moon_weights
    end
    'sun': begin
        restore, mask_dir+'sun_weights.sav'
        sun_list = [4,6,11,17,18,19]
        w = sun_weights
        w /= total(sun_weights[sun_list])
        ; set other weights to zero:
        whsun = intarr(20)
        tmp = indgen(20)
        for ii=0, n_elements(sun_list)-1 do whsun[where(tmp eq sun_list[ii])] = 1
        w *= whsun
    end
    else: begin
        restore, mask_dir+'field_weights.sav'
        w = field_weights
    end
endcase

;-------------------------
; check all fields
;-------------------------
if (n_elements(field_idx) eq 0) then begin
    list = indgen(20)
    if keyword_set(k11) then list = [0,1,3,4,5] ;;; 2008, 2009
    ;list = indgen(14) + 6 ;;; 2010, 2011

    ; not all fields get the sun jack
    if (stub eq 'sun') then list = [4,6,11,17,18,19]

; only check one field
endif else begin
    list = [field_idx]
endelse    
nlist = n_elements(list)

; Arrays to return
chisq  = dblarr(nlist)
dof    = dblarr(nlist)
pte    = dblarr(nlist)
pte2   = dblarr(nlist)
nsigma = dblarr(nlist)
npos   = intarr(nlist)
names  = strarr(nlist)

chisq_exp  = dblarr(nlist) ; do I want this here?
pte_exp    = dblarr(nlist) ; do I want this here?
nsigma_exp = dblarr(nlist)

; Arrays that are used
i0=1
i1=5
banddef = (1+findgen(6))*500.
dsfac=5.
nnbands=6
nbin=i1-i0+1
invkerns=dblarr(nbin,nbin,nlist)
diagcovs=fltarr(nbin,nlist)
dls=fltarr(nbin,nlist)


;-------------------------
; loop over fields
;-------------------------

print,' ' 
for i=0,nlist-1 do begin
;for i=1,nlist-1 do begin
    idx = list[i]
    field_name = f[idx].name
    jackfile = jackdir+'jack_150_'+field_name+'_'+stub+'_info.sav'

    name = strsplit(list[i],'/',/extract)
    names[i] = name[n_elements(name)-1]

    ;-------------------------
    ; Mode Coupling
    ;-------------------------

    use_mode_mixing = 1 ; do we want to take into account mode-mixing?
    if (use_mode_mixing) then begin

        ; restore end2end
        ;   gets 'transfer', 'beam_interp', 'kernel', 'ellkern', 
        restore, end_dir+'end_'+f[idx].name+'_'+run+'_kweight.sav'

        ; re-define things that get clobbered by end2end
        banddef = (1+findgen(6))*500.
        dsfac=5.
        nnbands=6
        nbin=i1-i0+1

        ; get the total calibration
        calib2use = keyword_set(use_e2e_calib) ? calib : mycalib

        ; re-bin to jackknife bins:
        tf_ks = interpol(transfer, ellkern, ellkern)
        beam_ks = interpol(beam_interp, ellkern, ellkern)

        use_tf = 1
        if use_tf then begin
            kern = rebin_coupling_matrix(kernel,ellkern, banddef, $
                                         transfer=tf_ks, beam=beam_ks)
        endif else begin
            kern = rebin_coupling_matrix(kernel,ellkern, banddef)
        endelse
        
        ; invert the coupling kernel
        invkerns[*,*,i]=invert(kern[i0:i1,i0:i1])

        ; diagonalize if requested
        if keyword_set(diagkernel) then begin  
            kk=kern[i0:i1,i0:i1]*0.0
            kk[indgen(i1-i0+1),indgen(i1-i0+1)]=$
              kern[indgen(i1-i0+1)+i0,indgen(i1-i0+1)+i0]
            invkerns[*,*,i]=kk
        endif

    ; otherwise ignore mode-mixing
    endif else begin
        for ii=0, nbin-1 do  $
          invkerns[ii,ii, i] = 1.0
    endelse


    ;-------------------------
    ; Retrieve the jackknife spectra, cov
    ;-------------------------

    ;print, 'Restoring jackfile: ', jackfile
    restore,jackfile
    ;undefine, dl
    ;if (size(dl))[0] eq 0 then dl=cl ; deal with old jack files that were miss-labled
    neff = n_elements(setdef)
    npair = n_elements(allspectra[0,*])


    ;-------------------------
    ; Calculate chisq using the stddev method
    ;-------------------------

    corspec = reform(invkerns[*,*,i]) ## transpose(allspectra[i0:i1,*]) * calib2use

    ; correct Dl for mode coupling
    fac = 1. ; convert sims to K^2 ??
    dls[*,i]= reform(invkerns[*,*,i]) ## dl[i0:i1] * calib2use * fac

    ; get errors on bins by directly taking the var
    for ii=i0,i1 do $
      diagcovs[ii-i0,i]=variance(corspec[*,ii-i0])*((neff-3.)/(neff-1.)/(npair-1.))
    ;corcov =  reform(invkerns[*,*,i]) ## cov[i0:i1,i0:i1] ##transpose(reform(invkerns[*,*,i])) * calib2use^2
    ;dcorcov = corcov[indgen(i1-i0+1)+i0,indgen(i1-i0+1)+i0]

    ; calculate chisq, pte, and nsigma
    chisq[i] = total((dls[*,i]/sqrt(diagcovs[*,i]))^2)
    dof[i] = nbin
    pte[i] = mpchitest(chisq[i],dof[i])
    pte2[i] = 1-chisqr_pdf(chisq[i],dof[i])
    nsigma[i] = mpchitest(chisq[i],dof[i],/sigma)
    wh = where(dls[*,i] gt 0, nwh)
    npos[i] = nwh

    if keyword_set(expected) then begin
        chisq_exp[i] = total(( (dls[*,i] - expected_dls[*,idx])/sqrt(diagcovs[*,i]))^2)
        pte_exp[i] = mpchitest(chisq_exp[i],dof[i])
        nsigma_exp[i] = mpchitest(chisq_exp[i],dof[i],/sigma)
        sigma_str_exp = (nsigma_exp[i] lt 10e5) ? sigfig(nsigma_exp[i],3) : strtrim(string(nsigma_exp[i]),2)
    endif

    ; for printing:
    sigma_str = (nsigma[i] lt 10e5) ? sigfig(nsigma[i],3) : strtrim(string(nsigma[i]),2)


    ;-------------------------
    ; Print results
    ;-------------------------

    ; wikiprint
    if keyword_set(print_wiki) then begin
        if ((pte[i] lt 0.05) or (pte[i] gt 0.95)) then begin
            mystr="'''"
            PRINT, FORMAT='(%"%s", F5.2, ", ", F7.4, ",", I2, %"%s")', mystr, chisq_exp[i], pte_exp[i], npos[i], mystr
        endif else begin
            ;PRINT, FORMAT='(F5.2, ", ", F7.4, ",", I2)', chisq[i], pte[i], npos[i]
            PRINT, FORMAT='(F5.2, ", ", F7.4, ",", I2)', chisq_exp[i], pte_exp[i], npos[i]
        endelse

    ; normal print
    endif
;     endif else begin
;         print, ' '
;         print,field_name,'  ',chisq[i],' ',pte[i],' (',sigma_str,'-sigma)'
;         if ((pte[i] lt 0.05) or (pte[i] gt 0.95)) then print,dls[*,i],(dls[*,i]/sqrt(diagcovs[*,i]))^2
;     endelse

    if ~keyword_set(print_wiki) then begin
        print, ' '              ; DEBUGGING
        print, '******* ', field_name, ', ', idx ; debugging
        print,'chisq, pte: ',chisq[i],' ',pte[i],' (',sigma_str,'-sigma)'
        print,'dls / err =  ', dls[*,i]/sqrt(diagcovs[*,i])
        print,'dls       = ', dls[*,i]*1d12
        print,'err       = ', sqrt(diagcovs[*,i])*1d12

        if keyword_set(expected) then begin 
            print, 'exp dls   = ', expected_dls[*,idx]*1d12
            print,'zero: ',chisq[i],' ',pte[i],' (',sigma_str,'-sigma)'
            print,'exp:  ', chisq_exp[i],' ',pte_exp[i],' (',sigma_str_exp,'-sigma)'
        endif
    endif

;stop
endfor ; loop over list

;-------------------------
; Combine jacks
;-------------------------

mndl = fltarr(nbin)
estcov=mndl*0.0

for j=0,nlist-1 do mndl += reform(w[j] * dls[*,j] )
mndlsq = mndl^2
                                ;and the covariance:                                                                                                                                                       
for j=0,nlist-1 do estcov += diagcovs[*,j]*w[j]^2

chisq_tot = total(mndlsq / estcov)

; print results
dof_tot = nbin
pte_tot    = mpchitest(chisq_tot,dof_tot)
nsigma_tot = mpchitest(chisq_tot,dof_tot,/sigma)
print, '-------------------'
print,'combined results:     ',chisq_tot,' ',pte_tot,' (',strtrim(string(nsigma_tot),2),'-sigma)'

; compare against expected?
if keyword_set(expected) then begin 
    mnexp = fltarr(nbin)
    for j=0,nlist-1 do mnexp += reform(w[j] * expected_dls[*,j] )
    chisq_tot_exp  = total( (mndl - mnexp)^2 / estcov)
    pte_tot_exp    = mpchitest(chisq_tot_exp,dof_tot)
    nsigma_tot_exp = mpchitest(chisq_tot_exp,dof_tot,/sigma)
    print,'combined exp results: ',chisq_tot_exp,' ',pte_tot_exp,' (',strtrim(string(nsigma_tot_exp),2),'-sigma)'
endif

print, 'npos points (out of 5): ', n_elements( where(mndl gt 0) )
print, 'dl  = ', mndl*1d12
print, 'err = ', sqrt(estcov)*1d12

;-------------------------
; plot combined jacks
;-------------------------
if keyword_set(plot_combined) then begin

    chars=2
    thick=2
    charthick=2

    window, 1
    ddl = sqrt(estcov)
    plot,banddef[i0:i1],mndl/ddl,$
      title='combined '+stub+' jack, '+strtrim(string(chisq_tot),2)+' chisq',$
      chars=chars, thick=thick, charthick=charthick, $
      xr=[min(banddef[i0:i1])-100,max(banddef[i0:i1])+100],/xst,xtitle='!12l!X!N',$
      ytitle='JK Power ('+textoidl('\sigma')+')' ;,yr=[-4,4],/yst
    errplot,banddef[i0:i1],(mndl-ddl)/ddl,(mndl+ddl)/ddl, thick=thick
    
    window, 3
    plot,banddef[i0:i1],mndl*1d12,$
      title='combined '+stub+' jack, '+strtrim(string(chisq_tot),2)+' chisq',$
      chars=chars, thick=thick, charthick=charthick, $
      xr=[min(banddef[i0:i1])-100,max(banddef[i0:i1])+100],/xst,xtitle='!12l!X!N',$
      ytitle='Dl ('+textoidl('\mu K^2')+')'
    errplot,banddef[i0:i1],(mndl-ddl)*1d12,(mndl+ddl)*1d12, thick=thick

    ; plot all jack ontop of each other
    window, 2, xsize=1000, ysize=800
    legend_str = strarr(20)
    color_arr = kscolor_array()
    yr = minmax(dls*1e12)
    plot,banddef[i0:i1],mndl*1e12,$
      title='combined '+stub+' jack, '+strtrim(string(chisq_tot),2)+' chisq',$
      chars=chars, thick=thick, charthick=charthick, $
      xr=[min(banddef[i0:i1])-100,max(banddef[i0:i1])+100],/xst,xtitle='!12l!X!N',$
      ytitle='Dl ('+textoidl('\mu K^2')+')', yr=yr ;,/yst
    
    ; make array of colors
    for ii=0, nlist-1 do begin
        fdl  = dls[*,ii]*1e12 &   fddl = sqrt(diagcovs[*,ii])*1e12
        oplot, banddef[i0:i1], fdl, color=color_arr[ii], linestyle=3, thick=thick
        errplot,banddef[i0:i1],(fdl-fddl),(fdl+fddl), color=color_arr[ii], linestyle=3, thick=thick
        legend_str[ii] = f[ii].name
    endfor
    ;legend, legend_str[0:(nlist-1)], linestyle=3, thick=thick, colors=color_arr[0:(nlist-1)], position=[2300, ((yr[1] - yr[0])*0.95 + yr[0])]
    legend, legend_str[0:(nlist-1)], linestyle=3, thick=thick, colors=color_arr[0:(nlist-1)], position=[2300, 14]
    oplot,banddef[i0:i1],mndl*1e12,thick=3
    
    ; plot all BAD jacks ontop of each other
    window, 4, xsize=1000, ysize=800
    legend_str = strarr(20)
    color_arr = kscolor_array()
    yr = minmax(dls*1e12)
    plot,banddef[i0:i1],mndl*1e12,$
      title='combined '+stub+' jack, '+strtrim(string(chisq_tot),2)+' chisq',$
      chars=chars, thick=thick, charthick=charthick, $
      xr=[min(banddef[i0:i1])-100,max(banddef[i0:i1])+100],/xst,xtitle='!12l!X!N',$
      ytitle='Dl ('+textoidl('\mu K^2')+')', yr=yr ;,/yst
    bad_list = [12,13,15,18]      ; azrms, run_05
    ;bad_list = [9, 14, 15]      ; azrms_95, run_05
    ;bad_list = [1, 6, 8, 10] ; lr
    for ilist=0, n_elements(bad_list)-1 do begin
        ii = bad_list[ilist]
        fdl  = dls[*,ii]*1e12 &   fddl = sqrt(diagcovs[*,ii])*1e12
        oplot, banddef[i0:i1], fdl, color=color_arr[ilist], linestyle=3, thick=thick
        errplot,banddef[i0:i1],(fdl-fddl),(fdl+fddl), color=color_arr[ilist], linestyle=3, thick=thick
        legend_str[ilist] = f[ii].name
    endfor
    legend, legend_str[0:3], linestyle=3, thick=thick, colors=color_arr[0:3], position=[2300, ((yr[1] - yr[0])*0.95 + yr[0])]
    
    oplot,banddef[i0:i1],mndl*1e12,thick=3
    errplot,banddef[i0:i1],(mndl-ddl)*1d12,(mndl+ddl)*1d12,thick=3
    
    ; save plots
    if keyword_set(save_plots) then begin
        plot_name = '/home/kstory/public_html/notebook/spt_lps12/jack_0514_'+stub
        err= tvread(/png, filename=plot_name, /nodialog)
    endif
    
    stop
    
    ; temporary sav file
    save, banddef, mndl, estcov, ddl, chisq_tot, pte_tot, nsigma_tot, $
      dls, diagcovs, w, chisq, pte, nsigma, $
      filename = '/home/kstory/lps12/scripts/sav_files/lr_0410.sav'
    
    print,' '
    
endif ; plot_combined

;; debugging
if keyword_set(stopit) then stop

RETURN, {dl:mndl*1d12,err:sqrt(estcov)*1d12,chisq:chisq_tot_exp,pte:pte_tot_exp}
END

