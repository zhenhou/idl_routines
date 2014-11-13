;;;
; NAME: read_jk_sim.pro
; PURPOSE:
;   read out simulated lr jackknife
;
; INPUTS:
;   stub,             name of jackknife, i.e. 'lr', '12', 'azrms', 'tweight', 'moon', 'sun'
;   field_idx,        use this for checking individual fields
;   use_e2e_calib,    use the 'calib' factor from end2end sav file.
;
; NOTES:
; 1) Adds expected Dls to sav file jackdir+'expected_dls_lrjack.sav'
;    dls[nbin=5, nfields=20]
;
; MODIFICATION HISTORY:
;  03/14/2012: (KTS) Created
;  03/26/2012: CR: adapted to include mode-mixing
;  04/01/2012: (KTS) Fix mode mixing index error
;  05/10/2012: (KTS) Copied from read_jk.pro
;  05/13/2012: (KTS) Add ability to read all types of jacks
;  05/31/2012: (KTS) Add run
;  06/10/2012: (KTS) Clean out un-needed command-line arguments, add
;                    use_e2e_calib keyword
;;;


;...................................................................
; Check overlap in apod masks between fields
; from read_jk.pro, RK
PRO read_jk_sim, stub, run, $
                 field_idx = field_idx, $
                 k11=k11, $
                 use_e2e_calib=use_e2e_calib, $
                 stopit=stopit

compile_opt IDL2, HIDDEN

; get the calibration
mycalib = 1.0
if keyword_set(use_e2e_calib) then begin
    print, 'READ_JK_SIM: using calib from end2end'
endif else begin
    case run of
        '01' : mycalib = 0.760
        '02' : mycalib = 0.760
        '03' : mycalib = 0.760
        '04' : mycalib = 0.760
        '05' : mycalib = 0.680
        '06' : use_e2e_calib = 1   ; calib = 0.825
        '07' : use_e2e_calib = 1   ; calib = 0.825
        '08' : mycalib = 0.825
        else : begin
            print, 'un-recognised run: '+run+', quitting.'
            RETURN
        endelse
    endcase
    print, 'READ_JK_SIM: using mycalib = ', mycalib
endelse


;-------------------------
; Setup
;-------------------------
f = lps12_fieldstruct()
list = ['']
;stub = 'lr_sim'
n_sim_coadds = 100.
reso = 1.0 ; arcmin per pixel
reso_rad = reso/60.*!dtor
nbig = 4320L
l = make_fft_grid(reso/60.*!dtor/2./!pi,nbig,nbig,fx=lx,fy=ly)

; directories
end_dir  = '/home/kstory/lps12/end2end/'
jackdir = '/home/kstory/lps12/jacks/sims/run_'+run+'/'
mask_dir = '/data/kstory/projects/lps12/masks/masks_50mJy/'

; Still be able to use run_03
if run eq '03' then begin
end_dir  = '/home/kstory/lps12/end2end/run_03/'
jackdir = '/home/kstory/lps12/jacks/sims/run_03/'
mask_dir = '/data/kstory/projects/lps12/masks/0522_run03/masks_50mJy/'

endif

; get the expected jack sav file
;expected_dls=fltarr(nbin,20) ; from sav file
savfile = jackdir+'expected_jack_dls_'+stub+'.sav'
if file_test(savfile) then begin 
    restore, savfile
endif else begin
    print, 'savfile does not exist!  Make new one: ' + savfile
    nbin = 5
    expected_dls=fltarr(nbin,20)
    save, expected_dls, filename=savfile
endelse

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
;chisq  = dblarr(nlist)
dof    = dblarr(nlist)
pte    = dblarr(nlist)
pte2   = dblarr(nlist)
nsigma = dblarr(nlist)
npos   = intarr(nlist)
names  = strarr(nlist)

; Arrays that are used
i0=1
i1=5
banddef = (1+findgen(6))*500.
nbin=i1-i0+1
invkerns=dblarr(nbin,nbin,nlist)
;diagcovs=fltarr(nbin,nlist)


;-------------------------
; loop over fields
;-------------------------

print,' ' 
for ilist=0,nlist-1 do begin
;for i=1,nlist-1 do begin
    idx = list[ilist]
    info = get_lps12_fieldinfo(idx)
    field_name = f[idx].name
    jackfile = jackdir+'jack_150_'+field_name+'_'+stub+'_sim_info.sav'

    name = strsplit(list[ilist],'/',/extract)
    names[ilist] = name[n_elements(name)-1]

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
        invkerns[*,*,ilist]=invert(kern[i0:i1,i0:i1])

        ; diagonalize if requested
        if keyword_set(diagkernel) then begin  
            kk=kern[i0:i1,i0:i1]*0.0
            kk[indgen(i1-i0+1),indgen(i1-i0+1)]=$
              kern[indgen(i1-i0+1)+i0,indgen(i1-i0+1)+i0]
            invkerns[*,*,ilist]=kk
        endif

    ; otherwise ignore mode-mixing
    endif else begin
        for ii=0, nbin-1 do  $
          invkerns[ii,ii, ilist] = 1.0
    endelse


    ;-------------------------
    ; Retrieve the jackknife spectra, cov
    ;-------------------------
    print, 'READ_JK_SIM: reading jack file: ', jackfile
    restore, jackfile

    ;-------------------------
    ; Calculate chisq using the stddev method
    ;-------------------------

    corspec = reform(invkerns[*,*,ilist]) ## transpose(allspectra[i0:i1,*]) * calib2use

    ; correct Dl for mode coupling
    ;dl /=1.e12 ; convert from uK^2 to K^2 (only for old lr jacks
    expected_dls[*,idx]= reform(invkerns[*,*,ilist]) ## dl[i0:i1] * calib2use

    ;-------------------------
    ; Print results
    ;-------------------------
    print, 'field '+strtrim(string(idx),2)+': '+field_name
    print,'dls [uK^2] = ', expected_dls[*,idx] * 1d12

endfor ; loop over list

;-------------------------
; Combine jacks
;-------------------------
; keep
mnexp = fltarr(nbin)
for j=0,nlist-1 do mnexp += reform(w[j] * expected_dls[*,j] )

; print results
print, '-------------------'
print,'SIM combined results:  Dl [uK^2]'
print, 'npos points (out of 5): ', n_elements( where(mnexp gt 0) );, mean(npos[[0,1,3,4,5]])
print, mnexp*1d12

; save output
save, expected_dls, calib2use, filename=savfile

;; debugging
if keyword_set(stopit) then stop

END

