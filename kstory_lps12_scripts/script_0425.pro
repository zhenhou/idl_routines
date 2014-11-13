;;;
; NAME: script_0425
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) make apod masks
; 2) make kernels
; 3) run end2end
; 4) run lr jacks
; 5) pad_mask, add padded masks to mask .sav files
; 6) plot sim power from RK and KTS
;
; MODIFICATION HISTORY:
;  04/25/2012: (KTS) Created
;;;

;...................................................................
; make coadds and masks
PRO mask_6
coadd_maps_0420, 6
make_mask_lps12, field_idx=6
END

PRO mask_11
coadd_maps_0420, 11
make_mask_lps12, field_idx=11
END

PRO mask_16
coadd_maps_0420, 16
make_mask_lps12, field_idx=16
END

;...................................................................
; make coupling kernel
PRO kern_6
make_coupling_kernel, 6
END

PRO kern_11
make_coupling_kernel, 11
END

PRO kern_16
make_coupling_kernel, 16
END

;...................................................................
; make coupling kernel
PRO end_0
try_end2end, 0, run='0425', /resume
END

PRO end_1345
try_end2end, 1, run='0425', /resume
try_end2end, 3, run='0425', /resume
try_end2end, 4, run='0425', /resume
try_end2end, 5, run='0425', /resume
END

;...................................................................
; Run jackknives on fields 10, 11, 16
PRO jack_lr
lps12_jack, 10, 'lr'
lps12_jack, 11, 'lr'
lps12_jack, 16, 'lr'
END

;...................................................................
; Add padded masks to mask files
PRO pad_mask, idx
f = lps12_fieldstruct()
fst = f[idx]
field_name = fst.name
info = get_lps12_fieldinfo(idx) & nbig = info.nbig

mask_file = '/home/kstory/lps12/masks/masks_50mJy/mask_'+field_name+'.sav'
restore, mask_file

mask_padded = dblarr(nbig, nbig)
nx = (size(mask))[1]
ny = (size(mask))[2]
mask_padded[0:nx-1, 0:ny-1] = mask

save, mask, mask_padded, filename=mask_file
END


;...................................................................
; Plot RK's simulated power
PRO plot_rkpower, doall=doall
field='ra5h30dec-55'
simdir = '/data/rkeisler/low_ell_sims/output/'


; re-make everything ?
if keyword_set(doall) then begin
nbig = 2160
spawn,'ls '+simdir+'coadd_a_'+field+'*.sav',list
nlist=n_elements(list)

; get mask
restore,'/data/rkeisler/ps09/mask_'+field+'_20101009_055240.sav'
; get bigmask (zero-padded mask)
bigmask = fltarr(nbig,nbig)
bigmask[0:n_elements(mask[*,0])-1, 0:n_elements(mask[0,*])-1] = mask

mask_factor = mean(mask^2.)
bigmask_factor = mean(bigmask^2.)

l = make_fft_grid(1.0/60.*!dtor/2./!pi,960,960,fx=lx,fy=ly) ; debugging
for j=0,nlist-1 do begin
    print,strtrim(j,2),'/',strtrim(nlist-1,2)
    restore,list[j]
    if j eq 0 then begin
        nx = n_elements(coadd[*,0])
        ny = n_elements(coadd[0,*])
        power = dblarr(nx,ny)
        bigpower = dblarr(nbig,nbig)
    endif
    bigcoadd = fltarr(nbig,nbig)
    bigcoadd[0:nx-1, 0:ny-1] = coadd
    
    power[*,*] += (abs(fft(coadd*mask))^2.)
    bigpower[*,*] += (abs(fft(bigcoadd*bigmask))^2.)

    tv_spt_map, coadd*mask, winnum=2
    tv_spt_map, shift((abs(fft(coadd*mask))^2.), 480, 480), winnum=3,reso=(l[1]-l[0]),title='RK, 1sim'
    stop
    
endfor
power *= (1./nlist)
bigpower *= (1./nlist)
save, power, bigpower, nx, ny, nbig, bigcoadd, list, mask, bigmask, filename='rk_power_tmp_0425.sav'

; otherwise restore
endif else begin
    restore, 'rk_power_tmp_0425.sav'
endelse

; restore kstory stuff
;restore, 'make_twodim_tfs.sav'
restore,'/home/kstory/lps12/masks/masks_50mJy/mask_ra5h30dec-55_2008.sav'
mcfiles = '/home/kstory/lps12/sims/coaddsim_ra5h30dec-55_2008.dat'
maps = read_dat_map(mcfiles, nx, ny, 100)
map_padded = pad_array(maps[*,*,0], 4320)
power_k = (abs(fft(map_padded*mask_padded))^2.)

reso=1.0
l = make_fft_grid(reso/60.*!dtor/2./!pi,960,960,fx=lx,fy=ly)
tv_spt_map,shift(power,480,480),reso=(l[1]-l[0]),winnum=5,xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='RK, all sims, ra5h30dec-55'

l = make_fft_grid(reso/60.*!dtor/2./!pi,4320,4320,fx=lx,fy=ly)
tv_spt_map,shift(power_k,2160,2160),reso=(l[1]-l[0]),winnum=6,xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title='KTS, coadd, ra5h30dec-55'


stop
END
