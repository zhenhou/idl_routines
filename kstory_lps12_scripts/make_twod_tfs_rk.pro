; use option /hundred to make TF for only 100 sims
pro make_twod_tfs_rk, doall=doall, plotit=plotit, dosave=dosave, hundred=hundred

;idx_list=[-1, 1];,1,3,4,5] ; Specify which files to run this on.
field_list = indgen(20)
nfields = n_elements(field_list)

; intermediate sav file for RK data
tf_idf = keyword_set(hundred) ? '/home/kstory/lps12/twod_tfs/tf_rk/create_twod_tfs_100s.sav' : $; DEBUGGING
  '/home/kstory/lps12/twod_tfs/tf_rk/create_twod_tfs.sav'
tf_dir = '/home/kstory/lps12/twod_tfs/tf_rk/'

restore, tf_idf

; set default values for structures going into sav file
for ii=0, nfields-1 do ex=execute('s'+strtrim(floor(ii),2)+' = -1')

if keyword_set(doall) then begin

;f = lps12_fieldstruct()
info = get_lps12_fieldinfo(0)
nbig = info.nbig

reso = 1.0; arcmin per pixel
reso_rad = reso/60.*!dtor

fields = field_vec()
nfields = n_elements(fields)

simdir = '/data/rkeisler/low_ell_sims/output/'

nbig = 2160

for i=0,nfields-1 do begin
;for i=0,0 do begin ;DEBUGGING
    field = fields[i]
    print, 'analize field '+field

; get list
    spawn,'ls '+simdir+'coadd_a_'+field+'*.sav',list
    nlist = keyword_set(hundred) ? 100 : n_elements(list)

; get mask
    restore,'/data/rkeisler/ps09/mask_'+field+'_20101009_055240.sav'
; get bigmask (zero-padded mask)
    bigmask = fltarr(nbig,nbig)
    bigmask[0:n_elements(mask[*,0])-1, 0:n_elements(mask[0,*])-1] = mask

    mask_factor = mean(mask^2.)
    bigmask_factor = mean(bigmask^2.)
    
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

    endfor
    power *= (1./nlist)
    bigpower *= (1./nlist)

    ex=execute('s'+strtrim(floor(i),2)+'={mask_factor:mask_factor, bigmask_factor:bigmask_factor, field:field, nx:nx, ny:ny, power:power, bigpower:bigpower}')
    
endfor

;save,s0,s1,s2,s3,s4,fields,filename='create_twod_tfs_rk.sav'
save,s0,s1,s2,s3,s4,fields,filename=tf_idf
endif else restore,tf_idf

;readcol,'/data/rkeisler/low_ell_sims/input/spectrum_file.txt',lth,dlth
readcol,'/data/rkeisler/low_ell_sims/input/dl_input_true_20101221_061503.txt',lth,dlth

reso = 1.0; arcmin per pixel
reso_rad = reso/60.*!dtor
for i=0,4 do begin
    ex=execute('s=s'+strtrim(i,2))
    if n_elements(nbig) eq 0 then nbig = n_elements((s.bigpower)[*,0])
    field = s.field
    power = s.power*1d12/s.mask_factor*(s.nx*s.ny)*(reso_rad^2.)
    bigpower = s.bigpower*1d12/s.bigmask_factor*(nbig*nbig)*(reso_rad^2.)

    l = make_fft_grid(reso/60.*!dtor/2./!pi,s.nx,s.ny,fx=lx,fy=ly)
    big_l = make_fft_grid(reso/60.*!dtor/2./!pi,nbig,nbig,fx=big_lx,fy=big_ly)

    dlth_interp = interpol(dlth,lth,l)
    big_dlth_interp = interpol(dlth,lth,big_l)

    clth_interp = dlth_interp/l/(l+1.)*2.*!pi
    big_clth_interp = big_dlth_interp/big_l/(big_l+1.)*2.*!pi
    
; get the beam (get the exact beam used by the sims)
    beamdir = '/data/rkeisler/low_ell_sims/input/rev2.0_copy_20Sep2010/'
    case field of 
        'ra5h30dec-55': beamfile=beamdir+'blgrid_2008_150.txt'
        'ra23h30dec-55': beamfile=beamdir+'blgrid_2008_150.txt'
        'ra3h30dec-60': beamfile=beamdir+'blgrid_2009_150.txt'
        'ra21hdec-60': beamfile=beamdir+'blgrid_2009_150.txt'
        'ra21hdec-50': beamfile=beamdir+'blgrid_2009_150.txt'
    endcase
    readcol,beamfile,l_beam,bl_beam,format='d,d'


    bl = interpol(bl_beam,l_beam,l)
    wh_high = where(l gt max(l_beam), n_high)
; replace the part of the beam that is undefined with the beam value
; at the highest ell that is defined.
    last_bl_beam = (bl_beam[n_elements(bl_beam)-1])[0]
    if n_high gt 0 then bl[wh_high]=last_bl_beam
    
; repeat for big array
    big_bl = interpol(bl_beam,l_beam,big_l)
    wh_high = where(big_l gt max(l_beam), n_high)
; replace the part of the beam that is undefined with the beam value
; at the highest ell that is defined.
    last_bl_beam = (bl_beam[n_elements(bl_beam)-1])[0]
    if n_high gt 0 then big_bl[wh_high]=last_bl_beam


; get the effective change in the TF from the mode-coupling
; correction.
    restore,'/home/rkeisler/ps09/save_mode_coupling.sav'
    mode_coupling = mode_couplings[i,*]
    mode_coupling[where(l_mode_coupling gt 5e3)]=1.
    corr_mc = interpol(mode_coupling,l_mode_coupling,l)
    wh_high = where(l gt max(l_mode_coupling), n_high)
    if n_high gt 0 then corr_mc[wh_high] = $
      mode_coupling[n_elements(mode_coupling)-1]
    big_corr_mc = interpol(mode_coupling,l_mode_coupling,big_l)
    wh_high = where(big_l gt max(l_mode_coupling), n_high)
    if n_high gt 0 then big_corr_mc[wh_high] = $
      mode_coupling[n_elements(mode_coupling)-1]


    tf = power/clth_interp/bl/bl*corr_mc
    big_tf = bigpower/big_clth_interp/big_bl/big_bl*big_corr_mc


    tf_w_beam = power/clth_interp*corr_mc
    big_tf_w_beam = bigpower/big_clth_interp*big_corr_mc

if keyword_set(plotit) then tv_spt_map,shift(tf,s.nx/2.,s.ny/2.),min=0,max=1,reso=(l[1]-l[0]),xtitle='!12l!X!N!Dx!X!N',ytitle='!12l!X!N!Dx!X!N',chars=1.3,title=s.field

;stop

README = 'All TFs are in POWER units, as opposed to TEMPERATURE units.'
if keyword_set(dosave) then begin
    save, $
      l,lx,ly,tf,tf_w_beam,README, $
      filename=tf_dir+'tf_'+field+'.sav'

    save, $
      big_l,big_lx,big_ly,big_tf,big_tf_w_beam,README, $
      filename=tf_dir+'tf_pad_'+field+'.sav'
endif



endfor

end

