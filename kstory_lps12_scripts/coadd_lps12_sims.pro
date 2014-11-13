;;;
; NAME: coadd_lps12_sims.pro
; PURPOSE:
;   coadd lps12 sims
;
; INPUTS:
;   stub,           extra stub to output file name.  Also for coadding randcut lists.
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  06/04/2012: (KTS) Created from /home/cr/code/spt/lowell/coadd_lowell_sims.pro
;;;

pro coadd_lps12_sims,field,jacklr=jacklr, lmax4500=lmax4500, stub=stub
print,field
bdir='/data23/hou/lps12/sim_mapmaking/output/'

bstub='/map'
bend='.bin'
;temporarily using jacklr due to write problem
ostub='coaddsim_lmax8000_'
if keyword_set(lmax4500) then ostub='coaddsim_'
if keyword_set(jacklr) then begin
    bstub='/jacklr'
    ostub='jacklr/coaddsim_lr_lmax8000_'
    if keyword_set(lmax4500) then ostub='jacklr/coaddsim_lr_'
    bend='.dat'
endif
if n_elements(stub) ne 0 then ostub = ostub+stub+'_'
rend='.txt'

; get the runlist
;rstub='/data/kstory/projects/lps12/runlists/xspec2/runlist_xspec2_lps12_'
;fobslist = rstub + field + rend
fobslist = '/home/kstory/lps12/runlists/runlist_lps12_'+field+'.txt'

; randcut jacks work from preaz runlists:
if n_elements(stub) ne 0 then begin
    if (stub eq 'randcut_95_seed1_azrms') then fobslist = '/home/kstory/lps12/runlists/runlist_randcut_95_seed1_'+field+'.txt'
    if (stub eq 'randcut_95_seed2_azrms') then fobslist = '/home/kstory/lps12/runlists/runlist_randcut_95_seed2_'+field+'.txt'
    if (stub eq 'randcut_95_seed3_azrms') then fobslist = '/home/kstory/lps12/runlists/runlist_randcut_95_seed3_'+field+'.txt'
endif

if not file_test(fobslist,/read) then begin
    print,'skipping '+field+' because missing runlist'
    return
endif
obslist = (read_ascii(fobslist)).field1
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
        readcol,fobslist,a,format='a'
        obslist[0,*]=a
    end
    2: begin
        readcol,fobslist,a,b,format='a,a'
        obslist[0,*]=a
        obslist[1,*]=b
    end
    
    4: begin
        readcol,fobslist,a,b,c,d,format='a,a,a,a'
        obslist[0,*]=a
        obslist[1,*]=b
        obslist[2,*]=c
        obslist[3,*]=d
    end
    else: stop
endcase
;stop
;IF NROW GT 1 AND NCOL GT 1 THEN STOP




files = file_search(bdir+field+bstub+'_150_lmax8000_*'+bend)
;stop
prestub=bdir+field+bstub+'_150_lmax8000_'
if keyword_set(lmax4500) then begin
    files = file_search(bdir+field+bstub+'_150_*'+bend)
    prestub=bdir+field+bstub+'_150_'
endif

;estub='.bin'
if n_elements(files) lt n_elements(obslist) then begin
    print,'stopping because seem to have file number mismatch (is the field done?)',field,n_elements(files) , n_elements(obslist)
    return
endif


ofile='/data/kstory/projects/lps12/sims/'+ostub+field+'.dat'
if file_test(ofile,/read) then begin
    print,'returning from ',field,' because output file exists: ',ofile
    return
endif


hdr = spt_read_binary_header(files[0])

;stop
npix = hdr.nx*long64(hdr.ny)
omap = fltarr(npix,hdr.nmap)
owt = fltarr(npix)
tmap = fltarr(npix,hdr.nmap)
twt = fltarr(npix)
nmappix = hdr.nmappix
if (hdr.isdbl or hdr.version lt 4) then stop ; not dealing with this

nsim=long64(hdr.nmap)
nf = n_elements(files)
offset = long64(4*8)

z64 = long64(0)
get_lun,lun
print, 'Coadding sims.  First file: ' + prestub+obslist[0,0]+bend;print, 'Coadding sims.  First file: ' + prestub+obslist[0,0]+estub
for i=0,nrow-1 do begin
    print,i,nrow, ',  '+prestub+obslist[0,i]+bend ;print,i,nrow, ',  '+prestub+obslist[0,i]+estub
    tmap*=0
    twt*=0
    for k=0,ncol-1 do begin
        infile=prestub+obslist[k,i]+bend ;infile=prestub+obslist[k,i]+estub
        pting = spt_read_binary_pting(infile)
        wt = spt_read_binary_wt(infile)
        twt[pting]+=wt
        npp=n_elements(pting)
        scr = fltarr(npp)
        openr,lun,infile
        point_lun,lun,offset
        for j=0,nsim-1 do begin
            readu,lun,scr
            tmap[pting,j]+=scr*wt
        endfor
        close,lun        
    endfor
    tind = where(twt gt 0)
    owt[tind]++
    for j=0,nsim-1 do begin
        omap[tind,j] +=tmap[tind,j]/twt[tind]
    endfor
endfor

ind = where(owt gt 0,ni)
for i=0,nsim-1 do $
  omap(ind,i) /= owt(ind)
omap /= 1e6


print, 'Write output to: ' + ofile
openw,lun,ofile
writeu,lun,omap
free_lun,lun
end
