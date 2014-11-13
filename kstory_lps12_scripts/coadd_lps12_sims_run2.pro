;;;
; NAME: coadd_lps12_sims_run2.pro
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
;  08/07/2012: (KTS) Coadd 400 sims
;;;

pro coadd_lps12_sims_run2,field;,jacklr=jacklr, lmax4500=lmax4500, lmax8000=lmax8000, stub=stub
print,field
bdir='/data23/hou/lps12/sim_run2/output/'

bstub='/map'
bend='.bin'
ostub='coaddsim_run2_'

; Previous sim runs
;if keyword_set(lmax4500) then ostub='coaddsim_'
; if keyword_set(lmax8000) then begin
;     bdir='/data23/hou/lps12/sim_mapmaking/output/'
;     ostub='coaddsim_lmax8000_'
; endif

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
; if n_elements(stub) ne 0 then begin
;     if (stub eq 'randcut_95_seed1_azrms') then fobslist = '/home/kstory/lps12/runlists/runlist_randcut_95_seed1_'+field+'.txt'
;     if (stub eq 'randcut_95_seed2_azrms') then fobslist = '/home/kstory/lps12/runlists/runlist_randcut_95_seed2_'+field+'.txt'
;     if (stub eq 'randcut_95_seed3_azrms') then fobslist = '/home/kstory/lps12/runlists/runlist_randcut_95_seed3_'+field+'.txt'
; endif

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

f0 = file_search(bdir+field+bstub+'_0sims99_150_lmax8000_*'+bend)
f1 = file_search(bdir+field+bstub+'_100sims199_150_lmax8000_*'+bend)
f2 = file_search(bdir+field+bstub+'_200sims299_150_lmax8000_*'+bend)
f3 = file_search(bdir+field+bstub+'_300sims399_150_lmax8000_*'+bend)
;files = [f0, f1, f2, f3]
ppre = strarr(4)
ppre[0] = bdir+field+bstub+'_0sims99_150_lmax8000_'
ppre[1] = bdir+field+bstub+'_100sims199_150_lmax8000_'
ppre[2] = bdir+field+bstub+'_200sims299_150_lmax8000_'
ppre[3] = bdir+field+bstub+'_300sims399_150_lmax8000_'

; if keyword_set(lmax8000) then begin
;     files = file_search(bdir+field+bstub+'_150_lmax8000_*'+bend)
;     prestub=bdir+field+bstub+'_150_lmax8000_'
; endif

; if keyword_set(lmax4500) then begin
;     files = file_search(bdir+field+bstub+'_150_*'+bend)
;     prestub=bdir+field+bstub+'_150_'
; endif

;estub='.bin'
if n_elements(f0) lt n_elements(obslist) then begin
    print,'stopping because seem to have file number mismatch (is the field done?)',field,n_elements(f0) , n_elements(obslist)
    stop
    return
endif

if ((n_elements(f0) + n_elements(f1) + n_elements(f2) + n_elements(f3))/4. ne n_elements(f0) )then begin
    print,'stopping because f0, f1, f2 and f3 do not have the same number of files.',field,n_elements(f0) , n_elements(obslist)
    return
endif


ofile='/data/kstory/projects/lps12/sims/'+ostub+field+'.dat'
if file_test(ofile,/read) then begin
    print,'returning from ',field,' because output file exists: ',ofile
    return
endif


hdr = spt_read_binary_header(f0[0])

;stop
nmap = long64(hdr.nmap)
nsim = nmap * 4
npix = hdr.nx*long64(hdr.ny)
omap = fltarr(npix,nsim)
owt = fltarr(npix)
tmap = fltarr(npix,nsim)
twt = fltarr(npix)
nmappix = hdr.nmappix
if (hdr.isdbl or hdr.version lt 4) then stop ; not dealing with this

nfiles = n_elements(f0)
offset = long64(4*8)

z64 = long64(0)
get_lun,lun
print, 'Coadding sims.  First file: ' + ppre[0]+obslist[0,0]+bend;print, 'Coadding sims.  First file: ' + prestub+obslist[0,0]+estub

for isim=0, 4-1 do begin        ; loop over 100-sim files
    prestub = ppre[isim]
    
;    for i=0,nrow-1 do begin     ; loop over observations
    for i=0,5 do begin     ; loop over observations
        print,i,nrow, ',  '+prestub+obslist[0,i]+bend
        tmap*=0
        twt*=0
        
        for k=0,ncol-1 do begin ; loop over LT pairs
            infile=prestub+obslist[k,i]+bend 
            pting = spt_read_binary_pting(infile)
            wt = spt_read_binary_wt(infile)
            twt[pting]+=wt
            npp=n_elements(pting)
            scr = fltarr(npp)
            openr,lun,infile
            point_lun,lun,offset
            for j=isim*nmap, (isim+1)*nmap-1 do begin
                print, "i, isim, k, j = ", i, isim, k, j
                readu,lun,scr
                tmap[pting,j]+=scr*wt
            endfor
            close,lun        
        endfor
        tind = where(twt gt 0)
        owt[tind]++
        for j=isim*nmap, (isim+1)*nmap-1 do begin
            omap[tind,j] +=tmap[tind,j]/twt[tind]
        endfor
    endfor
    stop
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
