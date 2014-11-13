;;;
; NAME: coadd_lps12_sims_pipetest2_ltpairs.pro
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
;  07/27/2012: (KTS) Modified to coadd sims from pipeline test2
;  08/17/2012: (KTS) Fix for the fact that these sim files have 2
;                    different spectra in them
;;;

pro coadd_lps12_sims_pipetest2_ltpairs;,field,jacklr=jacklr, lmax4500=lmax4500, stub=stub
field = 'ra3h30dec-60'
print,field
bdir='/data23/cr/lps12/ra3h30dec-60/maps/'

bstub='dmap';'/map'
bend='.dat';'.bin'
ostub='coaddsim_pipetest2_'

; get the runlist
fobslist = '/home/kstory/lps12/runlists/runlist_lps12_'+field+'.txt'

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



;--------------------
; Get the input sim maps; hard-code for ra3h30dec-60, 2 columns
prestub = bdir+bstub+'_'+field+'_150_'
poststub = '_isims-0-199'+bend
files  = strarr(ncol,nrow) ; list of input files

; output files
odir1 = '/data/kstory/projects/lps12/sims/pipetest2/spec1_xspec_maps/'
odir2 = '/data/kstory/projects/lps12/sims/pipetest2/spec2_xspec_maps/'
ofiles1 = strarr(nrow)
ofiles2 = strarr(nrow)

for i=0, nrow-1 do begin 
    for j=0, ncol -1 do begin
        files[j,i] = file_search(prestub + obslist[j,i] + poststub)
    endfor

    ofiles1[i] = odir1+'sim_xspec_pipe2_'+field+'_'+obslist[0,i]+'.dat'
    ofiles2[i] = odir2+'sim_xspec_pipe2_'+field+'_'+obslist[0,i]+'.dat'

endfor


; Checks
if n_elements(files) lt n_elements(obslist) then begin
    print,'stopping because seem to have file number mismatch (is the field done?)',field,n_elements(files) , n_elements(obslist)
    return
endif


if file_test(ofiles1[0],/read) or file_test(ofiles2[0],/read) then begin
    print,'returning from ',field,' because output file exists: ',ofiles1[0], ofiles2,[0]
    return
endif


hdr = spt_read_binary_header(files[0])

;stop
npix = hdr.nx*long64(hdr.ny)
nmappix = hdr.nmappix
; omap = fltarr(npix)
; owt = fltarr(npix)
; tmap = fltarr(npix)
; twt = fltarr(npix)
omap = fltarr(npix,hdr.nmap)
owt = fltarr(npix)
tmap = fltarr(npix,hdr.nmap)
twt = fltarr(npix)
if (hdr.isdbl or hdr.version lt 4) then stop ; not dealing with this

nsim=long64(hdr.nmap)
offset = long64(4*8)

z64 = long64(0)
get_lun,lun
print, 'Coadding sims.  First file: ' + prestub+obslist[0,0]+poststub;bend;print, 'Coadding sims.  First file: ' + prestub+obslist[0,0]+estub
for i=0,nrow-1 do begin ; loop over observations
    print,i,nrow, ',  '+prestub+obslist[0,i]+poststub;bend ;print,i,nrow, ',  '+prestub+obslist[0,i]+estub
    tmap*=0
    twt*=0
    omap *=0
    owt  *=0
    for k=0,ncol-1 do begin ; loop over LT pairs
        pting = spt_read_binary_pting(files[k,i])
        wt = spt_read_binary_wt(files[k,i])
        twt[pting]+=wt
        npp=n_elements(pting)
        scr = fltarr(npp)
        openr,lun,files[k,i]
        point_lun,lun,offset
        for j=0,nsim-1 do begin ; loop overs sims
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

    ; write out coadded LT pair
    ind = where(owt gt 0,ni)
    omap[ind] /= owt[ind]
    omap /= 1e6

    ;; split into two different sims
    omap1 = omap[*,0:99]
    omap2 = omap[*,100:199]

    print, 'Write output to: ' + ofiles1[i]
    openw,lun,ofiles1[i]
    writeu,lun,omap1
    free_lun,lun

    print, 'Write output to: ' + ofiles2[i]
    openw,lun,ofiles2[i]
    writeu,lun,omap2
    free_lun,lun
endfor

end
