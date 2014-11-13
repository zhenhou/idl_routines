;;;
; NAME: coadd_lps12_sims_pipetest2.pro
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
;  08/16/2012: (KTS) re-write because pipeline test2 sims have 200
;     sims in them, 0-99 are of one type and 100-199 are of the other type
;  08/21/2012: (KTS) Fix bug to write out correct omap
;  08/29/2012: (KTS) Modify to write out second pipeline test
;  09/20/2012: (KTS) Modify to write out third pipeline test
;;;

pro coadd_lps12_sims_pipetest3
field = 'ra3h30dec-60'
print,field
bdir='/data23/cr/lps12/ra3h30dec-60/maps_heal/'

;bstub='map';'/map'
bend='.dat';'.bin'
ostub='coaddsim_pipetest3_'
rend='.txt'

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



prestub = bdir+'map_'+field+'_150_'
poststub = '_isims-0-143'+bend
files = file_search(prestub+'*'+poststub)
;stop

;estub='.bin'
if n_elements(files) lt n_elements(obslist) then begin
    print,'stopping because seem to have file number mismatch (is the field done?)',field,n_elements(files) , n_elements(obslist)
    return
endif


ofile0='/data/kstory/projects/lps12/sims/pipetest3/'+ostub+'spec0_'+field+'.dat'
ofile1='/data/kstory/projects/lps12/sims/pipetest3/'+ostub+'spec1_'+field+'.dat'
ofile2='/data/kstory/projects/lps12/sims/pipetest3/'+ostub+'spec2_'+field+'.dat'
if (file_test(ofile1,/read) or file_test(ofile2,/read)) then begin
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
print, 'Coadding sims.  First file: ' + prestub+obslist[0,0]+poststub;bend;print, 'Coadding sims.  First file: ' + prestub+obslist[0,0]+estub
for i=0,nrow-1 do begin ; observations
    print,i,nrow, ',  '+prestub+obslist[0,i]+poststub;bend ;print,i,nrow, ',  '+prestub+obslist[0,i]+estub
    tmap*=0
    twt*=0
    for k=0,ncol-1 do begin ; LT pairs
        infile=prestub+obslist[k,i]+poststub;bend ;infile=prestub+obslist[k,i]+estub
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

;;split into two different sims
omap1 = omap[*,0:47]
omap2 = omap[*,48:95]
omap3 = omap[*,96:143]

print, 'Write output to: ' + ofile0
openw,lun,ofile0
writeu,lun,omap0
free_lun,lun

print, 'Write output to: ' + ofile1
openw,lun,ofile1
writeu,lun,omap1
free_lun,lun

print, 'Write output to: ' + ofile2
openw,lun,ofile2
writeu,lun,omap2
free_lun,lun
end
