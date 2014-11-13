;;;
; NAME: coadd_lps12_jack_sims.pro
; PURPOSE:
;   coadd lps12 sims
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  05/11/2012: (KTS) Created from /home/cr/code/spt/lowell/coadd_lowell_sims.pro
;  06/05/2012: (KTS) Modified for lmax8000 sims
;  06/16/2012: (KTS) Add run-specific output, i.e. lps12/sims/run_07/
;;;



PRO coadd_lps12_jack_sims,idx,stub, run=run, out_dir=out_dir

f = lps12_fieldstruct()
field = f[idx].name
print,field

; directories
bdir='/data23/hou/lps12/sim_mapmaking/output/'
if n_elements(out_dir) eq 0 then out_dir = '/data/kstory/projects/lps12/sims/run_'+run+'/'

bstub='/map'
bend='.bin'
;temporarily using jacklr due to write problem
ostub='coaddsim_lmax8000_'+stub+'_'
if stub eq 'lr' then begin
    bstub='/jacklr'
    bend='.dat'
endif
prestub=bdir+field+bstub+'_150_lmax8000_'

; get the runlist
rend='.txt'
;rstub='/data/kstory/projects/lps12/runlists/xspec2/runlist_xspec2_lps12_'
;fobslist = rstub + field + rend
fobslist = f[idx].runlist

; randcut jacks work from preaz runlists:
if ( (stub eq 'randcut_95_seed1_azrms') or (stub eq 'randcut_95_seed2_azrms') or (stub eq 'randcut_95_seed3_azrms') ) then begin
    fobslist = '/home/kstory/lps12/runlists/runlist_preaz_'+f[idx].name+'.txt'
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
;IF NROW GT 1 AND NCOL GT 1 THEN STOP



; get the file list
files = strarr(nrow)
for i=0, nrow-1 do begin
    files[i] = file_search(prestub+obslist[0,i]+bend)
endfor

if n_elements(files) lt n_elements(obslist[0,*]) then begin
    print,'stopping because seem to have file number mismatch (is the field done?)',field,n_elements(files) , n_elements(obslist)
    return
endif


ofile=out_dir+ostub+field+'.dat'
if file_test(ofile,/read) then begin
    print,'returning from ',field,' because output file exists: ',ofile
    return
endif


hdr = spt_read_binary_header(files[0])

npix = hdr.nx*long64(hdr.ny)
omap = fltarr(npix,hdr.nmap)
owt = fltarr(npix)
tmap = fltarr(npix,hdr.nmap)
twt = fltarr(npix)
nmappix = hdr.nmappix
if (hdr.isdbl or hdr.version lt 4) then stop ; not dealing with this

nsim=long64(hdr.nmap)
;nf = n_elements(files)
offset = long64(4*8)

; get the jack setdef
jackdef=get_lps12_jack_defs(idx,stub,files, /sim)
s = size(jackdef.setdef)

nsets = s[0] ; either 1 or 2
setsize = s[1] ; number of maps per set

; Coadd the maps
z64 = long64(0)
get_lun,lun
print, 'Coadding sims.  First file: ' + prestub+obslist[0,jackdef.setdef[0, 0]]+bend
for i=0,setsize-1 do begin
    print,i,setsize
    for j=0, nsets-1 do begin

        tmap*=0
        twt*=0

        ; make the second set negative
        jack_fac = j ? -1. : 1.

        for k=0,ncol-1 do begin
            ;sidx = jackdef.setdef[i, 0]
            infile = prestub+obslist[k,jackdef.setdef[i, j]]+bend
            print, infile
            pting = spt_read_binary_pting(infile)
            wt = spt_read_binary_wt(infile)
            twt[pting]+=wt
            npp=n_elements(pting)
            scr = fltarr(npp)
            openr,lun,infile
            point_lun,lun,offset
            for jsim=0,nsim-1 do begin
                readu,lun,scr
                tmap[pting,jsim]+=scr*wt*jack_fac
            endfor
            close,lun        
        endfor
        tind = where(twt gt 0)
        owt[tind]++
        for jsim=0,nsim-1 do begin
            omap[tind,jsim] +=tmap[tind,jsim]/twt[tind]
        endfor
    endfor
endfor

ind = where(owt gt 0,ni)
for i=0,nsim-1 do $
  omap(ind,i) /= owt(ind)
omap /= 1e6 ; convert to Kelvin


openw,lun,ofile
writeu,lun,omap
free_lun,lun
end
