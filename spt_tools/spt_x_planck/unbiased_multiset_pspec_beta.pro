;+
; NAME:
;
;  unbiased_multiset_pspec
;
; PURPOSE:
;
;  From a collection of maps, in one or more sets (ie frequencies,
;  or jacknife sets) estimate the temperature power spectrum.  If the
;  data is divided into multiple sets then the cross covariances are
;  also computed. 
;
; CATEGORY:
;
;
;
; Calling Sequence:
;
;
;
; INPUTS:
;
;  mapfile: The name of a file containing maps as dblarr's.  The size
;  of the maps is assumed to be the same as the size of the data
;  window function.   The data is assumed to be unwindowed.
;
;  If mapfile is an array, then it is interpreted as a list of files
;  from which to extract the maps.  
;  
;  The input files may be a .bin ("lotsomaps" binary format) a .fits, 
;  or a .sav file. 
;  If a save file is input,
;  it is assumed that the desired maps have the same name 
;  in every file.  The mapname input must also be specified for a 
;  file list to work.  The syntax for specifying the mapname is 
;  shown in these examples:
;  
;  For run1a-style .sav files:  mapname='mapstruct.map'
;  For Christian's simulation .sav files:
;    mapname='mapstruct.map[*,*,8]' if the desired simulation is number 8
;  For standard spt-style fits files:  mapname='map[8].map' for the
;    8th simulation.  Or, simply mapname='map.map' if it's a data file.
;
;  If the array 2 dimensional then each row in the array corresponds
;  to a different _set_ of maps.  
;
;  window:  A dblarr definition of the window function
;
;  reso_arcmin: The resolution of the map in arcminutes
;
; OPTIONAL INPUTS:
;
;  banddef:  (INPUT/OUTPUT)
;            An array definition of the bands to be used in creating the
;            power spectrum.   One element per band.  Each element
;            defines the upper ell-limit of that band. 
;            Default: bins are spaced at half the minimum spacing
;            and extend up to ell of 10000 or the maximum ell
;            permitted by the resolution, whichever is smaller
;            If undefined the default bands are stored in the variable
;            pointed to by banddef
;
;  persistdir:  (INPUT/OUTPUT)
;            The directory in which to store the preliminary data products
;            If this directory already exists then proceed from there.
;            Otherwise create this directory, and initialize all the 
;            subdirectories.  If undefined, create a new directory in 
;            'basedir'
;            If undefined then a new persistdir is created and 
;            the name of the dir is stored in the persistdir variable 
;
;  neverwrite: Do not overwrite. The code will crash if the fft's
;  aren't pre-existing. If set, resume is ignored and defaults to on.
;
;            
;  basedir:  Define a place to keep the temporary output of this
;            routine.  Default: Current working directory
;
;
;  mapname:  set to the name of the map within the files being input,
;            if mapname is an array of fits or sav files.  For example, 
;            'mapstruct.map', or 'mapstruct.map[*,*,10]', or 
;            'map[10].map'
;
;   WARNING! This is case-sensitive for Fits files.
;
;
;  kmask:  A fourier-space kmask to apply.  This should be an array of
;         the same size as the input maps.  This can also be used to 
;         weight the fourier bins.  The kmask is applied in an 
;         expression that looks like fft(i)*conj(fft(j))*kmask, and the 
;         effect on the number of modes per output bin is accounted for.
;  
;
; KEYWORD PARAMETERS:
;
;  resume: When used with the persistdir keyword this keyword enables
;          the use of previous fft's and/or cross spectra to resume a 
;          previous calculation.   If a persistdir is specified and
;          the resume keyword is not set then an warning message will 
;          come up before deleting any previous intermediate data
;          products.   Explicitly setting resume=0 will get rid of
;          this message.
;  
;  cmbweighting:  when set to 1, the l(l+1)/2pi weighting is applied
;
;  winfact:  when set to 1, the fsky factor associated with the
;            apodization window is accounted for 
;
;  npix; This sets the size of the lotsomaps format. Defaults to be
;  the same size as the kmask.
;
;  wt_by_unmasked_modes: divide by total(/double,kmask) within a bin instead
;  of nmodes in the bin.  
; *Special case* if total(/double,kmask) eq 0 then the output is set to zero
;
;
; OUTPUTS:
;
;  spectrum:  The band power spectrum
;  covariance:  The covariance matrix for the band power spectum.
;               (Contains no estimates of contributions from
;               sample variance, contains instrumental noise only) 
;  status:  A structure detailing the status of the analysis
;
; COMMON BLOCKS:
;
;  None yet
;
; SIDE EFFECTS:
;
;  Temporary cross spectra data products are stored in a temporary
;  'persist' directory.  These cross spectra can be used to resume a previous
;  (interupted) analysis by using the resume keyword.  It's a good
;  idea to keep track of these persist directories and delete any that 
;  are unneeded. 
; 
;
; RESTRICTIONS:
;
;  None yet known
;
; PROCEDURE:
;
;  Step 0:  Setup the temp directory and the status structure
;  Step 1:  Read in all the maps, window them, take their Fourier
;           Transform and save the transforms to disk.
;  Step 2:  Generate cross spectra from map ffts. 
;           Save Cross spectra to file. 
;  Step 3:  Calculate sample covariance matrix based on available
;           cross spectra.
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
; Feb 26, 2009, MVL
;  Created from unbiased_pspec.pro
; 2010/06/29, ES
;  Included Martin's convolkern keyword
;       for pre-whitening studies.
; 2010/07/16, ES
;       Automatically handle pixel vector format input.
;
; 2012/06/28, CR
;       Changed autospectrum covaraince defn (instead of
;       /nrealizations -> /nrealizations-1)
;
;-


;;;;;
;;
;;  Start with some Miscelaneous macros
;;
;;;;;

function unbiased_multiset_pspec_name_tempdir, basedir
while 1 do begin 
    rand=floor(randomu(seed, 1)*1e6)
    result=basedir+'/unbiased_multiset_pspec-'$
      +strcompress(string(rand, format='(I06)'))
    if file_test(result) eq 0 then break
endwhile
return, result
end

function unbiased_multiset_pspec_make_stat_struct
return, {UNBIASED_MULTISET_PSPEC_STATUS, mapsperset:0L, $
         nsets:0L, ffts_done:0L, ffts_coadded:0L, $
         xpecs_done:0L, auto_specs_done:0L}
end

pro load_map_from_file, mapfiles, filenumber, mapnumber, $
                        dest=fltbuffer, $
                        mapname=mapname, npix=npix, $
                        isfits=isfits, isbin=isbin, issav=issav
s=size(fltbuffer)
if((s[0] ne 2) or (s[1] ne npix[0]) or (s[2] ne npix[1])) then begin
    message, 'destinstion buffer does not equal map size'
endif
test=file_test(mapfiles[filenumber])
if test eq 0 then begin
    message, 'file '+ mapfiles[filenumber]+ ' does not exist' 
endif

if keyword_set(isbin) then begin
    openr, lun, mapfiles[filenumber], /get_lun
    point_lun, lun, (ulong64(npix[0])*ulong64(npix[1])*mapnumber)*4
    readu, lun, fltbuffer
    free_lun, lun
endif else begin
    if (mapnumber ne 0) then begin
        message, $
          "Multiple maps only supported for bin files"
    endif

    if(n_elements(mapname) eq 0) then begin
        
    endif

    if keyword_set(isfits) then begin
        mapnamestub = (strsplit(mapname,'.',/extract))[0]
;;        data=read_spt_fits(mapfiles[filenumber],include=mapnamestub)
        ;;modifed by ES 2010/07/16 to handle pixel vector format maps
        data=read_spt_fits(mapfiles[filenumber],include=[mapnamestub,'PIXELS'])

        if find_matching_tag(data,'PIXELS') eq 'PIXELS' then $
           data=spt_expand_smallmap(data,nsidex=npix[0],nsidey=npix[1])
; check for the variable
        teststring=strcompress('testn=n_elements(data.'+mapname+')',/rem)
        test=execute(teststring)
        mapstring=strcompress('fltbuffer[*]=data.'+mapname,/rem)
    endif
    if keyword_set(issav) then begin
        restore, file=mapfiles[filenumber]
        teststring=strcompress('testn=n_elements('+mapname+')',/rem)
        test=execute(teststring)
        mapstring=strcompress('fltbuffer[*]='+mapname,/rem)
    endif

    if testn eq ulong64(npix[0])*ulong64(npix[1]) then begin
        ok=execute(mapstring)
        if ok ne 1 then begin
            message, 'Could not read map from file: ' +mapfiles[filenumber]
        endif                 ; object wrong size
    endif else begin
        message, 'Desired map not found in file, or is the wrong size'
    endelse                     ; map not in file or is wrong size
endelse ; not flat file
end

pro take_and_reformat_ffts, mapfiles, fftfile_out, npix, nmapsperfile, $
                            nffts_max, window=window, kmask=kmask, $
                            winfact=winfact, cmbweighting=cmbweighting, $
                            ellgrid=ellgrid, ellsidx=ellsidx, $
                            reso=reso, ramlimit=ramlimit, mapname=mapname,$
                            isfits=isfits, isbin=isbin, issav=issav, $
                            convolkern=convolkern


; default ram limit: 16GB
if n_elements(ramlimit) eq 0 then ramlimit=16*(ulong64(2)^30)

winsize=(size(window))[1]

nfiles=n_elements(mapfiles)

;; number of bytes in a Dcomplex: 16
;; number of arrays we need to make to do this efficiently: 6 or less
;; number of pixels in an fft: winsize^2
comp_unit=16*6*ulong64(winsize)^2

parallelism=(fix(ramlimit/comp_unit) > 1)

bufferin=dblarr(winsize, winsize, parallelism)
smallbuffer=fltarr(npix[0], npix[1])
ellfact=ellgrid*(ellgrid+1)/2/!PI

scalingfact=kmask/winfact

if keyword_set(cmbweighting) then begin
    scalingfact*=ellfact
endif

scalingfact=sqrt(abs(scalingfact))

fileidx=0
mapidx=0

scaling3D=rebin(scalingfact, winsize, winsize, parallelism) 
window3D=rebin(window, winsize, winsize, parallelism)
; make a sorting-index array to sort multiple ffts by ell at once
ellsidx3D=rebin(ellsidx, ulong64(winsize)^2, parallelism)
; as it stands ellsidx maps all of the sorted arrays to the
; first array.  We need to add an offset to each row
; (ie add dindgen(parallelism)[i] to the ith row) of this array
ellsidx3D+=transpose(rebin(ul64indgen(parallelism), parallelism, ulong64(winsize)^2))$
  *ulong64(winsize)^2

code=reverse_linefeed_code()
print,''

openw, lun, fftfile_out, /get_lun

totalmaps=nfiles*nmapsperfile
idx=0


  
while ((fileidx lt nfiles) and (mapidx lt nmapsperfile)) do begin

    i=fileidx*nmapsperfile+mapidx

    ;print, code+'processing file: '+strtrim(string(fileidx), 2)+' of '+$
    ;  strtrim(string(nfiles), 2)+', map '+strtrim(string(mapidx), 2)+' of '+$
    ;  strtrim(string(nmapsperfile), 2)+'        '
    
    count=(parallelism < (totalmaps-i))
    if count lt parallelism then begin
        ;; this will be the last pass
        ;; shrink all of our work spaces accordingly 
        bufferin=dcomplexarr(winsize, winsize, count)
        scaling3d=rebin(scalingfact, winsize, winsize, count)
        window3d=rebin(window, winsize, winsize, count)
        ellsidx3d=ellsidx3d[*, 0:count-1]
    endif
    ;; load many maps into a buffer
    for j=0, count-1 do begin 
        load_map_from_file, mapfiles, fileidx, mapidx, $
          dest=smallbuffer, npix=npix, mapname=mapname,$
          isfits=isfits, issav=issav, isbin=isbin
        
        if(n_elements(convolkern) ne 0) then begin
            smallbuffer=convol(smallbuffer, convolkern, /edge_zero)
        endif
        bufferin[0:npix[0]-1, 0:npix[1]-1, j]=smallbuffer
        mapidx++
        if mapidx ge nmapsperfile then begin
            mapidx-=nmapsperfile
            fileidx++
        endif
    endfor
; take the fft of multiple data, applying windows, and
    ; kspace masks all in the same step.
    ; Not making new variables with each step saves ram. 
    thefft=(fft(fft(bufferin*window3D, dim=1, /double), dim=2, /double)*reso^2*winsize^2*scaling3D)
    ; sort the ffts by ell, delete the
    ; unsorted ffts (temporary command) to
    ; save ram 

    thefftsorted=(temporary(thefft))[ellsidx3d]
    writeu, lun, thefftsorted
endwhile
free_lun, lun
end

function load_cross_spectra_data_from_disk,tmpfftfile, winsize, $
                                           nffts, $
                                           start=start,$
                                           stop=stop
nmodes=stop-start+1
result=dcomplexarr(nmodes, nffts)
tmpresult=dcomplexarr(nmodes)
openr, lun, tmpfftfile, /get_lun
for i=0L, nffts-1 do begin
    point_lun, lun, (i*ulong64(winsize)^2+start)*16
;    point_lun, -1*lun, pos
;    print, string(pos/ulong64(3120)^2)
    readu, lun, tmpresult
    result[*, i]=tmpresult
endfor
free_lun, lun
return, result
end

pro generate_jackknife_ffts, processed_fftfile, jackknife_fftfile, winsize, $
                             setdef=setdef
buffera=dcomplexarr(winsize, winsize)
bufferb=dcomplexarr(winsize, winsize)
openr, lun_in, processed_fftfile, /get_lun
openw, lun_out, jackknife_fftfile, /get_lun
nsets=(size(setdef))[2]
setsize=(size(setdef))[1]
for i=0, setsize-1 do begin
    ; read the two relevant files
    point_lun, lun_in, setdef[i, 0]*ulong64(winsize)^2 * 16
    readu, lun_in, buffera
    point_lun, lun_in, setdef[i, 1]*ulong64(winsize)^2 * 16
    readu, lun_in, bufferb
    mapout=(buffera-bufferb)/2
    writeu, lun_out, mapout
endfor
free_lun, lun_in
free_lun, lun_out
setdef=reform(lindgen(setsize), setsize, 1)
end

function take_all_cross_spectra, fftfile, winsize, setdef, banddef, $
                                 bandstartidx=bandstartidx, $
                                 nmodes=nmodes_out, $
                                 wt_by_unmasked_modes=wt_by_unmasked_modes, $
                                 ramlimit=ramlimit, auto=auto, reso=reso, $
                                 no_cross_set=no_cross_set

code=reverse_linefeed_code()

; default ram limit: 16 GB
if n_elements(ramlimit) eq 0 then ramlimit=ulong64(16)*(ulong64(2)^30)

nsets=(size(setdef))[2]
setsize=(size(setdef))[1]

if (keyword_set(no_cross_set)) then nspectra=nsets $
else nspectra=(nsets*(nsets+1))/2

if not keyword_set(auto) then begin
    nrealizations=(setsize*(setsize-1))/2
endif else begin
    nrealizations=(setsize)
endelse
nbands=n_elements(banddef)-1
nffts=max(setdef)+1

allspectra_out=dblarr(nbands, nspectra, nrealizations)

nmodes_out=lonarr(nbands)

;; Step 1, copy all of the fft files and apply scalings masks etc


;; now average all the bands

print, ''

tmpresult=dblarr(setsize, setsize)

i=0
while i lt nbands do begin 
    
    maxnmodes=ramlimit/nffts/16/2

    istop=max(where((bandstartidx-bandstartidx[i]) lt maxnmodes))

    if istop le i then begin
        message, 'Insufficient ram for processing even a single bin'
    endif

    print, code+'loading bands '+strtrim(string(i), 2)$
      +' through '+strtrim(string(istop-1), 2)+'              '
    ; technical: delete the last iteration of banddata_big first
    banddata_big=0
    ;get data for as many bins as will fit in our ramlimit
    banddata_big=load_cross_spectra_data_from_disk(fftfile, winsize, $
                                                   nffts, $
                                                   start=bandstartidx[i],$
                                                   stop=bandstartidx[istop]-1)
    ;process this data
    for iprime=i, istop-1 do begin
        ;print, code+'processing band '+strtrim(string(iprime), 2)+'                 '
        nmodes=(bandstartidx[iprime+1]-bandstartidx[iprime])
        nmodes_out[iprime]=nmodes
        aidx=bandstartidx[iprime]-bandstartidx[i]
        banddata=banddata_big[aidx:(aidx+nmodes-1), *]

        spectrum_idx=0
        for j=0, nsets-1 do begin

            if (keyword_set(no_cross_set)) then kc=j else kc=nsets-1
            for k=j, kc do begin
                if not keyword_set(auto) then begin
                    tmpresult=real_part(banddata[*, setdef[*, j]]##transpose(conj(banddata[*, setdef[*, k]])))/(nmodes*reso^2*winsize^2)
                    ;; symmetrize the result
                    tmpresult/=2
                    tmpresult+=transpose(tmpresult)
                    a=0
                    for l=0, setsize-2 do begin
                        rowlength=setsize-l-1
                        allspectra_out[iprime, spectrum_idx, a:(a+rowlength-1)]=$
                          tmpresult[l, l+1:setsize-1]
                        a+=rowlength
                    endfor
                endif else begin
                    idx=lindgen(setsize)
                    tmpresult=total(/double,real_part(banddata[*, setdef[*, j]]*conj(banddata[*, setdef[*, k]])), 1)/(nmodes*reso^2*winsize^2)
                    allspectra_out[iprime, spectrum_idx, *]=tmpresult
                endelse
                spectrum_idx++
            endfor
        endfor
    endfor
    i=istop
endwhile

return, allspectra_out
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  FUNCTION UNBIASED_MULTISET_PSPEC starts here
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro unbiased_multiset_pspec_test, mapfile, window, reso_arcmin, $
                             ramlimit=ramlimit,$
                             setdef=setdef, $
                             spectrum=spectrum, covariance=cov, $
                             resume=do_resume,$
                             intfile_ident=intfile_ident, $ 
                             status=status, banddef=banddef0, $
                             persistdir=tempdir, basedir=basedir, $
                             cmbweighting=cmbweighting,  $
                             nmodes=nmodes, mapname=mapname, winfact=winfact,$
                             kmask=kmask, nmapsperfile=nmapsperfile,$
                             maxell=maxell, jackknife=jackknife, $
                             verbose=verbose, no_cross_set=no_cross_set, $
                             npix=npix, auto=auto, chisq=chisq,$
                             fftrdonly=fftrdonly, allspectra=allspectra, $
                             meansub=meansub,$
                             wt_by_unmasked_modes=wt_by_unmasked_modes, $
                             est1_cov=cov1, est2_cov=cov2, convol_kernel=convolkern, $
                             delete_intfile=delete_intfile

if keyword_set(wt_by_unmasked_modes) then begin 
    message, 'Sorry wt_by_unmasked_modes has not been re-implemented since the big overhaul'
endif


;;;; 
;;
;;  Step 0: do a lot of housekeeping related to temporary storage of
;;  ffts, cross spectra etc. 
;;  
;;;;

;
; First check to see that the window is square
; 

winsize=size(window)

if(winsize[0] ne 2) then begin 
    message, 'Window should be a 2-D array'
endif 

if(winsize[1] ne winsize[2]) then begin 
    message, 'Window must be square'
endif 

winsize=winsize[1]
mapsizebytes=long(winsize)^2*long(4)

if n_elements(npix) eq 1 then begin
   if npix gt winsize then begin
      message, 'Window must be at least as large as maps'
      return    
   endif
   npix = [npix,npix]
endif else if  n_elements(npix) eq 2 then begin
    if npix[0] gt winsize or npix[1] gt winsize then begin
      message, 'Window must be at least as large as maps'
      return    
   endif
   
endif else npix = [winsize,winsize]
;npix = winsize
inmapsizebytes=long(npix[0])*long(npix[1])*long(4)


; Next check to see 
; 1. If the input is actually a list of files.  
;    If it is an array and any of files are suffixed with .fits and .sav, 
;    we require that the mapname field also be filled.  
; 2. Does this list of files represent a multiple sets of files?
;    This is complicated by the fact the fits and save files can only
;    contain one map, whereas lotsomaps files are always
;    multiplefiles.  The number of maps per file is assumed to be 
;    
;   
;    Any array of lotsomaps files respresents multiple sets.

isfitsfile=stregex(mapfile, '\.fits$', /boolean)
issavfile=stregex(mapfile, '\.sav$', /boolean)
isbinfile=stregex(mapfile, '\.bin$', /boolean) or $
  stregex(mapfile, '\.dat$', /boolean) 
tmp=where(isfitsfile, cnt)
fitsfile=(cnt ne 0)
tmp=where(issavfile, cnt)
savfile=(cnt ne 0)
tmp=where(isbinfile, cnt)
binfile=(cnt ne 0)

tmp=where((isbinfile+isfitsfile+issavfile eq 0), cnt)
otherfile=(cnt ne 0)

mixture=((fitsfile+savfile+binfile) ne 1)

if (otherfile or mixture) then begin
    message, 'Files must be .fits, .sav, or .bin and must be all the same format'
endif

;
;  Next step: Find out how many files there are and how many maps per
;  file there are.
;

if n_elements(nmapsperfile) eq 0 then begin
    if (fitsfile or savfile) then begin
        nmapsperfile=1
    endif else begin
        ;; for binary files:
        ;; first check that all the files are the same size
        fsize=(file_info(mapfile)).size
        uidx=uniq(fsize, sort(fsize))
        if n_elements(uidx) gt 1 then begin
            message, 'Binary files must all be the same size!'
        endif
        fsize=fsize[0]
        if (fsize mod inmapsizebytes) ne 0 then begin 
            message, 'Binary file does not contain an integer number of maps'
        endif
        nmapsperfile=fsize/inmapsizebytes
    endelse 
endif 

nfiles=n_elements(mapfile)

;
;  Next step: Determine how these files are to be divided into sets.
;
;  Default behaviour:  If setdef is not set assume the following:
;  1.) Single map files can be listed in a 1-D array, are treated as a
;      single set
;  2.) Single map files can be listed in 2-D arrays to arrange them as
;      sets.
;  3.) A lone multi-map file can represent a single set.
;  4.) Multi-map files represent a list of sets if listed in a 1-D
;      array.  In this case aach file is a set.
;  5.) Multi-map files can also be listed as 2-D arrays and  each
;       row of the array specifies a set of size (nfiles-in-row)*nmapsperfile
;
;  If setdef is defined: 
;    Then the details of the file arrangement has no
;    bearing on the number of maps per set.  The maps are assigned an
;    index, based on the number of maps per file and order of the files.
;    2-D arrays of mapfiles are flattened to a 1-D array using 
;    mapfile=reform(mapfile, n_elements(mapfile)).   
;   
;    The maps in the first array are numbered 0, 1, 2... nmapsperfile-1
;    The maps in the second array are numbered nmapsperfile,
;    nmapsperfile+1, ... 2*nmapsperfile-1.  and so forth      
;  
;  Setdef should be a 1-D or 2-D, listing which maps (by index) fall into which
;  set.  1-D set definitions define a single set to be analyzed
;  2-D set definitions break the files into multiple sets (one set per row).


s=size(mapfile)

if n_elements(setdef) eq 0 then begin
    if(nmapsperfile eq 1) then begin
        if(s[0] eq 1) then begin
            setdef=reform(dindgen(s[1]), s[1], 1)
        endif 
        if(s[0] eq 2) then begin
            setdef=reform(dindgen(s[1]*s[2]), s[1], s[2])
        endif
    endif else begin
        if(s[0] eq 0) then begin
            setdef=reform(dindgen(nmapsperfile), nmapsperfile, 1)
        endif
        if(s[0] eq 1) then begin
            setdef=reform(dindgen(nmapsperfile*s[1]), nmapsperfile, s[1])
        endif
        if(s[0] eq 2) then begin 
            setdef=reform(dindgen(nmapsperfile*s[1]*s[2]), $
                          nmapsperfile*s[1], s[2])
        endif
    endelse
endif

;; check to see if setdef is defined NOW, if not then it means that
;; the user did not explicitely define a set definition, and the 
;; implicit set definition did not work out because either:
;; 1. Only one file, with a single map was specified, or 
;; 2. The array defining the map list had more then 2 dimensions 
;;
if n_elements(setdef) eq 0 then begin 
  message, 'Sorry I could not understand your set definition.  Map files should be arranged in 1-D or 2-D arrays'
endif 

;; for convenience flatten the mapfile.  We have already 
;; recorded the set definitions implicit in its original format.
mapfile=reform(mapfile, n_elements(mapfile))

s=size(setdef)
;; recast 1-D set definitions as 2-D arrays
if ( s[0] eq 1) then setdef=reform(setdef, s[1], 1)

s=size(setdef)
nsets=s[2]
setsize=s[1]

if keyword_set(jackknife) then begin
    if nsets ne 2 then begin
        message, $
          'Direct jackknife comparision of maps requires that the maps be broken up into exactly two data sets'
    endif
    nsets=1 ; make only one jacknife spectrum                                                              
    no_cross_set=1              ;  Do not do make cross-set-spectra.  
;   setsize=setsize/2
endif

; check to see if the tempdir is defined, if not define a new one
s=size(basedir, /type)
if s ne 7 then cd, '.', c=basedir

s=size(tempdir, /type)
if s ne 7 then begin 
                                ; check to see of the basedir is
                                ; defined, if not used the current
                                ; directory    
    tempdir=unbiased_multiset_pspec_name_tempdir(basedir)
    file_mkdir, tempdir
endif

; check to make sure that the directory exists and is
; a directory
if file_test(tempdir, /directory) eq 0 then begin
    ; does something else exist here that is not a directory?
    if file_test(tempdir) ne 0 then begin 
        print, 'WARNING: persistdir "'+tempdir+'" is not a directory'
        tempdir=unbiased_multiset_pspec_name_tempdir(basedir)
        print, 'WARNING: using directory "', tempdir, '" instead'
    endif
    file_mkdir, tempdir
endif 

if (keyword_set(intfile_ident)) then begin
    statfile=tempdir+'/status_'+intfile_ident+'.'+mapname+'.sav'
    fftfile=tempdir+'/ffts_'+intfile_ident+'.'+mapname+'.bin'
    processedfftfile=tempdir+'/ffts_processed_'+intfile_ident+'.'+mapname+'.bin'
    sumfftfile=tempdir+'/sumffts_'+intfile_ident+'.'+mapname+'.bin'
    xpecfile=tempdir+'/xspecs_'+intfile_ident+'.'+mapname+'.bin'
endif else begin
    statfile=tempdir+'/status.'+mapname+'.sav'
    fftfile=tempdir+'/ffts.'+mapname+'.bin'
    processedfftfile=tempdir+'/ffts_processed.'+mapname+'.bin'
    sumfftfile=tempdir+'/sumffts.'+mapname+'.bin'
    xpecfile=tempdir+'/xspecs.'+mapname+'.bin'
endelse

if not keyword_set(do_resume) then begin 
    if file_test(statfile) then begin
        if n_elements(do_resume) eq 0 then begin
            ; resume keyword was not EXPLICITLY SET to 0
            print, 'WARNING: Persistdir "'+tempdir+'" already contains xspec data '
            print, '         from a previous run of unbiased_multiset_pspec'
            print, 'WARNING: This data can be used to speed up the calculation of'
            print, '         the power spectrum'
            print, 'WARNING: To use this data, hit Ctrl-C now and use the /resume ', $
              'keyword'
            print, 'WARNING: Otherwise, this data will be deleted in 15 seconds' 
            print, 'NOTICE : To skip this message in the future, use resume=0'
            wait, 15.
        endif
        print, 'deleting old files'
        file_delete, statfile
        if file_test(fftfile) then file_delete, fftfile
        if file_test(xpecfile) then file_delete, xpecfile
    endif
endif

; delete the contents of the variable status
status='bogus_value'
if(file_test(statfile) eq 1) then begin
    restore, statfile
endif 

; make sure that the status variable is valid
s=size(status, /type)
;stop
if s ne 8 then begin   
    status=unbiased_multiset_pspec_make_stat_struct()
endif

n=tag_names(status, /str)
if n ne "UNBIASED_MULTISET_PSPEC_STATUS" then begin
    print, 'WARNING: status structure is not correctly defined, starting over'
    status=unbiased_multiset_pspec_make_stat_struct()    
endif

status.nsets=nsets
status.mapsperset=setsize

save, status, file=statfile

rev_linefeed=reverse_linefeed_code()

;; blank print statement
print

nfiles=n_elements(mapfile)
nffts_max=nmapsperfile*nfiles

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Generate bookkeeping in the to take xspectra of all available ffts
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

reso=double(reso_arcmin*!dtor/60)
deltau=1./reso/(winsize)

if n_elements(banddef0) eq 0 then begin
    if(n_elements(maxell) eq 0) then maxell=deltau*winsize*!PI
    maxell=min([maxell, deltau*winsize*!PI])
    print, 'constructing default bands.  deltau:', deltau, ', maxell:', maxell
    nbands=floor(maxell/!PI/deltau/2)
    banddef=dindgen(nbands+1)*!PI*deltau*2
    ; write the bands to bandef0, where they can be returned
    ; to the user
    banddef0=banddef[1:*]
endif else begin
    banddef=[0, banddef0]
endelse

nbands=n_elements(banddef)-1

bandcenters=(banddef[0:nbands-1]+banddef[1:nbands])/2

bandellfact=bandcenters*(bandcenters+1)/2/!PI


ell=shift(dindgen(winsize)-winsize/2, winsize/2)*deltau*2*!PI
ellgridx=(ell#replicate(1.0, winsize))
ellgrid=ellgridx^2+transpose(ellgridx)^2
ellgrid=sqrt(ellgrid)

;; while we're thinking about weights and factors, also calculate the 
;; factor that accounts for the effect of the window on fsky
if keyword_set(winfact) then begin
    winfactor=total(/double,window^2)/(double(winsize)^2)
endif else begin
    winfactor=1.0d
endelse

;; final potential weighting:  set up fourier space kmask:

;kmask  (it might be fine not to bother with doubles for this variable
;but just to play it safe I'm leaving everything in doubles that I can)

if n_elements(kmask) eq 0 then  begin
    kkmask = dblarr(winsize,winsize)+1.
endif else begin
    kkmask = double(kmask)
endelse

;; We are going to binning up MANY MANY cross spectra since it takes a
;; while to figure out which fft elements go into which bin, we can
;; buy a lot of efficiency by mapping out which fft bins are in which
;; band once and for all.
;;
;; We start by sorting all of the fft bins by |u|
;; then we iterate through them, noting the index in the sorted list
;; where each band begins (and the previous one ends)


ellsidx=sort(ellgrid, /L64)

bandstartidx=lonarr(nbands+1)
i=0

bandstartidx[0]=1
for h=1L, n_elements(ellgrid) do begin
                                ; Then banddef is DEFINED
                                ; to be the upper edge of each bin, 
                                ; and is included in that bin. 
    if ellgrid[ellsidx[h]] gt banddef[i+1] then begin
        bandstartidx[i+1]=h
        i++
    endif
    if i ge nbands then break
endfor

info=file_info(processedfftfile)
desired_size=nffts_max*16*ulong64(winsize)^2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  This is where we start by computing the FFTs for all files
;;  Before we put the ffts to disk, we mask them
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

if((not keyword_set(do_resume)) or (info.size lt desired_size)) then begin
    take_and_reformat_ffts, mapfile, processedfftfile, npix, nmapsperfile, $
      nffts_max, window=window, kmask=kkmask, winfact=winfactor, reso=reso, $
      cmbweighting=cmbweighting, ellgrid=ellgrid, ellsidx=ellsidx,$
      mapname=mapname, isfits=fitsfile, issav=savfile, isbin=binfile, $
      ramlimit=ramlimit, convolkern=convolkern
endif

if keyword_set(jackknife) then begin
    jackknifefftfile=processedfftfile+'.jackknife'
    generate_jackknife_ffts, processedfftfile, jackknifefftfile, winsize, $
      setdef=setdef
    ; from now on we will use the jackknife file for computation
    processedfftfile=jackknifefftfile

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  This is where we actually compute the xspectra
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 

allspectra=take_all_cross_spectra(processedfftfile, winsize, setdef, banddef, $
                                  bandstartidx=bandstartidx, $
                                  nmodes=nmodes, $
                                  wt_by_unmasked_modes=wt_by_unmasked_modes, $
                                  auto=auto, reso=reso, no_cross_set=no_cross_set)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;;  Generate means and sample variances 
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

print, "Correlating Cross Spectra"

if (keyword_set(no_cross_set)) then nspectra=nsets $
else nspectra=(nsets*(nsets+1))/2

spectrumreformed=dblarr(nbands,nspectra)
if not keyword_set(auto) then begin
    nrealizations=(1l*setsize*(setsize-1))/2
endif else begin
    nrealizations=(setsize)
endelse
allspectra=reform(allspectra, nbands*nspectra, nrealizations)
cov=dblarr(nbands*nspectra, nbands*nspectra)

spectrum=total(/double,allspectra, 2)/nrealizations
spectrum_2d=rebin(spectrum, nbands*nspectra, nrealizations)
cov1=(transpose(allspectra-spectrum_2d)##(allspectra-spectrum_2d))/nrealizations/(nrealizations-1)

if not keyword_set(auto) then begin

    realization_to_complement=dblarr(nrealizations, setsize)
    
    for i=0, setsize-1 do begin 
        realization_idx=0L 
        for j=0, setsize-1 do begin 
            for k=j+1, setsize-1 do begin 
                if ((i eq j) or (i eq k)) then begin
                    realization_to_complement[realization_idx, i]=1./(setsize-1) 
                endif 
                realization_idx++ 
            endfor 
        endfor 
    endfor 
    
    allcomplementspectra=realization_to_complement##allspectra
    spectrum_2d=rebin(spectrum, nbands*nspectra, setsize)
    
    cov2=transpose(allcomplementspectra-spectrum_2d)##(allcomplementspectra-spectrum_2d)
    cov2/=(long(setsize)^2/2)
    cov=2*cov2-cov1
endif else begin
    cov=cov1*(nrealizations) ;CR - June 28, 2012
endelse

if (keyword_set(delete_intfile)) then begin
    ;statfile=tempdir+'/status_'+intfile_ident+'.'+mapname+'.sav'
    ;fftfile=tempdir+'/ffts_'+intfile_ident+'.'+mapname+'.bin'
    ;processedfftfile=tempdir+'/ffts_processed_'+intfile_ident+'.'+mapname+'.bin'
    ;sumfftfile=tempdir+'/sumffts_'+intfile_ident+'.'+mapname+'.bin'
    ;xpecfile=tempdir+'/xspecs_'+intfile_ident+'.'+mapname+'.bin'
    file_delete, [statfile, processedfftfile]
endif

end

