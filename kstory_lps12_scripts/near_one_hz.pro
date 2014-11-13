;;;
; NAME: near_one_hz
; PURPOSE:
;   Find bolometers that have excess power near 1hz
;
; CALLING SEQUENCE: near_one_hz, 17
;
; INPUTS: 
;   field_idx,           index in lps12_fieldstruct()
;
; OUTPUTS:
;   sav file at '/data/kstory/projects/lps12/1hzBolo/near_one_hz_'+field_name+'.sav'
;
; NOTES: (from RK)
;   1) Default is to use 20 observations, which are chosen randomly
;   from the runlist
;
; MODIFICATION HISTORY:
;   02/09/2012: (KTS) Created from /data/rkeisler/lines/near_one_hz.pro
;;;

pro near_one_hz, field_idx, nobs=nobs, dosave=dosave, everyobs=everyobs
compile_opt IDL2, HIDDEN

savecount = -1
if n_elements(nobs) eq 0 then nobs=20

; Get the IDF directory
field_arr = lps12_fieldstruct()
fitsdir = field_arr[field_idx].idf_lpfds_dirs[0]
;fitsdir = field_arr[field_idx].idf_dirs[0]
field_name = field_arr[field_idx].name

print, "Analyzing field, ", field_name

; get the low ell run list
get_runlist_by_field, field_idx, list=lowlist
list = fitsdir + 'field_scan_150_'+lowlist+'.fits'
list = scramble(list) ; scramble the list order

nlist = n_elements(list)
date = extract_date_from_filename(list)
mjd = date_string_to_mjd(date)

if keyword_set(everyobs) then nobs=nlist

; npts Is the number of TOD data points do use
;npts0 = 2.^floor(alog(7e4)/alog(2.))
npts0 = 2.^ceil(alog(30000.)/alog(2.))
npts = npts0

sample_rate = (100./6.)
this_hz = make_fft_grid(1./sample_rate,npts)
whnear1 = where(abs(abs(this_hz)-1.0) lt 0.1,nnear1)
hz_fit = abs(this_hz[whnear1])

hann = hanning(npts)

; Initialize data arrays
nnotch = fltarr(nobs)
notch_freq = fltarr(nobs,5)
notch_width = fltarr(nobs,5)

nnotch2 = fltarr(nobs)
notch_freq2 = fltarr(nobs,5)
notch_width2 = fltarr(nobs,5)


hz = fltarr(nobs,npts0)
p = fltarr(nobs,npts0)

d=read_spt_fits(list[0])
nbolostot = n_elements(d.bolodata.adc[0,*])
boloind = d.observation.BOLO_READOUT_IDX

bolowts = fltarr(nobs,nbolostot)

fitamp = fltarr(nobs,nbolostot)
fitfreq = fltarr(nobs,nbolostot)
fitfwhm = fltarr(nobs,nbolostot)

;-----------------------
; Loop over observations
;-----------------------
list_idx = -1 ; start a 0
for i=0,nobs-1 do begin
    list_idx++
    ; Initialize
    pp = dblarr(npts) ; total power

    ; read in data
    found_file = 0
    while ~found_file do begin
        filename = list[list_idx]
        print,'=== processing obs ',strtrim(i,2),'/',strtrim(nobs-1,2),'==='
        print,'    file: ' + filename
        timea=systime(/sec)
        d=read_spt_fits(list[list_idx])
        CATCH, Error_status
        if Error_status eq 0 then begin
            found_file=1
        endif else begin ; file not found, go to next file
            print, " ---- caught error reading in file, ", filename
            list_idx++
        endelse
    endwhile

    timeb = systime(/sec)
    print,'...reading in fits file took ',timeb-timea,' seconds.'

    if n_elements(d.bolodata.adc[*,0]) lt npts0 then continue

    bf = d.observation.bolo_flags
    whgb = where(bf eq 0, ngb)
    
    ; get bolo weights
    bolowts[i,*] = (d.observation.BOLO_PSD_1TO3*d.observation.bolo_relcal)^(-2.)
    
    ; calibrate timestreams
    ; tt = systime(/sec)
    for j=0,ngb-1 do $
      d.bolodata.adc[0:npts-1,whgb[j]] *= (d.observation.bolo_relcal)[whgb[j]]
    ; print,systime(/sec)-tt
    
    ;-----------------------
    ; Loop over Bolometers
    ;-----------------------
    if 0 then begin
    ; since we're just interested in line centers and widths, let's try
    ; coadding all (good) timestreams and then taking a single fft, rather
    ; than taking coadding the fft's.
        yy = dblarr(npts)
        for j=0,ngb-1 do begin
;            print,strtrim(j,2),'/',strtrim(ngb-1,2)
            k = whgb[j]
            y = reform(d.bolodata.adc[0:npts-1,k])
            yy += y
;            this_p = abs(fft(y*hann))^2.
;            pp += this_p
;        plot,this_hz,this_p,/yl,xr=[0,6] & wait,.5
        endfor
        yy *= (1./ngb)
        pp = abs(fft(yy*hann))^2.
        pp *= npts ;normalize fft
    endif else begin
        for j=0,ngb-1 do begin
            ;print,strtrim(j,2),'/',strtrim(ngb-1,2)
            k = whgb[j]
            y = reform(d.bolodata.adc[0:npts-1,k])
            this_p = abs(fft(y*hann))^2.
            this_p *= npts ;normalize fft

            ; fit for line near 1 Hz
            p_fit = this_p[whnear1]
            this_result = GAUSSFIT(hz_fit, p_fit, aa, nterms=4)
            this_amp = aa[0]
            this_freq = aa[1]
            this_sigma = aa[2]
            this_fwhm = this_sigma*2.355

            fitamp[i,k] = this_amp
            fitfreq[i,k] = this_freq
            fitfwhm[i,k] = this_fwhm

            pp += this_p

            ;; Watch plots of ffts
            ;plot,this_hz,this_p,/yl,xr=[0,6] & wait,.5

        endfor
    endelse
    
    ; search for notched lines in 0 < f < 5 Hz.
    ;wh15 = where(this_hz ge 1.0 and this_hz le 5.)
    wh15 = where(this_hz ge 0.5 and this_hz le 5.)
    medp = median(pp[wh15])
    wh_notch = intersect(wh15,where(pp le medp/1e2))
    nwh_notch = wh_notch[0] eq -1 ? 0 : n_elements(wh_notch)
    
    this_nnotch = 0
    
    if nwh_notch-1 gt 0 then begin
        
        start_freq = this_hz[wh_notch[0]]
        last_freq = this_hz[wh_notch[0]]
        
        for j=1L,nwh_notch-1 do begin
            this_freq = this_hz[wh_notch[j]]
            if (this_freq-last_freq) gt 0.1 or j eq nwh_notch-1 then begin
                
                this_nnotch++
                if j eq nwh_notch-1 then last_freq = this_freq
                
                this_width = last_freq-start_freq
                this_center = 0.5*(last_freq+start_freq)
                
                notch_width2[i,this_nnotch-1] = this_width
                notch_freq2[i,this_nnotch-1] = this_center
                
                start_freq = this_freq
                if this_nnotch eq 5 then goto,done_w_search
            endif
            last_freq = this_freq
        endfor
        
    endif
done_w_search:
    nnotch2[i] = this_nnotch
    
    ;stop
    
    ntmp = min([npts0,npts])
    hz[i,0:ntmp-1] = this_hz[0:ntmp-1]
    p[i,0:ntmp-1] = pp[0:ntmp-1]
    
    timec = systime(/sec)
    print,'...processing took ',timec-timeb,' seconds.'
    
    if keyword_set(dosave) and (i mod 10) eq 0 then begin
        savecount++
        evenodd = is_even(savecount) ? '_even' : '_odd'
        save,hz,p,nobs,list,mjd,nlist,nnotch,notch_freq,notch_width,nnotch2,notch_freq2,notch_width2,fitamp,fitfwhm,fitfreq,boloind,bolowts,$
          filename='/data/kstory/projects/lps12/1hzBolo/near_one_hz_'+field_name+evenodd+'.sav'
        timed = systime(/sec)
        print,'...saving took ',timed-timec,' seconds.'

    endif
    
endfor


; Save the output
if keyword_set(dosave) then begin

    save,hz,p,nobs,list,mjd,nlist,nnotch,notch_freq,notch_width,nnotch2,notch_freq2,notch_width2,fitamp,fitfwhm,fitfreq,boloind,bolowts,$
      filename='/data/kstory/projects/lps12/1hzBolo/near_one_hz_'+field_name+'.sav'

endif

end




