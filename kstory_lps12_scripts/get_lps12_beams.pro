;;;
; NAME: get_lps12_beam.pro
; PURPOSE:
;   Return a 2d vector of the beam
;
; CALLING SEQUENCE:
;      IDL> make_noise_psd_lps12
;
; INPUTS: year [int]
;    l_beam,         empty, returned vector of l-values
;    bl_beam,        empty, returned vector of beam values
;
; OPTIONAL INPUTS:
;    beamfiles,      return the beam file if requested
;
; OUTPUTS: 
;    beamfiles,      The string path to the beam files used
;
; NOTES:
; 1) beams for 2010, 2011 are preliminary (as of 0427)
; 2) sim_run_03 => 2009 beams for 09, 10, 11.  Using ellcut=4500 sims
;
; MODIFICATION HISTORY:
;  04/27/2012: (KTS) Created with preliminary beams for 2010, 2011
;  05/18/2012: (KTS) Add option 'sim_run_03'
;  06/05/2012: (KTS) Add b11/rev3.0 beams
;  06/28/2012: (KTS) Add b11/rev3.1 beams (run_08)
;  07/16/2012: (KTS) Add change directory to /data/sptdat/beams/rev3.1
;  08/27/2012: (KTS) Fix sim beams, copy beams to my directory
;;;


;...................................................................
; Return the beam
FUNCTION get_lps12_beams, year, l_beam, bl_beam, sim_lmax4500=sim_lmax4500, sim_lmax8000=sim_lmax8000, $
                          d08=d08, s08=s08, $
                          d09=d09, s09=s09


;------------------------
; run_09 beams
;------------------------
; data beams
if keyword_set(d09) then begin
    syear = strtrim(string(year),2)
    beamfiles = '/data/kstory/projects/lps12/beams/run_09/data_09/blgrid_'+syear+'_150.txt'
    if ~file_test(beamfiles) then begin
        print, 'GET_LPS12_BEAM: year, (', year, ') failed.  Returning -1'
        RETURN, -1
    endif
endif    

; sim beams
if keyword_set(s09) then begin
    syear = strtrim(string(year),2)
    beamfiles = '/data/kstory/projects/lps12/beams/run_09/sims_09/blgrid_'+syear+'_150.txt'
    if ~file_test(beamfiles) then begin
        print, 'GET_LPS12_BEAM: year, (', year, ') failed.  Returning -1'
        RETURN, -1
    endif
endif    


;------------------------
; Older, obsolete beams
;------------------------

; Used in run_08
if keyword_set(d08) then begin
    case year of
        2008 : beamfiles = '/data/sptdat/beams/rev3.1/blgrid_2008_150.txt'
        2009 : beamfiles = '/data/sptdat/beams/rev3.1/blgrid_2009_150.txt'
        2010 : beamfiles = '/data/sptdat/beams/rev3.1/blgrid_2010_150.txt'
        2011 : beamfiles = '/data/sptdat/beams/rev3.1/blgrid_2011_150.txt'
        else : begin
            print, 'GET_LPS12_BEAM: year, (', year, ') not recognized.  Returning -1'
            RETURN, -1
        endelse
    endcase
endif

; Used in run_05 to run_07
; case year of
;     2008 : beamfiles = '/home/rkeisler/b11/rev3.0/blgrid_2008_150.txt'
;     2009 : beamfiles = '/home/rkeisler/b11/rev3.0/blgrid_2009_150.txt'
;     2010 : beamfiles = '/home/rkeisler/b11/rev3.0/blgrid_2010_150.txt'
;     2011 : beamfiles = '/home/rkeisler/b11/rev3.0/blgrid_2011_150.txt'
;     else : begin
;         print, 'GET_LPS12_BEAM: year, (', year, ') not recognized.  Returning -1'
;         RETURN, -1
;     endelse
; endcase

;------------------------
; Deal with special cases
;------------------------

; Simulations used prelimitary beams
if keyword_set(sim_lmax8000) then begin
    case year of
        2008 : beamfiles = '/data/rkeisler/beams2009/prelim_release_19Sep2010/blgrid_2008_150.txt'
        2009 : beamfiles = '/data/rkeisler/beams2009/prelim_release_19Sep2010/blgrid_2009_150.txt'
        2010 : beamfiles = '/home/rkeisler/b9/blgrid_2011_150_early27Apr2012.txt'
        2011 : beamfiles = '/home/rkeisler/b9/blgrid_2011_150_early27Apr2012.txt'
        else : begin
            print, 'GET_LPS12_BEAM: year, (', year, ') not recognized.  Returning -1'
            RETURN, -1
        endelse
    endcase
endif


; sim_lmax4500 used 2009 beams for 2009, 2010 and 2011
if keyword_set(sim_lmax4500) then begin
    case year of
        2008 : beamfiles = '/data/rkeisler/beams2009/prelim_release_19Sep2010/blgrid_2008_150.txt'
        2009 : beamfiles = '/data/rkeisler/beams2009/prelim_release_19Sep2010/blgrid_2009_150.txt'
        2010 : beamfiles = '/data/rkeisler/beams2009/prelim_release_19Sep2010/blgrid_2009_150.txt'
        2011 : beamfiles = '/data/rkeisler/beams2009/prelim_release_19Sep2010/blgrid_2009_150.txt'
        else : begin
            print, 'GET_LPS12_BEAM: year, (', year, ') not recognized.  Returning -1'
            RETURN, -1
        endelse
    endcase
endif

; if nothing matched:
if n_elements(beamfiles) eq 0 then begin
    print, 'GET_LPS12_BEAM: something failed.  Returning -1'
    RETURN, -1
endif

readcol,beamfiles,l_beam,bl_beam,format='d,d'
RETURN, beamfiles
END
