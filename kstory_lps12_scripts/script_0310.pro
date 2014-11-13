;;;
; NAME: script_0310
; PURPOSE:
;   General script for today
;
; NOTES:
; 1) 
; 2) make nmiss arrays
;
; MODIFICATION HISTORY:
;  03/10/2012: (KTS) Created
;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; apod masks
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;...................................................................
; Check overlap in apod masks between fields
pro make_all_apodmasks_0310
compile_opt IDL2, HIDDEN



end


;...................................................................
; Check overlap in apod masks between fields
pro plot_fullmap_overlap
compile_opt IDL2, HIDDEN

restore, '/home/kstory/lps12/scripts/apod_test1_0305.sav'

;------------------------
; Setup
;------------------------
field_arr = lps12_fieldstruct()
reso_arcmin = 1.0


;------------------------
; Make big map
;------------------------

radec0 = [15, -55]
npixels = [floor(960*5), floor(960*4)]
bigmap = fltarr(npixels[0], npixels[1])

for ii=0, 19 do begin

    fst = field_arr[ii] & field_name = fst.name & npixels = get_lps12_map_npix(ii) & radec0 = [ fst.ra0, fst.dec0 ] ; Get the apod masks
    restore, '/home/kstory/lps12/masks/apod_0305/apod_'+field_name+'_60_0.0500_40.sav'

    pix2ang_proj5, npixels, radec0, reso_arcmin, ra, dec
    
    ; insert apod into bigmap
    print, "loop ", ii
    npix_sm = long(npixels[0]) * long(npixels[1])
    for ii=0, npix_sm-1 do begin
        if (apod[ii] ne 0) then begin
            ang2pix_proj5, ra[ii], dec[ii], npixels, radec0, 1.0, ipix
            bigmap[ipix] += apod[ii]
        endif
    endfor

endfor


;------------------------
; Plots
;------------------------

; Plot 1: show overlap greater than 1
ss = bigmap*0
wh = where(bigmap gt 1)
ss[wh] = 1
tv_spt_map, ss, /norms, scale=0.2, /forcesize, winnum=2, title='show where above 1'

; Plot 2: show where equal to 0
; ss0 = bigmap*0
; wh = where(bigmap le 0)
; ss0[wh] = 1
; tv_spt_map, ss0[*,150:*], /norms, scale=0.6, /forcesize, winnum=2, title='show where equal to 0'

; Plot 3: apod masks
tv_spt_map, bigmap, /norms, scale=0.2, /forcesize, winnum=1, title='5h-7h, dec-65 to -50'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_0.5overlap_0306'
;figname = '/home/kstory/lps12/masks/figs_0305/apod_Nooverlap_0306'
figname = '/home/kstory/lps12/masks/figs_0305/apod_1deg_0306'
;err= tvread(/png, filename=figname, /nodialog)

; save output for easy recycling
save, apod1, apod2, apod3, bigmap, ss, filename='apod_test1_0305.sav'
stop
end

