;;;
; NAME: get_lps12_kweight.pro
; PURPOSE:
;   Return a 2-D kweight
;
; INPUTS: field_idx,      index in lps12_fieldstruct()
;
; OUTPUTS:
;   kweight,         2-D weight mask
;
; NOTES:
;   1) give option to use box cut of all kx < ell=500
;
; MODIFICATION HISTORY:
;  03/07/2012: (KTS) Created
;  05/08/2012: (KTS) use actual twod_kweight
;  06/05/2012: (KTS) Re-arrange setup into box-mask section.
;;;


;...................................................................
; Return the kweight
function get_lps14_kweight_mdw, field_idx, band, box_kmask=box_kmask, mask_spots=mask_spots

f = lps12_fieldstruct()
home = getenv('HOME')

; Get the kweight
if ~keyword_set(box_kmask) then begin
    ;restore, '/home/kstory/lps12/twod_kweights/weight_2d_'+f[field_idx].name+'.sav'
    restore, home+'/data/spt_data/kweights/weight_2d_'+f[field_idx].name+'_'+band+'.sav'

    ;restore, '/home/kstory/lps12/twod_kweights/weight_2d_'+f[field_idx].name+'.sav'
    kweight = weight_2d

    if keyword_set(mask_spots) then begin
        reso_arcmin=1.0
        reso=reso_arcmin/60.*!dtor
    
        info = get_lps12_fieldinfo(field_idx)
        npix = info.npix
        nbig = info.nbig
        nbigx=npix[0]
        nbigy=npix[1]

        if field_idx eq 8 then begin

            xellmin = 720       ;550
            xellmax = 930     ;1250
            yellmin = 420
            yellmax = 690       ;800
            ;ellcut = 340

        endif else begin
            if field_idx eq 7 then begin
            
            xellmin = 670       ;550
            xellmax = 1080      ;1250
            yellmin = 400
            yellmax = 770       ;800
            ;ellcut = 340

            endif else begin

                if field_idx eq 13 then begin
                    
                    xellmin = 800 ;550
                    xellmax = 850 ;1250
                    yellmin = 480
                    yellmax = 630
 ;800
                    ;ellcut = 340

                endif else begin
                    if field_idx eq 18 then begin
                           xellmin = 730       ;550
                           xellmax = 1000 ;1250
                           yellmin = 390
                           yellmax = 720 ;800
                           ;ellcut = 340
                    endif else begin

                xellmin = 600       ;550
                xellmax = 600   ;1250
                yellmin = 600
                yellmax = 600   ;800
                ;ellcut = 340

            endelse
        endelse
        endelse
        endelse

        dell  = 2*!pi/(nbig*reso)
        ;delly = 2*!pi/(nbigy*reso)
        jjxmin = fix(xellmin/dell)
        jjxmax = fix(xellmax/dell)
        jjymin = fix(yellmin/dell)
        jjymax = fix(yellmax/dell)
	;jj = fix(ellcut/dell)

        print,jjxmin,jjxmax

        ;kweight[*,*]=1.

        kweight[jjxmin:jjxmax,jjymin:jjymax]=0.
        kweight[nbig-jjxmax:nbig-jjxmin,nbig-jjymax:nbig-jjymin]=0.

        kweight[jjxmin:jjxmax,nbig-jjymax:nbig-jjymin]=0.
        kweight[nbig-jjxmax:nbig-jjxmin,jjymin:jjymax]=0. 
    	
	;kweight[0:jj,*]=0
    	;kweight[nbig-jj:nbig-1,*]=0

        ;tv_spt_map,shift(kweight,nbig/2,nbig/2),min=0.,max=1.,res=5.,scale=0.4/3.,/force

    endif


; if asked, just use a box mask
endif else begin
    reso_arcmin=1.0
    reso=reso_arcmin/60.*!dtor
    
    info = get_lps12_fieldinfo(field_idx)
    npix = info.npix
    nbig = info.nbig
    
    ell_cut = 300.

    dell    = 2*!pi/(nbig*reso)
    jj      = fix(ell_cut/dell)
    kweight = fltarr(nbig,nbig)+1
    
    kweight[0:jj,*]=0
    kweight[nbig-jj:nbig-1,*]=0
endelse

return, kweight
end

