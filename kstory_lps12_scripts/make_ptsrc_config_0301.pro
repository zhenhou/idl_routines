;;;
; NAME: make_ptsrc_config_0301
; PURPOSE:
;   make point-source config files
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/17/2012: (KTS) Created
;  02/20/2012: (KTS) added procedure add_rk_ptsrc_to_config
;  03/01/2012: (KTS) call add_rk_ptsrc_to_config automatically
;;;

;--------------------------------
; function to call make_ptsrc_config_single_field
; on all fields
pro make_ptsrc_config_0301

for ii=0, 20 do begin
    make_ptsrc_config_single_field, ii
endfor
end


;--------------------------------
; Make the ptsrc_config file for a single field
; indexed by lps12_fieldstruct
;
pro make_ptsrc_config_single_field, field_idx

; get the field properties
field_arr = lps12_fieldstruct()

field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name

tmp = read_spt_fields()
wh = where( strcmp(tmp.name, field_dir_name) eq 1, nwh)
if(nwh eq 0) then begin
    print, "field name miss-match"
    return
endif

fstruct = tmp[wh]

;------------
; Set the input values to compile_source_mask_list_from_tcfiles
;------------
radec0 = [fstruct.ra0, fstruct.dec0]
pixel_mask = 0.
ra_min  = radec0[0] - (fstruct.dx / 2.)
ra_max  = radec0[0] + (fstruct.dx / 2.)
dec_min = radec0[1] - (fstruct.dy / 2.)
dec_max = radec0[1] + (fstruct.dy / 2.)
proj = 5
mytime = '20120301';date_toolkit(systime(/julian), 'file')


;------------
; Mask 1: mask all sources above 50mJy
maskrad_arcmin_50mJy  = 5
maskthresh_50mJy      = 50
thresh4maskrad_50mJy  = 50
configfile_out_50mJy  = '/home/kstory/lps12/ptsrc_lists/ptsrc_config_50mJy_'+field_name+'_'+mytime+'.txt'

;------------
; Mask 2: mask all sources with flux in [SNR=5, 50mJy] with 2 arcmin mask,
;         and all sources above 50mJy with 5 arcmin mask

; Get the mJy cut for SNR = 5:
restore, '/data/sptdat/point_source_lists/current/'+field_dir_name+'_list_0p5_loosepm.sav'
mJy_at_5snr = median( mjy_src_150 / sn_src_150 ) * 5.
print, 'mJy_at_5snr = ', mJy_at_5snr

maskrad_arcmin_5snr  = [5, 2]
maskthresh_5snr      = [50, mJy_at_5snr]
thresh4maskrad_5snr  = [50, mJy_at_5snr]
configfile_out_5snr  = '/home/kstory/lps12/ptsrc_lists/ptsrc_config_5snr_'+field_name+'_'+mytime+'.txt'


;;; print out arguments:
print, '******** field ', field_name, ', ', field_idx
; print, '               ', fstruct.name
; print, 'radec0 = ', radec0
; print, 'pixel_mask = ', pixel_mask
; print, 'maskthresh_50mJy = ', maskthresh_50mJy
; print, 'ra_min = ', ra_min
; print, 'ra_max = ', ra_max
; print, 'dec_min = ', dec_min
; print, 'dec_max = ', dec_max
; print, 'maskrad_arcmin = ', maskrad_arcmin
; print, 'proj = ', proj
print, ' - configfile_out_50mJy = ',configfile_out_50mJy
print, ' - configfile_out_5snr  = ',configfile_out_5snr
print, ''
;stop

;----------------------
; PS cut at 50 mJy
print, '--- Make config for 50 mJy cut ---'

compile_source_mask_list_from_tcfiles, radec0, pixel_mask, maskthresh_50mJy, output_struct, $
  ra_min=ra_min, ra_max=ra_max, $
  dec_min=dec_min, dec_max=dec_max, $
  maskrad_arcmin=maskrad_arcmin_50mJy, $
  bands4thresh=150, $
  thresh4maskrad=thresh4maskrad_50mJy, $
  proj=proj, $
  /flux_select, $
  configfile_out=configfile_out_50mJy

;----------------------
; PS cut at 5 SNR
print, '--- Make config for split cut, [2-50, 50+] ---'
compile_source_mask_list_from_tcfiles, radec0, pixel_mask, maskthresh_5snr, output_struct, $
  ra_min=ra_min, ra_max=ra_max, $
  dec_min=dec_min, dec_max=dec_max, $
  maskrad_arcmin=maskrad_arcmin_5snr, $
  bands4thresh=[150, 150], $
  thresh4maskrad=thresh4maskrad_5snr, $
  proj=proj, $
  /flux_select, $
  configfile_out=configfile_out_5snr


;----------------------
; Add RK sources to the fields that need it
add_rk_ptsrc_to_config, field_idx

end


;--------------------------------
; Add sources to specific fields
;   - extra point sources from
;     ~rkeisler/ps09/ptsrc/make_ptsrc_configs_from_flux_cut.pro
; To be run after compile_source_mask_list_from_ctfiles.pro
;
pro add_rk_ptsrc_to_config, field_idx

; Input files to be fixed
in_dir = '/home/kstory/lps12/ptsrc_lists/'
need_fix = 1

; get the field properties
field_arr = lps12_fieldstruct()

field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name


;-------------------
; Get point sources specific to individual fields
;-------------------

case field_dir_name of

    ; if fields is ra23h30dec-55, make sure these two really bright 
    ; guys on the bottom edge are masked:
    'ra23h30dec-55': begin ; index 1
        index_50  = [15]
        index_5snr = [119]
        ra = [357.352386]
        dec = [-49.540691]
        rad = [5.]
        
        index_50  =  [index_50, 16]
        index_5snr = [index_5snr, 120]
        ra = [ra, 356.928894]
        dec = [dec, -49.774326]
        rad = [rad, 5.]
    end
    
    ; if fields is ra21hdec-60, make sure this really bright guy on the
    ; bottom right edge is masked:
    'ra21hdec-60': begin ; index 3
        index_50  = [17]
        index_5snr = [289]
        ra = [330.223541]
        dec = [-55.334740]
        rad = [5.]
    end

    ; if fields is ra3h30dec-60, make sure this really bright guy on the
    ; bottom edge is masked:
    'ra3h30dec-60': begin ; index 4
        index_50   = [30]
        index_5snr = [349]
        ra = [43.367134]
        dec = [-54.696003]
        rad = [5.]
    end

    ; all other fields need to be checked (as of Feb 20, 2012)
    else: need_fix=0
endcase

; Return now if there are no sources to add
if ~need_fix then begin
    print, "No sources to add from RK-by-hand for field: ", field_name
    return
endif

;-------------------
; Write out the extra sources
;-------------------

; Set file names
spawn, 'ls ' + in_dir + 'ptsrc_config_50mJy_'+field_name+'_*.txt ', list 
fname_50mJy = list[0]
spawn, 'ls ' + in_dir + 'ptsrc_config_5snr_'+field_name+'_*.txt ', list
fname_5snr  = list[0]

; copy file to temparary location
spawn, 'less ' + fname_50mJy + ' >> tmp_50mJy.txt'
spawn, 'less ' + fname_5snr +  ' >> tmp_5snr.txt'

; Over-write ptsrc_config
spawn, 'rm -f ' + fname_50mJy
spawn, 'rm -f ' + fname_5snr

; loop over the two cuts
for icut=0, 1 do begin
    tmp_fname = '/home/kstory/lps12/ptsrc_lists/tmp.txt'
    get_lun, lun1
    openw, lun1, tmp_fname

    ; Loop over sources to add
    print, "*** cut ", icut,", Add ", n_elements(ra), " point sources to config file"
    for ii=0, n_elements(ra)-1 do begin

        case icut of
            0: indstr = string(index_50[ii],form='(i4)')
            1: indstr = string(index_5snr[ii],form='(i4)')
        endcase

        rastr = string(ra[ii],form='(f11.5)')
        decstr = string(dec[ii],form='(f11.5)')
        radstr = string(rad[ii]/60.,form='(f8.4)')
        printf,lun1,indstr+'  '+rastr+'   '+decstr+'     '+radstr
    endfor
    close,lun1

    ; write out new file
    case icut of
        0: begin
            spawn, 'less tmp_50mJy.txt ' + tmp_fname + ' >> ' + fname_50mJy
            spawn, 'rm -f tmp_50mJy.txt'
        end
        1: begin
            spawn, 'less tmp_5snr.txt ' + tmp_fname + ' >> ' + fname_5snr
            spawn, 'rm -f tmp_5snr.txt'
        end
    endcase

    ; clean up
    spawn, 'rm -f ' + tmp_fname
endfor


end

