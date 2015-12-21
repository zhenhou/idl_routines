;;;
; NAME: lps12_fieldstruct
; PURPOSE:
;   Make structures for each field that hold useful information
;
; CALLING SEQUENCE: make_lps_runlists
;
; INPUTS: none (edit script)
;
; OUTPUTS:
;   idl sav file with an array of structs
;
; NOTES:
;   1) access data in the following way:
;      field_arr = lps12_fieldstruct()
;      (*field_arr[0]).name
;
; MODIFICATION HISTORY:
;   
;   09/11/2011: (KTS) Created
;   09/11/2011: (KTS) fix dirs for fields that have 2010 obs in 2011 dirs
;   02/23/2012: (KTS) add ret_read_spt_fields() function
;   03/01/2012: (KTS) add obs_map_dir, xspec_map_dir
;   04/23/2012: (KTS) add year
;   05/03/2012: (KTS) fix ra0, dec0 in lps12_fieldstruct.  In the previous commit, these numbers were being truncaced to an int.
;   05/23/2012: (KTS) Fix bug in f[ii].dy.  Previously, this was being set to fst.dx.  CRAP!
;   06/12/2012: (KTS) Change for new runlists (azrms_95 cut).
;   06/28/2012: (KTS) Change mapdirs to 0620
;   03/04/2013: (KTS) Only run read_spt_fields() once
;;;


;...................................................................
; Return additional information struct from read_spt_fields()
;
FUNCTION ret_read_spt_fields, field_info, field_dir_name 
wh = where( strcmp(field_info.name, field_dir_name) eq 1, nwh)
fstruct = field_info[wh]
return, fstruct
end


;...................................................................
; Main function
;
FUNCTION spt_fieldsinfo
compile_opt IDL2, HIDDEN

;;;;;;;;;;;;;;;;;;;;;
; Make the field_arr
;;;;;;;;;;;;;;;;;;;;;

nfields = 21

; initialize the array of structures
tmp = {name:'xx', $
       dir_name:'xx', $
       clt_name:'xx', $
       lead_trail:-1, $
       obs_map_dir:'xx', $
       xspec_map_dir:'xx', $
       autoprocessed_map_dirs:'xx', $
       idf_dirs:'xx', $         ; full idf's
       idf_lpfds_dirs:'xx', $ ; low-passed, filtered idf's from RK, Jan 2012
       $;obs_runlist:'xx', $
       $;xspec_runlist:'xx', $
       runlist:'xx', $
       $
       $ ; structure members from ret_read_spt_fields.pro
       NX:-1, NY:-1, $
       RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
       RA0_CLT:-1.0d0, DEC0_CLT:-1.0d0, NX_CLT:-1., NY_CLT:-1.}

f = replicate(tmp, nfields)
;field_arr_ = ptrarr(nfields, /allocate_heap)

;; 2008
f[0] = {name:'ra5h30dec-55_2008', dir_name:'ra5h30dec-55', clt_name:'ra5h30dec-55_2year', lead_trail:0, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra5h30dec-55_2008/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra5h30dec-55_2008/', $
        autoprocessed_map_dirs:'/data/sptdat/run2b/ra5h30dec-55/maps/', $
        idf_dirs:'/data/sptdat/run2b/ra5h30dec-55/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra5h30dec-55_2008/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra5h30dec-55_2008.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra5h30dec-55_2008.txt', $
        NX:960, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:82.7026d0, DEC0_CLT:-55.0005d0, NX_CLT:3120., NY_CLT:3120.}

f[1] = {name:'ra23h30dec-55_2008', dir_name:'ra23h30dec-55', clt_name:'ra23h30dec-55_2year', lead_trail:1, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra23h30dec-55_2008/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra23h30dec-55_2008_lt/', $
        autoprocessed_map_dirs:'/data/sptdat/run2b/ra23h30dec-55/maps/', $
        idf_dirs:'/data/sptdat/run2b/ra23h30dec-55/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra23h30dec-55_2008/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra23h30dec-55_2008.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra23h30dec-55_2008.txt', $
        NX:960, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:352.502d0, DEC0_CLT:-55.0004d0, NX_CLT:3360., NY_CLT:3360.}

f[2] = {name:'ra23h30dec-55_2010', dir_name:'ra23h30dec-55', clt_name:'ra23h30dec-55_2010', lead_trail:0, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra23h30dec-55_2010/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra23h30dec-55_2010/', $
        autoprocessed_map_dirs:'/data/sptdat/run2_2010/ra23h30dec-55_2010/maps/', $
        idf_dirs:'/data/sptdat/run1_2010/ra23h30dec-55/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra23h30dec-55_2010/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra23h30dec-55_2010.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra23h30dec-55_2010.txt', $
        NX:960, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:352.502d0, DEC0_CLT:-55.0004d0, NX_CLT:3360., NY_CLT:3360.}

;; 2009
f[3] = {name:'ra21hdec-60', dir_name:'ra21hdec-60', clt_name:'ra21hdec-60', lead_trail:1, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra21hdec-60/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra21hdec-60_lt/', $
        autoprocessed_map_dirs:'/data/sptdat/run1_2009/ra21hdec-60/maps/', $
        idf_dirs:'/data/sptdat/run1_2009/ra21hdec-60/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra21hdec-60/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra21hdec-60.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra21hdec-60.txt', $
        NX:1280, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:315.004d0, DEC0_CLT:-60.0009d0, NX_CLT:4600., NY_CLT:3000.}

f[4] = {name:'ra3h30dec-60', dir_name:'ra3h30dec-60', clt_name:'ra3h30dec-60', lead_trail:1, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra3h30dec-60/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra3h30dec-60_lt/', $
        autoprocessed_map_dirs:'/data/sptdat/run1_2009/ra3h30dec-60/maps/', $
        idf_dirs:'/data/sptdat/run1_2009/ra3h30dec-60/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra3h30dec-60/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra3h30dec-60.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra3h30dec-60.txt', $
        NX:2048, NY:1024, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:52.5044d0, DEC0_CLT:-60.0010d0, NX_CLT:6700., NY_CLT:3000.}

f[5] = {name:'ra21hdec-50', dir_name:'ra21hdec-50', clt_name:'ra21hdec-50', lead_trail:1, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra21hdec-50/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra21hdec-50_lt/', $
        autoprocessed_map_dirs:'/data/sptdat/run1_2009/ra21hdec-50/maps/', $
        idf_dirs:'/data/sptdat/run1_2009/ra21hdec-50/fits/', $ 
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra21hdec-50/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra21hdec-50.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra21hdec-50.txt', $
        NX:1536, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:315.001d0, DEC0_CLT:-50.0012d0, NX_CLT:5800., NY_CLT:3600.}

;; 2010 only
f[6] = {name:'ra4h10dec-50', dir_name:'ra4h10dec-50', clt_name:'ra4h10dec-50', lead_trail:0, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra4h10dec-50/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra4h10dec-50/', $
        autoprocessed_map_dirs:'/data/sptdat/run2_2010/ra4h10dec-50/maps/', $
        idf_dirs:'/data/sptdat/run1_2010/ra4h10dec-50/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra4h10dec-50/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra4h10dec-50.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra4h10dec-50.txt', $
        NX:1536, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:62.5009d0, DEC0_CLT:-50.0018d0, NX_CLT:4700., NY_CLT:3000.}

f[7] = {name:'ra0h50dec-50', dir_name:'ra0h50dec-50', clt_name:'ra0h50dec-50', lead_trail:0, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra0h50dec-50/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra0h50dec-50/', $
        autoprocessed_map_dirs:'/data/sptdat/run2_2010/ra0h50dec-50/maps/', $
        idf_dirs:'/data/sptdat/run1_2010/ra0h50dec-50/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra0h50dec-50/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra0h50dec-50.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra0h50dec-50.txt', $
        NX:1536, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:12.5016d0, DEC0_CLT:-50.0010d0, NX_CLT:4700., NY_CLT:3000.}

f[8] = {name:'ra2h30dec-50', dir_name:'ra2h30dec-50', clt_name:'ra2h30dec-50', lead_trail:0, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra2h30dec-50/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra2h30dec-50/', $
        autoprocessed_map_dirs:'/data/sptdat/run2_2010/ra2h30dec-50/maps/', $
        idf_dirs:'/data/sptdat/run1_2010/ra2h30dec-50/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra2h30dec-50/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra2h30dec-50.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra2h30dec-50.txt', $
        NX:1536, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:37.5015d0, DEC0_CLT:-50.0006d0, NX_CLT:4700., NY_CLT:3000.}

f[9] = {name:'ra1hdec-60', dir_name:'ra1hdec-60', clt_name:'ra1hdec-60', lead_trail:0, $
        obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra1hdec-60/', $
        xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra1hdec-60/', $
        autoprocessed_map_dirs:'/data/sptdat/run2_2010/ra1hdec-60/maps/', $
        idf_dirs:'/data/sptdat/run1_2010/ra1hdec-60/fits/', $
        idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra1hdec-60/', $
        $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra1hdec-60.txt', $
        runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra1hdec-60.txt', $
        NX:1280, NY:960, $
        RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:15.0025d0, DEC0_CLT:-60.0009d0, NX_CLT:4600., NY_CLT:3000.}

f[10] = {name:'ra5h30dec-45', dir_name:'ra5h30dec-45', clt_name:'ra5h30dec-45', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra5h30dec-45/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra5h30dec-45/', $
         autoprocessed_map_dirs:'/data/sptdat/run2_2010/ra5h30dec-45/maps/', $
         idf_dirs:'/data/sptdat/run1_2010/ra5h30dec-45/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra5h30dec-45/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra5h30dec-45.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra5h30dec-45.txt', $
         NX:1024, NY:1024, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:82.5013d0, DEC0_CLT:-44.9999d0, NX_CLT:3200., NY_CLT:3200.}

;; 2010 & 2011
f[11] = {name:'ra6h30dec-55', dir_name:'ra6h30dec-55', clt_name:'ra6h30dec-55', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra6h30dec-55/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra6h30dec-55/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra6h30dec-55/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra6h30dec-55/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra6h30dec-55/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra6h30dec-55.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra6h30dec-55.txt', $
         NX:960, NY:960, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
         RA0_CLT:97.5028d0, DEC0_CLT:-55.0001d0, NX_CLT:2700., NY_CLT:3000.}
         
f[12] = {name:'ra23hdec-62.5', dir_name:'ra23hdec-62.5', clt_name:'ra23hdec-62.5', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra23hdec-62.5/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra23hdec-62.5/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra23hdec-62.5/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra23hdec-62.5/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra23hdec-62.5/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra23hdec-62.5.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra23hdec-62.5.txt', $
         NX:1152, NY:640, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:344.999d0, DEC0_CLT:-62.4984d0, NX_CLT:4100., NY_CLT:1800.}

f[13] = {name:'ra21hdec-42.5', dir_name:'ra21hdec-42.5', clt_name:'ra21hdec-42.5', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra21hdec-42.5/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra21hdec-42.5/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra21hdec-42.5/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra21hdec-42.5/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra21hdec-42.5/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra21hdec-42.5.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra21hdec-42.5.txt', $
         NX:1728, NY:640, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:315.005d0, DEC0_CLT:-42.4995d0, NX_CLT:5900., NY_CLT:1800.}

f[14] = {name:'ra22h30dec-55', dir_name:'ra22h30dec-55', clt_name:'ra22h30dec-55', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra22h30dec-55/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra22h30dec-55/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra22h30dec-55/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra22h30dec-55/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra22h30dec-55/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra22h30dec-55.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra22h30dec-55.txt', $
         NX:960, NY:960, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:337.503d0, DEC0_CLT:-54.9963d0, NX_CLT:2700., NY_CLT:3000.}
         
f[15] = {name:'ra23hdec-45', dir_name:'ra23hdec-45', clt_name:'ra23hdec-45', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra23hdec-45/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra23hdec-45/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra23hdec-45/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra23hdec-45/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra23hdec-45/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra23hdec-45.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra23hdec-45.txt', $
         NX:1728, NY:960, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:345.005d0, DEC0_CLT:-45.0001d0, NX_CLT:5900., NY_CLT:3000.}

f[16] = {name:'ra6hdec-62.5', dir_name:'ra6hdec-62.5', clt_name:'ra6hdec-62.5', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra6hdec-62.5/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra6hdec-62.5/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra6hdec-62.5/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra6hdec-62.5/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra6hdec-62.5/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra6hdec-62.5.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra6hdec-62.5.txt', $
         NX:1152, NY:640, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
         RA0_CLT:89.9994d0, DEC0_CLT:-62.5008d0, NX_CLT:4000., NY_CLT:1800.}

f[17] = {name:'ra3h30dec-42.5', dir_name:'ra3h30dec-42.5', clt_name:'ra3h30dec-42.5', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra3h30dec-42.5/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra3h30dec-42.5/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra3h30dec-42.5/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra3h30dec-42.5/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra3h30dec-42.5/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra3h30dec-42.5.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra3h30dec-42.5.txt', $
         NX:2160, NY:960, $
         RA0:-1.1, DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
         RA0_CLT:52.5055d0, DEC0_CLT:-42.5019d0, NX_CLT:8700., NY_CLT:1800.}

f[18] = {name:'ra1hdec-42.5', dir_name:'ra1hdec-42.5', clt_name:'ra1hdec-42.5', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra1hdec-42.5/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra1hdec-42.5/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra1hdec-42.5/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra1hdec-42.5/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra1hdec-42.5/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra1hdec-42.5.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra1hdec-42.5.txt', $
         NX:1728, NY:640, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:15.0055d0, DEC0_CLT:-42.4998d0, NX_CLT:5900., NY_CLT:1800.}

f[19] = {name:'ra6h30dec-45', dir_name:'ra6h30dec-45', clt_name:'ra6h30dec-45', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra6h30dec-45/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra6h30dec-45/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra6h30dec-45/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra6h30dec-45/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra6h30dec-45/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra6h30dec-45.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra6h30dec-45.txt', $
         NX:1024, NY:1024, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:97.5054d0, DEC0_CLT:-45.0003d0, NX_CLT:3200., NY_CLT:3200.}

f[20] = {name:'ra5h30dec-55_2011', dir_name:'ra5h30dec-55', clt_name:'ra5h30dec-55_2011', lead_trail:0, $
         obs_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra5h30dec-55/', $
         xspec_map_dir:'/data/kstory/projects/lps12/maps/20120620/ra5h30dec-55/', $
         autoprocessed_map_dirs:'/data17/sptdat/run2_2011/ra5h30dec-55/maps/', $
         idf_dirs:'/data/sptdat/run1_2011/ra5h30dec-55/fits/', $
         idf_lpfds_dirs:'/home/rkeisler/lowellfits/ra5h30dec-55_2011/', $
         $;obs_runlist:'/data/kstory/projects/lps12/runlists/runlist_dateOnly_lps12_ra5h30dec-55_2011.txt', $
         runlist:'/data/kstory/projects/lps12/runlists/runlist_lps12_ra5h30dec-55_2011.txt', $
         NX:960, NY:960, $
         RA0:-1., DEC0:-1., DX:-1., DY:-1., AREA:-1., URL:'xx', OBSERVED:-1., YEAR:-1, COMMENT:'xx', $
        RA0_CLT:-1.0d0, DEC0_CLT:-1.0d0, NX_CLT:-1., NY_CLT:-1.}

;;; Add information from read_spt_fields()
field_info = read_spt_fields()
for ii=0, nfields-1 do begin
    field_dir_name = (f[ii]).dir_name
    fst = ret_read_spt_fields(field_info, field_dir_name)
    f[ii].RA0=fst.RA0 & f[ii].DEC0=fst.DEC0 & f[ii].DX=fst.DX & f[ii].DY=fst.DY & f[ii].AREA=fst.AREA & f[ii].URL=fst.URL & f[ii].YEAR=fst.YEAR
endfor

; special case, fix year for ra23h30dec-55_2010, ra5h30dec-55_2011
f[2].year = 2010L
f[20].year = 2011L


RETURN, f
END

