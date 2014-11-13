;;;
; NAME: get_lps12_jack_defs.pro
; PURPOSE:
;   Return jackknife definitions
;
; INPUTS
;   field_idx,            index from lps12_fieldstruct.pro
;   stub,                 name of the jackknife
;   files150,             list of files, i.e. files150 = get_lps12_runlist(field_idx, /xspec_maps)
;
; NOTES:
; 1) reference: ~cr/code/spt/jackknives/spt_lowell_jack.pro
;
; MODIFICATION HISTORY:
;  03/07/2012: (KTS) Created from ~cr/code/spt/jackknives/spt_get_lowell_jack_defs.pro
;  03/20/2012: (KTS) Change to using field_idx from lps12_fieldstruct,
;                    fix azrms jack
;  05/16/2012: (KTS) Add azrms_70 and azrms_90 jacks
;  08/05/2012: (KTS) Fix LT naming convention to 'lt_map_...'
;;;

;...................................................................
; Return jackknife definitions
function get_lps12_jack_defs,field_idx,stub,files150, sim=sim

;------------------------
; Setup
;------------------------
field_arr_ = lps12_fieldstruct()
fst = field_arr_[field_idx]
field_name = fst.name
field_dir_name = fst.dir_name

bdir = keyword_set(sim) ? '/data23/hou/lps12/sim_mapmaking/output/' : fst.xspec_map_dir
bfreq = '150'

; Error checking: make sure we properly identified sims or not
if ( (strsplit(files150[0], '/', /extract))[0] ne (strsplit(bdir, '/', /extract))[0] ) then begin
    print, 'GET_LPS12_JACK_DEFS: Input files do not match the bdir.  Did you [not] use the sim flag?'
    print, 'bdir = ', bdir
    print, 'files150[0] = ', files150[0]
    print, 'returning -1'
    RETURN, -1
endif

nfiles=n_elements(files150)
jackflag=1
dataname='MAP.MAP'
postprocess=0

case stub of 
    'lr': begin
        jackflag=0
        i1=lindgen(nfiles)
        i2=0
        dataname='DMAP.MAP'
        setdef=[[i1]]
    end
; OBSOLETE, sim is taken care of by passing in the right files in files150
;     'lr_sim': begin
;         jackflag=0
;         i1=lindgen(nfiles)
;         i2=0
;         dataname='DMAP.MAP'
;         setdef=[[i1]]
;     end
    '12': begin
        en=nfiles/2
        i1=lindgen(en)
        i2=lindgen(en)+en
        setdef=[[i1],[i2]]
    end
    'azrms': begin
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_highrms.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'azrms_50': begin ; keep the 50% of maps with lowest azrms
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms_50.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_highrms_50.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'azrms_70': begin ; keep the 70% of maps with lowest azrms
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms_70.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_highrms_70.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'azrms_90': begin ; keep the 90% of maps with lowest azrms
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms_90.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_highrms_90.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'azrms_95': begin ; keep the 95% of maps with lowest azrms
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_lowrms_95.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_az_highrms_95.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'randcut_95_seed1_azrms': begin ; azrms jack on list with randcut_95_seed1
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_randcut_95_seed1_az_lowrms.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_randcut_95_seed1_az_highrms.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'randcut_95_seed2_azrms': begin ; azrms jack on list with randcut_95_seed2
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_randcut_95_seed2_az_lowrms.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_randcut_95_seed2_az_highrms.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'randcut_95_seed3_azrms': begin ; azrms jack on list with randcut_95_seed3
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_randcut_95_seed3_az_lowrms.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_randcut_95_seed3_az_highrms.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'tweight': begin
        jfile  = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_tweight_low.txt'
        jfile2 = '/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_tweight_high.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'

        postprocess=1
        
    end
    'moon': begin
        jfile='/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_moon_up.txt'
        jfile2='/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_moon_down.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'
        
        postprocess=1
        
    end
    'sun': begin
        jfile='/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_sun_up.txt'
        jfile2='/data/kstory/projects/lps12/jacks/jackgoodfiles/goodfiles_'+field_name+'_sun_down.txt'
        readcol,jfile,g1,format='a'
        readcol,jfile2,g2,format='a'
        
        postprocess=1

        ; check that this field has a reasonable amount of sun-up data:
        restore, '/data/kstory/projects/lps12/masks/masks_50mJy/sun_weights.sav'
        if (p_used[field_idx] lt 0.25) then begin
            print, 'GET_LPS12_JACK_DEFS: this field does not have many sun-up obs, p_used = ', p_used[field_idx]
        endif
        
    end
    'azdsl': begin
        stop                    ; this needs work

    end
    else: stop
endcase

if postprocess eq 1 then begin
    ij = fltarr(nfiles)+1
    ng = n_elements(g1)
    
    ; special case of lt fields
    name_stub=''
    if fst.lead_trail then name_stub = 'lt_'

    fbase = bdir+name_stub+'map_'+field_dir_name+'_'+ bfreq+'_'
    fend  = '.fits'

    if keyword_set(sim) then begin
        fbase = bdir+field_name+'/map_150_lmax8000_'
        fend  = '.bin'
    endif

    ; loop over set 1
    for i=0,ng-1 do begin
        fg=fbase+g1[i]+fend
        ind = where(fg eq files150,nfound)
        ij[ind[0]] = 0
    endfor
    i1 = where(ij eq 0,ni1)

    ; loop over set 2
    ij = fltarr(nfiles)+1
    ng = n_elements(g2)
    for i=0,ng-1 do begin
        ;fg=bdir+name_stub+'map_'+field_dir_name+'_'+bfreq+'_'+g2[i]+'.fits'
        fg=fbase+g2[i]+fend
        ind = where(fg eq files150,nfound)
        ij[ind[0]] = 0
    endfor
    i2 = where(ij eq 0,ni2)
    ni = min([ni1,ni2])
    
    setdef=[[i1[0:ni-1]],[i2[0:ni-1]]]
endif

RETURN,{dataname:dataname,setdef:setdef,jackflag:jackflag}
END
