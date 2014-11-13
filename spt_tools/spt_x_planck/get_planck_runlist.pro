;;;
; NAME: get_planck_runlist.pro
; PURPOSE:
;   Get a list of files of projected Planck map on a single SPT field
;
; INPUTS:
;   field_idx,      field index from lps12_fieldstruct()
;   freqs,          Planck HFI frequency(ies) to be loaded 

function get_planck_runlist, field_idx, freqs=freqs, type=type

    f = lps12_fieldstruct()
    fst = f[field_idx]
    field_name = fst.name

    runlist_dir = '/data23/hou/projects/spt_x_planck/reproj/'

    if (not keyword_set(freqs)) then freqs = 143
    num_freqs = n_elements(freqs)

    if (not keyword_set(type)) then type = 'nominal_ringfull'
    
    filenames = 'hfi_SkyMap_'+type+'_'

    for i_freq=0, num_freqs-1 do begin
        tmp = file_search(runlist_dir+field_name,'hfi_SkyMap_'+strcompress(string(freqs[i_freq]),/remove)+'_'+type+'*_maxOrder4_prj.fits')
        if (i_freq eq 0) then files = tmp else files = [files,tmp]
    endfor

    return, files
end
