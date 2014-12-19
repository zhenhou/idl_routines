pro read_sims_spec, ifield, savfile_path, savfile_root, nsims, ils, dls_sims
    
    f = lps12_fieldstruct()
    field_name = f[ifield].name
    sidx = strcompress(string(nsims),/remove)
    
    savfile = savfile_path+'/'+field_name+'/sim_dls/'+savfile_root+'_0sims'+sidx+'.sav'

    restore, savfile

    ils = bandcenters

    return
end
