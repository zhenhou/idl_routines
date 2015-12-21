pro get_lps14_runlist_mdw, ifield, obslist, listfile=listfile

    runlist_path = '/home/zhenhou/scratch-data/projects/sptsz_lowl/runlists/'
    f = lps12_fieldstruct()
    field_name = f[ifield].name

    fobslist = runlist_path+'runlist_lps14_'+field_name+'_0602.txt'
    listfile = fobslist

    obslist = (read_ascii(fobslist)).field1
    nrow=n_elements(obslist[0,*])
    ncol=n_elements(obslist[*,0])
    if nrow eq 1 and ncol gt 1 then begin
        tmp=nrow
        nrow=ncol
        ncol=tmp
    endif
    obslist=strarr(ncol,nrow)
    case ncol of
        1: begin
            readcol,fobslist,a,format='a'
            obslist[0,*]=a
        end
        2: begin
            readcol,fobslist,a,b,format='a,a'
            obslist[0,*]=a
            obslist[1,*]=b
        end
    
        4: begin
            readcol,fobslist,a,b,c,d,format='a,a,a,a'
            obslist[0,*]=a
            obslist[1,*]=b
            obslist[2,*]=c
            obslist[3,*]=d
        end
    else: stop
    endcase

    return
end

