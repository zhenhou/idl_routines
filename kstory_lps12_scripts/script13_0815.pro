;;;
; NAME: script13_0815
;
; NOTES:
;  1) Make runlist_lps12_ra5h30dec-55_2011.txt
;;;

PRO make_azrms_runlist
skeep = '95'
keep = 0.95

f = lps12_fieldstruct()
goodfile_dir ='/home/kstory/lps12/jacks/jackgoodfiles/'
runlist_dir  = '/home/kstory/lps12/runlists/'

idx = 20 ; ra5h30dec-55_2011

    new_runlist = runlist_dir+'runlist_lps12_'+f[idx].name+'.txt'
    f1 = goodfile_dir + 'goodfiles_preaz_'+f[idx].name+'_az_highrms.txt'
    f2 = goodfile_dir + 'goodfiles_preaz_'+f[idx].name+'_az_lowrms.txt'
    
    
    ; read jack files
    print, 'read f1: ' + f1
    readcol,/silent,f1,date_1,format='a'
    
    print, 'read f2: ' + f2
    readcol,/silent,f2,date_2,format='a'
    
    dates = [date_1, date_2]
    dates = dates[sort(dates)]
    nmaps = n_elements(dates)

    nmaps = floor(nmaps * keep)

    ; write new runlist
    print, 'write out new runlist: ' + new_runlist
    get_lun, lun1
    openw, lun1, new_runlist
    for ii=0, nmaps-1 do begin
        printf, lun1, dates[ii]
    endfor
    close, lun1
    free_lun,lun1

stop
END
