;;;
; NAME: make_azrms_runlists.pro
; PURPOSE:
;   make azrms_95 runlists
;
; NOTES:
;
; MODIFICATION HISTORY:
;  06/10/2012: (KTS) Created from script_0531.pro
;;;

PRO make_azrms_runlists, skeep
;skeep = '90'

; setup
f = lps12_fieldstruct()
goodfile_dir ='/home/kstory/lps12/jacks/jackgoodfiles/'
runlist_dir  = '/home/kstory/lps12/runlists/'

; loop over fields
for idx=0, 19 do begin
    
    new_runlist = runlist_dir+'runlist_azrms_'+skeep+'_'+f[idx].name+'.txt'
    f1 = goodfile_dir + 'goodfiles_'+f[idx].name+'_az_highrms_'+skeep+'.txt'
    f2 = goodfile_dir + 'goodfiles_'+f[idx].name+'_az_lowrms_'+skeep+'.txt'
    
    
    ; read jack files
    print, 'read f1: ' + f1
    readcol,/silent,f1,date_1,format='a'
    
    print, 'read f2: ' + f2
    readcol,/silent,f2,date_2,format='a'
    
    dates = [date_1, date_2]
    dates = dates[sort(dates)]
    nmaps = n_elements(dates)

    ; write new runlist
    print, 'write out new runlist: ' + new_runlist
    get_lun, lun1
    openw, lun1, new_runlist
    for ii=0, nmaps-1 do begin
        printf, lun1, dates[ii]
    endfor
    close, lun1
    free_lun,lun1
    
endfor


stop
END
