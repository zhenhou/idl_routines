;;;
; NAME: make_runlist_0612
; PURPOSE:
;   Make new runlists which incorporate the azrms_95 cut.
;
; NOTES:
; 1) make_runlist_rawmap:  copies runlist_lps12_field.txt to runlist_rawmap_field.txt
; 2) make_runlist_raw:     copies runlist_dateOnly_lps12_field.txt to runlist_raw_field.txt
; 3) make_runlist_preaz:   makes runlist with matched lt pairs, no azrms_95 cut
; 4) make_runlist_lps12:   makes runlists with azrms_95 cut
;
; 5) Some test functions to make sure it worked.
;
; MODIFICATION HISTORY:
;  06/12/2012: (KTS) Created
;;;

;-----------------------------
; Make runlist_rawmap_field.txt
PRO make_runlist_rawmap
rdir_old = '/home/kstory/lps12/runlists_0612/'
rdir_new = '/home/kstory/lps12/runlists/'
f = lps12_fieldstruct()

for i=0, 20 do begin
    command = 'cp ' + rdir_old+'runlist_lps12_'+f[i].name+'.txt  ' + $
      rdir_new+'runlist_rawmap_'+f[i].name+'.txt'
    print, command
    spawn, command
endfor
END


;-----------------------------
; Make runlist_raw_field.txt
PRO make_runlist_raw
rdir_old = '/home/kstory/lps12/runlists_0612/'
rdir_new = '/home/kstory/lps12/runlists/'
f = lps12_fieldstruct()

for i=0, 20 do begin
    command =  'cp ' + rdir_old+'runlist_dateOnly_lps12_'+f[i].name+'.txt  ' + $
      rdir_new+'runlist_raw_'+f[i].name+'.txt'
    print, command
    spawn, command
endfor
END



;-----------------------------
; Make runlist_preaz_field.txt
PRO make_runlist_preaz
rdir_old = '/home/kstory/lps12/runlists_0612/'
rdir_new = '/home/kstory/lps12/runlists/'
f = lps12_fieldstruct()

for i=0, 20 do begin
    ;-----------------------------
    ; lead-trail
    ;-----------------------------
    if f[i].lead_trail then begin

        ; get ncol, nrow
        fobslist = rdir_old+'lt/runlist_lt_'+f[i].name+'.txt'
        obslist = (read_ascii(fobslist)).field1
        nrow=n_elements(obslist[0,*])
        ncol=n_elements(obslist[*,0])
        if nrow eq 1 and ncol gt 1 then begin
            tmp=nrow
            nrow=ncol
            ncol=tmp
        endif
        obslist=strarr(ncol,nrow)

        ; read obs list
        case ncol of
            4: begin
                readcol,fobslist,a,b,c,d,format='a,a,a,a'
                obslist[0,*]=a
                obslist[1,*]=b
                obslist[2,*]=c
                obslist[3,*]=d
            end
             2: begin
                readcol,fobslist,a,b,format='a,a'
                obslist[0,*]=a
                obslist[1,*]=b
            end
        endcase

        ; sort
        sind = sort(obslist[0,*])
        sobslist = strarr(ncol,nrow)
        for jcol=0, ncol-1 do begin
            sobslist[jcol,*] = obslist[jcol,sind]
        endfor


        ; write new runlist
        new_runlist = rdir_new+'runlist_preaz_'+f[i].name+'.txt'
        get_lun, lun1
        openw, lun1, new_runlist
        for irow=0, nrow-1 do begin
            case ncol of 
                4: printf, lun1, sobslist[0,irow] +' '+ sobslist[1,irow] +' '+ sobslist[2,irow] +' '+ sobslist[3,irow]
                2: printf, lun1, sobslist[0,irow] +' '+ sobslist[1,irow]
            endcase
        endfor
        close, lun1
        free_lun,lun1

    ;-----------------------------
    ; not lead trail
    ;-----------------------------
    endif else begin
        new_runlist = rdir_new+'runlist_preaz_'+f[i].name+'.txt'
        ;print, 'new runlist: ' + new_runlist
        command =  'cp ' + rdir_old+'runlist_xspec_lps12_'+f[i].name+'.txt  ' + new_runlist
        print, command
        spawn, command
    endelse
endfor
stop
END



;-----------------------------
; Make runlist_preaz_field.txt
PRO make_runlist_lps12

; setup
f = lps12_fieldstruct()
goodfile_dir ='/home/kstory/lps12/jacks/jackgoodfiles/'
rdir_old  = '/home/kstory/lps12/runlists_0612/'
rdir_new  = '/home/kstory/lps12/runlists/'

; loop over fields
for idx=0, 19 do begin
    new_runlist = rdir_new+'runlist_lps12_'+f[idx].name+'.txt'


    ;-----------------------------
    ; Deal with normal maps first
    ;-----------------------------
    if (not f[idx].lead_trail) then begin
        command =  'cp ' + rdir_old+'runlist_azrms_95_'+f[idx].name+'.txt  ' + new_runlist
        print, command
        spawn, command
    endif else begin


    ;-----------------------------
    ; Lead Trail maps, oh what a pain
    ;-----------------------------

        f1 = goodfile_dir + 'goodfiles_'+f[idx].name+'_az_highrms_95.txt'
        f2 = goodfile_dir + 'goodfiles_'+f[idx].name+'_az_lowrms_95.txt'
        
    
        ;----------------
        ; read jack files to get the azrms_95 cut
        ;----------------
        print, 'read f1: ' + f1
        readcol,/silent,f1,date_1,format='a'
        
        print, 'read f2: ' + f2
        readcol,/silent,f2,date_2,format='a'
        
        dates = [date_1, date_2]
        dates = dates[sort(dates)]
        nmaps = n_elements(dates)

        ;----------------
        ; read the preaz runlist to get all observations
        ;----------------
        preaz_runlist = rdir_new+'runlist_preaz_'+f[idx].name+'.txt'

        ; get ncol, nrow
        obslist = (read_ascii(preaz_runlist)).field1
        nrow=n_elements(obslist[0,*])
        ncol=n_elements(obslist[*,0])
        if nrow eq 1 and ncol gt 1 then begin
            tmp=nrow
            nrow=ncol
            ncol=tmp
        endif
        obslist=strarr(ncol,nrow)

        ; read obs list
        case ncol of
            4: begin
                readcol,preaz_runlist,a,b,c,d,format='a,a,a,a'
                obslist[0,*]=a
                obslist[1,*]=b
                obslist[2,*]=c
                obslist[3,*]=d
            end
             2: begin
                readcol,preaz_runlist,a,b,format='a,a'
                obslist[0,*]=a
                obslist[1,*]=b
            end
        endcase

    
        ;----------------
        ; Make the azrms_95 cut
        ;----------------
        whnew = intarr(nmaps)
        for ii=0, nmaps-1 do begin
            whnew[ii] = where(dates[ii] eq obslist[0,*])
        endfor

        new_obslist = strarr(ncol, nmaps)
        new_obslist = obslist[*,whnew]

        ;----------------
        ; write new runlist
        ;----------------
        print, 'write out new runlist: ' + new_runlist
        get_lun, lun1
        openw, lun1, new_runlist
        for ii=0, nmaps-1 do begin
            case ncol of 
                4: printf, lun1, new_obslist[0,ii] +' '+ new_obslist[1,ii] +' '+ new_obslist[2,ii] +' '+ new_obslist[3,ii]
                2: printf, lun1, new_obslist[0,ii] +' '+ new_obslist[1,ii]
            endcase
        endfor
        close, lun1
        free_lun,lun1

    endelse        

endfor

stop
END





;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check preaz runlists
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO check_preaz_field1
ncol = 4
nrow = 272

fold = '/home/kstory/lps12/runlists_0612/lt/runlist_lt_ra23h30dec-55_2008.txt'
readcol,fold,a,b,c,d,format='a,a,a,a'

fnew = '/home/kstory/lps12/runlists/runlist_preaz_ra23h30dec-55_2008.txt'
readcol,fnew,a1,b1,c1,d1,format='a,a,a,a'

whidx_a = intarr(nrow)
whfail_a = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_a[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_a[ir] = 1
endfor

whidx_b = intarr(nrow)
whfail_b = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_b[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_b[ir] = 1
endfor

whidx_c = intarr(nrow)
whfail_c = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_c[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_c[ir] = 1
endfor

whidx_d = intarr(nrow)
whfail_d = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_d[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_d[ir] = 1
endfor


stop
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
PRO check_preaz_field3
ncol = 2
nrow = 270

fold = '/home/kstory/lps12/runlists_0612/lt/runlist_lt_ra23h30dec-55_2008.txt'
readcol,fold,a,b,format='a,a'

fnew = '/home/kstory/lps12/runlists/runlist_preaz_ra23h30dec-55_2008.txt'
readcol,fnew,a1,b1,format='a,a'

whidx_a = intarr(nrow)
whfail_a = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_a[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_a[ir] = 1
endfor

whidx_b = intarr(nrow)
whfail_b = intarr(nrow)
for ir=0, nrow-1 do begin
    whidx_b[ir] = where(a[ir] eq a1, nwh)
    if nwh ne 1 then whfail_b[ir] = 1
endfor

stop
END
