;;;
; NAME: make_runlist_dateOnly_lps12
; PURPOSE:
;   Convert a runlist made with make_runlist_for_lps12 into a runlist
;   with dates only
;
; INPUTS:
;   field_idx,       field index from lps12_fieldstruct()
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/23/2012: (KTS) Created
;  03/05/2012: (KTS) use extract_date_from_filename.pro
;;;

pro run_all
for ii=0, 19 do begin
    make_runlist_dateOnly_lps12, ii
endfor
end

;...................................................................
;
pro make_runlist_dateOnly_lps12, field_idx

field_arr = lps12_fieldstruct()
field_name = field_arr[field_idx].name

runlist_dir = '/data/kstory/projects/lps12/runlists/'
old_runlist = runlist_dir + 'runlist_lps12_'+field_name+'.txt'
new_runlist = runlist_dir + 'runlist_dateOnly_lps12_'+field_name+'.txt'

; Get the date list
field_arr = lps12_fieldstruct()
fst = field_arr[field_idx]

field_name = fst.name
;runlist    = fst.runlist

;maps = get_lps12_runlist(field_idx, /idfs)
readcol,/silent,old_runlist,c1,maps,format='a,a'
nmaps = n_elements(maps)

; Write the new runlist
print, "make new runlis: "
print, "   " + new_runlist
get_lun, lun1
openw, lun1, new_runlist
for ii=0, nmaps-1 do begin
    date = extract_date_from_filename(maps[ii])
    printf, lun1, date
endfor
close, lun1
free_lun,lun1
end
