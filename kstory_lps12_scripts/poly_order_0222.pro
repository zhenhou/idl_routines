;;;
; NAME: poly_order_0222
; PURPOSE:
;   Calculate the poly order for all fields
;
; INPUTS:
;
; NOTES:
; 1) 
;
; MODIFICATION HISTORY:
;  02/22/2012: (KTS) Created
;  02/23/2012: (KTS) Change order of print satements
;;;


;------------------------
; test noise psds
pro poly_order
for ii=0, 19 do begin
    poly_order_single_field, ii
endfor
end 

;------------------------
; test noise psds
pro poly_order_single_field, field_idx

; Get the field information
field_arr = lps12_fieldstruct()

field_name = field_arr[field_idx].name
field_dir_name = field_arr[field_idx].dir_name
is_lt = field_arr[field_idx].lead_trail

tmp = read_spt_fields()
wh = where( strcmp(tmp.name, field_dir_name) eq 1, nwh)
if(nwh eq 0) then begin
    print, "field name miss-match"
    return
endif

; Get the field width
dx = tmp[wh].dx
dec0 = tmp[wh].dec0

print, 9 * ( (45*cos(60*!dtor)) / (15*cos(55*!dtor)) )

ltfac = 1.
if is_lt then ltfac = 2.

; Calculate the poly order
print, "-------- Poly order for field[",field_idx,"]: ", field_name
print, '  dx = ', dx, ', dec0 = ', dec0, ', ltfac = ', ltfac

; ra5h30dec-55
poly1 = 15 * ( ( (dx/ltfac)*cos(dec0*!dtor)) / (15*cos(55*!dtor)) )
    
; ra23h30dec-55
poly2 = 9 * ( ( (dx/ltfac)*cos(dec0*!dtor)) / ( (15/2.)*cos(55*!dtor)) )

; ra21hdec-60
poly3 = 12 * ( ( (dx/ltfac)*cos(dec0*!dtor)) / ( (30/2.)*cos(60*!dtor)) )

; ra3h30dec-60
poly4 = 18 * ( ( (dx/ltfac)*cos(dec0*!dtor)) / ( (45/2.)*cos(60*!dtor)) )

; ra21hdec-50
poly5 = 15 * ( ( (dx/ltfac)*cos(dec0*!dtor)) / ( (30/2.)*cos(50*!dtor)) )

; print extra info for each field, if you want to
if 0 then begin
    print, "ra5h30dec-55: ", poly1
    print, "ra23h30dec-55:", poly2
    print, "ra21hdec-60:", poly3
    print, "ra3h30dec-60:", poly4
    print, "ra21hdec-50:", poly5
endif

; The actual poly order we plan to use: average of  2009 fields divided by 3
print, ' *** Use poly order = ', mean( [poly3, poly4, poly5] ) / 3.

end
