;;;
; NAME: script_0210
; PURPOSE:
;   Run near_one_hz on 2011 fields
;
; NOTES:
; 1) Now using the lpf-ds IDF's from Ryan, Jan3 2012
;
; MODIFICATION HISTORY:
;  02/10/2012: (KTS) Created
;;;

;;;
pro run_11to15
ind = [11, 12, 13, 14, 15]
for ii=0, n_elements(ind)-1 do begin
    near_one_hz, ind[ii], nobs=100, /dosave
endfor
end

;;;
pro run_16to19
ind = [16, 17, 18, 19]
for ii=0, n_elements(ind)-1 do begin
    near_one_hz, ind[ii], nobs=100, /dosave
endfor
end



;--------------------
; Check if flagged bolometers live on the same squids
pro check_squids

; 2008
bb08 = [22,  58,  66,  99, 116, 127, 129, 137, 326, 341, 355, 356, 362, 370, 426, 472, 502, 529, 543, 639]
rk08 = [22,  58,  66,  93,  99, 116, 127, 326, 341, 344, 355, 356, 362, 370, 412, 426, 435, 472, 502, 543, 639]

ac08 = read_array_config(starttime='10-Mar-2008:06:24:22', endtime='10-Mar-2008:16:54:26')

sq08 = ac08.squid_addr[bb08]
rk_sq08 = ac08.squid_addr[rk08]

uniq_sq08 = uniqval( union(sq08, rk_sq08) )

print, '=== 2008 ==='
tab = STRING(9B)
print, 'squid'+tab+tab+' kst'+tab+'     rk'
for ii=0, n_elements(uniq_sq08)-1 do begin
    wh = where(sq08 eq uniq_sq08[ii], nwh)
    wh = where(rk_sq08 eq uniq_sq08[ii], rk_nwh)
    print, uniq_sq08[ii], nwh, rk_nwh
endfor

; 2009
bb09 = [58, 158, 166, 171, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 220, 235, 268, 362, 389, 407, 465, 644, 768]
rk09 = [128, 166, 171, 176, 181, 187, 192, 193, 194, 195, 196, 203, 204, 207, 208, 210, 219, 223, 235, 255, 268, 362, 381, 389, 397, 407, 465, 644, 768]

ac09 = read_array_config(starttime='25-Mar-2009:11:18:49', endtime='25-Mar-2009:12:21:38')

sq09 = ac09.squid_addr[bb09]
rk_sq09 = ac09.squid_addr[rk09]

uniq_sq09 = uniqval( union(sq09, rk_sq09) )

print, '=== 2009 ==='
tab = STRING(9B)
print, 'squid'+tab+tab+' kts'+tab+'     rk'
for ii=0, n_elements(uniq_sq09)-1 do begin
    wh = where(sq09 eq uniq_sq09[ii], nwh)
    wh = where(rk_sq09 eq uniq_sq09[ii], rk_nwh)
    print, uniq_sq09[ii], nwh, rk_nwh
endfor
stop
end

