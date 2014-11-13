;;;
; NOTES:
;  1) Play with new lrange A_L chains
;;;

;--------------------
; degeneracy between r and tau
;--------------------
PRO ss
m15_l = [0.131316063367, 0.140709878613, 0,0] & m15_l[2]=mean([m15_l[0],m15_l[1]]) & m15_l[3]=(m15_l[1]-m15_l[0])/2.
m30_l = [0.126373494498, 0.134652222239, 0,0] & m30_l[2]=mean([m30_l[0],m30_l[1]]) & m30_l[3]=(m30_l[1]-m30_l[0])/2.
m15_a = [0.13013681663, 0.13930806664, 0.1346169147, 0]   & m15_a[3] = (m15_a[1]-m15_a[0])/2.
m30_a = [0.12888981261, 0.13848147852, 0.133760267792, 0] & m30_a[3] = (m30_a[1]-m30_a[0])/2.

a15_a = [1.16379330998, 1.8274715411, 1.4855121717, 0]      & a15_a[3] = (a15_a[1]-a15_a[0])/2.
a30_a = [0.687654369648, 0.985415467341, 0.830802459096, 0] & a30_a[3] = (a30_a[1]-a30_a[0])/2.


print, "lowell"
xx = m15_l & print, "omh2 (lcdm ) = ", xx[2]," +- ", xx[3]
xx = m15_a & print, "omh2 (afree) = ", xx[2]," +- ", xx[3]
print, "% Change: ", m15_a[3] / m15_l[3]

print, "highell"
xx = m30_l & print, "omh2 (lcdm ) = ", xx[2]," +- ", xx[3]
xx = m30_a & print, "omh2 (afree) = ", xx[2]," +- ", xx[3]
print, "% Change: ", m30_a[3] / m30_l[3]

print, ""
print,"Alens"
xx = a15_a & print, "Alens (lowell ) = ", xx[2]," +- ", xx[3]
xx = a30_a & print, "Alens (highell) = ", xx[2]," +- ", xx[3]
sig_a = (a15_a[2] - a30_a[2]) / sqrt( a15_a[3]^2 + a30_a[3]^2)

stop
END
