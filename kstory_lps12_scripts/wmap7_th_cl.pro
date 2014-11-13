function wmap7_th_cl, l_in, nu, nolensing=nolensing, dl=dl

if keyword_set(nolensing) then $
  readcol,'camb_89736026_scalcls.dat',l,v1,v2,v3,v4,v5 $
  else $
  readcol,'camb_89736026_lensedcls.dat',l,v1,v2,v3,v4

dl_uk2 = interpol(v1,l,l_in) ; = DL = CL*L*(L+1)/(2PI)
cl_uk2 = dl_uk2/l_in/(l_in+1.)*2.*!pi

dl=dl_uk2

return,cl_uk2
end



