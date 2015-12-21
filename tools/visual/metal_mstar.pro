function myfunc,x,p
  return, p[0]+p[1]*x+p[2]*x*x
end
;----------------
; Sample: sdss+hi
; X-axis: M* from MPA
; Y-axis: metallicity from MPA
; color-coding: f(gas) = M_gas/(M_gas + M_star), here M_gas=1.4*(M_HI+M_H2)
pro look

  d = mrdfits('SDSS_HI_catalog.fits.gz',1)

  ii= where(d.oh gt 0.0)
  d = d[ii]

  class= d.class
  otoh = d.oh
  mour = d.log_m_our[3]
  mmpa = d.log_m_mpa

  mhi  = 10.^(d.log_m_hi)
  t    = d.t_fukugita
  h2tohi=3.7-0.8*t+0.043*t*t             ; M(H2)-to-M(HI) ratio
  mh2  = h2tohi*mhi                      ; M(H2)
  mgas = 1.4d0*mhi*(1+h2tohi)            ; M(gas)

  mstar= mour

  mstar= 10.^mstar
  mtot = alog10(mgas+mstar)              ; total mass
  fgas = alog10(mgas)-mtot               ; f(gas)
  mstar= alog10(mstar)
  mgas = alog10(mgas)
  mhi  = alog10(mhi)
  mh2  = alog10(mh2)

  xmin = 6.0
  xmax =12.0
  ymin = 7.5
  ymax = 9.5

  zmin = -1.0
  zmax =  0.0

  x = mstar
  y = otoh
  z = fgas

  xx=x
  yy=y
  zz=z

  dx= 0.30
  dy= dx * ((ymax-ymin)/(xmax-xmin))
  nx= long((xmax-xmin)/dx+1.e-3)
  ny= long((ymax-ymin)/dy+1.e-3)
  help, nx, ny

  dist=fltarr(nx,ny)

  xbin = findgen(nx+1)*dx + xmin
  ybin = findgen(ny+1)*dy + ymin

  distmin=1.e6
  distmax=-distmin
  for i=0,nx-1 do begin
      xl=xbin[i]
      xu=xbin[i+1]
      for j=0,ny-1 do begin
          yl=ybin[j]
          yu=ybin[j+1]

          ii=where(x ge xl and x lt xu and $
                   y ge yl and y lt yu, nii)
          this = -1000.
          if(nii gt 0)then begin
             this = avg(z[ii])
             if(distmin gt this)then distmin=this
             if(distmax lt this)then distmax=this
          endif
          dist[i,j] = this
      endfor
  endfor

  help, distmin, distmax

  for i=0,nx-1 do dist[i,where(dist[i,*] eq -1000.)]=-1000.

;goto,hehe
  mydevice = !D.name
  myfont   = !P.font
  !P.font = 0
  set_plot,'PS'
  device,filename='metal_mstar.ps', /color,bits_per_pixel=8, $
    set_font='Times Italic', /tt_font, $
    xsize=16.0,ysize=16.0,xoffset=0.0,yoffset=0.0

hehe:
  loadct, 5
  tvim, dist, noaxis=0,                   $
   xrange=[xmin,xmax],yrange=[ymin,ymax], $
   xtitle=textoidl('log(M_\ast/M_{sun})'),$
   ytitle=textoidl('12+log(O/H)'),        $
   stitle=textoidl('log(f_g)'),           $
   range=[zmin,zmax,0.2],                 $
   pcharsize=1.2, /scale, /rct

; xyouts, 6.0, 9.55, alignment=0.0, textoidl('High S/N SF gals from SDSS DR4 + H I sample'), charsize=1.2, charthick=2.0

  j=-1
  for i=0,nx-1 do begin
      xl=xbin[i]
      xu=xbin[i+1]
      ii=where(xx ge xl and xx lt xu, nii)
      if(nii gt 0)then begin
         tax = avg(xx[ii])
         tay = avg(yy[ii])
         j=j+1
         if(j eq 0)then begin
            ax = tax
            ay = tay
         endif else begin
            ax = [ax,tax]
            ay = [ay,tay]
         endelse
      endif
  endfor 

  hogg_usersym, 30
; oplot, ax, ay, psym=8, symsize=2.0, color=djs_icolor('green')

  p0=[-1.492,1.847,-0.08026]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Tremonti et al. (2004)
;
  tx = findgen(61)*0.1+6.0
  ty = p0[0]+p0[1]*tx+p0[2]*tx*tx
; oplot, tx,ty,color=djs_icolor('default'),thick=3.0
  print,p0
; oplot, [6.2,6.8], [9.35,9.35], color=djs_icolor('default'),thick=3.0
; xyouts, 7.0, 9.35, 'Tremonti et al.', color=djs_icolor('default'), charthick=2.0, charsize=1.5

  ee   = xx
  ee[*]=1.0
  p =mpfitfun('myfunc',xx,yy,ee,p0,yfit=yfit,/quiet)
  print,p
  ty = p[0]+p[1]*tx+p[2]*tx*tx
  oplot, tx,ty,color=djs_icolor('green'),thick=3.0
; oplot, [6.2,6.8], [9.20,9.20], color=djs_icolor('green'),thick=3.0
; xyouts, 7.0, 9.20, 'This work', color=djs_icolor('green'), charthick=2.0, charsize=1.5
  xyouts, 11.8,7.60, textoidl('Y='+strmid(strtrim(string(p[0]),2),0,4)+'+'+strmid(strtrim(string(p[1]),2),0,4)+'\cdotX'+strmid(strtrim(string(p[2]),2),0,6)+'\cdotX^2'), $
    align=1.0,color=djs_icolor('green'), charthick=2.0, charsize=1.5

;return
  device,/close
  set_plot,mydevice
  !P.font= myfont

end

