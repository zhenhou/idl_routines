pro ps_open,psfile,color=color,landscape=landscape,encap=encap,xsize=xsize,ysize=ysize,inches=inches,portrait=portrait
	set_plot,'ps'
	if keyword_set(color) then begin
                device,filename=psfile+'.ps',/color,bits=8,landscape=landscape,encap=encap,xsize=xsize,ysize=ysize,inches=inches,portrait=portrait	
        endif else device,filename=psfile+'.ps',landscape=landscape,encap=encap,xsize=xsize,ysize=ysize,inches=inches,portrait=portrait,bits=8
end

pro ps_close
	device,/close
	set_plot,'x'

end
