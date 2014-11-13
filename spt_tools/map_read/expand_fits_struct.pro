function expand_fits_struct,m

if find_matching_tag(m,'PIXELS') ne 'PIXELS' then return,m

pixels = convert_pixel_array(m)
map = fltarr(m.mapinfo.nsidex,m.mapinfo.nsidey)
weight=map

map(pixels) = m.map.map
weight(pixels) = m.weight.map

if find_matching_tag(m,'DMAP') eq 'DMAP' then begin
    dmap = fltarr(m.mapinfo.nsidex,m.mapinfo.nsidey)
    dmap(pixels) = m.dmap.map
    return,{map:{map:map},weight:{map:weight},mapinfo:m.mapinfo,processing:m.processing,dmap:{map:dmap}}
endif else $
  return,{map:{map:map},weight:{map:weight},mapinfo:m.mapinfo,processing:m.processing}

end
