function convert_pixel_array,struct
nn=n_elements(struct.map.map)
pixels = lonarr(nn)
npair = n_elements(struct.pixels.map)/2
k=0l
for i=0l,npair-1 do for j=struct.pixels.map[2l*i],struct.pixels.map[2l*i+1l] do pixels[k++]=j
if (k ne nn) then stop
return,pixels
end
