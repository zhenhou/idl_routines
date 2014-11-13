pro print_ml_chain_all,files,skip=skip,nopsclus=nopsclus
del=1 ; the column of some variables changes depending on whether we have 15 params (with PS cluster amp) or only 14
if keyword_set(nopsclus) then del=0
if n_elements(skip) ne 1 then skip = 2000

a = read_ascii(files[0])
a=a.field01[*,*]

nparam = n_elements(a[*,0])
nf = n_elements(files)
final = fltarr(nparam, 1)

; read all files
for i=1, nf-1 do begin
    print, 'Read file ', i
    b = read_ascii(files[i])
    b = b.field01[*,*]
    final = [ [final], [b] ]
endfor
final = final[*,1:*]

like = final[1,*]
ind = (where(like eq min(like)))[0]
vector = final[*,ind]
print,vector[2:17]
print,vector[0:1]
print,vector
end
