pro print_ml_chain_all,file,skip=skip,nopsclus=nopsclus
del=1 ; the column of some variables changes depending on whether we have 15 params (with PS cluster amp) or only 14
if keyword_set(nopsclus) then del=0
if n_elements(skip) ne 1 then skip = 2000

a = read_ascii(file)
a=a.field01[*,*]
like = a[1,*]
ind = (where(like eq min(like)))[0]
vector = a[*,ind]
print,vector[2:17]
print,vector[0:1]
print,vector
end
