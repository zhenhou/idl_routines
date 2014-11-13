
function scale_ps,c0,f0,f1,alpha
c1 =c0* (f1/f0)^(alpha*2.0)/ dbdt_ratio(f1,f0)^2

return,c1


end
