function dbdt_ratio,nu,nu0

x0 = nu0/56.78d
dBdT0 = x0^4 * exp(x0) / (exp(x0)-1)^2

x = nu/56.78d
return, x^4 * exp(x) / (exp(x)-1)^2 / dbdT0

end
