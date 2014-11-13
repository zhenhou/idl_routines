pro array2poly, xarr, y1, y2, xpoly, ypoly

    nx = n_elements(xarr)

    xpoly = dblarr(2*nx)
    xpoly[0:nx-1] = xarr
    for i=nx, 2*nx-1 do xpoly[i] = xarr[2*nx-1-i]

    ypoly = dblarr(2*nx)
    ypoly[0:nx-1] = y1
    for i=nx, 2*nx-1 do ypoly[i] = y2[2*nx-1-i]

    return
    
end
