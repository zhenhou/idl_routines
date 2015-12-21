pro calc_powspec, map_field, mask_field, il, cls

    ellmax = 4000.
    delta_ell = 5.
    nbins = floor(ellmax/delta_ell)
    ell_nom = findgen(nbins)*delta_ell + delta_ell/2.
    reso_arcmin = 1.
    ngrid_big = 4320
    ellg_big = make_ellgrid(ngrid_big,!dtor/60.)

    map_big = dblarr(ngrid_big,ngrid_big)
    map_big[0,0] = map_field

    mask_padded = dblarr(ngrid_big,ngrid_big)
    mask_padded[0,0] = mask_field
    maskfac = mean(mask_padded^2, /double)

    cltmp = double((xspec_flatsky(map_big*mask_padded, map_big*mask_padded, reso_arcmin, delta_ell=delta_ell)/maskfac))

    il = ell_nom
    cls = cltmp

    return

end
