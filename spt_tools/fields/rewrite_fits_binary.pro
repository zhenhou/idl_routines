pro rewrite_fits_binary, fits_file, bin_file

    res = read_spt_fits(fits_file)
    openw, 5, bin_file
    writeu,5, res.map.map
    close, 5

end
