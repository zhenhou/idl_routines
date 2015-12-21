pro inline_print, i, n
    i_str = strcompress(string(i),/remove)
    n_str = strcompress(string(n),/remove)

    print, format='(%"\33[1M processing %s of %s \33[1A")',i_str, n_str
end
