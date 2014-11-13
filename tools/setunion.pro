FUNCTION SetUnion, a, b
    if a[0] lt 0 then return, b    ;A union NULL = a
    if b[0] lt 0 then return, a    ;B union NULL = b
    return, where(histogram([a,b], OMIN = omin)) + omin ;Return combined set
end
