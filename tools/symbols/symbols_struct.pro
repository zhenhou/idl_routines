function symbols_struct

    greek = create_struct('phi',    '!9'+String("146B)+'!X', $
                          'kappa',  '!9'+string("153B)+'!X', $
                          'mu',     '!9'+String("155B)+'!X', $
                          'DDelta', '!9'+String("104B)+'!X')

    math  = create_struct('grad',   '!9'+String("321B)+'!X')

    sym = create_struct('greek',greek, 'math',math)

    return, sym
end
