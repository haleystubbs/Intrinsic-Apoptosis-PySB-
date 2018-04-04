
monomer('A', ['state'], {'state': ['A', 'M']})
monomer('B', ['state'], {'state': ['A', 'M']})
monomer('C')

equilibrate(A(state='A'), A(state='M'), [apoptosome_formation, apoptosome_degradation])
equilibrate(B(state='A'), B(state='M'), [apoptosome_formation, apoptosome_degradation])

Parameter('apoptosome_formation', 1e-1)
Parameter('apoptosome_degradation', )