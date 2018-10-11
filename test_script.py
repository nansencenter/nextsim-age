from pynextsim.nextsim_bin import NextsimBin
nb = NextsimBin('/output/test00/mesh/field_20061201T000000Z.bin')
nb.plot_var('Concentration')
