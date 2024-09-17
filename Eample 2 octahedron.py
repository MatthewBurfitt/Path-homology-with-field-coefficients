import PathHomology as ph
import matplotlib.pyplot as plt

#define digraph
G = ph.Digraph(edges =
               [('u','v1'),('u','v2'),('u','v3'),('u','v4'),
                ('v1','v3'),('v1','v4'),('v2','v3'),('v2','v4'),
                ('w','v1'),('w','v2'),('w','v3'),('w','v4')])

#obtain path homology of the digraph
homology, _, _, _ = G.path_homology(max_dim = 3, coefficients = 0)
ph.display_homology(homology)

#display the digraph
G.plot()
plt.show()

