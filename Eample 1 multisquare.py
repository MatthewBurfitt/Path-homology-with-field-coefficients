import PathHomology as ph
import matplotlib.pyplot as plt

#define digraph
G = ph.Digraph(edges =
               [('u','v1'),('u','v2'),('u','v3'),
                ('v1','w'),('v2','w'),('v3','w')])

#obtain basis of the path chain complex of the digraph G
null_spaces, paths = G.path_chain_basis(max_dim = 3, coefficients = 0)
ph.display_path_chain_basis(null_spaces, paths, dim = 2)

#display the digraph
G.plot()
plt.show()

#add an additional edge to the digraph between vertices 'u' and 'w'
G.add_edge(('u','w'))

#obtain basis of the path chain complex of the digraph G
null_spaces, paths = G.path_chain_basis(max_dim = 3, coefficients = 0)
ph.display_path_chain_basis(null_spaces, paths, dim = 2)

#display the new digraph
G.plot()
plt.show()
