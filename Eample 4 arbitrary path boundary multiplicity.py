import PathHomology as ph
import matplotlib.pyplot as plt

#generate digraph edges and a set of vertex positions
def M_t(t = 2):
    edges = [('T','u1A'),('T','u2A'),('T','uB'),
             ('u1A','v1A'),('u1A','v2A'),('u2A','v1A'),('u2A','v2A'),
             ('v1A','wA'),('v2A','wA'),
             ('wA','H')]
    vertex_positions = {'T':(0,0),
                        'u1A':(1,1),'u2A':(2,1),'uB':(-t,1),
                        'v1A':(1,2),'v2A':(2,2),
                        'wA':(1.5,3),
                        'H':(0,4)}
    for i in range(0,2*t):
        edges.append(('uB','v'+str(i)+'B'))
        edges.append(('v'+str(i)+'B','w'+str(i)+'B'))
        edges.append(('v'+str(i)+'B','w'+str(((i+1)%(2*t)))+'B'))
        edges.append(('w'+str(i)+'B','H'))
        vertex_positions['v'+str(i)+'B'] = (-i,2)
        vertex_positions['w'+str(i)+'B'] = (-i,3)
    for i in range(0,t):
        edges.append(('u1A','v'+str(((2*i)%(2*t))+1)+'B'))
        edges.append(('u2A','v'+str(((2*i)%(2*t)))+'B'))
        edges.append(('v1A','w'+str(((2*i)%(2*t)))+'B'))
        edges.append(('v2A','w'+str(((2*i)%(2*t))+1)+'B'))
    return edges, vertex_positions

#define digraph
edge, positions = M_t(t = 2)
G = ph.Digraph(edges = edge)

#obtain path homology and path boundary matrices of the digraph
print('\n')
homology, bnd_matrices, null_spaces, paths = G.path_homology(max_dim = 5, coefficients = 0)
ph.display_homology(homology)
print('\n')
print('boundary matrix from dimension 4 to 3:')
print(bnd_matrices[3])

#display the digraph
G.plot(positions = positions)
plt.show()

