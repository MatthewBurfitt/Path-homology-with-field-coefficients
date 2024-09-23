import PathHomology as ph
import matplotlib.pyplot as plt

#generate digraph edges and a set of vertex positions
def E_t(t = 3):
    edges = [('T','u1A'),('T','u2A'),
             ('u1A','vA'),('u2A','vA'),]
    vertex_positions = {'T':(-3*t+1,0),
                        'u1A':(1,1),'u2A':(2,1),
                        'vA':(1.5,2),
                        'H':(-3*t-0.5,4)}
    for i in range(0,2*t):
        edges.append(('T','u'+str(i)+'C'))
        edges.append(('T','v'+str(i)+'C'))
        #
        edges.append(('u1A','v'+str(i)+'B1'))
        edges.append(('u2A','v'+str(i)+'B2'))
        edges.append(('u'+str(i)+'C','v'+str(i)+'B1'))
        edges.append(('u'+str(i)+'C','v'+str(i)+'B2'))
        edges.append(('u'+str(i)+'C','v'+str(i)+'C'))
        #
        edges.append(('v'+str(i)+'B1','w'+str((i+1)%(2*t))))
        edges.append(('v'+str(i)+'B2','w'+str(i)))
        edges.append(('v'+str(i)+'C','w'+str(i)))
        edges.append(('v'+str(i)+'C','w'+str((i+1)%(2*t))))
        edges.append(('vA','w'+str(i)))
        #
        edges.append(('v'+str(i)+'B1','H'))
        edges.append(('v'+str(i)+'B2','H'))
        edges.append(('w'+str(i),'H'))
        #
        vertex_positions['u'+str(i)+'C'] = ((-3*i-1),1)
        vertex_positions['v'+str(i)+'C'] = ((-3*i-1),2)
        vertex_positions['v'+str(i)+'B1'] = ((-3*i-2),2)
        vertex_positions['v'+str(i)+'B2'] = ((-3*i),2)
        vertex_positions['w'+str(i)] = ((-3*i-1),3)
    return edges, vertex_positions

#define digraph
edge, positions = E_t(t = 3)
G = ph.Digraph(edges = edge)

#show chain rank vector with different coefficients
print('path chain rank vector over Q:   ', G.chain_rank_vector(max_dim = 5, coefficients = 0))
print('path chain rank vector over Z/2Z:', G.chain_rank_vector(max_dim = 5, coefficients = 2))
print('path chain rank vector over Z/3Z:', G.chain_rank_vector(max_dim = 5, coefficients = 3))
print('path chain rank vector over Z/5Z:', G.chain_rank_vector(max_dim = 5, coefficients = 5))

#display the digraph
G.plot(positions = positions)
plt.show()

