import PathHomology as ph
import matplotlib.pyplot as plt

#define digraph
G = ph.Digraph(edges =
               [('x0','x10'),('x0','x11'),('x0','x12'),
                 ('x10','x20'),('x11','x21'),('x12','x22'),('x10','x21'),('x11','x22'),('x12','x20'),
                 ('x20','x30'),('x21','x31'),('x22','x32'),('x20','x31'),('x21','x32'),('x22','x30'),
                 ('x30','x4'),('x31','x4'),('x32','x4'),
                 ('x10','x30'),('x11','x31'),('x12','x32'),('x10','x32'),('x11','x30'),('x12','x31')])

#show chain rank vector with different coefficients
print('path chain rank vector over Q:   ', G.chain_rank_vector(max_dim = 5, coefficients = 0))
print('path chain rank vector over Z/2Z:', G.chain_rank_vector(max_dim = 5, coefficients = 2))
print('path chain rank vector over Z/3Z:', G.chain_rank_vector(max_dim = 5, coefficients = 3))

#display the digraph
positions = {
            'x0':(0,0),
            'x10':(-1,1),'x11':(0,1),'x12':(1,1),
            'x20':(-1,2),'x21':(0,2),'x22':(1,2),
            'x30':(0,3),'x31':(1,3),'x32':(2,3),
            'x4':(1,4)
            }
G.plot(positions = positions)
plt.show()

