# Path-homology-with-field-coefficients

# Table of Contents

1. [Overview](#Overview)
2. [Installation](#Installation)  
3. [Usage](#Usage)
      1. [Examples](#Examples)
      3. [Reference manual](#Reference-manual)

# Overview

Computes (though not especially efficiently) the path homology, path homology boundary matrices and bases of the path chain complex of a digraph with respect to coefficients that are rational or a finite field.

# Installation

Required python packages: numpy (version 1.26.4), matplotlib (version 3.8.4), networkx (version 3.3) and sympy (version 1.12).

Place a copy of "PathHomology.py" in your python path, nothing else required other than Python 3 (version 3.10.4).

# Usage

Load the *NeighbourhoodBoundaryVolume* library.

```python
import PathHomology as ph
```

## Examples

Example computing each basis of the path chain complex of a multisquare digraph and gain after adding an additional edge across the multisquare.
In the first case a basis in dimensions $`2`$ consists of two elements obtained as the difference of two $`2`$-paths around distinct directed squares, in the second the basis consists of all $`2`$-paths.

```python
import PathHomology as ph
import matplotlib.pyplot as plt

#Define digraph
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
```

Example computing the path homology of a digraph $`G`$ whose edges lie on the one-skeleton of a tetrahedron.
In particulare, we see that the path homology groups of $`G`$ are
```math
H_n(G;\mathbb{Q})
\begin{cases} 
        \mathbb{Q}
        &
        \text{if} \; n=0 \: \text{or} \: n=1
        \\
        0
        &
        \text{otherwise.}
\end{cases}
```

```python
#Define digraph
G = ph.Digraph(edges =
               [('u','v1'),('u','v2'),('u','v3'),('u','v4'),
                ('v1','v3'),('v1','v4'),('v2','v3'),('v2','v4'),
                ('v1','w'),('v2','w'),('v3','w'),('v4','w')])

#obtain path homology of the digraph
homology, _, _, _ = G.path_homology(max_dim = 3, coefficients = 0)
ph.display_homology(homology)

#display the digraph
G.plot()
plt.show()
```


## Reference-manual
