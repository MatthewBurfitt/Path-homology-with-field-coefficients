# Path-homology-with-field-coefficients

# Table of Contents

1. [Overview](#Overview)
2. [Methods](#Methods)
3. [Installation](#Installation)  
4. [Usage](#Usage)
      1. [Examples](#Examples)
      3. [Reference manual](#Reference-manual)

# Overview

Computes (though not especially efficiently) the path homology, path homology boundary matrices and bases of the path chain complex of a digraph with respect to coefficients that are rational or a finite field.

This code accompanies the paper: Matthew Burfitt and Tyrone Cutler, [Inductive construction of path homology chains](???).
The original purpose of the code was to independently check the examples Section 6 of the paper, which we include as example scrips.
The algorithm is independent of the method developed in the paper above and based on the procedure outlaid in the paper: Alexander Grigor’yan, [Advances in path homology theory of digraphs](https://intlpress.com/site/pub/files/_fulltext/journals/iccm/2022/0010/0002/ICCM-2022-0010-0002-a007.pdf) Section 1.7.

# Methods

Let $`K`$ be the field of rational numbers $`\mathbb{Q}`$ or a finite field $`\mathbb{Z}/p\mathbb{Z}`$ for some prime number $`p`$.
The path homology of a digraph $`G`$ is computed in the fllowing steps.

1. Allowed paths are vertices of $`G`$ in dimension $`0`$, edges of $`G`$ in dimension $`1`$ and obtained in dimension $`n \geq 2`$ by concatenating any allowed paths in dimensions $`n-1`$ and an edge whose source vertex is the same as the final vertex in the $`(n-1)`$-path.

2. A basis for the path chains $`C_*(G;K)`$ is obtained as the null space of the magnitude homology differential by identifying allowed paths with the basis of the diagonal magnitude chains.

3. The boundary matrix of the path homology differential  $`M(\partial_n)`$ with respect to the computed basis of the path chains and is obtained directly form the basis of the path chains.

4. The rank of the homology in dimension $`n \geq 0`$ is given by the standard formula
```math
\text{rank}\left(H_n(G;K)\right) = \text{dim}\left(C_n(G;K)\right) - \text{rank}\left(M(\partial_n)\right) - \text{rank}\left(M(\partial_{n-1})\right)
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;where $`\text{rank}\left(M(\partial_{-1})\right) = 0`$.

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
H_n(G;\mathbb{Q}) =
\begin{cases} 
        \mathbb{Q}
        &
        \text{if} \; n=0 \: \text{or} \: n=2
        \\
        0
        &
        \text{otherwise}
\end{cases}
```
for each $`n \in \mathbb{N}`$.

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
