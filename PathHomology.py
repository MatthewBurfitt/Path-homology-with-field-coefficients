import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from sympy import Matrix, QQ, GF
from sympy.polys.matrices import DM
from sympy.matrices.normalforms import smith_normal_form


#display path chain basis
def display_path_chain_basis(null_spaces, paths, dim = 0):
    print('\n')
    print('Path chain basis in dimension ' + str(dim))
    for i in range(len(null_spaces[dim])):
        print('Element '+str(i+1)+':')
        first = True
        for j in range(len(paths[dim])):
            if null_spaces[dim][i][j] != 0:
                if null_spaces[dim][i][j] > 0:
                    if first:
                        print(' '+str(null_spaces[dim][i][j]), paths[dim][j])
                        first = False
                    else:
                        print('+'+str(null_spaces[dim][i][j]), paths[dim][j])
                else:
                    print(null_spaces[dim][i][j], paths[dim][j])
                    first = False


#display homology dictionary
def display_homology(homology):
    for key in homology:
        if key[0:10] == 'dimension ':
            if str(homology[key][1]) == '1':
                print('homology in '+ key+':', homology[key][0])
            else:
                print('homology in '+ key+':', homology[key][0]+'^'+homology[key][1])
    if homology['last dimension computed'] == None:
        print('and 0 otherwise')
    else:
        print('last dimension computed:', homology['last dimension computed'])


#class representing a digraph
class Digraph:
    def __init__(self, vertices = [], edges = []):
        self.vertices = []
        self.add_vertex(vertices)
        self.edges = []
        self.add_multiple_edges(edges)

    #retrive basic digraph infomation
    def vertex_number(self):
        return len(self.vertices)
    def edge_number(self):
        return len(self.edges)

    #extend digraph by adding additional vertice(s)
    def add_vertex(self, vertex = []):
        if type(vertex) == list:
            temp = len(self.vertices)
            self.vertices = list(set(self.vertices + vertex))
            if len(self.vertices) < temp + len(vertex):
                print("Some vertices already exist!")
        else:
            if vertex not in self.vertices:
                self.vertices.append(vertex)
            else:
                print("Vertex already exists!")

    #remove vertice(s) and all associated edges from the digraph
    def remove_vertex(self, vertex = []):
        if type(vertex) != list:
            vertex = [vertex]
        for i in range(len(self.edges)-1,-1,-1):
            if self.edges[i][0] in vertex or self.edges[i][1] in vertex:
                del self.edges[i]
        for i in range(len(self.vertices)-1,-1,-1):
            if self.vertices[i] in vertex:
                del self.vertices[i]

    #extend digraph by adding one additional edge
    def add_edge(self, edge):
        if len(edge) < 2:
            print("edge "+str(edge)+" not of the correct form and not added!")
        else:
            if edge[0] not in self.vertices:
                self.add_vertex(edge[0])
            if edge[1] not in self.vertices:
                self.add_vertex(edge[1])
            if edge[0] != edge[1]:
                    self.edges.append([edge[0], edge[1]])
            else:
                print("Edge at " + str(edge[0]) + " on vertex was not added!")

    #extend digraph by adding additional edges
    def add_multiple_edges(self, edges = []):
        temp = len(edges)
        edges = list(set(edges))
        if len(edges) < temp:
            print("multiple edges detected and removed!")
        for edge in edges:
            self.add_edge(edge)

    #remove an edge from the digraph
    def remove_edge(self, edge):
        for i in range(len(self.edges)-1,-1,-1):
            if self.edges[i][0] == edge[0] and self.edges[i][1] == edge[1]:
                del self.edges[i]
                
    #remove one or more edges from the digraph
    def remove_multiple_edges(self, edges = []):
        for edge in edges:
            self.remove_vertex(edge)

    #display digraph
    def plot(self, positions = None):
        fig, ax = plt.subplots()
        G = nx.DiGraph()
        G.add_edges_from(self.edges)
        if type(positions) == dict:
            pos = positions
        else:
            pos = nx.spring_layout(G,seed=5)
        nx.draw_networkx_nodes(G, pos, ax=ax)
        nx.draw_networkx_labels(G, pos, ax=ax)
        curved_edges = [edge for edge in G.edges() if reversed(edge) in G.edges()]
        straight_edges = list(set(G.edges()) - set(curved_edges))
        nx.draw_networkx_edges(G, pos, ax=ax, edgelist=straight_edges)
        arc_rad = 0.25
        nx.draw_networkx_edges(G, pos, ax=ax, edgelist=curved_edges, connectionstyle=f'arc3, rad = {arc_rad}')

    #allowed paths up to a given dimension
    def get_allowed_paths(self, max_dim = 3):
        allowed_paths = [[]]
        for v in self.vertices:
            allowed_paths[0].append([v])
        if max_dim >= 1:
            allowed_paths.append(self.edges)
        dim = 2
        for dim in range(2,max_dim+1):
            new_allowed_paths = []
            for path in allowed_paths[-1]:
                for edge in self.edges:
                    if edge[0] == path[-1]:
                        new_allowed_paths.append(list(path)+[edge[1]])
            allowed_paths.append(new_allowed_paths)
        return allowed_paths

    #magnitude differential matrices up to a given dimension
    def magnitude_differentials(self, max_dim = 3):
        allowed_paths = self.get_allowed_paths(max_dim = max_dim)
        differential_matrices = [np.zeros(len(self.vertices)).reshape(-1,1).astype(int),
                                 np.zeros(len(self.edges)).reshape(-1,1).astype(int)]
        path_images = []
        for dim in range(2, max_dim + 1):
            matrix = []
            all_rows = []
            img = []
            #find all images of allowed paths under the magnitude differntial and record the image indices
            for path in allowed_paths[dim]:
                row = []
                for i in range(1, dim):
                    position_removed_path = path[:i]+path[i+1:]
                    if position_removed_path not in img:
                        img.append(position_removed_path)
                    not_found_edge = True
                    if path[i-1] == path[i+1]:
                        row.append(-1)
                        not_found_edge = False
                    else:
                        for edge in self.edges:
                            if path[i-1] == edge[0] and path[i+1] == edge[1]:
                                row.append(-1)
                                not_found_edge = False
                                break
                    if not_found_edge:
                        row.append(img.index(position_removed_path))
                all_rows.append(row)
            #build matrix form row indices
            for r in range(len(all_rows)):
                matrix_row = [0]*len(img)
                for i in range(len(all_rows[r])):
                    if all_rows[r][i] != -1:
                        matrix_row[all_rows[r][i]] += (-1)**(i+1)
                matrix.append(matrix_row)
            path_images.append(img)
            differential_matrices.append(np.array(matrix))
        return differential_matrices, allowed_paths, path_images

    #path chain basis up to a given dimension with given coefficients
    def path_chain_basis(self, max_dim = 3, coefficients = 0):
        if type(coefficients) == int and coefficients > 0:
            domain = GF(coefficients)
        else:
            domain = QQ
        all_null_spaces = [np.eye(len(self.vertices)).astype(int), np.eye(len(self.edges)).astype(int)]
        all_required_allowed_paths = [self.vertices, self.edges]
        magnitude_differential_matrices, allowed_paths, _ = self.magnitude_differentials(max_dim = max_dim)
        for dim in range(2, max_dim + 1):
            if len(magnitude_differential_matrices[dim]) > 0:
                #compute nulspaces with correct coefficients
                magnitude_differential_matrices[dim] = magnitude_differential_matrices[dim].astype(int)
                null_space = DM(magnitude_differential_matrices[dim].T.tolist(), domain).nullspace()
                null_space = np.array(null_space.to_Matrix())#clear domain first???
                if len(null_space) > 0:
                    null_space.reshape(len(null_space),-1)
                #reduce the collections of allowed paths to just those that appear in the null spaces
                for i in range(null_space.shape[1]-1,-1,-1):
                    if not np.any(null_space[:,i]):
                        null_space = np.delete(null_space, i, 1)
                        allowed_paths[dim] = allowed_paths[dim][:i] + allowed_paths[dim][i+1:]
                all_null_spaces.append(null_space)
            else:
                all_null_spaces.append(np.array([]))
        return all_null_spaces, allowed_paths

    #chain rank vector up to a given dimension with given coefficients
    def chain_rank_vector(self, max_dim = 3, coefficients = 0):
        null_spaces, _ = self.path_chain_basis(max_dim = max_dim, coefficients = coefficients)
        rank_vector = []
        for i in range(max_dim + 1):
            rank_vector.append(len(null_spaces[i]))
        return rank_vector

    #path homology differential matrices up to a given dimension with given coefficients
    def path_differentials(self, max_dim = 3, coefficients = 0):
        if type(coefficients) == int and coefficients > 0:
            domain = GF(coefficients)
        else:
            domain = QQ
        null_spaces, paths = self.path_chain_basis(max_dim = max_dim, coefficients = coefficients)
        differential_matrices = []
        for dim in range(1, max_dim + 1):
            matrix = []
            #find column space, its inverse and indices of the columns to later obtain 
            #linear combinations of chain basis expressions form path basis expressions
            null_spaces[dim-1] = null_spaces[dim-1].astype(int)
            null_space_matrix = DM(null_spaces[dim-1].tolist(), domain)
            null_space_coll_space = null_space_matrix.columnspace()
            null_space_coll_space_index = []
            for i in range(null_space_coll_space.shape[1]):
                for j in range(null_space_matrix.shape[1]):
                    if null_space_coll_space[:,i] == null_space_matrix[:,j]:
                           null_space_coll_space_index.append(j)
                           break
            null_space_coll_space_inv = null_space_coll_space.inv()
            #determine the image of the differential on a single element
            for element in null_spaces[dim]:
                diff_paths = []
                diff_path_coefs = []
                for j in range(len(element)):
                    if element[j] != 0:
                        for i in range(dim + 1):
                            image = paths[dim][j][:i]+paths[dim][j][i+1:]
                            if image in diff_paths:
                                diff_path_coefs[diff_paths.index(image)] += (-1)**i * element[j]
                            else:
                                diff_paths.append(image)
                                diff_path_coefs.append((-1)**i * element[j])
                #correct for coefficients in finite field
                if type(coefficients) == int and coefficients > 0:
                    for i in range(len(diff_path_coefs)):
                        diff_path_coefs[i] = diff_path_coefs[i] % coefficients
                #first compute row corresponding to alowed paths
                path_row = [0]*len(paths[dim - 1])
                for i in range(len(diff_paths)):
                    if diff_path_coefs[i] != 0:
                        path_row[paths[dim-1].index(diff_paths[i])] = diff_path_coefs[i]
                #then get matrix row by determining it at a linear combination of null space vectors
                path_row_reduced = [path_row[i] for i in null_space_coll_space_index]
                path_row_reduced = np.array(path_row_reduced).astype(int).reshape(1,-1).tolist()
                path_row_reduced = np.array(path_row_reduced).astype(int).tolist()
                row = DM(path_row_reduced, domain).matmul(null_space_coll_space_inv)
                row = np.array(row.to_Matrix()).astype(int).reshape([-1]).tolist()
                matrix.append(row)
            #correct the coefficients in the finite field case
            if type(coefficients) == int and coefficients > 0:
                matrix = np.mod(np.array(matrix), coefficients).astype(int)
            else:
                matrix = np.array(matrix)
            differential_matrices.append(matrix)
        return differential_matrices, null_spaces, paths


    #path homology and homology rank vector up to a given dimension with given coefficients
    def path_homology(self, max_dim = 3, coefficients = 0, as_vector = False):
        homology_vector = np.zeros(max_dim+1).astype(int)
        homology = {}
        homology['coefficients'] = coefficients
        homology['last dimension computed'] = max_dim
        bnd_matrices, null_spaces, paths = self.path_differentials(max_dim = max_dim+2, coefficients = coefficients)
        last_rank = 0
        for dim in range(max_dim+1):
            if type(coefficients) == int and coefficients > 1:
                homology_vec = np.abs(np.array(smith_normal_form(Matrix(bnd_matrices[dim].astype(int)),
                                                                 domain = GF(coefficients))).diagonal())
            else:
                homology_vec = np.abs(np.array(smith_normal_form(Matrix(bnd_matrices[dim]),
                                                                 domain = QQ)).diagonal())
            rank = np.sum(homology_vec != 0)
            free_rank = len(null_spaces[dim]) - last_rank - rank
            if free_rank > 0:
                homology_vector[dim] = free_rank
                if type(coefficients) == int and coefficients > 1:
                        homology['dimension '+str(dim)] = ('(Z/'+str(coefficients)+'Z)',str(free_rank))
                else:
                        homology['dimension '+str(dim)] = ('Q',str(free_rank))
            last_rank = rank
        if as_vector:
            return homology_vector, bnd_matrices, null_spaces, paths
        else:
            return homology, bnd_matrices, null_spaces, paths


    
