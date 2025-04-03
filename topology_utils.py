import numpy as np
import dionysus as d

# Convert a dionysus diagram object to a numpy array

def diagram_to_array(diagram, dimension=0) -> np.ndarray:
    """ Converts a persistence diagram to a numpy array.
        
        Takes in a dionysus diagram object and returns a m x 2 numpy array with the birth and death times
        of each persistence feature of a given dimension in the diagram. Of the form
                   feature_1 feature_2 ... feature_n 
        birth time     .         .     ...     .
        death time     .         .     ...     .
        To be called upon by compute_persistence.run_persistence().

        Parameters
        ----------    
        diagram : List[d.Diagram] 
            List of dionysus diagrams, one for each dimension.
        dimension : int, default=0
            Which homology dimension to take the persistent features from. 


        Returns
        ----------
        np.ndarray
            Array containing birth and death features for each feature in the specified dimension of the diagram. 
            Dimensions 2 x n_features.
            
    """
    
    res = np.zeros(( len(diagram[dimension]), 2 ))
    
    for i, pt in enumerate(diagram[dimension]):
        res[i,0] = pt.birth
        res[i,1] = pt.death
        
    return res 


# Compute the p-norm of a diagram 

def p_norm(diagram: np.ndarray, p: int =2) -> float:
    """ Computes the p-norm of a diagram given in array form. To be called upon by compute_persistence.run_persistence().

        Parameters
        ----------
        diagram : np.ndarray
            Array containing birth and death times of features from a persistence diagram. Of the form
                       feature_1 feature_2 ... feature_n 
            birth time     .         .     ...     .
            death time     .         .     ...     .
        p : int, default=2:
            Which norm to compute. 

        Returns 
        ----------
        float
            p-norm of the given diagram.
            
    """
    
    # We first remove features with no death time 
    # In particular, for a 0-dim persistence module, there will always be one feature
    # (the single cc that is the point cloud) with infinite lifetime
    features = diagram[diagram[:,1]!=np.inf]
    
    # Compute the liftime of each feature
    lifetimes = np.abs(features[:,1] - features[:,0])
    
    # Compute the p-norm of the vector of lifetimes
    if len(lifetimes)==0:
        res = 0
    elif p==np.inf:
        res = lifetimes.max()
    else:
        res = (np.sum(lifetimes**p))**(1/p)
    
    return res


# Construct the filtration for a specific set of vertices, edges, and values defined on vertices

def function_filtration(values: np.ndarray, edges: np.ndarray) -> d.Simplex:
    """ Computes the upper star filtration given am adjacency structure and values at each vertex.
        
        Takes in an adjacency structure (which vertices are adjacent to each other) a list of function values at a set of vertices, and
        computes the upper star filtration for this function on the given network structure. In the resulting filtration,     
            - index of each vertex = value of the function at that vertex
            - index of an edge = value of the function on the lowest vertex of the edge
        Higher dimensional faces not included as they do not affect the 0D PH. 
        To be called upon by compute_persistence.run_persistence().

        Parameters
        ----------
        values : np.ndarray
            Values of the function at each vertex. Dimensions n_vertices x 1
        edges: np.ndarray
            A 2d matrix specfying which vertices are adjacent. Each row [a,b]  specfies that a and b are adjacent,
            and adds the edge [a,b] to the simplex. Dimensions n_edges x 2.
            
        Returns
        ----------
        d.Simplex
            Simplex representing the function filtration.
        
    """
    
    # Created filtration using list comprehension
    num_vertices = values.size
    f = d.Filtration( [ ([i], values[i]) for i in range(num_vertices) ] +
                      [ ([edges[i,0], edges[i,1]], np.min((values[edges[i,0]], values[edges[i,1]]))) for i in range(edges.shape[0]) ] )
    
    # Sort the simplices in descending order
    f.sort(reverse=True)
    
    return f