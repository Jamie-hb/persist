import numpy as np

# Compute the dtm from a point to the weighted empirical measure

def distance_to_measure_point(weights: np.ndarray, distances: np.ndarray, m: float) -> float:
    """ Computes the dtm at a single point.
    
        Given a sorted vector of probability masses and squared distances to the neighbours of a single point,  
        computes the distance to measure at that point. To be called on by smoothed_expression.distance_to_measure_weighted().

        Parameters
        ----------
        weights : np.ndarray
            Probability masses at neighbour points, pre sorted by increasing distance from central vertex. 
        distances : np.ndarray
            Squared distances to neighbour vertices, in increasing order.
        m : float
            Threshold probability mass to compute distance to. Should lie in (0,1).

        Returns
        ----------
        float
            Distance to measure at input point.
            
        """
    
    # Create indexes to sort via distance (distances will come in order of weight so this will break ties)
    indices = distances.argsort()
    # Sort the weights via distance (ascending), ties will end up broken via weight (descending)
    weights = weights[indices]
    
    # Set up counters
    mass = 0 # current mass
    nn = 0 # number of nodes added
    
    # Compute dtm
    while mass < m:
        mass += weights[nn]
        nn += 1
        
    res = distances[indices[:nn]].sum()
    
    return res


def distance_to_measure_weighted(weights: np.ndarray, dmat: np.ndarray, m: float) -> np.ndarray:
    """ Computes the distance to measure at each point.
    
        Given a vector of expression weights, a well-well distance matrix, and an ordered list of 
        distances from the central vertex in a radial network, computes the
        distance to measure at each vertex. Designed to be called upon by compute_persistence.run_persistence()

        Parameters
        ----------
        weights : np.ndarray
            Expression in each well, should be normalised to sum to 1. Dimensions n_wells x 1
        dmat : np.ndarray
            Well-well distance matrix where dmat[i,j] = d(well_i, well_j). Dimensions n_wells x n_wells
        network_distances : np.ndarray
             Sorted array of distances from a central well to its nearest neighbour wells, 
             distances[i] = d(central vertex, ith nearest other vertex). Dimensions n_wells x 1.
        m : float
            Probability mass threshold to use in distance to measure calculations. Should lie in (0,1).

        Returns
        ----------
        np.ndarray
            Distance to measure at each well. Dimensions n_wells x 1.
    
    """
     
    # Create vector to store dtms
    res = np.zeros(dmat.shape[0])
    
    for i in range(dmat.shape[0]):        
        # Compute dtm for well i
        res[i] = distance_to_measure_point(weights, dmat[i,:], m)
    
    # Invert dtm before returning
    res = res.max() - res
    
    return res