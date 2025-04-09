import kneed
import numpy as np
import pandas as pd

def find_knee(v: pd.Series) -> float:
    """ Computes a threshold score for spatial variability.
    
        Takes in a Series v of spatial structrue scores and returns a cutoff point to use for deciding 
        which genes are spatially variable. To be called by compute_persistence.run_persistence().

        Parameters
        ----------
        v : pd.Series
            Series of CoSS's.
            
        Returns
        ----------
        float
            Threshold for declaring a gene to be Spatially Variable.
        
    """
    
    elbow = kneed.KneeLocator(x=np.linspace(1, v.size, v.size), 
                              y=v.sort_values(ascending=False), 
                              curve="convex", direction="decreasing", S=1)
    
    return elbow.knee