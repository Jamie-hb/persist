import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from typing import List

def plot_many_genes(df: pd.DataFrame, genes: List[str], s: float =25, figwidth: float =4, figheight:float =8, numcols: int =5, 
                    title: str =None, title_fontsize: float =20, dpi: float =None, show_zero_wells: bool =False) -> None:
    """ Produces a figure with subplots showing the spatial expression of genes from the same tissue sample.

        Produces a suplot for the expression of a specified list of genes, whose expression is contained in 
        the columns of the given data frame.

        Parameters
        ----------
        df : pd.DataFrame
            pandas DataFrame of the form
                   x , y , gene1 , gene2, ... , geneN
                   .   .     .       .            .    
                   .   .     .       .            .    
            where (x, y) are the co-ordinates of each well, and genei is the expression of gene i in each well.
        genes : List[str]
            List of genes to make expression plots of. Each entry should be the name of a column in df.
        s : float, default=25
            Size of spots at each well. 's' parameter in sns.scatterplot().
        figwidth : float, default=4
            Width of each subplot. The width of the overall figure will then be this multiplied by the number of columns, passed as the 
            first element of the figsize argument in plt.subplots().
        figheight : float, default=8
            Height of each subplot. The height of the overall figure will then be this multiplied by the number of rows, passed as the 
            second element of the figsize argument in plt.subplots().
        numcols : int, default=5
            The number of genes plotted in each row of the overal figure.
        title : float, default=None
            Title for the overall figure. Passed to plt.suptitle().
        title_fontsize : float, default=20
            Font size of figure title. Passed to the fontsize argument in plt.suptitle()
        dpi : float, default=None
            Dots per inch of the overall figure, passed to the dpi parameter in plt.subplots().
        show_zero_wells : bool, default=False
            Specify whether to plot wells with zero expression. 

        Returns
        ----------
        None

    """
    
    n = len(genes)
    
    num_rows = ( (n-1)//numcols ) + 1
    num_cols = np.min([n,numcols])
    
    df = np.array(df.loc[:,["x_position", "y_position"] + genes])
    
    if not dpi==None:
        plt.subplots(num_rows, num_cols, figsize=(figwidth*num_cols, figheight*num_rows), dpi=dpi)
    else:
        plt.subplots(num_rows, num_cols, figsize=(figwidth*num_cols, figheight*num_rows))
    
    for i in range(num_rows * num_cols):
        
        plt.subplot(num_rows, num_cols, i+1)
        
        # If gene[i], make the plot
        if i < n:
            if not show_zero_wells:
                sample = df[df[:,i+2]>0]
            
            sns.scatterplot(x=sample[:,0], y=sample[:,1], hue=sample[:,i+2], s=s)
            plt.axis("equal")
            plt.title(genes[i])
        
        # We want to remove ticks and frame regardless
        plt.xticks([])
        plt.yticks([])
        plt.box(False)
    
    if not title==None:
        plt.suptitle(title, fontsize=title_fontsize)
        
    plt.show()