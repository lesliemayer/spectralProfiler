import numpy as np
import pandas as pd
import scipy.stats as ss


def continuum_correct(spectrum, nodes=None, method='linear'):
    """
    Apply a continuum correction to a given spectrum

    Parameters
    ==========
    spectrum : pd.Series
               A pandas series or Spectrum object

    nodes: list
           A list of the nodes between which piecewise continuum
           will be fit

    method : {'linear', 'regresison', 'cubic'}
             The type of regression to be fit, where 'linear' is a piecewise
             linear fit, 'regression' is an Ordinary Least Squares fit, and 
             'cubic' is a 2nd order polynomial fit.

    Returns
    =======
     : pd.Series
       The continuum corrected Spectrum
     
     : pd.Series
       The continuum line
    """
    x = spectrum.index
    y = spectrum


    if not nodes:
        nodes = [x[0], x[-1]]

    if not nodes[-1] == x[-1]:
        nodes + [x[-1]]

    nlist = list(zip(nodes, nodes[1:]))

    if len(nlist) == 1:
        # Single continuum
        n = nlist[0]
        ny = y.loc[n[0]:n[1]]
        nx = ny.index
        continuum = correction_methods[method](ny, ex=x.values)
        corrected = y / continuum
    else:
        # Piecewise continuum
        corrected = pd.Series(index=x)
        continuum = pd.Series(index=x)
        for i, n in enumerate(nlist):
            if i == 0:
                # Start
                ny = y.loc[n[0]:n[1]] # Values for continuum
                nx = y.loc[:n[1]] # Full length
                ey = y.loc[:n[1]] # Full length
                c = correction_methods[method](ny, ex=nx.index.values)
                continuum.loc[nx.index] = c
                corrected.loc[nx.index] = ey / c
            elif i == len(nlist) - 1:
                #Stop
                ny = y.loc[n[0]:n[1]] # Values for continuum
                nx = y.loc[n[0]:] # Full length
                ey = y.loc[n[0]:] # Full length
                c = correction_methods[method](ny, ex=nx.index.values)
                continuum.loc[nx.index] = c
                corrected.loc[nx.index] = ey / c
            else:
                #Mid
                ny = y.loc[n[0]:n[1]]
                c = correction_methods[method](ny, ex=ny.index.values)
                continuum.loc[ny.index] = c
                corrected.loc[ny.index] = ny / c

    return pd.Series(corrected, index=x), pd.Series(continuum, index=x)

def regression(ny, ex=None):
    """
    Parameters
    ==========
    specturm : pd.series
               Pandas Series object

    nodes : list
            of nodes to be used for the continuum

    Returns
    =======
    corrected : array
                Continuum corrected array

    continuum : array
                The continuum used to correct the data

    x : array
        The potentially truncated x values
    """

    m, b, r_value, p_value, stderr = ss.linregress(ny.index.values, ny.values)
    if not isinstance(ex, (np.ndarray, pd.Series)):
        ex = ny.index
    c = m * ex + b
    return c


def linear(ny, ex=None):
    """
    Compute a linear continuum between nx[0] and nx[-1] for the given ny
    values.  If ex is supplied, compute the slope and intercept using nx,
    ny and compute the continuum over the entire extent of ex.


    Parameters
    ----------
    ny : ndarray
         The input y values for the fit area

    ex : ndarray
         The y values for the total extent of the continuum fit area

    Returns
    -------
    c : ndarray
        The continuum
    """
    y1 = ny.iloc[0]
    y2 = ny.iloc[-1]

    wv1 = ny.index[0]
    wv2 = ny.index[-1]

    m = (y2-y1) / (wv2-wv1)
    b = y1 - (m * wv1)
    if not isinstance(ex, (np.ndarray, pd.Series)):
        ex = ny.index

    y = m * ex + b
    return y


def cubic(spectrum, nodes):
    raise(NotImplemented)

correction_methods = {'linear':linear,
                      'regression':regression,
                      'cubic': cubic}