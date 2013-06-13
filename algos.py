# -*- coding: utf-8 -*-
"""Les algorithmes

Created in MATLAB by E. Grandchamp in 2000
Adapted to PYTHON by F. Magimel in 2013
Minor modifications by J. Gergaud and D. Ruiz

Institution : INPT - ENSEEIHT
Auteurs: E. Grandchamp, F. Magimel, J. Gergaud et D. Ruiz
Copyright (C) 2013 — Tous droits réservés.
      Ce programme est un logiciel libre ; vous pouvez le redistribuer ou le
      modifier suivant les termes de la “GNU General Public License” telle que
      publiée par la Free Software Foundation : soit la version 3 de cette
      licence, soit (à votre gré) toute version ultérieure.

      Ce programme est distribué dans l’espoir qu’il vous sera utile, mais SANS
      AUCUNE GARANTIE : sans même la garantie implicite de COMMERCIALISABILITÉ
      ni d’ADÉQUATION À UN OBJECTIF PARTICULIER. Consultez la Licence Générale
      Publique GNU pour plus de détails.
                        
      Vous devriez avoir reçu une copie de la Licence Générale Publique GNU avec
      ce programme ; si ce n’est pas le cas, consultez :
      <http://www.gnu.org/licenses/>.
"""

from numpy import linalg as la
import numpy as np


## MdcNE
def MdcNE(A, b):
    """
    Résolution du problème des moindres carrés linéaire :
                Min_{alpha} || b - A*alpha ||^2
    par factorisation de Cholesky du système des équations normales.

    Parameters
    ----------
    A : np.array ou np.matrix
    b : np.array ou np.matrix

    Returns
    -------
    alpha : np.array (dans tous les cas)
        solution du problème des moindres carrés linéaire
    """
    S = np.matrix(A).T * np.matrix(A)

    # Vérification au préalable du conditionnement du système et de la stabilité
    # numérique de la résolution qui va suivre
    c = la.cond(S)
    if c > 1e16:
        print('Attention : le conditionnement de la matrice des équations')
        print('            normales est très grand ---> %0.5g' % c)

    # Factorisation suivi de la résolution
    L = la.cholesky(S)           # matrice triangulaire inférieure
    m = b.size
    bvect = np.matrix( b.reshape(m,1) )
    y = A.T * bvect
    z = la.solve(L, y)
    alpha = np.array( la.solve(L.T, z) )

    return alpha



## normrnd
def normrnd(mu, sigma, m, n):
    """
    Matrice de valeurs pseudo-aléatoires à distribution normale N(mu, sigma^2)

    Parameters
    ----------
    mu : int
        moyenne
    sigma : int
        écart-type
    m : int
        nombre de lignes
    n : int
        nombre de colonnes
    """
    #np.random.seed(0)
    return sigma * np.random.randn(m, n) + mu

    #Alternative avec le module STAT de SCIPY
    #import scipy.stats as stats
    #return stats.norm.rvs(mu, sigma, (m, n))
