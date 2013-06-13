# -*- coding: utf-8 -*-
"""
Ce programme simule une constellation de satellites GPS et le déplacement
d'un mobile (avion, par exemple) sur une trajectoire 3D.
Le mobile reçoit (de manière silulée) à des intervalles réguliers diverses
informations provenant des satellites en visibilité à l'instant T.
Ces paramètres permettent de mettent en oeuvre le calcul de triangulation
et de résoudre le problème d'estimation de paramètres pour le calcul
de positionnement par GPS.

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

import numpy as np
import matplotlib.pyplot as plt


def walker(T, P, F, a, e, i, w, RE, v):
    """
    Calcul des paramètres de la constellation.

    :param a: en mètres
    :param i: en degrés
    :param v: indique s'il faut afficher la sortie ou non

    :type v: bool

    :rtype: ndarray
    """
    CDR = np.pi / 180

    CONST = np.zeros((T, 12))
    ## modification des colonnes
    CONST[:, 0] = a + RE  # a (m)
    CONST[:, 1] = e
    CONST[:, 2] = i * CDR  # i (rad)
    CONST[:, 3] = w

    #for j in xrange(P):
    for j in range(P):
        DM = j * 2 * np.pi * F / T
        ab = j * 2 * np.pi / P
        CONST[int(j*T/P):int((j+1)*T/P), 4] = ab
        for k in range(int(T/P)):
            _or = (DM + k * 2 * np.pi * P / T) % (2 * np.pi)
            CONST[int((k+1)+j*T/P)-1, 5] = _or

    if v:
        ## affichage
        plt.plot(CONST[:, 4], CONST[:, 5], 'r.')  # red points
        # axes
        plt.xlim(0, 2*np.pi)
        plt.ylim(0, 2*np.pi)

        plt.suptitle("Position")
        plt.show()

    return CONST
