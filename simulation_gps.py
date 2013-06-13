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
import time

import calculs
import affichage
from params import POSINIT, DEBUT, NBIT, PAS, c  # simulation_gps
from params import NBSAT, w, incl, We, a  # propagation
from params import OMEGA, M0  # position_ut
from positionnement import posGaussNewton2 as posFUN


## propagation
def propagation(tgps):
    """Propagation des satellites"""
    LONS = np.zeros((NBSAT, 1))
    LATS = np.zeros((NBSAT, 1))
    COORDSAT = np.zeros((3, NBSAT))

    for sat in range(NBSAT):
        # Propagation des satellites
        D = w * tgps + M0[sat]
        D1 = D % (2 * np.pi)
        if (D1 >= 0) and (D1 <= np.pi/2):
            D2 = 0
        elif (D1 <= np.pi) and (D1 > np.pi/2):
            D2 = np.pi
        elif (D1 > np.pi) and (D1 <= 3*np.pi/2):
            D2 = -np.pi
        else:
            D2 = 0

        # Calcul longitude, latitude
        LAT = np.arcsin(np.sin(incl) * np.sin(D))
        LON = OMEGA[sat] - We * tgps + np.arctan(np.tan(D) * np.cos(incl)) + D2
        LON = LON + 2 * np.pi
        LON = LON % (2 * np.pi)

        # Coord (X, Y, Z) du satellite en repère géocentrique
        r_sat = a * np.cos(LAT)
        CSAT = np.zeros(3)  # init CSAT
        CSAT[0] = r_sat * np.cos(LON)  # X
        CSAT[1] = r_sat * np.sin(LON)  # Y
        CSAT[2] = a * np.sin(LAT)  # Z

        LONS[sat] = LON
        LATS[sat] = LAT
        COORDSAT[:, sat] = CSAT

    return COORDSAT, LONS, LATS


## simulation_gps
def main(DLAT_VIS, DLON, contour):
    P0 = POSINIT

    ########################################################
    # Simulation de l'évolution de la constellation
    # Simulation de la triangulation
    ########################################################
    ## initialisation des plots
    point_trajet, = plt.plot([], [], marker='*', color='k',
                             markerfacecolor='g', markersize=12)
    cercleg, = plt.plot([], [], 'r-')
    cercled, = plt.plot([], [], 'r-')
    sat_plt = []
    for k in range(NBSAT):
        sat_plt.append(plt.plot([], [], marker='d')[0])

    # Boucle sur le temps
    for it in range(NBIT+1):
        ################################################
        # Propagation des satellites
        ################################################
        t = DEBUT + PAS * it
        coordsat, lons, lats = propagation(t)

        ################################################
        # Position utilisateur
        ################################################
        lon_ut, lat_ut, h_ut, Pexacte = calculs.calc_pos_ut(it)

        ################################################
        # Determination des satellites visibles
        ################################################
        nbvisi, indvisi = calculs.sat_visi(coordsat, Pexacte)

        ################################################
        # Calcul de position
        ################################################
        if nbvisi >= 4:
            # Acquisition des signaux
            TGps, TRecep, X, Y, Z = calculs.reception(coordsat, Pexacte,
                                                      indvisi, nbvisi, t)

            # Triangulation
            PSOL = posFUN(P0, TGps, TRecep, X, Y, Z, c, Pexacte)

            #print('PSOL = \t\t\tPexacte = ')
            #print(np.array_str(np.concatenate((PSOL, Pexacte), 1)))
            time.sleep(0.25)

        ################################################
        # Affichages
        ################################################
        #try:
        fig1 = plt.figure(1)
        affichage.evo_constellation(P0, indvisi, lon_ut, lat_ut,
                                    DLAT_VIS, DLON,
                                    lons, lats,
                                    point_trajet, contour,
                                    cercleg, cercled, sat_plt, fig1)
        time.sleep(0.25)
        affichage.zoom_depart(P0, PSOL, lon_ut, lat_ut)
        #except Exception as e:
        #    print("Fin prématurée : %s - %s !" % (type(e).__name__, e))
        #    break

        # Mémorisation position précédente
        P0 = PSOL
    print("Fin de la simulation\a")
