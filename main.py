#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Programme principal

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

import configparser
import matplotlib.pyplot as plt
import numpy as np
import os.path

import affichage
from calculs import calc_cercle_visi
from params import (LON_UT_DEP, LON_UT_AR, LON_UT, LON_INI,
                    LAT_UT_DEP, LAT_UT_AR, LAT_UT, LAT_INI,
                    H_UT_DEP, H_UT_AR, PSOL, P0)
#from params import *
import simulation_gps


plt.ion()

if __name__ == '__main__':
    ## reinitialisation du random generateur
    np.random.seed(0)

    ## affichage carte + trajet
    affichage.worldmap(0)

    inputs = input("Saisir les coordonnées du trajet ? [o/N]\n")
    while inputs not in ['o', 'n', 'N', '']:
        print("Réponses possibles : o, n, N, <ENTREE>")
        inputs = input("Saisir les coordonnées du trajet ? [o/N]\n")

    if inputs == 'o' or not os.path.exists("trajet.cfg"):
        if inputs == 'o':
            from saisie_pos_user import *

        ## écrit les paramètres dans un fichier
        config = configparser.ConfigParser()
        config['depart'] = {
                            'lon': LON_UT_DEP,
                            'lat': LAT_UT_DEP,
                            'h': H_UT_DEP,
                           }
        config['arrivee'] = {
                            'lon': LON_UT_AR,
                            'lat': LAT_UT_AR,
                            'h': H_UT_AR,
                           }
        with open('trajet.cfg', 'w') as configfile:
            config.write(configfile)

    ## Trajet
    coord_lon_ut = LON_UT_DEP, LON_UT_AR, LON_UT, LON_INI
    coord_lat_ut = LAT_UT_DEP, LAT_UT_AR, LAT_UT, LAT_INI
    contour = affichage.trajet_ut(coord_lon_ut, coord_lat_ut)

    ## zoom
    affichage.zoom_depart(P0, PSOL, LON_UT, LAT_UT)

    ## calcul du cercle de visibilité
    DLAT_VIS, DLON = calc_cercle_visi()

    ## début de la simulation
    simulation_gps.main(DLAT_VIS, DLON, contour)

plt.ioff()
input("\aAppuyez sur Entrée pour quitter")
#plt.show()
