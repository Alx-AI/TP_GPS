# -*- coding: utf-8 -*-

"""Diverses fonctions de calculs

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
import algos


## calcul_cercle_visi
def calc_cercle_visi():
    """Calcul du cercle de visibilité"""
    from params import ELEVMIN, h, Re  # , a, mu

    # Rayon du cercle de visibilité (theta)
    phi = np.arcsin(np.cos(ELEVMIN) / (1 + h / Re))
    theta = np.pi / 2 - phi - ELEVMIN

    ## XXX : pas utilisés ici, déjà dans params.py
    #T = 2 * np.pi * np.sqrt((a**3) / mu) # Période orbitale (s)
    #w = np.sqrt(mu / (a**3)) # vitesse angulaire du satellite (rad.s^(-1))

    # Calcul des points extremes du cercle de visibilite
    pasR = 0.01
    DLAT_VIS = np.arange(-theta, theta, pasR)  # np.array([..])
    DLON = []
    #DLON = np.zeros(len(DLAT_VIS))
    for lat_vis in DLAT_VIS:
    #for i, lat_vis in enumerate(DLAT_VIS):
        dlat = abs(lat_vis)
        st2 = np.sin(theta)**2
        sdlat = np.sin(dlat)**2  # pas utilisé
        DLON.append(np.arcsin(np.sqrt(abs((st2 - sdlat) / (1 - sdlat)))))
        #DLON[i] = np.arcsin(np.sqrt(abs((st2 - np.sin(dlat)**2) /
        #    np.cos(dlat)**2)))

    DLON = np.array(DLON)
    return DLAT_VIS, DLON


## visi : pourrait être directement dans simulation_gps...
def sat_visi(COORDSAT, Pexacte):
    """Détermination des satellites visibles"""
    from params import NBSAT, ELEVMIN

    SatVisi = np.zeros((NBSAT, 1))
    nbvisi = 0
    indvisi = []
    for sat in range(NBSAT):
        CSAT = COORDSAT[:, sat]
        # Calcul de visibilité
        Vutsat = CSAT - Pexacte[:3, 0]  # Vecteur Ut-Sat
        Nutsat = np.linalg.norm(Vutsat, 2)  # Norme du vecteur Ut-Sat
        Nut = np.linalg.norm(Pexacte[:3, 0], 2)  # Norme vecteur position Ut
        ProdScal = sum(Vutsat * Pexacte[:3, 0])  # Produit scalaire
        elevation = np.arcsin(ProdScal / (Nutsat * Nut))  # Calcul elevation
        # Test de visibilite
        if elevation > ELEVMIN:
            indvisi.append(sat)
            SatVisi[sat] = 1
            #colorsat[sat] = 'r'
            nbvisi += 1
        #else
        #    colorsat[sat] = 'b'

        #indvisi = np.array(indvisi)
        #indvisi = np.nonzero(SatVisi == 1)[0]
        #indvisi = indvisi.tolist()

    indvisi = np.array(indvisi)
    return nbvisi, indvisi


## position_ut
def calc_pos_ut(it):
    """Calcul la position courante de l'utilisateur"""
    from params import (LON_UT_DEP, PAS_LON,
                        LAT_UT_DEP, PAS_LAT,
                        H_UT_DEP, PAS_H, Re, BH)

    ## affectations
    #lon_ut_dep, pas_lon = longitude
    #lat_ut_dep, pas_lat = latitude
    #h_ut_dep, pas_ut = hauteur

    lon_ut = LON_UT_DEP + it * PAS_LON
    lat_ut = LAT_UT_DEP + it * PAS_LAT
    h_ut = H_UT_DEP + it * PAS_H

    # vars locales
    dut = h_ut + Re
    r_ut = dut * np.cos(lat_ut)

    Pexacte = np.zeros((4, 1))
    Pexacte[0] = r_ut * np.cos(lon_ut)  # x
    Pexacte[1] = r_ut * np.sin(lon_ut)  # y
    Pexacte[2] = dut * np.sin(lat_ut)   # z
    Pexacte[3] = BH

    return lon_ut, lat_ut, h_ut, Pexacte


## reception
def reception(coordSat, Pexacte, indvisi, nbvisi, t):
    """
    Reception du signal
    function [TGps,TRecep,X,Y,Z] = mesure

    Parameters
    ----------
    coordSat : array_like
        Coordonnées des satellites
    Pexacte : array_like
        Position exacte de l'utilisateur
    indvisi : array_like
        Indices des satellites visibles
    nbvisi : int
        Nombre de satellites visibles (len(indvisi)...)
    t : int

    Return
    ------
    TGps : vecteur colonne de taille le nombre de satellites visibles
        indiquant le temps GPS à l'émission du signal
    TRecep : vecteur colonne de taille le nombre de satellites visibles
        indiquant le temps Recepteur à la réception du signal
    X : vecteur colonne de taille le nombre de satellites visibles
        indiquant la coordonnée X des satellites visibles
    Y : vecteur colonne de taille le nombre de satellites visibles
        indiquant la coordonnée Y des satellites visibles
    Z : vecteur colonne de taille le nombre de satellites visibles
        indiquant la coordonnée Z des satellites visibles
    """
    from params import c, sigma

    # Distance reelle ut-sat (sans erreur de mesure, sans biais d'horloge)
    sat = coordSat[:, indvisi]

    TGps = t * np.ones((nbvisi, 1))
    TRecep = np.zeros((nbvisi, 1))  # init
    for i in range(nbvisi):
        D = sat[:, i:i+1] - Pexacte[:3]
        N = np.linalg.norm(D, 2)
        TRecep[i] = TGps[i] + N / c + Pexacte[3]

    EPS = np.zeros((nbvisi, 1))
    EPS = algos.normrnd(0, sigma, nbvisi, 1)
    X = np.array([coordSat[0, indvisi]]).T + EPS
    EPS = algos.normrnd(0, sigma, nbvisi, 1)
    Y = np.array([coordSat[1, indvisi]]).T + EPS
    EPS = algos.normrnd(0, sigma, nbvisi, 1)
    Z = np.array([coordSat[2, indvisi]]).T + EPS

    return TGps, TRecep, X, Y, Z
