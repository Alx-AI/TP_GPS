# -*- coding: utf-8 -*-

"""
Initialisation des paramètres :

* constantes physiques
* paramètres de la simulation
* de la figure

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

from scipy import constants
import numpy as np
import configparser

from walker import walker
import algos


## init_params
########################################################
## Définition des constantes
########################################################
CDR = np.pi / 180  # Conversion degrés -> radians
CRD = 180 / np.pi  # Conversion radians -> degrés
Re = 6378136  # Rayon équatorial (m)
mu = 3.986005e14  # Constante gravitationnelle (m^3/s^2)
We = 7.2921151467e-5  # Vitesse angulaire de la terre (rad.s^(-1))
c = constants.c  # 2.99792458e8  # Vitesse de la lumière

########################################################
## Paramètres de la simulation
########################################################
NBIT = 20
PAS = 300  # Pas de temps de la simulation (s)
DEBUT = 0  # Temps début simulation (s)
FIN = PAS * NBIT  # Temps fin simulation (s)

########################################################
## Paramètres de la loi sur les erreurs de mesure
########################################################
# Les moyennes sont nulles
BH = 100  # écart type sur le biais d'horloge
sigma = 50  # écart type sur la mesure de distance
# Temps GPS
tgps = DEBUT
trecep = tgps + BH

########################################################
## Paramètres de la constellation
########################################################
ELEVMIN = 5 * CDR  # Elevation minimale
NBSAT = 24  # Nombre de satellites
incl = 55 * CDR  # Inclinaison du plan orbital (rad)
h = 20183614  # Altitude des satellites
a = h + Re  # demi grand axe de l'orbite (m)
## fct
CONST = walker(24, 6, 1, a, 0, incl, 0, Re, 0)  # Calcul params constellation
M0 = CONST[:, 5].T  # Anomalie moyenne initiale (rad)
OMEGA = CONST[:, 4].T  # Ascension droite du noeud ascendant (rad)

## cercle de visibilité (calcul_cercle_visi)
T = 2 * np.pi * np.sqrt((a**3) / mu)  # Periode orbitale (s)
w = np.sqrt(mu / (a**3))  # vitesse angulaire du satellite (rad.s^(-1))


########################################################
## Affichage
########################################################
# Taille du zoom en degrés
zoom = 6


## saisie_ut
########################################################
## Saisie utilisateur...
########################################################

try:
    config = configparser.ConfigParser()
    config.read("trajet.cfg")
    LON_UT_DEP, LAT_UT_DEP, H_UT_DEP = map(int, config['depart'].values())
    LON_UT_AR, LAT_UT_AR, H_UT_AR = map(int, config['arrivee'].values())
except:
    # Coordonnées par défaut
    LON_UT_DEP = 100
    LAT_UT_DEP = 45
    H_UT_DEP = 350

    LON_UT_AR = 280
    LAT_UT_AR = -25
    H_UT_AR = 8000


LON_UT = LON_UT_DEP
LAT_UT = LAT_UT_DEP
H_UT = H_UT_DEP


LON_UT = LON_UT * CDR
LAT_UT = LAT_UT * CDR
LON_UT_DEP = LON_UT_DEP * CDR
LAT_UT_DEP = LAT_UT_DEP * CDR
LON_UT_AR = LON_UT_AR * CDR
LAT_UT_AR = LAT_UT_AR * CDR


PAS_LON = (LON_UT_AR - LON_UT_DEP) / NBIT
PAS_LAT = (LAT_UT_AR - LAT_UT_DEP) / NBIT
PAS_H = (H_UT_AR - H_UT_DEP) / NBIT

DUT_DEP = H_UT_DEP + Re
DUT_AR = H_UT_AR + Re

IT = 0

## position_ut
########################################################
# Calcul de la position de l'utilisateur
########################################################

### XXX XXX
# calculs.calc_pos_user(it)
LON_UT = LON_UT_DEP + IT * PAS_LON
LAT_UT = LAT_UT_DEP + IT * PAS_LAT
H_UT = H_UT_DEP + IT * PAS_H

# vars locales
DUT = H_UT + Re
r_ut = DUT * np.cos(LAT_UT)

Pexacte = np.zeros((4, 1))
Pexacte[0] = r_ut * np.cos(LON_UT)  # X
Pexacte[1] = r_ut * np.sin(LON_UT)  # Y
Pexacte[2] = DUT * np.sin(LAT_UT)   # Z
Pexacte[3] = BH
### XXX XXX

## init_pos
# Tirage d'un point initial proche de la solution
P0 = Pexacte
EPS = np.zeros((4, 1))
EPS[:3] = algos.normrnd(0, 200000, 3, 1)  # fct
P0 = P0 + EPS
P0[3] = 0
LAT_INI = np.arcsin(P0[2] / np.linalg.norm(P0[:3], 2))
LON_INI = np.arccos(P0[0] / np.linalg.norm(P0[:2], 2))
if P0[1] < 0:
    LON_INI = 2 * np.pi - LON_INI

# XXX
POSINIT = P0
LON_SOL = []
LAT_SOL = []
PSOL = []
