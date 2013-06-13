# -*- coding: utf-8 -*-

"""Saisie position utilisateur

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

from params import *

## saisie_ut
try:
    # Coordonnées utilisateur
    LON_UT_DEP = int(input('Longitude de départ (degre): '))
    LAT_UT_DEP = int(input('Latitude de départ (degre): '))
    H_UT_DEP = int(input('Altitude de départ (m): '))
    LON_UT_AR = int(input("Longitude d'arrivée (degre): "))
    LAT_UT_AR = int(input("Latitude d'arrivée (degre): "))
    H_UT_AR = int(input("Altitude d'arrivée (m): "))
except ValueError:
    print("Coordonnées par défaut utilisées")


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
