# -*- coding: utf-8 -*-
"""Gestion de l'affichage

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


## init_param
def worldmap(show):
    """Afficher la carte du monde"""
    # intialisation de la figure
    fig = plt.figure(1)
    fig.suptitle("Position")
    plt.axis([0, 360, -90, 90])
    plt.xticks(np.arange(0, 360, 50))
    plt.yticks(np.arange(-80, 100, 20))

    # chargement de la carte
    worldmap = np.loadtxt("worldmap.bin")
    worldmap[:, 0] = worldmap[:, 0] + 180
    plt.plot(worldmap[:, 0], worldmap[:, 1], 'k,')  # noir + pixels

    if show:
        plt.show()


## saisie_ut
def trajet_ut(coord_ut_lon, coord_ut_lat):
    """Affichage du trajet de l'utilisateur"""
    from params import CRD as crd
    from params import zoom

    ## coord utilisateur départ - arrivée
    lon_ut_dep, lon_ut_ar, lon_ut, lon_ini = coord_ut_lon
    lat_ut_dep, lat_ut_ar, lat_ut, lat_ini = coord_ut_lat

    plt.figure(1)
    plt.axis([0, 360, -90, 90])
    plt.xticks(np.arange(0, 360, 50))
    plt.yticks(np.arange(-80, 100, 20))

    lon_ut_ = np.array([lon_ut_dep, lon_ut_ar]) * crd
    lat_ut_ = np.array([lat_ut_dep, lat_ut_ar]) * crd
    plt.plot(lon_ut_, lat_ut_, 'g-')  # trajet
    plt.plot(lon_ut * crd, lat_ut * crd, 'g*')  # départ

    lon = [k * crd for k in [lon_ut, lon_ini]]
    lat = [k * crd for k in [lat_ut, lat_ini]]
    lons = list(map(int, sorted(lon)))  # floor
    lats = list(map(int, sorted(lat)))  # floor
    lons[0] = lons[0] - zoom
    lons[1] = lons[1] + zoom
    lats[0] = lats[0] - zoom
    lats[1] = lats[1] + zoom

    lon_ = [lons[int(k/2) % 2] for k in range(5)]  # 0 0 1 1 0
    lat_ = [lats[int(k/2) % 2] for k in range(1, 6)]  # 0 1 1 0 0
    contour, = plt.plot(lon_, lat_, 'k-.')  # contour

    #plt.show()
    plt.draw()

    return contour


## affichage_zoom
def zoom_depart(P0, psol, lon_ut, lat_ut):
    """Zoom sur la position de départ de l'utilisateur"""
    from params import CRD as crd
    from params import zoom

    psol = np.array(psol)  # force le type
    lon_sol = np.array([])
    lat_sol = np.array([])
    if psol.size:
        lat_sol = np.arcsin(psol[2] / np.linalg.norm(psol[:3], 2))
        lon_sol = np.arccos(psol[0] / np.linalg.norm(psol[:2], 2))
        if psol[1] < 0:
            lon_sol = 2 * np.pi - lon_sol

    lat_ini = np.arcsin(P0[2] / np.linalg.norm(P0[:3], 2))
    lon_ini = np.arccos(P0[0] / np.linalg.norm(P0[:2], 2))
    if P0[1] < 0:
        lon_ini = 2 * np.pi - lon_ini

    lon = [k * crd for k in [lon_ut, lon_ini]]
    lat = [k * crd for k in [lat_ut, lat_ini]]
    lons = list(map(int, sorted(lon)))  # floor
    lats = list(map(int, sorted(lat)))  # floor

    plt.figure(2)
    plt.clf()  # TODO : simplifier avec des animations
    plt.axis([lons[0]-zoom, lons[1]+zoom, lats[0]-zoom, lats[1]+zoom])

    # chargement de la carte XXX
    worldmap = np.loadtxt("worldmap.bin")
    worldmap[:, 0] = worldmap[:, 0] + 180
    plt.plot(worldmap[:, 0], worldmap[:, 1], 'k,', markersize=5)

    plt.plot(lon_ut*crd, lat_ut*crd, marker='p', color='k',
             markerfacecolor='g', markersize=16)
    plt.plot(lon_ini*crd, lat_ini*crd, marker='v', color='c',
             markerfacecolor='c', markersize=8)
    if lon_sol.size:
        plt.plot(lon_sol*crd, lat_sol*crd, marker='*', color='m',
                 markerfacecolor='m', markersize=8)

    plt.show()
    #plt.draw()
    plt.figure(1)


## aff_cercle_visi
def aff_cercle_visi(lon, lat, dlon, dlat, col, fig):
    """
    Affichage des cercles de visibilité

    lon, lat : position utilisateur (float)
    dlon, dlat : vecteurs
    col : param affichage des cercles
    """
    from params import CRD
    from matplotlib.patches import Ellipse

    ## Ellipse
    wi = dlat.max()
    he = dlon.max()
    ell = Ellipse(xy=(lon*CRD, lat*CRD), height=he*CRD,
                  width=wi*CRD, facecolor='none',
                  edgecolor='r')
    ax = fig.add_subplot(111)
    ax.add_patch(ell)
    print(ell)

    s = len(dlon)
    lon_vis = np.zeros((2, s))
    lon_vis[0, :] = lon + dlon
    lon_vis[1, :] = lon - dlon
    #if min(lon_vis < 0) or max(lon_vis > 2 * np.pi):

    dlat2 = dlat[(s-1)::-1]  # np.array (250,)

    #lon_vis_array = np.array([lon_vis[0], lon_vis[1], lon_vis[0, 0]]) * CRD
    #dlat_array = (np.array([dlat, dlat2, dlat[0]]) + lat) * CRD
    #plt.plot(lon_vis_array, dlat_array, col)

    latv = dlat + lat

    ## bordures
    ## haut du demi-cercle
    indlat = np.nonzero(latv > np.pi / 2)[0]  # 1 dimension
    latv[indlat] = np.pi - latv[indlat]
    lon_vis[:, indlat] = lon_vis[:, indlat] + np.pi
    #plt.plot(lon_vis[0, indlat]*CRD, latv[indlat]*CRD, 'g-')
    #toto = np.nonzero(latv >= np.pi)[0].size
    #if toto:
    #    print("*"*5, toto)

    ## bas du demi-cercle
    indlat = np.nonzero(latv < -np.pi / 2)[0]  # 1 dimension
    latv[indlat] = -np.pi - latv[indlat]
    lon_vis[:, indlat] = lon_vis[:, indlat] + np.pi
    lon_vis[1, :] = lon_vis[1, (s-1)::-1]
    #plt.plot(lon_vis[1, indlat]*CRD, latv[indlat]*CRD, 'b-')
    #toto = np.nonzero(latv <= -np.pi)[0].size
    #if toto:
    #    print(toto, "*"*5)

    ## côtés du demi-cercle
    lon_vis = lon_vis + 2 * np.pi
    lon_vis = lon_vis % (2 * np.pi)

    latv2 = latv[(s-1)::-1]

    #plt.plot(lon_vis[0, :]*CRD, latv*CRD, col)
    #plt.plot(lon_vis[1, :]*CRD, latv2*CRD, col)
    cercled = lon_vis[0, :]*CRD, latv*CRD
    cercleg = lon_vis[1, :]*CRD, latv2*CRD

    #plt.plot(np.array([lon_vis[0, 0], lon_vis[1, 0]])*CRD,
    #         np.array([latv2[0], latv2[0]])*CRD, col)

    #plt.show()
    return cercleg, cercled


## affichage
def evo_constellation(P0, indvisi, lon_ut, lat_ut, dlat_vis, dlon, lons, lats,
                      point_trajet, contour, cercleg, cercled, sat_plt, fig):
    """Affichage de la simulation GPS"""
    from params import NBSAT, CRD, zoom
    import configparser

    ## coord utilisateur départ - arrivée
    config = configparser.ConfigParser()
    config.read("trajet.cfg")
    try:
        lon_ut_dep, lat_ut_dep, h_ut_dep = map(eval, config['depart'].values())
        lon_ut_ar, lat_ut_ar, h_ut_ar = map(eval, config['arrivee'].values())
    except KeyError as e:
        print("Trajet : coordonnées %s manquantes" % e)
        raise
    except Exception as e:
        print(type(e).__name__, e)
        raise

    lat_ini = np.arcsin(P0[2] / np.linalg.norm(P0[:3], 2))
    lon_ini = np.arccos(P0[0] / np.linalg.norm(P0[:2], 2))

    if P0[1] < 0:
        lon_ini = 2 * np.pi - lon_ini

    bl = 'b'
    colorsat = np.array([bl] * NBSAT)
    colorsat[indvisi] = 'r'

    #point_trajet, = plt.plot(lon_ut*CRD, lat_ut*CRD, marker='*', color='k',
    #                         markerfacecolor='g', markersize=12)
    point_trajet.set_data(lon_ut*CRD, lat_ut*CRD)
    dcercle = aff_cercle_visi(lon_ut, lat_ut, dlon, dlat_vis, 'r-', fig)
    cercleg.set_data(dcercle[0])
    cercled.set_data(dcercle[1])

    # Affichage des satellites
    for sat in range(NBSAT):
        #plt.plot(lons[sat]*CRD, lats[sat]*CRD, marker='d', color=colorsat[sat],
        #         markerfacecolor=colorsat[sat])
        psat_plt = sat_plt[sat]
        psat_plt.set_data(lons[sat]*CRD, lats[sat]*CRD)
        psat_plt.set_color(colorsat[sat])
        psat_plt.set_markerfacecolor(colorsat[sat])
        # XXX : il y a plusieurs cercles !
        #colorsat_sat = colorsat[sat] + '.'
        #dcercle = aff_cercle_visi(lons[sat], lats[sat], dlon, dlat_vis, colorsat_sat)
        #cercleg.set_data(dcercle[0])
        #cercled.set_data(dcercle[1])

    ## cf. L48 (trajet_ut)
    lon = [k * CRD for k in [lon_ut, lon_ini]]
    lat = [k * CRD for k in [lat_ut, lat_ini]]
    lons = list(map(int, sorted(lon)))  # floor
    lats = list(map(int, sorted(lat)))  # floor
    lons[0] = lons[0] - zoom
    lons[1] = lons[1] + zoom
    lats[0] = lats[0] - zoom
    lats[1] = lats[1] + zoom

    lon_ = [lons[int(k/2) % 2] for k in range(5)]  # 0 0 1 1 0
    lat_ = [lats[int(k/2) % 2] for k in range(1, 6)]  # 0 1 1 0 0
    #plt.plot(lon_, lat_, 'k-.')  # contour
    contour.set_data(lon_, lat_)

    #plt.show()
    plt.draw()


