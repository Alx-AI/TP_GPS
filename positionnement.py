# -*- coding: utf-8 -*-

"""
Calcul de la position courante de l'utilisateur par résolution d'un problème
de moindres carrés non linéaire.

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



def positionnement( M0, TGps, TRecep, X, Y, Z, c, Pexacte ):
    """
    Calcul de la position de l'utilisateur
    function M = positionnement(M0,TGps,TRecep,X,Y,Z,c)

    Entrees :
        M0 : point initial M0 = [X0,Y0,Z0,BH0]' (4x1)
        TGps : vecteur colonne de taille le nombre de satellites visibles
            indiquant le temps GPS à l'émission du signal
        TRecep : vecteur colonne de taille le nombre de satellites visibles
            indiquant le temps Utilisateur à la réception du signal
        X : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée X des satellites visibles
        Y : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée Y des satellites visibles
        Z : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée Z des satellites visibles
        c : vitesse de propagation du signal
        Pexacte : vecteur colonne (utile pour les tests de convergence uniquement)

    Sortie :
        M : point solution (4x1)

    fonctions utiles :
        test_convergence
    """

    # Initialisations / récupération des données
    coordSat = np.concatenate((X.T, Y.T, Z.T))
    nbSat = len(coordSat[1])
    Mn = M0.copy()

    nbIter = 0
    critere = 0

    #####################  A COMPLETER #######################
    # Iterations de Gauss-Newton :
    #
    #####################  A COMPLETER (fin) ################

    # Affichage du nombre d'iterations et du critere d'arret ainsi
    # que des erreurs sur la position calculee et le biais d'horloge:
    ##test_convergence(Mn[:3], Mn[3], nbIter, critere, Pexacte)

    # Affectation de la solution retournee:
    M = Mn

    return M




######################################################################
#   test_convergence
######################################################################
def test_convergence(pos, BH, nbIter, critere, Pexacte):
    """
    Teste la distance entre M et la solution exacte
    Réalise une pause
    function test_convergence(Pos, BH, NbIter, critere, Pexacte);

    Parameter
    ---------
    pos : solution courante (3x1)
    BH : Biais horloge courant (1x1)
    nbIter : Nombre d'iterations effectuees
    critere : valeur du critere d'arret courant
    Pexacte : array_like (4x1)
    """
    ecart = (pos - Pexacte[:3]).T
    biais = BH - Pexacte[3]

    # Affichage du nombre d'iterations et du critere d'arret
    print()
    print("nombre d'itérations et critère d'arrêt : %3d %0.5e" % (nbIter, critere))

    # Affichage des erreurs sur la position calculee et le biais d'horloge
    print("écart par rapport à la position exacte :",)
    print(np.array2string(ecart,
                          formatter={'float_kind':lambda x: "%.4g" % x}))
    print("écart sur le biais d'horloge exact : %0.5e" % biais)

    input("Appuyez sur une touche pour continuer")




######################################################################
#   Diverses variantes pour la résolution du problème non-linéaire   #
######################################################################

def posGaussNewton1( M0, TGps, TRecep, X, Y, Z, c, Pexacte ):
    """
    Calcul de la position de l'utilisateur
    function M = positionnement(M0,TGps,TRecep,X,Y,Z,c,Pexacte)

    Entrees :
        M0 : point initial M0 = [X0,Y0,Z0,BH0]' (4x1)
        TGps : vecteur colonne de taille le nombre de satellites visibles
            indiquant le temps GPS à l'émission du signal
        TRecep : vecteur colonne de taille le nombre de satellites visibles
            indiquant le temps Utilisateur à la réception du signal
        X : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée X des satellites visibles
        Y : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée Y des satellites visibles
        Z : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée Z des satellites visibles
        c : vitesse de propagation du signal
        Pexacte : vecteur colonne (utile pour les tests de convergence uniquement)

    Sortie :
        M : point solution (4x1)

    fonctions utiles :
        test_convergence
    """

    # Initialisations / récupération des données
    coordSat = np.concatenate((X.T, Y.T, Z.T))
    nbSat = len(coordSat[1])
    Mn = M0.copy()

    # Calcul de la pseudo-distance
    R = c * (TRecep - TGps)

    # Itération de Gauss-Newton :
    nbIter = 0
    conv = 0

    while not conv:  # == 0
        # Constitution de la matrice A
        A = c * np.ones((nbSat, 4))
        b = (c*Mn[3]) * np.ones((nbSat, 1)) - R[:nbSat]  # vecteur colonne

        for i in range(nbSat):
            D = Mn[:3] - coordSat[:, i:i+1]          # vecteur colonne
            N = np.linalg.norm(D, 2)
            A[i, :3] = D.T / N
            b[i] += N

#       dM = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, b))
        dM = algos.MdcNE( A, b )

        # Calcul d'erreur relative et test de convergence
        E = abs(Mn)
        diffmax = max(abs(dM) / (E+1e-15))
        if diffmax <= 1e-8:
            conv = 1

	# Mise à jour de la position utilisateur
        Mn = Mn - dM
        nbIter += 1

#       test_convergence(Mn[:3], Mn[3]/c, nbIter, diffmax, Pexacte)

    # Affichage du nombre d'iterations et du critere d'arret ainsi
    # que des erreurs sur la position calculee et le biais d'horloge:
    test_convergence(Mn[:3], Mn[3]/c, nbIter, diffmax, Pexacte)

    # Affectation de la solution retournee
    M = Mn

    return M




def posGaussNewton2( M0, TGps, TRecep, X, Y, Z, c, Pexacte ):
    """
    Calcul de la position de l'utilisateur
    function M = positionnement(M0,TGps,TRecep,X,Y,Z,c,Pexacte)

    Entrees :
        M0 : point initial M0 = [X0,Y0,Z0,BH0]' (4x1)
        TGps : vecteur colonne de taille le nombre de satellites visibles
            indiquant le temps GPS à l'émission du signal
        TRecep : vecteur colonne de taille le nombre de satellites visibles
            indiquant le temps Utilisateur à la réception du signal
        X : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée X des satellites visibles
        Y : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée Y des satellites visibles
        Z : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée Z des satellites visibles
        c : vitesse de propagation du signal
        Pexacte : vecteur colonne (utile pour les tests de convergence uniquement)

    Sortie :
        M : point solution (4x1)

    fonctions utiles :
        test_convergence
    """

    # Initialisations / récupération des données
    coordSat = np.concatenate((X.T, Y.T, Z.T))
    nbSat = len(coordSat[1])
    Mn = M0.copy()

    # On travaille par rapport à c*tau pour éviter les pbs numériques
    Mn[3] = c * Mn[3]

    # Calcul de la pseudo-distance
    R = c * (TRecep - TGps)

    # Itération de Gauss-Newton :
    nbIter = 0
    conv = 0

    while not conv:  # == 0
        # Constitution de la matrice A
        A = np.ones((nbSat, 4))
        b = Mn[3] * np.ones((nbSat, 1)) - R[:nbSat]  # vecteur colonne

        for i in range(nbSat):
            D = Mn[:3] - coordSat[:, i:i+1]          # vecteur colonne
            N = np.linalg.norm(D, 2)
            A[i, :3] = D.T / N
            b[i] += N

#       dM = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, b))
        dM = algos.MdcNE( A, b )

        # Calcul d'erreur relative et test de convergence
        E = abs(Mn)
        diffmax = max(abs(dM) / (E+1e-15))
        if diffmax <= 1e-8:
            conv = 1

	# Mise à jour de la position utilisateur
        Mn = Mn - dM
        nbIter += 1

#       test_convergence(Mn[:3], Mn[3]/c, nbIter, diffmax, Pexacte)

    # Affichage du nombre d'iterations et du critere d'arret ainsi
    # que des erreurs sur la position calculee et le biais d'horloge:
    ##test_convergence(Mn[:3], Mn[3]/c, nbIter, diffmax, Pexacte)

    # Affectation de la solution retournee
    M = Mn
    M[3] = M[3] / c

    return M




def posNewton( M0, TGps, TRecep, X, Y, Z, c, Pexacte ):
    """
    Calcul de la position de l'utilisateur
    function M = positionnement(M0,TGps,TRecep,X,Y,Z,c,Pexacte)

    Entrees :
        M0 : point initial M0 = [X0,Y0,Z0,BH0]' (4x1)
        TGps : vecteur colonne de taille le nombre de satellites visibles
            indiquant le temps GPS à l'émission du signal
        TRecep : vecteur colonne de taille le nombre de satellites visibles
            indiquant le temps Utilisateur à la réception du signal
        X : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée X des satellites visibles
        Y : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée Y des satellites visibles
        Z : vecteur colonne de taille le nombre de satellites visibles
            indiquant la coordonnée Z des satellites visibles
        c : vitesse de propagation du signal
        Pexacte : vecteur colonne (utile pour les tests de convergence uniquement)

    Sortie :
        M : point solution (4x1)

    fonctions utiles :
        test_convergence
    """

    # Initialisations / récupération des données
    coordSat = np.concatenate((X.T, Y.T, Z.T))
    nbSat = len(coordSat[1])
    Mn = M0.copy()

    # On travaille par rapport à c*tau pour éviter les pbs numériques
    Mn[3] = c * Mn[3]

    # Calcul de la pseudo-distance
    R = c * (TRecep - TGps)

    # Itération de Newton :
    nbIter = 0
    conv = 0

    while not conv:  # == 0
        # Constitution de la matrice A
        A = np.ones((nbSat, 4))
        H = np.zeros((4, 4))
        b = Mn[3] * np.ones((nbSat, 1)) - R[:nbSat]  # vecteur colonne

        for i in range(nbSat):
            D = Mn[:3] - coordSat[:, i:i+1]          # vecteur colonne
            N = np.linalg.norm(D, 2)
            A[i, :3] = D.T / N
            b[i] += N
            H[:3, :3] += (b[i] / N) * (np.eye(3) - np.dot(D / N, (D / N).T))

        dM = np.linalg.solve(np.dot(A.T, A) + H, np.dot(A.T, b))

        # Calcul d'erreur relative et test de convergence
        E = abs(Mn)
        diffmax = max(abs(dM) / (E+1e-15))
        if diffmax <= 1e-8:
            conv = 1

	# Mise à jour de la position utilisateur
        Mn = Mn - dM
        nbIter += 1

        #test_convergence(Mn[:3], Mn[3]/c, nbIter, diffmax, Pexacte)

    # Affichage du nombre d'iterations et du critere d'arret ainsi
    # que des erreurs sur la position calculee et le biais d'horloge:
    test_convergence(Mn[:3], Mn[3]/c, nbIter, diffmax, Pexacte)

    # Affectation de la solution retournee
    M = Mn
    M[3] = M[3] / c

    return M

