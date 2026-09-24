import numpy as np
from scipy.linalg import hilbert
import time
import matplotlib.pyplot as plt

# Première version : Méthode classique
def gauss_seidel_classique(A, b, x, epsilon, max_iter):
    n = len(b)
    erreur = epsilon + 1  
    y = np.zeros(n)
    m = 0
    for k in range(max_iter):
        if erreur < epsilon:
            break
        else:
            m += 1
            for i in range(n):
                somme1 = np.dot(A[i, :i], y[:i])  
                somme2 = np.dot(A[i, i+1:], x[i+1:])  
                y[i] = (b[i] - somme1 - somme2) / A[i, i]  

            erreur = np.max(np.abs(y - x))
            x = np.copy(y)

    return x, erreur, m

# Deuxième version : Méthode matricielle
def gauss_seidel_matricielle(A, b, x, epsilon, max_iter):
    m = 0
    D = np.diag(np.diag(A))
    L = np.tril(A, -1)
    U = np.triu(A, 1)
    
    DL_inv = np.linalg.inv(D + L)  # Inversion de D+L
    F = -np.dot(DL_inv, U)         # F = -(D+L)^(-1) * U
    G = np.dot(DL_inv, b)          # G = (D+L)^(-1) * b
    
    erreur = epsilon + 1  # Initialiser l'erreur à une valeur élevée
    
    for k in range(max_iter):
        m += 1
        x_new = np.dot(F, x) + G
        erreur = np.linalg.norm(x_new - x, ord=np.inf)  # Norme infinie pour l'erreur
        x = x_new
        if erreur <= epsilon:
            break
    
    return x, erreur, m

# Partie test et comparaison des deux versions
def test_gauss_seidel(n, epsilon, max_iter, w):
    A = hilbert(n) + w * np.eye(n)
    b = np.sum(A, axis=1)
    x0 = np.zeros(n)
    
    # Test méthode classique
    start_time_1 = time.time()
    x1, erreur1, iter1 = gauss_seidel_classique(A, b, x0, epsilon, max_iter)
    time_1 = time.time() - start_time_1
    
    # Test méthode matricielle
    start_time_2 = time.time()
    x2, erreur2, iter2 = gauss_seidel_matricielle(A, b, x0, epsilon, max_iter)
    time_2 = time.time() - start_time_2

    return {
        "classique": {"temps": time_1, "iterations": iter1, "erreur": erreur1},
        "matricielle": {"temps": time_2, "iterations": iter2, "erreur": erreur2},
    }

# Comparaison sur plusieurs tailles de matrice
def comparer_performances(taille_max, epsilon, max_iter, w):
    tailles = range(2, taille_max + 1)
    temps_classique = []
    temps_matricielle = []
    iterations_classique = []
    iterations_matricielle = []
    
    for n in tailles:
        resultats = test_gauss_seidel(n, epsilon, max_iter, w)
        temps_classique.append(resultats["classique"]["temps"])
        temps_matricielle.append(resultats["matricielle"]["temps"])
        iterations_classique.append(resultats["classique"]["iterations"])
        iterations_matricielle.append(resultats["matricielle"]["iterations"])
    
    plt.figure(figsize=(12, 6))
    
    # Temps d'exécution
    plt.subplot(1, 2, 1)
    plt.plot(tailles, temps_classique, label="Classique", marker="o")
    plt.plot(tailles, temps_matricielle, label="Matricielle", marker="x")
    plt.xlabel("Taille de la matrice (n)")
    plt.ylabel("Temps d'exécution (s)")
    plt.title("Comparaison des temps d'exécution")
    plt.legend()
    plt.grid(True)

    # Nombre d'itérations
    plt.subplot(1, 2, 2)
    plt.plot(tailles, iterations_classique, label="Classique", marker="o")
    plt.plot(tailles, iterations_matricielle, label="Matricielle", marker="x")
    plt.xlabel("Taille de la matrice (n)")
    plt.ylabel("Nombre d'itérations")
    plt.title("Comparaison du nombre d'itérations")
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()

# Paramètres
taille_max = 10
epsilon = 1e-6
max_iter = 100
w = 1

comparer_performances(taille_max, epsilon, max_iter, w)

