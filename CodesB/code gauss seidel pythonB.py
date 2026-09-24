import numpy as np
from scipy.linalg import hilbert

# Initialisation de la matrice A, du vecteur b et du vecteur x
def initialisation(n):
    A = np.zeros((n, n))
    b = np.zeros(n)
    x = np.zeros(n)
    
    print("Saisie des coefficients de la matrice A :")
    for i in range(n):
        for j in range(n):
            A[i, j] = float(input(f"Entrer A[{i+1},{j+1}]: "))

    print("Saisie des coefficients du vecteur b :")
    for i in range(n):
        b[i] = float(input(f"Entrer b[{i+1}]: "))

    return A, b, x

# Vérification de la dominance diagonale
def verif_dominance_diagonale(A):
    n = A.shape[0]
    for i in range(n):
        somme_ligne = np.sum(np.abs(A[i, :])) - np.abs(A[i, i])
        if abs(A[i, i]) <= somme_ligne:
            return False
    return True

def verif_n(n, max_size):
    return 0 < n <= max_size

def verif_A_diagonale_non_nulle(A):
    return np.all(np.diag(A) != 0)

# Fonction pour calculer l'erreur (différence maximale entre y et x)
def calcul_erreur(x, y):
    return np.max(np.abs(y - x))

def gauss_seidel(A, b, x, epsilon, max_iter):
    n = len(b)
    erreur = epsilon + 1  
    y= np.zeros(n)
    m=0
    for k in range(max_iter):
        if erreur < epsilon:
            break
        else:
            m=m+1
            for i in range(n):
                
                somme1 = np.dot(A[i, :i], y[:i])  
                somme2 = np.dot(A[i, i+1:], x[i+1:])  
                y[i] = (b[i] - somme1 - somme2) / A[i, i]  

            erreur = calcul_erreur(x, y)
            x = np.copy(y)

    return x, erreur,m

max_size=100
epsilon = float(input("Entrer la tolérance (epsilon) : "))
max_iter = int(input("Entrer le nombre maximal d'itérations : "))
n = int(input("Entrer la taille de la matrice A (entier n) : "))
w=float(input("Entrer le omega"))
if not verif_n(n, max_size):
    print("Erreur : Taille de la matrice invalide.")
else:
    x = np.zeros(n)
    
    
    A = hilbert(n) + w * np.eye(n)  # Matrice de Hilbert avec w ajouté sur la diagonale
    b = np.sum(A, axis=1)

    if not verif_A_diagonale_non_nulle(A):
        print("Erreur : La matrice A a un élément diagonal nul, elle n'est pas inversible.")
    elif not verif_dominance_diagonale(A):
        print("Erreur : La matrice A n'est pas strictement dominante.")
    else:
        x, erreur,m = gauss_seidel(A, b, x, epsilon, max_iter)
        if erreur<=epsilon:
            print("\nLa solution x est :")
            print(x)
            print("le nombre d'itération est:")
            print(m)
            print("\nL'erreur finale est :")
            print(erreur)
        else:
            print("Il n'y a pas convergence")

    
