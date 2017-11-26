from numpy import pi, exp, log, sqrt, cos, sin, real, imag, conjugate, floor, ceil, \
    ndarray, array, arange, linspace, empty, shape, meshgrid, amax
import numpy.linalg as la
import matplotlib.pyplot as plt


### Definitions des parametres


# Lecture du fichier de configuration

a, k = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 0
for line in open('parametres.txt', 'r'):
    a[k] = line[:].rstrip('\n')
    k += 1
open('parametres.txt', 'r').close()
a[0:1] = []

# Initialisation des parametres externes

N = float(a[0])         # Ordre d'approximation des fenetres duales. Tests effectues a l'ordre 0 par manque de temps.
lambda0 = float(a[1])   # Longueur d'onde en metres
theta = float(a[2])     # Colatitude du vecteur ki en radians
phi = float(a[3])       # longitude du vecteur ki en radians
E = float(a[4])         # Intensité du champ E en volts par metres
Lx = float(a[5])        # Largeur de fenetre gaussienne en metres. Tenter 6, 8, 10, 20, 30 lambda0
nux = float(a[6])       # Coefficient d'echantillonnage. 0.25, 0.5, 0.375 ou 0.75
epsilon = float(a[7])   # Seuil de precision
Nx = float(a[8])        # Nombre de valeurs de x ou kx
Nx = int(Nx)

# Definition de l'onde incidente Ei

k0 = 2 * pi / lambda0                                                       # Nombre d'onde
ki = k0 * array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])     # Vecteur d'onde (coordonnees spheriques)


def Ei(x, y):
    """ Onde plane incidente initiale dans le plan z = 0 """
    return E * exp(1j * (x*ki[0] + y*ki[1]))


# Definition des frames: on choisit le meme pour les variables x et y

A = {0.25: 3.9703, 0.375: 2.5025, 0.5: 1.6693, 0.75: 0.5954}    # Limites du frame pour w gaussienne.
B = {0.25: 4.0299, 0.375: 2.8309, 0.5: 2.3607, 0.75: 2.0725}
Ox = 2 * pi / Lx        # Largeur spectrale
xb = sqrt(nux) * Lx     # Pas de translation spatiale
kxb = sqrt(nux) * Ox    # Pas de translation spectrale
pas = lambda0 / 10      # Pas de translation pour les valeurs de x pour E (xE)
kpas = 2 * pi / pas     # Pas de translation pour les valeurs de kx pour E (kxE)

# Valeurs de x pour calculer les wtd

xmax = 4 * Lx
xmin = - xmax

# Valeurs de x et y pour les champs incidents et reconstruits Ei et Er.
# shape(X) renvoit un tuple contenant la taille de X. Il faut que shape(xE)[0] > shape(yE)[0]

xmaxE = 8 * lambda0
xminE = - xmaxE
xE = arange(xminE, xmaxE+pas, pas)
IxE = arange(shape(xE)[0], dtype=int)    # Ensemble des indices de xE

ymaxE = 6 * lambda0
yminE = - ymaxE
yE = arange(yminE, ymaxE+pas, pas)
IyE = arange(shape(yE)[0], dtype=int)    # Ensemble des indices de yE

X, Y = meshgrid(xE, yE, indexing='ij')   # Pour tracer Ei (matplotlib).
# "Matrice" avec les valeurs de xE et yE pour lignes et colonnes respectivement.


# Valeurs de kx pour tracer les wtd

kxmax = 4 * Ox
kxmin = - kxmax
kx = linspace(kxmin, kxmax, Nx, endpoint=True)

# Valeurs de kx pour tracer les wtd

kxmaxE = 8 * k0
kxminE = - kxmaxE
kxE = arange(kxminE, kxmaxE + kpas, kpas)


### Definition des indices et intervalles de calculs et sommations


we = Lx * sqrt(-log(epsilon)/pi)    # Demi-largeur spatiale pour un seuil epsilon
oe = Ox * sqrt(-log(epsilon)/pi)    # Demi-largeur spectrale pour un seuil epsilon


# Indices de translation spatiales pour les fenetres duales

Mmin = floor((xmin - we)/xb)
Mmax = ceil((xmax + we)/xb)
IM = arange(Mmin, Mmax+1, dtype=int)

# Indices de translations spectrales

Nmin = floor((kxmin - oe)/kxb)
Nmax = ceil((kxmax + oe)/kxb)
IN = arange(Nmin, Nmax+1, dtype=int)

# Indices de translation spatiales pour le frame wmn decomposant Er

MminE = floor((xminE - we)/xb)
MmaxE = ceil((xmaxE + we)/xb)
IME = arange(MminE, MmaxE+1, dtype=int)

# Indices de translation spectrales pour wmn

NminE = floor((kxminE - oe)/kxb)
NmaxE = ceil((kxmaxE + oe)/kxb)
INE = arange(NminE, NmaxE+1, dtype=int)


# Indices de translation spatiales dans l'expression de alpha

def Mamin(r):
    """ Minimum des indices m tels de le produit de wmn et wrs soit non negligeable au seuil epsilon """
    return floor(r - 2*we/xb)


def Mamax(r):
    """ Maximum des indices m tels de le produit de wmn et wrs soit non negligeable au seuil epsilon """
    return ceil(r + 2*we/xb)


def IMa(r):
    """ Intervalle des indices m tels de le produit de wmn et wrs soit non negligeable au seuil epsilon """
    return arange(Mamin(r), Mamax(r)+1, dtype=int)


# Indices de translation spectrales dans l'expression de alpha

def Namin(s):
    """ Minimum des indices n tels de le produit de wtnm et wtsr soit non negligeable """
    return floor(s - 2*oe/kxb)


def Namax(s):
    """ Maximum des indices n tels de le produit de wtnm et wtsr soit non negligeable """
    return ceil(s + 2*oe/kxb)


def INa(s):
    """ Intervalle des indices n tels de le produit de wtnm et wtsr soit non negligeable """
    return arange(Namin(s), Namax(s)+1, dtype=int)


### Définition des fonctions utiles à la décomposition de l'onde


# Precalcul des constantes
a = sqrt(sqrt(2)/Lx)
b = sqrt(sqrt(2)*Lx)
pLx2 = pi / Lx * Lx
pOx2 = pi / (Ox * Ox)


def w(x):
    """ Fenetre gaussienne dans le domaine spatial"""
    return a * exp(-x*x*pLx2)


def wt(kx):
    """ Fenetre gaussienne dans le domaine spectral/transforme """
    return b * exp(-kx*kx*pOx2)


def wmn(m, n, x):
    """ Frame de Gabor de limites A et B """
    return w(x - m * xb) * exp(1j * n * x * kxb)


def wtnm(n, m, kx):
    """ Frame de Gabor dans le domaine spectral, limites 2piA et 2piB """
    return wt(kx - n * kxb) * exp(-1j * m * xb * kx)


#def ps(m, n, r, s):
#    """ Produit scalaire de wmn et wrs """
#    return exp(-pi*(((m-r)*(m-r)+(n-s)*(n-s))/2+1j*(s-n)*(m+r))*nux)


def pst(n, m, s, r):
    """ Produit scalaire de wtnm et wtsr """
    return 2 * pi * exp(-pi*(((m-r)*(m-r)+(n-s)*(n-s))/2+1j*(m-r)*(n+s))*nux)


#def alpha(r, s, N):
#    """ Coefficients de l'approximation de la fenetre duale dans le domaine spatial"""
#    if N == 0:
#        return 2 * (0 == r) * (0 == s) / (A[nux] + B[nux])
#
#    elif N == 1:
#        return 4 / (A[nux]+B[nux]) * ((0 == r)*(0 == s) - ps(0, 0, r, s)/(A[nux]+B[nux]))
#
#    else:
#        c = array([[alpha(i, j, N-1) * ps(i, j, r, s) for j in INa(s)] for i in IMa(r)])
#        return alpha(r, s, 0) + alpha(r, s, N-1) - 2/(A[nux]+B[nux]) * ndarray.sum(c)


def alphat(s, r, N):
    """ Coefficients de l'approximation de la fenetre duale dans le domaine transforme"""
    if N == 0:
        return 2 * (0 == r) * (0 == s) / (2 * pi * (A[nux] + B[nux]))

    elif N == 1:
        return 4 / (2 * pi * (A[nux] + B[nux])) * ((0 == r)*(0 == s) - pst(0, 0, s, r)/(2 * pi * (A[nux] + B[nux])))

    else:
        d = array([[alphat(j, i, N-1) * pst(j, i, s, r) for i in IMa(r)] for j in INa(s)])
        return alphat(s, r, 0) + alphat(s, r, N-1) - 2/(2 * pi * (A[nux] + B[nux])) * ndarray.sum(d)


#def wd0(x):
#    """ Fenêtre duale à  l'ordre 0 """
#    return alpha(0, 0, 0) * w(x)


def wtd0(kx):
    """ Fenêtre transformee duale approximee a l'ordre 0 """
    return alphat(0, 0, 0) * wt(kx)


#def wd(N, x):
#    """ Fenêtre duale à  l'ordre N """
#    return sum(sum(array([[alpha(m, n, N) * wmn(m, n, x) for n in IN] for m in IM])))


def wtd(N, kx):
    """ Fenêtre transformee duale a  l'ordre N """
    return sum(sum(array([[alphat(n, m, N) * wtnm(n, m, kx) for m in IM] for n in IN])))


#def wdmn(wd, N, m, n, x):
#    """ Frame dual de  wd à l'ordre N """
#    return wd(N, x - m * xb) * exp(1j * n * x * kxb)


def wtdnm(wtd, N, n, m, kx):
    """ Frame dual de  wtd approxime à l'ordre N """
    return wtd(N, kx - n * kxb) * exp(-1j * m * xb * kx)


### Décomposition de l'onde sur wmn


# Précalcul des valeurs de wmn(m, n, x) et wmn(p, q, y) pour le calcul ultérieur de Er(N

W = empty((shape(IME)[0], shape(INE)[0], shape(IxE)[0]), dtype=complex)
for r in IME:
    for s in INE:
        for z in IxE:
            W[r][s][z] = wmn(r, s, xE[z])

print('W = ', W)


# Précalcul des coefficients de la décomposition

AmnpqN = empty((shape(IME)[0], shape(INE)[0], shape(IME)[0], shape(INE)[0]), dtype=complex)
for m in IME:
    for n in INE:
        for p in IME:
            for q in INE:
                AmnpqN[m][n][p][q] = conjugate(wtdnm(wtd, N, n, m, ki[0]) * wtdnm(wtd, N, q, p, ki[1]))\
                                     * exp(-1j*(m*n + p*q)*xb*kxb)

print('A = ', AmnpqN)


# Décomposition finale

Er = E * sum(sum(sum(sum(array([[[[[[AmnpqN[m][n][p][q] * W[m][n][x] * W[p][q][y] for y in IyE] for x in IxE]
                 for q in INE] for p in IME] for n in INE] for m in IME])))))

print('Er = ', Er)


### Courbes et figures


fenetres = plt.figure()     # Fenetres duales approximees aux ordres 0, 1 et 2 dans le domaine spectral


plt.subplot(121)            # Parties réelles des $wtd^N$
plt.title('Parties réelles des $wtd^N$ pour \n'
          'lambda0 = {:.2f}, Lx = {:.1f}, nux = {:.3f}, epsilon = {:.3f}, Nx = {:.0f},\n'
          'E = {:.1f}, theta = {:.3f}, phi = {:.3f}\n'
          'Mmin = {:.0f}, Mmax = {:.0f}, Nmin = {:.0f}, Nmax = {:.0f}'
          .format(lambda0, Lx, nux, epsilon, Nx, E, theta, phi, Mmin, Mmax, Nmin, Nmax))
plt.xlabel('$k_x/O_x$')
plt.ylabel('$Re(wtd^N(k_x/O_x))$')
plt.plot(kx/Ox, real(wtd0(kx/Ox)), label='wtd0', color='grey')
plt.plot(kx/Ox, real(wtd(1, kx/Ox)), label='$wtd^1$', color='red')
plt.plot(kx/Ox, real(wtd(2, kx/Ox)), label='$wtd^2$', color='blue', linestyle='--')
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')


plt.subplot(122)            # Parties imaginaires des $wtd^N$
plt.title('Parties imaginaires des $wtd^N$ pour \n'
          'lambda0 = {:.2f}, Lx = {:.1f}, nux = {:.3f}, epsilon = {:.3f}, Nx = {:.0f},\n'
          'E = {:.1f}, theta = {:.3f}, phi = {:.3f}\n'
          'Mmin = {:.0f}, Mmax = {:.0f}, Nmin = {:.0f}, Nmax = {:.0f}'
          .format(lambda0, Lx, nux, epsilon, Nx, E, theta, phi, Mmin, Mmax, Nmin, Nmax))
plt.xlabel('$k_x/O_x$')
plt.ylabel('$Im(wtd^N(k_x/O_x))$')
plt.plot(kx/Ox, imag(wtd0(kx/Ox)), label='wtd0', color='grey')
plt.plot(kx/Ox, imag(wtd(1, kx/Ox)), label='$wtd^1$', color='red')
plt.plot(kx/Ox, imag(wtd(2, kx/Ox)), label='$wtd^2$', color='blue', linestyle='--')
plt.grid(True)
plt.legend(loc='upper right', fontsize='small')


champs = plt.figure()       # Champs incident et reconstruits


plt.subplot(121)
plt.title(('Partie réelle de Ei(x,y) pour\n'
            'N = {:.0f}, lambda0 = {:.2f}, Lx = {:.1f}, nux = {:.3f}, epsilon = {:.3f}, Nx = {:.0f},\n'
            'E = {:.1f}, theta = {:.3f}, phi = {:.3f}, Mmin = {:.0f}, Mmax = {:.0f}, Nmin = {:.0f}, Nmax = {:.0f}'
            .format(N, lambda0, Lx, nux, epsilon, Nx, E, theta, phi, Mmin, Mmax, Nmin, Nmax)))
plt.xlabel('$x_E$')
plt.ylabel('$y_E$')
plt.contourf(X, Y, Ei(X, Y), cmap=plt.cm.hot)
plt.grid(True)
plt.colorbar()


plt.subplot(122)
plt.title('Partie réelle de $Er(N)$ pour\n'
            'N = {:.0f}, lambda0 = {:.2f}, Lx = {:.1f}, nux = {:.3f}, epsilon = {:.3f}, Nx = {:.0f},\n'
            'E = {:.1f}, theta = {:.3f}, phi = {:.3f}, Mmin = {:.0f}, Mmax = {:.0f}, Nmin = {:.0f}, Nmax = {:.0f}'
            .format(N, lambda0, Lx, nux, epsilon, Nx, E, theta, phi, Mmin, Mmax, Nmin, Nmax))
plt.xlabel('$x_E$')
plt.ylabel('$y_E$')
plt.contourf(X, Y, Er, cmap=plt.cm.hot)
plt.grid(True)
plt.colorbar()


erreur = plt.figure()


plt.subplot(121)
plt.title(('Partie réelle de Ei(x,y) pour\n'
            'N = {:.0f}, lambda0 = {:.2f}, Lx = {:.1f}, nux = {:.3f}, epsilon = {:.3f}, Nx = {:.0f},\n'
            'E = {:.1f}, theta = {:.3f}, phi = {:.3f}, Mmin = {:.0f}, Mmax = {:.0f}, Nmin = {:.0f}, Nmax = {:.0f}'
            .format(N, lambda0, Lx, nux, epsilon, Nx, E, theta, phi, Mmin, Mmax, Nmin, Nmax)))
plt.xlabel('$x_E$')
plt.ylabel('$y_E$')
plt.contourf(X, Y, Ei(X, Y), cmap=plt.cm.hot)
plt.grid(True)
plt.colorbar()


plt.subplot(122)
plt.title('Partie réelle de $||Ei -Er(N)||/||Eimax||$ pour\n'
            'N = {:.0f}, lambda0 = {:.2f}, Lx = {:.1f}, nux = {:.3f}, epsilon = {:.3f}, Nx = {:.0f},\n'
            'E = {:.1f}, theta = {:.3f}, phi = {:.3f}, Mmin = {:.0f}, Mmax = {:.0f}, Nmin = {:.0f}, Nmax = {:.0f}'
            .format(N, lambda0, Lx, nux, epsilon, Nx, E, theta, phi, Mmin, Mmax, Nmin, Nmax))
plt.xlabel('$x_E$')
plt.ylabel('$y_E$')
plt.contourf(X, Y, la.norm(Ei(X, Y)-Er)/amax(la.norm(Ei(X, Y))), cmap=plt.cm.hot)
plt.grid(True)
plt.colorbar()


plt.show()
