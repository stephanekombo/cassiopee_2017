# Cassiopée 2017 - Intérêt des LFG pour les calculs de SER

## Notations, conventions de lecture

Le nom des variables ayant un sens physique a été choisi en cohérence avec le rapport, et se base
globalement sur la phonétique (en particulier pour les lettres grecques: epsilon, phi, etc - sauf pour omega, abrégé en o).

Les indices n'ont pas été signalés explicitement afin de ne pas alourdir les notations.
Par exemple, le nombre d'onde, noté $k_0$ en LaTeX, a simplement été noté k0 dans ce programme.

Pour les fonctions définies à partir de transformations d'autres définies précédemment,
les noms ont été créés par addition de lettres symbolisant la transformation en question.
En particulier, ft est l'analogue de la fonction spatiale f dans le domaine transformé.
De plus, les translations sont d'abord effectuées avant les modulations. 
Ainsi, si w est une fenêtre gaussienne, wt en est sa transformée, wd en est sa fenêtre duale,
et wtdnm est la fenêtre duale de la fenêtre transformée, translatée n fois et modulée m fois...

###### Note:

Ces conventions sont les mêmes que celles ayant servi à définir les labels et références dans
les fichiers .tex du rapport.

## Lecture et utilisation

### Structure du code

Le nombre de '#' (1 ou 3) dans les commentaires de parties sont là pour délimiter les parties
et sous parties cruiciales du script.

#### Lecture du fichier de configuration parametres.txt (lignes 10 à 17)

Pour simplifier la lecture du code, il a fallu accepter une structure extrêmement rigide
pour le fichier de configuration. Ainsi:

* Seule la première ligne (commentée) du fichier est modifiable sans incidence.
Il y est inscrite une référence à cette documentation, ainsi que l'instruction d'utiliser
le systême international pour les unités.

* L'ordre des paramètres n'est pas supposé être modifié. Toute modification du fichier de configuration
doit être accompagné de la correction adéquate des lignes 12 à 30 du code.
Notamment, l'ajout de tout paramètre dans ce fichier nécessite l'ajout d'un élément dans la liste a
(ligne 12) et du "cast" adéquat après la ligne 29.
L'ordre des paramètres à entrer est consultable entre les lignes 21 et 29 comprises du script.
Idem pour les unités et les conditions que doivent vérifier ces paramètres.

* La méthode de lecture du fichier ne sais convertir que certaines chaînes de caractères en flottants.
Pour éviter tout problème, indiquer les décimaux sous la forme 0.5 ou utiliser la notation scientifique
5e-3 si besoin (*a priori*, seulement pour epsilon).

#### Paramétrisation (lignes 19 à 247)

À partir de la lecture du fichier de configuration, et à partir des formules du problèmes, on définit
toutes les variables servant à définir (ou tracer) les fonctions du problème, puis les fonctions
en question dans un ordre thématique qui se veut le plus logique possible.

##### Précisions supplémentaires

###### À propos de certains commentaires

Des parties du script (lignes 184-6, 194~204, 220-2, 230-2, 240-2) sont présentes, mais commentées car
inutiles à l'exécution du programme en lui-même. Ces lignes sont là car montrent bien la symétrie du
problème, selon que l'on cherche à décomposer notre onde sur un frame du domaine spatial ou spectral.
Ces lignes aident aussi à comprendre les notations adoptées dans les parties analogues du script.

###### Précalculs

Les lignes 158~161, 255~261 et 267~272 sont là car il vaut mieux appeler des valeurs précalculées
que des répéter des calculs intervenant régulièrement lors du traitement. Plus précisément:

* Les lignes 158-61 sont utiles car les calculs de puissances, de racines carrées et les divisions
sont des opérations relativement longues. Il vaut mieux les effectuer avant tout traitement plutôt
que de les faire intervenir dans le calcul final.

* Les lignes 255-261 permettent de précalculer une seule fois les valeurs de fenêtres translatées,
dont on aura besoin pour le calcul final ligne 279. Ce stratagème n'est fonctionel que
si yE est strictement inclus dans xE, d'où le choix crucial des valeurs de pas, xmaxE, ymaxE
et des tailles de IxE et IyE.

###### Prints

Les lignes 261, 274 et 282 ne sont présentent que pour vérifier la bonne exécution du script,
et peuvent être commentées/supprimées. En effet, Python n'affiche pas tous les coefficients de tableaux d'aussi grandes
dimensions, et avec des choix pas forcément pertinents quant au nombre de décimales
à afficher pour les flottants.
De plus, la bonne reconstruction des fenêtres duales et du champ incident est difficile à juger avec la données
des coefficients individuels, d'où l'intérêt des tracés en fin de script, avoir une image globale et visuelle
de la reconstruction globale.

#### Décomposition de l'onde (lignes 250 à 280)

Tout ce qui précède est là pour simplifier au maximum cette partie, au moins du point de vue lecture/algorithmique.
Cette solution n'est donc optimale *a priori*.
Pour améliorer le traitement, il peut suffir de modifier les intervalles de sommation des indices
dans les boucles correspondantes. Attention à définir les tailles des matrices W (ligne 255)
et AmnpqN (ligne 266) en cohérence!

###### Note

shape(X) renvoit un tuple, il faut considérer son élément shape(X)[0] pour le traitement.

#### Tracé des fonctions (285 à fin)

Matplotlib s'inspire grandement de MATLAB au niveau de la syntaxe et de l'utilisation.
Lire la [documentation de matplotlib.pyplot](http://matplotlib.org/api/pyplot_api.html)
pour des détails supplémentaires.

### Utilisation

#### Prérequis

Ce script a été rédigé avec [Pycharm](https://www.jetbrains.com/pycharm) mais peut-être exécuté avec
n'importe quelle IDE ou console Python, dès lors que sont installés sur l'ordinateur utilisé:

* [Python 3.X](https://www.python.org/downloads/release/python-361) (la version 3.6.1 a été utilisée, or il n'y a
*a priori* pas rétro-compatibilité avec Python2),
avec un shell (ou une console Python pour les amateurs de Windows) 
* le [ScyPy Stack](https://www.scipy.org/install.html), avec notamment:
    * [Numpy](https://www.scipy.org/scipylib/download.html) (pour les habitués de MATLAB,
    il peut être utile de consulter cette [ressource](https://docs.scipy.org/doc/numpy-dev/user/numpy-for-matlab-users.html))
    * [Matplotlib](http://matplotlib.org/users/installing.html)

Il est recommandé de toujours choisir la version la plus récente de ces language et bibliothèques.

L'ensemble {cassiopee2017_script.py, paramteres.txt} doit être placé dans un seul et même dossier pour fonctionner.
Pour être prêt-à-l'emploi, il est placé dans le dossier cassiopee2017 avec les autres livrables.
Le plus simple reste de ne pas modifier l'arborescence des fichiers.

#### Instructions

Il faut au préalable avoir rempli correctement le fichier de configuration parametres.txt
(voir les consignes précédentes). De plus, il faut s'armer de patience... Il est recommandé de lancer le programme
le soir avant de se coucher, et de regarder les résultats le lendemain, pas forcément au réveil.
Pour la validation, il vaut mieux utiliser des paramètres abusément faibles
pour obtenir un résultat certes sans intérêt physique, mais en un temps fini.
Avec un processeur Intel i5 6200U par exemple, les résultats de tests intermédiaires ont été trouvés
en 2h environ avec les paramètres:
0
0.2
0.79
0.79
130
6
0.25
1e-0
2
.

1. Se placer dans le répertoire contenant les fichiers sources (cassiopee2017 *a priori*) avec la commande *cd*
(que ce soit dans un shell ou directement dans la console Python utilisée).

##### Avec un shell

2. Exécuter la commande 'python cassiopee2017_script.py'

##### Avec une console Python

2. Lancer la console

3. Exécuter la commande 'run cassiopee2017_script.py'

Si les résultats ne viennent pas après une ou deux demie.s journée.s de calcul (ou pour toute autre raison le justifiant),
le raccourci clavier 'Ctrl+C' stoppera l'exécution du script, que ce soit dans le shell ou la console.

___

###### Remarque

* Er n'admet pas d'arguments mais dépend quand meme de xE et yE par construction.
En effet, dans le temps imparti et compte tenu du reste du code, il était beaucoup plus simple de concevoir
le champ reconstruit comme un tableau de valeur selon des valeurs variantes de x et y pour les lignes et colonnes
(respectivement) que de rédiger une fonction ne faisant pas de calcul inutil.

* L'index de sommation pour AmnpqN est discutable... Plusieurs améliorations on été suggérées par Me Letrou,
mais n'ont pas donné de résultats fructueux à cause d'un débogage du script incomplet.
Cette solution a été gardée pour le rendu final car elle demeure la plus simple.

* Pour une intégration aux subroutines C++ préexistante, il peut être intéresant
de regarder l'interfacage en C et Python via [Cython](http://cython.org/).
