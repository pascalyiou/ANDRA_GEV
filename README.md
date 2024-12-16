
    Calcul des GEV sur des modèles CMIP6
    Laura Hasbini (LSCE)

Pour le calcul des gevs il faut lancer le code «gev_ns_script.R » en ligne de commande et il prend 4 arguments (Sur les machines Obelix ça tourner avec le module R/4.0.3):
Rscript gev_ns_script.R return_period variable method return_horizon

  return_period : nombre (par défaut 100 pour le niveau à 100ans)
  variable : ‘tasmax’, ‘tasmin’ (par défaut ‘tasmax’)
  method:  R2D2 ou CDFt (par défaut R2D2) -> Ici ça correspond à la méthode de correction de biais utilisée
  return_horizon : nombre (par défaut 2100)

 

A l’intérieur de ce script il y’a plusieurs méthode de calcul de la GEV qu’il faut commenter/décommenter au besoin :

    Calcul sans la période historique. Section “Computation without using the historical data” l151
    Calcul sans le point le plus extreme. A partir de “#### Computation without the most extreme points ####” l195 et tout ce qui suit

Le code peut tourner en laissant tout décommenté parce que les fichiers ont des noms différent mais c’est nettement plus long.

 

    Figures des résultats (R)

Tout est dans le fichier « Script_GEV.Rmd »
Le fichier "Script_GEV_plot.ipynb" est un peu nettoyé. Le fichier est en .ipynb mais le code est en R à l'intérieur. 
Il n'y a pas toutes les fonctions mais celles du rapport y sont. 

