<center> <h1> <br> Explication du programme analyse_bulk.py </br> </h1> </center>

## Description

Ce programme a pour but de réaliser l'analyse secondaire et tertiaire de bulk RNA-seq. Cela comprend donc l'expression différentielle dans un premier temps, puis le Gene Set Enrichment Analysis (GSEA) et l'Overrepresentation Analysis (ORA).

## Installation

### Dépendance Python

Toutes les prochaines actions seront réalisés dans le Terminal

#### Linux 

Assurez d'avoir Python et pip d'installés, si ce n'est pas le cas, utilisez la commande suivante :

`apt install python3-pip`

Une fois cela fait, vous devez installer **Miniconda** qui contient toutes les dépendances nécessaires pour éxécuter le programme. 

Il est préférable de l'installer à la racine, c'est pourquoi utilisez la commande : 

`cd /`

Et ensuite pour installer Miniconda

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
```

```
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

Fermez votre terminal et ouvrez le de nouvau ensuite. 

#### Windows

Si python n'est pas installé, ouvrez PowerShell, faites un clic droit et sélectionnez "Windows PowerShell (Admin)". Une fois cela fait, utilisez la commande :

`Invoke-WebRequest -Uri https://www.python.org/ftp/python/3.10.4/python-3.10.4-amd64.exe -OutFile python-3.10.4-amd64.exe`


Ensuite vouslancer l'éxécutable :

`.\python-3.10.4-amd64.exe`

Une fois cela fait, il vous faut installer **MiniConda** :

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe`

`Miniconda3-latest-Windows-x86_64.exe`

Fermez votre terminal et ouvrez-le de nouveau.

#### MacOS ####

Dans un premier temps, installez python3 si ce n'est pas fait :

`curl -O https://www.python.org/ftp/python/3.10.4/python-3.10.4-macosx10.9.pkg`

`sudo installer -pkg python-3.10.4-macosx10.9.pkg -target / `

Ensuite, il faut installer **Miniconda**, et donc pour cela :

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`

`bash Miniconda3-latest-MacOSX-x86_64.sh`

Fermez votre terminal et ouvrez-le de nouveau.

### Récupération du code

Une fois sur le lien : `https://gitlab.univ-nantes.fr/E183148Z/stage_crci2na`, télécharger le code source qui contient les programmes.


### Création du fichier de design

Le programme a besoin de deux choses : un fichier de comptage, qui est le résultat de l'analyse primaire du bulk RNA-seq. Ensuite il faut créer un fichier de design qui doit respecter certaines conditions :

- Le fichier doit être au format csv (à faire sur excel)

- La première colonne doit s'appeler "SampleName" et contient toutes les conditions qui sont dans le fichier de comptage.

- La deuxième colonne s'appelle "Condition" et tout doit être séparé par des tirets "-".

- Pour la condition vacciné, il faut mettre un "v", pour les contrôles "ctrl"

## Utilisation

### Se rendre dans le dossier ou est stocké le code
Ouvrez votre terminal, et utilisez la commande `cd` pour vous rendre dans le dossier contenant le code.
Identifiez sur votre ordinateur ou se situe le code.

Sur un environnement <b>Linux</b> vous pouvez faire un clic droit, et ouvrir le dossier dans le terminal, cela vous amènera directement à l'endroit ou les codes sont stockés et donc vous n'avez plus qu'a suivre le reste.

Pour les environnements <b>MacOS</b>, faites un clic droit sur le dossier ou sont stockés les codes, puis sélectionner <i>"Lire les informations" </i>. Vous trouverez en haut le chemin d'accès. Copiez ce chemin, et une fois dans le terminal faites : `cd "coller votre chemin" `et cela vous amènera vers le dossier où sont stockés les code et maintenant vous pouvez continuer les étapes.

### Lancer les programmes

Dans le dossier où le code se situe, activez votre environnement conda, pour ce faire :

Mettez à jour pip au cas où : `python -m pip install --upgrade pip`

`conda env create -f environnement.yaml`

Puis après pour l'activer : 

`conda activate bulkEnv`

Si vous voulez le désactivez 

`conda deactivate`

Ensuite vous devez lancer le premier programme en utilisant la commande :

`python3 bulk_analyse.py` 

Une interface graphique s'ouvre ensuite ou vous devrez renseigner le fichier de comptage, le fichier de design, le nom du dossier de sortie et enfin le nom du fichier pikcle. Le fichier pickle permet de sauvegarder la normalisation et le traitement des données du fichier de comptage. Enfin, vous devez sélectionner le mode, selon le type de donnéees que vous avez. Renseignez également l'organisme sur lequel vous travaillez.

De nombreuses sorties seront générées, qui seront détaillés plus tard.

Une fois fini, lancer le second programme :

`python3 boxplot.py`

Vous devez renseigner le fichier pickle, précédemment créé, renseigner le fichier de sortie, ensuite vous pouvez rentrer soit un nombre, afin de ressortir les gènes ayant le plus haut "Rank" (valeur qui prend en compte la p-value et la Log2Fold Change), soit une liste de gènes qui devra être écrite de cette sorte "x,y,b,c". Vous devrez une nouvelle fois renseigner le mode.

Une nouvelle interface s'ouvrira, où vous pourrez sélectionner des conditions, afin de créer un fichier csv qui contient les gènes en communs entre les conditions, ainsi que les valeurs de p-value et de Log2FoldChange.

## Output

Trois dossiers vont être générés : 

- DE
- GSEA
- ORA

### DE 
    - Tout d’abord, nous avons un fichier csv généré pour chaque condition comparée à la condition contrôle. 
      
    - Nous avons également deux fichiers csv, nommé « all_genes_significatifs » contenant les gènes significatifs pour chaque condition, et dans chaque lignée. 

    - Des boxplots sont crées, pour les 20 gènes significatifs ayant la plus forte valeur absolu de Log2FoldChange et ce, pour chaque condition. Les valeurs affichés de moyenne et de médiane, sont celles correspond à tous les gènes significativement exprimés. 
      
    - Un graphique nommé « expr_genes_all_conditions » est généré, pour observer le nombre de gènes différentiellement exprimés dans chaque condition par rapport au nombre total de gènes. 
      
    - Une APC est également générée selon la lignée d’un côté, et le temps de l’autre.
      
    - Deux UpSet plot sont générés, permettant d’observer le nombre de gène en commun entre les différentes conditions. Un cut-off de 8 a été utilisé pour rendre la figure plus compréhensible.

    - Un dossier "volcano_plot" contient les volcano plots pour chaque conditions, avec à chaque fois une échelle relative, et une autre arbitraire (-9/9)
    - Un dossier "Boxplot" qui correspond soit au nombre de gènes, soit à la liste de gènes fourni. On retrouve pour chaque gène un boxplot, avec en abscisse la valeur de baseMean pour tous les réplicats, et un fichier csv, qui donne les valeurs pour chacun des points.

### GSEA 

Le GSEA (Gene Set Enrichissement Analaysis) se base sur des ensembles de gènes, en fonction de la base de données fourni (GO,KEGG,…). Un score d’enrichissement va être calculé (ES) est calculé, qui reflète le degré de sur représentation d’un ensemble S aux extrêmes d’une liste classé (ici selon le Rank) L. La liste est parcouru, en augmentant une statistique de somme courante si on rencontre un gène dans S, et la diminuant si on ne trouve pas de gènes. L’ampleur de l’incrément dépend de la corrélation du gène avec le phénotype. Ensuite, la proportion de faux positifs est calculé (FDR) pour chaque Normalised Enrichissement Score (NES). Le FDR est la probabilité qu’un ensemble avec un NES donné représente un faux positif.

    - Des fichiers csv pour chaque condition sont générés, avec une version « filtered »  ne contenant que les pathways significatifs. 
      
    - Des graphiques d’enrichissement score sont générés montrant les 5 pathways les plus significatifs.
      
    - Des graph d’interactions entre les 22 pathways avec la p-value la plus faible sont créés pour chaque condition. 

    - Des volcanos-plots sont générés pour chaque condition, avec le nom des 5 pathways les plus sous et sur exprimés.

### ORA


L’analyse d’Enrichissement des Fonctions Biologiques ou l’Analyse par SurReprésentation, est une méthode pour déterminer si des fonctions biologiques sont sur représentées dans un ensemble de gènes étudié par rapport à ce à quoi on pourrait s’attendre par hasard.

Si pour un pathway donné, nous avons 40 gènes qui y corresponde dans l’ensemble de notre génome, et que pour un groupe  de  gène d’intérêt, on en retrouve 15, alors un test de Fisher est effectué afin de savoir s’il y a bien une sur représentation dans notre groupe d’intérêt. Est-ce le hasard qui fait qu’on retrouve 15 gènes, ou alors ce pathway est bien surreprésenté dans le groupe d’intérêt. 

    - Génération de csv pour chaque condition.
      
    • Des graphiques, contenant les pathways ayant la p-value la plus basse, montrant le pourcentage de gène dans le set, ainsi que le combined score.

Les bases de deux données utilisés pour la version "Time" sont KEGG et GO_Biological_Process, et pour le mode vaccine c'est le gene_set m8 de MSigDB ainsi que KEGG.
