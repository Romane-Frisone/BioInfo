# Bioinformatique et génomique des populations - Rapport Markdown
## Romane FRISONE et Aurore PUYOOU
### 19 Novembre 2025

##  2 — Quality control
**Goal:** Understand how to inspect raw sequencing reads using a standard quality-control tool.

**_1.	Run FastQC on the fastq file._**

```
cp SRR034310_10pc.fastq ./fastqc_v0.12.1/FastQC/ # copie le fichier de donnée dans le dossier de FastQC
cd fast/fastqc_v0.12.1/FastQC/ # on se place dans le dossier FastQC pour faire tourner la fonction 
./fastqc SRR034310_10pc.fastq # on utilise fastqc sur notre jeu de données 
```

**_2.	On your github, in markdown:_**  
  **_o	Include screenshots or saved HTML summaries_**  
  **_o	Describe in plain language:_**  
      **_-	sequence length distribution_**  
      **_-	quality drop at 5' and 3' ends of reads (if any)_**  
      **_- presence of adapters_**  
      **_-	overrepresented sequences_**  

![Tableau représentant les statistiques basiques de notre fichier fastq](/plots/FastQC_Basic-Statistics.png)

Notre fichier 'SRR034310_10pc.fastq' a été obtenu par séquença ge Sanger / Illumina. Il continent un total et 891 003 séquences et 32Mbp de bonne qualité (aucune des séquences n'a été marquée comme de mauvaise qualité). Les séquences ont une longueur de 36 nucléotide et leur % en GC est de 54. 

![Graphique représentant la distribution de la taille des séquences de notre fichier fastq](/plots/FastQC_Sequence-Length-Distribution.png)

Le graphique représentant la distribution de la taille des séquences contenues dans notre fichier 'SRR034310_10pc.fastq' montre que **toutes les séquences ont la même taille**, c'est-à-dire **36bp**.  (En effet, pour les valeurs 35 et 37, le nombre de séquences est égale à zéro et les autres tailles ne sont même par représentées. On est en présence d'un pic unique pour la valeur de 36bp).

![Graphique représentant la qualité du séquençage pour chaque base](/plots/FastQC_Per-Base-Sequence-quality.png)

Les boites à moustaches représentant la qualité de toutes les bases séquencées le long de la séquence montrent des **différences de qualité selon la position des bases**. 
En effet, pour les 4 premières bases (correspondant aux codes barres), les séquences ont la qualité maximale (40). On observe ensuite une diminution de qualité pour les bases 5 à 11 ce qui correspond certainement à un problème d'amorçage sur cette zone. Pour finir, on obersve une diminution progressive de la qualité lorsque l'on s'approche de la fin des séquences (fin des reads vers 25 à 36 bp). Ce dernier phénomène est classique avec ce type de séquançage. 

![Graphique représentant la présence d'adaptateurs dans nos séquences](/plots/FastQC_Adapter-Content.png)

Le graphique témoignant de la présence d'adaptateurs dans les séquences de notre fichier 'SRR034310_10pc.fastq' montre que **nous n'avons pas d'adaptateurs dans nos séquences** : toutes les valeurs sont à 0 et ce pour tous les types d'adaptateurs référencés ici (en haut à droite de la figure). Cela signifie donc que la baisse de qualité observée précedemment n'est pas liée à la présence de séquences d'adaptateurs dans les reads. 

![Tableau représentant les séquences surreprésentées](/plots/FastQC_Overrepresented-Sequences.png)

Le tableau montre que **nous avons des séquences surreprésentées dans notre fichier** 'SRR034310_10pc.fastq'. Ces séquences représentent une faible portions des séquences totales de notre fichier (approximativement 0.17%, 0.14% et 0.14%, soit au total 0.45% des séquences totales). Les séquences sont composées essentiellement de N, ce qui représente que le séquenceur n'a pas réussi à déterminer les bases présentes. De plus, comme aucune correspondance n'est établie avec ces séquences, il s'emblerait qu'il s'agisse simplement d'un problème technique (la suspiçion de contaminations peut être éliminée). Ces reads devraient être enlevés pour la suite des analyses. 


**_3.	which restriction enzyme was used to create these data?_**

L'enzyme utilisée pour obtenir ces données est **_SbfI_** du fait de son site de restriction retrouvé à la suite de tous les codes barres pour chacune de nos séquences (TGCAGG).

**_4.	what is the 4 nt sequence preceeding the enzyme overhang?_**

La séquence de 4 nucléotides précédant le site de restrcition de l'enzyme est la séquence du **code barre** ajouté en laboratoire par ligation afin de pouvoir identifier les individus des échantillons tout en les poolant tous pour effectuer le séquence de plusieurs individus en simultanné. Cela permet de les reconnaitre par la suite notamment pour les traitements bioinformatiques des données. 
