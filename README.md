# Bioinformatique et génomique des populations - Rapport Markdown
## Romane FRISONE et Aurore PUYOOU
### 19 Novembre 2025

##  Partie 1 — Bases de Linux
**Goal:** maitriser les commandes et opérations basiques.

**_1. Connexion au cluster, organisation de l'espace de travail, chargement des fichiers_**

**_o Accéder au cluster_**
```
ssh tpXXXXXX@core.cluster.france-bioinformatique.fr
```
**_o Demander un nœud de calcul_**
```
srun --pty bash
```
**_o Charger les modules nécessaires_**
```
module load sra-tools/3.1.1
module load seqtk/1.3
module load fastqc/0.12.1
module load bwa-mem2/2.2.1 
module load samtools/1.21
module load bcftools/1.16
module load vcftools/0.1.16
```
**_o Accéder au dossier partagé contenant les données de l’examen_**
```
cd /shared/projects/tp_2556_intro_pop_gen_182845
```
**_o Créer une direction "finalexam"_**

(Créer une direction "finalexam" nous permet d'y placer les fichiers, afin de ne pas travailler sur l'espace partagé, commun à tous les étudiants.)
```
mkdir finalexam # Création du dossier
ls # vérification que le dossier à bien été créé
cd finalexam # on se déplace dans le dossier que l'on vient de créer
```
**_o Copier les fichiers de l'examen dans "finalexam" et vérifier le succès de l'opération_**
```
cp Details_Barcode_Population_SRR034310.txt ~/finalexam
cp Reference_genome_chrI.fasta ~/finalexam
cp SRR034310_10pc.fastq ~/finalexam
cd ~/finalexam
ls
```

**_2. exploration de contenu des fichiers avec bash_**

**_o Combien y a-t-il de lignes dans le fichier fastq ?_**
```
 wc -l SRR034310_10pc.fastq
```
OUTPUT : 3564012 SRR034310_10pc.fastq  
Le fichier fastq contient  3564012 lignes.

**_o Quelle est la longueur de read ?_**
```
cat SRR034310_10pc.fastq | head
```
OUTPUT :  
<img width="268" height="134" alt="image" src="https://github.com/user-attachments/assets/fb618ff2-063d-4070-9e28-51a4e9529179" />  
La longueur de read est de 36 nucléotides.

**_o Combien de codes-barres y a-t-il dans le fichier Details_Barcode_Population_SRR034310 ?_**
```
cat Details_Barcode_Population_SRR034310.txt | head
```
--> Dans ce fichier, on a un barcode par ligne. Il suffit donc de connaître le nombre de lignes du fichier.
```
wc -l Details_Barcode_Population_SRR034310.txt
```
OUTPUT : 16 Details_Barcode_Population_SRR034310.txt  
Il y a donc 16 barcodes dans ce fichier.

## Partie 2 — Quality control
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

##  Partie 3 — Démultiplexage à l’aide des Barcodes
**Goal:** Utiliser des commandes Linux classiques pour séparer les reads du fichier SRR034310_10pc.fastq.

RAPPEL :  Le Barcoding  !!!!!!!!!!!!!!!!!!!!! redire mieux !!!!!!!!!!!!!!!!!!!
Le barcoding est une technique qui consiste à marquer des échantillons d’ADN avec de courtes séquences uniques, appelées codes-barres (souvent 4 à 12 nucléotides).

**_1. Création de 16 nouveaux fichiers FASTQ, un pour chaque échantillon, en utilisant les informations de Details_Barcode_Population_SRR034310._**
```
grep -B1 -A2 "^CCCC" SRR034310_10pc.fastq | sed '/^--$/d'>BearPaw1.fastq
grep -B1 -A2 "^CCAA" SRR034310_10pc.fastq | sed '/^--$/d'>BearPaw2.fastq
grep -B1 -A2 "^CCTT" SRR034310_10pc.fastq | sed '/^--$/d'>BearPaw3.fastq
grep -B1 -A2 "^CCGG" SRR034310_10pc.fastq | sed '/^--$/d'>BearPaw4.fastq
grep -B1 -A2 "^CACA" SRR034310_10pc.fastq | sed '/^--$/d'>BearPaw5.fastq
grep -B1 -A2 "^CAAC" SRR034310_10pc.fastq | sed '/^--$/d'>BearPaw6.fastq
grep -B1 -A2 "^CATG" SRR034310_10pc.fastq | sed '/^--$/d'>BearPaw7.fastq
grep -B1 -A2 "^CAGT" SRR034310_10pc.fastq | sed '/^--$/d'>BearPaw8.fastq
grep -B1 -A2 "^CTCT" SRR034310_10pc.fastq | sed '/^--$/d'>RabbitSlough1.fastq
grep -B1 -A2 "^CTAG" SRR034310_10pc.fastq | sed '/^--$/d'>RabbitSlough2.fastq
grep -B1 -A2 "^CTTC" SRR034310_10pc.fastq | sed '/^--$/d'>RabbitSlough3.fastq
grep -B1 -A2 "^CTGA" SRR034310_10pc.fastq | sed '/^--$/d'>RabbitSlough4.fastq
grep -B1 -A2 "^GGGG" SRR034310_10pc.fastq | sed '/^--$/d'>RabbitSlough5.fastq
grep -B1 -A2 "^GGAA" SRR034310_10pc.fastq | sed '/^--$/d'>RabbitSlough6.fastq
grep -B1 -A2 "^GGTT" SRR034310_10pc.fastq | sed '/^--$/d'>RabbitSlough7.fastq
grep -B1 -A2 "^GGCC" SRR034310_10pc.fastq | sed '/^--$/d'>RabbitSlough8.fastq
```
- B1 est utilisé pour pour garder la ligne au dessus de celle qui comporte CCCC (la ligne avec le nom de la séquence)
- A2 est utilisé pour garder les deux lignes qui suivent (+ et valeurs de qualité)
- sed est utilisé pour enlever les -- que l'on obtient entre chacune de nos séquencés

Avec la commande ls, on peut vérifier que nos 16 fichiers (un pour chaque individu) ont bien été créés.

**_2. Quel est le nombre de reads assignés à chaque échantillon_**

Pour obtenir cette information, on divise le nombre de lignes de chaque fichier par 4, puisque 4 lignes correspondent à un seul read.
```
nseq=$((`wc -l < BearPaw1.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < BearPaw2.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < BearPaw3.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < BearPaw4.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < BearPaw5.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < BearPaw6.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < BearPaw7.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < BearPaw8.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < RabbitSlough1.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < RabbitSlough2.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < RabbitSlough3.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < RabbitSlough4.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < RabbitSlough5.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < RabbitSlough6.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < RabbitSlough7.fastq` / 4))
echo $nseq
```
```
nseq=$((`wc -l < RabbitSlough8.fastq` / 4))
echo $nseq
```
Récapitulatif des OUTPUTs :  
<img width="305" height="122" alt="image" src="https://github.com/user-attachments/assets/125a4aaf-0e02-4230-be34-39596e52f7b7" />

**_3. Retirer les Barcodes qui ont servi au démultiplexage, à l'aide de la fonction seqtk_**
```
seqtk trimfq -b 4 BearPaw1.fastq > BearPaw1_trimmed.fastq
 seqtk trimfq -b 4 BearPaw2.fastq > BearPaw2_trimmed.fastq
 seqtk trimfq -b 4 BearPaw3.fastq > BearPaw3_trimmed.fastq
 seqtk trimfq -b 4 BearPaw4.fastq > BearPaw4_trimmed.fastq
 seqtk trimfq -b 4 BearPaw5.fastq > BearPaw5_trimmed.fastq
 seqtk trimfq -b 4 BearPaw6.fastq > BearPaw6_trimmed.fastq
 seqtk trimfq -b 4 BearPaw7.fastq > BearPaw7_trimmed.fastq
 seqtk trimfq -b 4 BearPaw8.fastq > BearPaw8_trimmed.fastq
 seqtk trimfq -b 4 RabbitSlough1.fastq > RabbitSlough1_trimmed.fastq
 seqtk trimfq -b 4 RabbitSlough2.fastq > RabbitSlough2_trimmed.fastq
 seqtk trimfq -b 4 RabbitSlough3.fastq > RabbitSlough3_trimmed.fastq
 seqtk trimfq -b 4 RabbitSlough4.fastq > RabbitSlough4_trimmed.fastq
 seqtk trimfq -b 4 RabbitSlough5.fastq > RabbitSlough5_trimmed.fastq
 seqtk trimfq -b 4 RabbitSlough6.fastq > RabbitSlough6_trimmed.fastq
 seqtk trimfq -b 4 RabbitSlough7.fastq > RabbitSlough7_trimmed.fastq
 seqtk trimfq -b 4 RabbitSlough8.fastq > RabbitSlough8_trimmed.fastq
```
Les 4 premiers nucléotides (correspondants aux Barcodes) sont retirés. En utilisant la fonction head sur nos fichiers "trimmed", on peut vérifier que nos reads ne comptent plus que 36-4 = 32 nucléotides.



##  Partie 4 — Alignement sur un génome de référence
**Goal:** Comprendre les notions de base de l’indexation et de l’alignement.

**_1._**
**_1._**
**_1._**
**_1._**


A FINIR :)       (partie 4 + rapper Barcoding dans partie 3) 
Bonne nuit !  



