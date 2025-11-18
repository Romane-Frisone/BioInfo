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

##  2 — Controle de qualité
**Objectif:** Comprendre comment inspecter les reads bruts issus du séquençage avec un outil de controle qualité.

**_1.	Utilise FastQC sur le fichier fastq_**

```
cp SRR034310_10pc.fastq ./fastqc_v0.12.1/FastQC/ # copie le fichier de donnée dans le dossier de FastQC
cd fast/fastqc_v0.12.1/FastQC/ # on se place dans le dossier FastQC pour faire tourner la fonction 
./fastqc SRR034310_10pc.fastq # on utilise fastqc sur notre jeu de données 
```    
**_2.	Sur le github, en markdown:_**  
  **_o Inclure des captures d’écran ou des résumés HTML enregistrés_**  
  **_o Décrire en termes simples :_**  
      **_-	la distribution de la longueur des séquences_**  
      **_-	la baisse de qualité aux extrémités 5' et 3' des lectures (le cas échéant)_**  
      **_- la présence d’adaptateurs_**  
      **_-	les séquences surreprésentées_**  

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


**_3. Quelle enzyme de restriction a été utilisée pour créer ces données ?_**

L'enzyme utilisée pour obtenir ces données est **_SbfI_** du fait de son site de restriction retrouvé à la suite de tous les codes barres pour chacune de nos séquences (TGCAGG).

**_4.	Quelle est la séquence de 4 nucléotides précédant le site de restriction de l'enzyme ?_**

La séquence de 4 nucléotides précédant le site de restrcition de l'enzyme est la séquence du **code barre** ajouté en laboratoire par ligation afin de pouvoir identifier les individus des échantillons tout en les poolant tous pour effectuer le séquence de plusieurs individus en simultanné. Cela permet de les reconnaitre par la suite notamment pour les traitements bioinformatiques des données. 


##  Partie 3 — Démultiplexage à l’aide des Barcodes
**Goal:** Utiliser des commandes Linux classiques pour séparer les reads du fichier SRR034310_10pc.fastq.

RAPPEL :  Le Barcoding  
Le barcoding est une technique qui consiste à marquer des échantillons d’ADN avec de courtes séquences uniques, appelées Barcodes (souvent 4 à 12 nucléotides).


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

**_1. Indexer la référence_**

On souhaite préparer le génome de référence pour permettre un alignement rapide et efficace des reads.
```
bwa-mem2 index Reference_genome_chrI.fasta 
```
**_2. Mapper chaque fichier FASTQ démultiplexé_**  
Autrement dit, on soihaite aligner chaque FASTQ démultiplexé sur le génome de référence.
```
bwa-mem2 mem Reference_genome_chrI.fasta BearPaw1_trimmed.fastq > BearPaw1.sam
bwa-mem2 mem Reference_genome_chrI.fasta BearPaw2_trimmed.fastq > BearPaw2.sam
bwa-mem2 mem Reference_genome_chrI.fasta BearPaw3_trimmed.fastq > BearPaw3.sam
bwa-mem2 mem Reference_genome_chrI.fasta BearPaw4_trimmed.fastq > BearPaw4.sam
bwa-mem2 mem Reference_genome_chrI.fasta BearPaw5_trimmed.fastq > BearPaw5.sam
bwa-mem2 mem Reference_genome_chrI.fasta BearPaw6_trimmed.fastq > BearPaw6.sam
bwa-mem2 mem Reference_genome_chrI.fasta BearPaw7_trimmed.fastq > BearPaw7.sam
bwa-mem2 mem Reference_genome_chrI.fasta BearPaw8_trimmed.fastq > BearPaw8.sam
bwa-mem2 mem Reference_genome_chrI.fasta RabbitSlough1_trimmed.fastq > RabbitSlough1.sam
bwa-mem2 mem Reference_genome_chrI.fasta RabbitSlough2_trimmed.fastq > RabbitSlough2.sam
bwa-mem2 mem Reference_genome_chrI.fasta RabbitSlough3_trimmed.fastq > RabbitSlough3.sam
bwa-mem2 mem Reference_genome_chrI.fasta RabbitSlough4_trimmed.fastq > RabbitSlough4.sam
bwa-mem2 mem Reference_genome_chrI.fasta RabbitSlough5_trimmed.fastq > RabbitSlough5.sam
bwa-mem2 mem Reference_genome_chrI.fasta RabbitSlough6_trimmed.fastq > RabbitSlough6.sam
bwa-mem2 mem Reference_genome_chrI.fasta RabbitSlough7_trimmed.fastq > RabbitSlough7.sam
bwa-mem2 mem Reference_genome_chrI.fasta RabbitSlough8_trimmed.fastq > RabbitSlough8.sam
```
**_3. Convertir le fichier d’alignement SAM en BAM, trier le fichier BAM, puis créer un index du BAM trié._**
```
samtools view -bS BearPaw1.sam > BearPaw1.bam
samtools sort BearPaw1.bam -o BearPaw1_sorted.bam
samtools index BearPaw1_sorted.bam

samtools view -bS BearPaw2.sam > BearPaw2.bam
samtools sort BearPaw2.bam -o BearPaw2_sorted.bam
samtools index BearPaw2_sorted.bam

samtools view -bS BearPaw3.sam > BearPaw3.bam
samtools sort BearPaw3.bam -o BearPaw3_sorted.bam
samtools index BearPaw3_sorted.bam

samtools view -bS BearPaw4.sam > BearPaw4.bam
samtools sort BearPaw4.bam -o BearPaw4_sorted.bam
samtools index BearPaw4_sorted.bam

samtools view -bS BearPaw5.sam > BearPaw5.bam
samtools sort BearPaw5.bam -o BearPaw5_sorted.bam
samtools index BearPaw5_sorted.bam

samtools view -bS BearPaw6.sam > BearPaw6.bam
samtools sort BearPaw6.bam -o BearPaw6_sorted.bam
samtools index BearPaw6_sorted.bam

samtools view -bS BearPaw7.sam > BearPaw7.bam
samtools sort BearPaw7.bam -o BearPaw7_sorted.bam
samtools index BearPaw7_sorted.bam

samtools view -bS BearPaw8.sam > BearPaw8.bam
samtools sort BearPaw8.bam -o BearPaw8_sorted.bam
samtools index BearPaw8_sorted.bam

samtools view -bS RabbitSlough1.sam > RabbitSlough1.bam
samtools sort RabbitSlough1.bam -o RabbitSlough1_sorted.bam
samtools index RabbitSlough1_sorted.bam

samtools view -bS RabbitSlough2.sam > RabbitSlough2.bam
samtools sort RabbitSlough2.bam -o RabbitSlough2_sorted.bam
samtools index RabbitSlough2_sorted.bam

samtools view -bS RabbitSlough3.sam > RabbitSlough3.bam
samtools sort RabbitSlough3.bam -o RabbitSlough3_sorted.bam
samtools index RabbitSlough3_sorted.bam

samtools view -bS RabbitSlough4.sam > RabbitSlough4.bam
samtools sort RabbitSlough4.bam -o RabbitSlough4_sorted.bam
samtools index RabbitSlough4_sorted.bam

samtools view -bS RabbitSlough5.sam > RabbitSlough5.bam
samtools sort RabbitSlough5.bam -o RabbitSlough5_sorted.bam
samtools index RabbitSlough5_sorted.bam

samtools view -bS RabbitSlough6.sam > RabbitSlough6.bam
samtools sort RabbitSlough6.bam -o RabbitSlough6_sorted.bam
samtools index RabbitSlough6_sorted.bam

samtools view -bS RabbitSlough7.sam > RabbitSlough7.bam
samtools sort RabbitSlough7.bam -o RabbitSlough7_sorted.bam
samtools index RabbitSlough7_sorted.bam

samtools view -bS RabbitSlough8.sam > RabbitSlough8.bam
samtools sort RabbitSlough8.bam -o RabbitSlough8_sorted.bam
samtools index RabbitSlough8_sorted.bam
```
view -bS :  transforme le fichier SAM en BAM (la version binaire moins lourde et plus pratique pour les calculs mais pas lisible par l’homme)  
sort : trie les reads par position sur le génome  
index : crée un index pour accéder rapidement aux reads (fichier BAI pour BAM Index)

**_4. Calculer les statistiques de mapping_**
```
samtools flagstat BearPaw1_sorted.bam
samtools flagstat BearPaw2_sorted.bam
samtools flagstat BearPaw3_sorted.bam
samtools flagstat BearPaw4_sorted.bam
samtools flagstat BearPaw5_sorted.bam
samtools flagstat BearPaw6_sorted.bam
samtools flagstat BearPaw7_sorted.bam
samtools flagstat BearPaw8_sorted.bam
samtools flagstat RabbitSlough1_sorted.bam
samtools flagstat RabbitSlough2_sorted.bam
samtools flagstat RabbitSlough3_sorted.bam
samtools flagstat RabbitSlough4_sorted.bam
samtools flagstat RabbitSlough5_sorted.bam
samtools flagstat RabbitSlough6_sorted.bam
samtools flagstat RabbitSlough7_sorted.bam
samtools flagstat RabbitSlough8_sorted.bam
```
On s'attend a environ 80% de reads qui matchent. Donc si on a 10000 reads, on s'attend à 8000 matchs. Et sachant qu'on a un seul chromosome, on s'attend à ce qu'il n'y ait plus que 10% qui matchent donc environ 800.  
C'est l'ordre de grandeur que l'on obtient !  
Consulter les statistiques de mapping ici :

[Mapping statistics](Mapping%20statistics)




## Partie 5 — Appel des SNP / Variant calling

**_1. Indexez votre référence pour_** 

```samtools faidx Reference_genome_chrI.fasta```

Cette ligne de commande permet d'obtenir l'index du fichier dans le fichier 'Reference_genome_chrI.fasta.fai'. Cela permet de savoir exactement où est ce que l'on cherche pour les analyses. 

**_2. Créez un dossier ```vcf``` et placez y vous._**

```
mkdir vcf # création du dossieur vcf
cd vcf #on se place dans le dossier vcf
```

**_3. Variant calling a partir des fichiers BAM_**

```
bcftools mpileup -Ou -f ../Reference_genome_chrI.fasta ../BearPaw1_sorted.bam ../BearPaw2_sorted.bam ../BearPaw3_sorted.bam ../BearPaw4_sorted.bam ../BearPaw5_sorted.bam ../BearPaw6_sorted.bam ../BearPaw7_sorted.bam ../BearPaw8_sorted.bam ../RabbitSlough1_sorted.bam ../RabbitSlough2_sorted.bam ../RabbitSlough3_sorted.bam ../RabbitSlough4_sorted.bam ../RabbitSlough5_sorted.bam ../RabbitSlough6_sorted.bam ../RabbitSlough7_sorted.bam ../RabbitSlough8_sorted.bam | bcftools call -mv -Ov -o raw_variants.vcf
```
Cette ligne permet de créer un fichier 'raw_variants.vcf' contenant toutes les variations obtenues par rapport au génome de référence et à partir des alignements bam triées de tous nos individus. Ce fichier indique les variants qui correspondent à chaque séquence et leur position par rapport à la séquence du génome de référence.  

**_4. Combien de SNP avez-vous analysés ?_**  

```vcftools --vcf raw_variants.vcf --freq```

La sortie du terminal indique : **After filtering, kept 127 out of a possible 127 Sites**.  
Nous avons donc 127 SNPs (ou variants) dans nos données. 

De plus, cette fonction permet de connaitre la fréquence de chaque allèle. 

**_5. En utilisant ```vfctools```, filtre les SNPs, en utilisant la commande suivante :_**   
**_```vcftools --vcf raw_variants.vcf --minDP 5 --max-missing 1 --min-alleles 2 --max-alleles 2 --recode --out filtered```_**  
**_Expliquez le rôle de cette commande. Combien de SNP restent-ils ?_**

```vcftools``` est l'outil qui permet de manipuler les fichiers vcf  
```--vcf raw_variants.vcf``` permet d'appeler notre fichier vcf qui est le fichier contenant toutes les variations de nos séquences  et que l'on va filtrer
```--minDP 5```  permet de filtrer les données : la profondeur de couverture (donc le nombre qui couvrent le SNP ici) doit être à minima égale à 5. Cela veut dire que les SNPs présentant trop peu de lectures sont enlevés car ils sont peu fiables (pourraient être issus d'erreurs de séquençage).     
```--max-missing 1``` permet de filtrer les données : il ne peut pas y avoir de donnée manquante sur la séquence, si il y en a la séquence est retirée pour faire les analyses  
```--min-alleles 2``` permet de filtrer les données : il doit y avoir a minima deux  allèles différents, ce qui veut dire que les sites ne présentant aucun polymorphismes ne sont pas gardés    
```--max-alleles 2``` permet de filtrer les données : il ne peut y avoir plus de 2 allèles sur chaque site   
```--recode``` permet de créer un nouveau fichier   
```--out filtered``` permet de définir le nom du fichier de sortie   'filtered.recode.vcf'   

Cette fonction permet donc d'appliquer une quantité de filtres correspondant à divers critères afin de n'avoir qu'un certain type de polymorphisme (ici bi-allélique).
La sortie du terminal indique : **After filtering, kept 4 out of a possible 127 Sites**.  
Nous n'avons donc plus que 4 SNPs (ou variants) dans nos données qui correspondent à tous ces critères. 

**Tableau récapitulatif des sites variables :**
|  | Avant filtration  | Après filtration |
| --- | --- | --- |
| Nombre de sites | 127 | 4 |


**_6. En utilisant ```vfctools```, calcule_**  
**_o les fréquences des allèles ```allele frequences``` (en utilisant les options ```--freq --out allele_freqs```)_**  
**_o les FST par sites entre les deux populations (créer des fichiers de popualtions basés sur ```Details_Barcode_Population_SRR034310```)(utilise l'option ```--weir-fst-pop```)_**  

```
vcftools --vcf raw_variants.vcf --freq --out allele_freqs 
```
Cette commande permet d'obtenir un fichier allele_freqs.frq donnant la fréquence de chacun des allèles. 

```
touch popbear.txt # Création d'un fichier texte contenant les noms appartenant à la population Bear
cat <<EOF >> popbear.txt # le EOF (end of file) permet de définir quand on aura fini de définir le texte
../BearPaw1_sorted.bam
../BearPaw2_sorted.bam
../BearPaw3_sorted.bam
../BearPaw4_sorted.bam
../BearPaw5_sorted.bam
../BearPaw6_sorted.bam
../BearPaw7_sorted.bam
../BearPaw8_sorted.bam
EOF # permet de dire que c'est la fin du fichier et donc qu'on a fini de rentrer tout les éléments que l'on voulait à l'intérieur

touch poprabbit.txt # de la même façon que pour la population Bear, on crée un fichier pour la population Rabbit 
cat <<EOF >>poprabbit.txt
../RabbitSlough1_sorted.bam
../RabbitSlough2_sorted.bam
../RabbitSlough3_sorted.bam
../RabbitSlough4_sorted.bam
../RabbitSlough5_sorted.bam
../RabbitSlough6_sorted.bam
../RabbitSlough7_sorted.bam
../RabbitSlough8_sorted.bam
EOF

 vcftools --vcf raw_variants.vcf  --weir-fst-pop popbear.txt --weir-fst-pop poprabbit.txt --out fst_pop
```
Ces commandes ont permis de créer les fichiers représentant le nom des individus de deux populations (Bear et Rabbit) afin de calculer par la suite les FST de Weir et Cockerham pour chaque locus entre ces deux populations. Le fichier 'fst_pop.fst', récapitule les différentes valeurs de FST ainsi que la position des locus considérés sur le chromosome I considéré ici. 


**_7. Visualise (dans R) les fréquences des allèles et les valeurs de FST._**

Pour répondre à cette question, nous avons importé les fichiers obtenus dans le cluster et les avons chargés avec R dans RStudio en local. 
Concernant la premmière partie de cette question, nous n'avons pas réussi à importer les données correctement dans R et de fait n'avons pas pu représenter les fréquences alléliques. 

/!\ Ces codes sont donc écrits en R : 

```r
FST_pop <- read.table("fst_pop.weir.fst", header = TRUE) #Chargement des données

FST_pop$WEIR_AND_COCKERHAM_FST <- as.numeric(FST_pop$WEIR_AND_COCKERHAM_FST) #transforme les valeurs de Fst en valeurs numériques
FST_pop_clean <- FST_pop[!is.nan(FST_pop$WEIR_AND_COCKERHAM_FST), ] #enlève les lignes/locus pour lequels les calculs de Fst n'ont pas donné de valeur mais des Nan

plot(FST_pop_clean$POS, FST_pop_clean$WEIR_AND_COCKERHAM_FST, 
     pch = 20, col = "blue",
     xlab = "Position",
     ylab = "FST",
     main = "Distribution des valeurs de FST") #trace le graphique qui suit avec les valeurs de Fst obtenues entre nos deux populations pour chaque locus représentées le long du chromosome I
```

![Graphique représentant la distribution des valeurs de FST entre les populations Rabbit et Bear sur le chromosome I](/plots/R_Distribution-Valeurs-FST.png)

Ce graphe permet de montrer les valeurs de Fst obtenues pour chaque locus entre nos deux populations Rabbit et Bear en les représentant le long du chromosome I. 



## Part 6 — Interpretation des concepts

   **_1. Quelle est la différence entre la couverture et la profondeur de couverture?_**

   
   **_2. Pourquoi les jeux de données RADseq contiennent-ils de nombreux loci avec des données manquantes ?_**

    
   **_3. Pourquoi est-il important d'appliquer un filtrage SNP avant les analyses de génomique des populations ?_**

   Il important et même nécessaire d'appliquer des filtres sur les SNPs pour les analyses en génomique afin de :  
   - n'avoir que les séquences des individus que l'on souhaite étudier (et non pas des contaminants, fragments chimériques dûs à la PCR, séquences d'adapatateurs uniquement par exemples)
   - avoir des séquences représentatives de la diversité réelle au sein de la population ou entre populations (éviter la surestimation du polymorphisme du fait d'erreurs de PCRs, meilleure estimation des fréquences alléliques)  
 - garder des informations de qualité (enlever les séquences de moins bonnes quallités ne permettant pas d'être sûrs permet de réduire les temps de calculs). 
   
   **_4. Qu'est-ce qu'un outlier locus ? Donnez une définition tirée du CM1._**

Un locus outlier est un locus qui se comporte différemment de la majorité du génome, il se démarque notamment par ses valeurs de Fst qui sont supérieures aux autres lors de comparaisons entre populations (du fait potentiellement d'une sélection positive liée à de  la sélection naturelle ou d'une adaptation locale par exemples). 

   **_5. Donnez un exemple d'organisme marin chez lequel des analyses génomiques ont détecté des îlots de différenciation._**

Un ilot de différenciation a été détecté chez le Bar en comparant des populations provenant de l'Atlantique et de la Mediterranée. Cette région du génome présentait un fort Fst témoignage d'un signal probable de sélection positive. 

![Diapositive du CM1 montrant un exmeple d'ilot de différenciation avec les test outliers](/plots/CM1_Exemple-ilot-differenciation.png)
