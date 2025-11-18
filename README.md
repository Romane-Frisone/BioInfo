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
















## Part 5 — Calling SNPs

**_1. Index your reference for_** 

```samtools faidx Reference_genome_chrI.fasta```

Cette ligne de commande permet d'obtenir l'index du fichier dans le fichier 'Reference_genome_chrI.fasta.fai'. Cela permet de savoir exactement où est ce que l'on cherche pour les analyses. 

**_2. Create a ```vcf``` folder and enter it._**

```
mkdir vcf # création du dossieur vcf
cd vcf #on se place dans le dossier vcf
```

**_3. Call variants from the two BAM files_**

```
bcftools mpileup -Ou -f ../Reference_genome_chrI.fasta ../BearPaw1_sorted.bam ../BearPaw2_sorted.bam ../BearPaw3_sorted.bam ../BearPaw4_sorted.bam ../BearPaw5_sorted.bam ../BearPaw6_sorted.bam ../BearPaw7_sorted.bam ../BearPaw8_sorted.bam ../RabbitSlough1_sorted.bam ../RabbitSlough2_sorted.bam ../RabbitSlough3_sorted.bam ../RabbitSlough4_sorted.bam ../RabbitSlough5_sorted.bam ../RabbitSlough6_sorted.bam ../RabbitSlough7_sorted.bam ../RabbitSlough8_sorted.bam | bcftools call -mv -Ov -o raw_variants.vcf
```
Cette ligne permet de créer un fichier 'raw_variants.vcf' contenant toutes les variations obtenues par rapport au génome de référence et à partir des alignements bam triées de tous nos individus. Ce fichier indique les variants qui correspondent à chaque séquence et leur position par rapport à la séquence du génome de référence. 

**_4. How many SNPs did you call ?_** 

```vcftools --vcf raw_variants.vcf --freq```

La sortie du terminal indique : **After filtering, kept 127 out of a possible 127 Sites**.  
Nous avons donc 127 SNPs (ou variants) dans nos données. 

De plus, cette fonction permet de connaitre la fréquence de chaque allèle. 

**_5. Using ```vfctools```, filter SNPs, using :_**   
```vcftools --vcf raw_variants.vcf --minDP 5 --max-missing 1 --min-alleles 2 --max-alleles 2 --recode --out filtered```  
**_Explain what this command is doing. How many SNPs remain ?_**

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

**_6. Using ```vfctools```, compute_**  
**_o ```allele frequencies``` (use the ```--freq --out allele_freqs option```)_**  
**_o per-site FST between two populations (create pop files based on ```Details_Barcode_Population_SRR034310```)(use the ```--weir-fst-pop option```)_**  

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

**_7. Visualize (in R) allele frequencies and FST values._**

Pour répondre à cette question, nous avons importé les fichiers obtenus dans le cluster et les avons chargés avec R dans RStudio en local.   
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


