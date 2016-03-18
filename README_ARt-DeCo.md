ARt-DeCo (Ancestral Reconstruction through DeCo)
Contact: severine.berard@umontpellier.fr / yoann.anselmetti@openmailbox.org
Last modification: 2016/03/017
Author: Yoann Anselmetti

################
### Overview ###
################

ARt-DeCo is a program for reconstructing gene adjacencies in ancestral and
extant genome. ART-DeCo is based on a previous software named DeCo
(Detection of Coevolution) (Bérard et al., 2012) that permit to determine
gene order in ancestral genome by reconstructing gene adjacency evolution.


#################
### Reference ###
#################

Y. Anselmetti, V. Berry, C. Chauve, A. Chateau, E. Tannier and S. Bérard:
Ancestral gene synteny reconstruction improves extant species scaffolding.
Accepted to RECOMB-CG 2015.


###########################
### System Requirements ###
###########################

ART-DeCo has been run on Linux system on a X86_64 architecture.
The program is encoded in C++ language using Bio++ library 
(Dutheil et al., 2006) and use a perl script to reduce the
species tree to the list of selected species.
Library requirement:
   - perl
   - gcc/g++
   - Bio++ 2.1.0: -lbpp-phyl -lbpp-seq -lbpp-core
     (See http://biopp.univ-montp2.fr/wiki/index.php/Installation)


####################
### Installation ###
####################

Unpack the ARt-DeCo package:
   unzip ARt-DeCo.gz   
For compiling ARt-DeCo:
   cd ARt-DeCo/BIN;
   make clean;   # To be sure there isn't previous compiled version
   make;

!!!
If it doesn't work, then edit the makefile file and go to the line
CFLAGS1 = -pg -g -Wall -O2 -I/usr/include/Bpp/ -L/usr/lib
Change the two "/usr" in the options to the directory where Bio++ is
installed on your system (it can be /usr/local or somewhere in your
home directory for example).
!!!

This step should create 6 binary files:
   Step0        => Use to create INPUT data files for ARt-DeCo from Ensembl database
                   (Need config file to fill the different INPUT files, reference species tree
                   and selected species file for which we have to get data).
   Step1        => Reconciliation of gene trees and clear gene in gene file not present
                   in gene trees.
   Step2        => Creation of Adjacencies classes
   Step3        => DeCo algorithm for reconstruction of ancestral adjacencies and adjacency
                   evolution.
   Step3_proba  => ARt-DeCo algorithm for reconstruction of ancestral adjacencies, adjacency
                   evolution AND prediction of new extant adjacencies.
   Step4        => Statistics on ARt-DeCo results.


###################
### Config File ###
###################

You have to prepare a configuration file. It is a set of lines where each line has
two words separated by a space. The first four lines are compulsory and tell where
your input files are:

trees_file YOUR_GENE_TREE_FILE
genes_file YOUR_GENE_FILE
species_file YOUR_SPECIES_TREE_FILE
adjacencies_file YOUR_ADJACENCY_FILE

exp_name NAME          # Will tag the output files of your run
directory DIR_NAME     # Will place the output file in directory DIR_NAME
ReconcilDone REC       # "REC=false" if gene trees are not reconciled. "REC=true" if they are
OUTPUT_FORMAT FORMAT   # FORMAT = 0 or 1, for Newick or NHX
sep CHARACTER          # Choose a character which doesn't appear in your dataset (by default |)
Adj_percentage NUMBER  # Between 0 and 1, corresponds to the percentage of time ARt-DeCo will
                         choose C1 when C0=C1
Crea COST              # Cost of an adjacency gain, by default 3
Break COST             # Cost of an adjacency breakage, by default 1

genomes_file SELECTED_SPECIES_FILE    # List of species to analyze with chromosome number expected
species_tree_input REF_SPECIES_TREE   # Reference species tree (Default species tree of Ensembl)
embl_dir DIRECTORY_EMBL_FILE          # Directory for store EMBL files if you generate data from
                                        Ensembl database (keep default if not generate data).
emf_dir DIRECTORY_EMF_FILE            # Directory for store EMF files if you generate data from
                                        Ensembl database (keep default if not generate data).
coeff_b FACT_MULT                     # Multiplicative factor for base log coefficient
                                        (by default 1.05).


###################
### INPUT files ###
###################

Here, we describe the different INPUT files for ARt-DeCo program.
If you generate your INPUT data files from Ensembl database by using "Step0".
You need to provide only 2 files and the config file with path for the different
directories/files to fill:
   - 1 file with the list of species and the chromosome number in each one of these species
      Structure:
         Species_name	Chr_number
      Example:
         Ailuropoda_melanoleuca	21
         Bos_taurus	31
         Callithrix_jacchus	24
   - 1 file containing the reference species tree on Newick or NHX format.

If you use your own data, you should provide the following files:
   - 1 file containing gene adjacencies of species that you want analyze:
      Structure:
         Gene1_ID	Gene2_ID
      Example:
         ENSMODP00000016865	ENSMODP00000008552
         ENSMODP00000008552	ENSMODP00000008623
         ENSMODP00000008623	ENSMODP00000015313
   - 1 file for Species - Gene association:
      Structure:
         Species_name	Gene_ID
      Example:
         Monodelphis_domestica	ENSMODP00000016865
         Monodelphis_domestica	ENSMODP00000008552
         Monodelphis_domestica	ENSMODP00000008623
   - 1 file for ROOTED gene trees on Newick or NHX format (Reconciled or not => "ReconcilDone"
     boolean in config file).
   - 1 file for species tree on Newick or NHX format


#############
### Usage ###
#############

To use the ART-DeCo binaries just type:
   ./Step0 YOUR_CONFIG_FILE

A bash script is available for running all steps in one command line:
   ARt-DeCo/BIN/script/deco_run.sh YOUR_CONFIG_FILE 0/1 0/1
   => Explanation of script in the head of the script code

An example of config file with corresponding files is available in ARt-DeCo/CONF_FILE.
To run ART-DeCo on this dataset, just type:
   ARt-DeCo/BIN/script/deco_run.sh CONF_FILE/DeCo.conf_79_11mammals_prot 0 1


####################
### OUTPUT files ###
####################

Art-DeCo, while running, outputs several informations on the screen.
At the end of execution, ARt-DeCo generates six output files:
   - DIR_NAME/NAME_OUTPUT_adjacencies: all the ancestral adjacencies computed by DeCo/ART-DeCo.
     It's a 3 columns file: species gene1 gene2. The genes are coded by the number of the tree
     they belong to, followed by the separator, followed by their id in this tree (the numbers
     of trees and id in trees are defined by DeCo/ARt-DeCo)
      Structure:
         Gene1_IDnb = GeneTree1_IDnb|Node_nb_Gene1
         Gene2_IDnb = GeneTree2_IDnb|Node_nb_Gene2
         Species_IDnb	Gene1_IDnb	Gene2_IDnb
      Exemple:
         70	1000|52	2288|52
         68	1000|50	2288|50
         67	1000|2	2288|2
         64	1000|49	2288|49
         54	1000|47	2288|47

   - DIR_NAME/NAME_OUTPUT_adjacencies_classes: all the ancestral adjacencies computed by DeCo.
     It's a 3 columns file: species gene1 gene2. The genes are coded by the number of the tree they
     belong to, followed by the separator, followed by their id in this tree (the numbers of trees
     and ID in trees are defined by DeCo/ARt-DeCo)
      Structure:
         Adj_IDnb = line number of the adjacency in DIR_NAME/NAME_OUTPUT_adjacencies
         AdjClass_ID = GeneTree1_IDnb|GeneTree2_IDnb|Node_nb_LCA_Gene
         AdjClass_ID	Adj_nb_in_AdjClass	List(Adj_IDnb)
      Exemple:
         0|0|1169	2	233566	233567
         0|0|1210	6	301935	301936	301937	516869	516870	516871
         0|0|1231	3	301941	301943	516873
         0|0|1239	2	301934	516868

   - DIR_NAME/NAME_OUTPUT_adjacencies_trees: all the adjacencies trees computed by DeCo/ARt-DeCo
     presented by classes.
      Example:
         Adjacency tree(s) solution of class 0|0|1110

         Tree 1:
         (((((((((ENSOARP00000009108-ENSOARP00000011840|Ovis_aries[&&NHX:D=?:Ev=Extant:S=65:ND=9],Breakage|933-1095[&&NHX:D=?:Ev=Break:S=64:ND=10])[&&NHX:D=F:Ev=Spec:S=66:ND=8])[&&NHX:D=T:Ev=GDup:S=66:ND=7],ADJLoss|940-1098[&&NHX:D=?:Ev=ALos:S=67:ND=11])[&&NHX:D=F:Ev=Spec:S=68:ND=6],ADJLoss|942-1100[&&NHX:D=?:Ev=ALos:S=69:ND=12])[&&NHX:D=F:Ev=Spec:S=70:ND=5],Loss|950-1102[&&NHX:D=?:Ev=GLos:S=71:ND=13])[&&NHX:D=F:Ev=Spec:S=72:ND=4],Loss|952-1094[&&NHX:D=?:Ev=GLos:S=63:ND=14])[&&NHX:D=F:Ev=Spec:S=73:ND=3])[&&NHX:D=T:Ev=GDup:S=73:ND=2],Loss|957-1105[&&NHX:D=?:Ev=GLos:S=50:ND=15])[&&NHX:D=F:Ev=Spec:S=74:ND=1],((((((((((Loss|959-1065[&&NHX:D=?:Ev=GLos:S=30:ND=26],Breakage|960-1064[&&NHX:D=?:Ev=Break:S=31:ND=27])[&&NHX:D=F:Ev=Spec:S=32:ND=25],ADJLoss|962-1067[&&NHX:D=?:Ev=ALos:S=33:ND=28])[&&NHX:D=F:Ev=Spec:S=34:ND=24],ADJLoss|964-1069[&&NHX:D=?:Ev=ALos:S=35:ND=29])[&&NHX:D=F:Ev=Spec:S=36:ND=23],ADJLoss|966-1071[&&NHX:D=?:Ev=ALos:S=37:ND=30])[&&NHX:D=F:Ev=Spec:S=38:ND=22],Loss|970-1073[&&NHX:D=?:Ev=GLos:S=29:ND=31])[&&NHX:D=F:Ev=Spec:S=39:ND=21],ENSCJAP00000016860-ENSCJAP00000046977|Callithrix_jacchus[&&NHX:D=?:Ev=Extant:S=40:ND=32])[&&NHX:D=F:Ev=Spec:S=41:ND=20],ADJLoss|974-1076[&&NHX:D=?:Ev=ALos:S=42:ND=33])[&&NHX:D=F:Ev=Spec:S=43:ND=19],Loss|976-1062[&&NHX:D=?:Ev=GLos:S=24:ND=34])[&&NHX:D=F:Ev=Spec:S=44:ND=18],ADJLoss|978-1079[&&NHX:D=?:Ev=ALos:S=45:ND=35])[&&NHX:D=F:Ev=Spec:S=46:ND=17],((Loss|980-1054[&&NHX:D=?:Ev=GLos:S=19:ND=38],Loss|981-1053[&&NHX:D=?:Ev=GLos:S=18:ND=39])[&&NHX:D=F:Ev=Spec:S=20:ND=37],Loss|983-1058[&&NHX:D=?:Ev=GLos:S=11:ND=40])[&&NHX:D=F:Ev=Spec:S=21:ND=36])[&&NHX:D=F:Ev=Spec:S=47:ND=16])[&&NHX:D=F:Ev=Spec:S=75:ND=0];

         Tree 2:
         ((((((ENSMUSP00000107144-ENSMUSP00000099669|Mus_musculus[&&NHX:D=?:Ev=Extant:S=12:ND=47],Loss|1007-1049[&&NHX:D=?:Ev=GLos:S=13:ND=48])[&&NHX:D=F:Ev=Spec:S=14:ND=46],ADJLoss|1009-1051[&&NHX:D=?:Ev=ALos:S=15:ND=49])[&&NHX:D=F:Ev=Spec:S=16:ND=45],ENSSTOP00000017187-ENSSTOP00000023409|Ictidomys_tridecemlineatus[&&NHX:D=?:Ev=Extant:S=17:ND=50])[&&NHX:D=F:Ev=Spec:S=18:ND=44])[&&NHX:D=T:Ev=GDup:S=18:ND=43],ADJLoss|1014-1054[&&NHX:D=?:Ev=ALos:S=19:ND=51])[&&NHX:D=F:Ev=Spec:S=20:ND=42],(ENSOCUP00000004565-ENSOCUP00000007021|Oryctolagus_cuniculus[&&NHX:D=?:Ev=Extant:S=9:ND=53],ADJLoss|1017-1057[&&NHX:D=?:Ev=ALos:S=10:ND=54])[&&NHX:D=F:Ev=Spec:S=11:ND=52])[&&NHX:D=F:Ev=Spec:S=21:ND=41];

   - DIR_NAME/NAME_OUTPUT_duplicated_gene_pairs: the genes found by DeCo/ART-DeCo to have been
     duplicated at the same time (i.e. an Adjacency duplication). It's a 3 columns file.
      Structure:
         Gene1_IDnb = GeneTree1_IDnb|Node_nb_Gene1
         Gene2_IDnb = GeneTree2_IDnb|Node_nb_Gene2
         Species_IDnb	Gene1_IDnb	Gene2_IDnb
      Example:
         76	0|835	0|1110
         57	0|767	0|1087
         14	0|1200	0|1209
         14	0|1221	0|1230
         14	0|1182	0|1251

   - DIR_NAME/NAME_OUTPUT_extant_adjacencies_file: same file as adjacencies INPUT file but with
     new adjacencies predicted by ARt-DeCo.
      Structure:
         Gene1_ID	Gene2_ID
      Example:
         ENSMODP00000016865	ENSMODP00000008552
         ENSMODP00000008552	ENSMODP00000008623
         ENSMODP00000008623	ENSMODP00000015313

   - DIR_NAME/NAME_OUTPUT_genes: the correspondence between ancestral genes names defined by
     DeCo/ARt-DeCo and extant genes. It's a multi-column file: species gene_code followed by
     all its descendant extant genes
      Structure:
         Gene_IDnb = GeneTree_IDnb|Node_nb_Gene
         Species_IDnb	Gene_IDnb	list(EXTANT_Gene_ID)
      Example:
         18	0|10	ENSMUSP00000097372	ENSRNOP00000064162	ENSMUSP00000088202
         16	0|8	ENSMUSP00000097372	ENSRNOP00000064162	ENSMUSP00000088202
         14	0|2	ENSMUSP00000097372	ENSRNOP00000064162

   - DIR_NAME/NAME_OUTPUT_new_adjacencies: the list of new extant adjacencies predicted by
     ARt-DeCo. It's a 4 columns file with the ID of the 2 genes involved in the new adjacency
     followed by the number of adjacencies of the two genes in the INPUT files.
      Structure:
         EXTANT_Gene1_ID	EXTANT_Gene2_ID	Adj_nb_Gene1	Adj_nb_Gene2
      Example:
         ENSPCAP00000016046	ENSPCAP00000016116	0	0
         ENSVPAP00000004291	ENSVPAP00000006001	1	0
         ENSAMEP00000002369	ENSAMEP00000002427	0	1
         ENSLAFP00000020729	ENSLAFP00000022796	1	0

   - DIR_NAME/NAME_OUTPUT_species: the correspondence between ancestral species names defined by
     DeCo/ARt-DeCo and extant species. It's a multi-column file: species followed by all
     its descendant extant species.
      Structure:
         Species_IDnb	list(EXTANT_Species_name)
      Example:
         4	Procavia_capensis	Loxodonta_africana	Echinops_telfairi
         2	Procavia_capensis	Loxodonta_africana
         0	Procavia_capensis

   - DIR_NAME/NAME_OUTPUT_species_stats: Statistics summary file on gene content of extant species
     analyze by ARt-DeCo.
      Structure:
         Species_name	#G_tot	#G_0A	#G_1A	#G_2A	#G_>2A	#A_tot	#Chr	#Contig	P(Adj)	c0(0,0)	c0(1,0/0,1)	c0(1,1)	c1(0,0)	c1(1,0/0,1)	c1(1,1)	Base_log
      Example:
         Species_name	#G_tot	#G_0A	#G_1A	#G_2A	#G_>2A	#A_tot	#Chr	#Contig	P(Adj)	c0(0,0)	c0(1,0/0,1)	c0(1,1)	c1(0,0)	c1(1,0/0,1)	c1(1,1)	Base_log
         Procavia_capensis	15793	7680	5556	2557	0	5335	28	10458	4.76868e-05	1.90772e-05	9.53817e-06	4.76897e-06	0.856487	0.925804	0.995121	22018.5
         Loxodonta_africana	19896	228	708	18960	0	19314	29	582	0.000817704	0.000457756	0.000228691	0.000114298	0.799593	0.896442	0.99329	1283.1
         Echinops_telfairi	16098	9835	4922	1341	0	3802	21	12296	4.05975e-05	1.59836e-05	7.99149e-06	3.99566e-06	0.85876	0.926979	0.995198	25863.6

   - DIR_NAME/NAME_OUTPUT_stats_DECO: Intermediate file to bring statistics from Step3(_proba) to
     Step4.

   - DIR_NAME/NAME_OUTPUT_stats_human_readable: Summary statistics on GENE, SPECIES, ADJACENCIES,
     ADJACENCIES CLASSES and ADJACENCIES TREES on human readable format.

   - DIR_NAME/NAME_OUTPUT_stats_machine: Summary statistics on GENE, SPECIES, ADJACENCIES,
     ADJACENCIES CLASSES and ADJACENCIES TREES on machine format.

   - DIR_NAME/NAME_reconcil: contains the reconciled trees generated by DeCo/ARt-DeCo from your
     gene trees file on Newick or NHX format (Depending of your choice for "OUTPUT_FORMAT" in
     config file)
