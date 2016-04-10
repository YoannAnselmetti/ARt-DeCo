/*========================================================================
This file is part of DeCo.

Copyright (c) 2012, UM2, UL1, INRIA

DeCo is free software. This software is governed by the CeCILL license
under French law and abiding by the rules of distribution of free
software

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided
only with a limited warranty and the software's author, the holder of
the economic rights, and the successive licensors have only limited
liability. See the file license_CeCILL_V2-en.txt for more details.

File: Step1_reconciliation.cpp                  Last modified on: 2016/04/08
Created by: Yoann Anselmetti/Sèverine Bérard    Created on: 2014/02/28
--------------------------------------------------------------------------
Specification:
File use to make reconciliation of genes trees with species tree.
If you produce Input files with Step0, you can pass directly to Step2.
=========================================================================*/
#include "Step1_reconciliation.h"

//Les pré-requis

// 1) G, l'arbre de gène à réconcilier, et S, l'arbre des espèces,
// doivent être deux arbres binaires

// 2) Toutes les feuilles de S doivent avoir un nom (considéré comme
// nom d'espèce) ne contenant pas le séparateur utilisé dans G (par
// défaut | stocké dans la var sep). Si les noeud internes sont nommés
// ce nom est récupéré pour la visualisation

// 3) Toutes les feuilles de G doivent avoir un nom composé du nom de
// gène quelconque suivi d'un séparateur (par défaut | stocké dans la
// var sep) suivi d'un nom d'espèce ne contenant pas ce séparateur
// mais n'appartenant pas forcément à S


void RewriteSpeciesGeneFile(map<string,string> &gene_species_EXT){
   
   cout<<"Rewrite Species-Gene association file after Genes clean..."<<flush;
   if (file_log)
      Offile_log<<"FI\tFrom Step1_reconciliation.cpp (RewriteSpeciesGeneFile): Rewrite Species-Gene association file after Genes clean..."<<flush;
   
   ifstream IffilSpeGene(fic_gene.c_str(),ios::in);
   if (!IffilSpeGene){
      cout<<"\nFrom Step1_reconciliation.cpp (RewriteSpeciesGeneFile): ERROR while opening file "<<fic_gene<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step1_reconciliation.cpp (RewriteSpeciesGeneFile): ERROR while opening file "<<fic_gene<<endl;
      exit(EXIT_FAILURE);
   }
   
   //Maintenant on crée un nouveau fichier (dont le nom est déclaré en
   //global) pour ré-écrire les gènes-espèces après nettoyage
   ofstream OffilSpeGene(fic_gene_clean.c_str(), ios::out|ios::trunc);
   if (!OffilSpeGene){
      cout<<"\nFrom Step1_reconciliation.cpp (RewriteSpeciesGeneFile): ERROR while opening file "<<fic_gene_clean<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step1_reconciliation.cpp (RewriteSpeciesGeneFile): ERROR while opening file  "<<fic_gene_clean<<endl;
      exit(EXIT_FAILURE);
   }

   string buffer="";
   string gene="";
   string line;
   while(getline(IffilSpeGene, line)){
      istringstream bufferline(line);
      bufferline>>buffer;
      bufferline>>gene;
      if (gene_species_EXT.find(gene)!=gene_species_EXT.end())
         OffilSpeGene<<line<<endl;
   }
   IffilSpeGene.close();
   OffilSpeGene.close();
   
   cout<<" DONE"<<endl<<endl;
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
}


void RewriteAdjFile(map<string,string> &gene_species_EXT){
   cout<<"Rewrite adjacencies file after Genes clean..."<<flush;
   if (file_log)
      Offile_log<<"FI\tFrom Step1_reconciliation.cpp (RewriteAdjFile): Rewrite adjacencies file after Genes clean..."<<flush;
   
   ifstream IffilAdj (fic_adjacence.c_str(), ios::in);
   if (!IffilAdj){
      cout<<"\nFrom Step1_reconciliation.cpp (RewriteAdjFile): ERROR while opening file "<<fic_adjacence<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step1_reconciliation.cpp (RewriteAdjFile): ERROR while opening file "<<fic_adjacence<<endl;
      exit(EXIT_FAILURE);
   }
   
   vector<adjacence> adjacencies;
   adjacence adj;
   string stored_G1="";
   string stored_G2="";
   string gene1="";
   string gene2="";
   while (!IffilAdj.eof()){
      //on lit la ligne en entier;
      IffilAdj>>gene1;      
      //Pour éviter un dernier tour de trop à cause d'un blanc en fin de fichier
      if (!IffilAdj.eof()){
         IffilAdj>>gene2;
         // If 
         if (gene1==stored_G2){
            if (gene_species_EXT.find(gene2)!=gene_species_EXT.end()){
		         if (gene_species_EXT.find(stored_G1)!=gene_species_EXT.end()){
		            adj.gene1=stored_G1;
		            adj.gene2=gene2;
		            adjacencies.push_back(adj);
						string stored_G1="";
						string stored_G2="";
		         }
				}
            else{
               // If gene2 is not in gene trees (gene_species_EXT), replace current stored_G2 by gene2.
               stored_G2=gene2;
            }
         }
         else{
		      //on vérifie si les deux gènes de l'adjacence sont dans gene_species_EXT (c-à-d dans nos arbres)
		      if (gene_species_EXT.find(gene1)!=gene_species_EXT.end()){
		         if (gene_species_EXT.find(gene2)!=gene_species_EXT.end()){
		            adj.gene1=gene1;
		            adj.gene2=gene2;
		            adjacencies.push_back(adj);
		            //Réinitialisation pour le cas où on change d'espèce ou qu'on arrive au début d'1 contig/chr.
						string stored_G1="";
						string stored_G2="";
		         }
               // IF gene2 is not in gene_species_EXT (i.e in gene trees store). Store gene1 and gene2 in stored_G1 and stored_G2
		         else{
		            stored_G1=gene1;
		            stored_G2=gene2;
		         }
		      }
         }
      }
   }
   IffilAdj.close();

   cout<<"Adjacencies of "<<fic_adjacence<<" are read, now we write them in "<<fic_adjacence_clean<<endl;
   ofstream OffilAdj(fic_adjacence_clean.c_str(), ios::out|ios::trunc);
   if (!OffilAdj){
      cout<<"\nFrom Step1_reconciliation.cpp (RewriteAdjFile): ERROR while opening file "<<fic_adjacence_clean<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step1_reconciliation.cpp (RewriteAdjFile): ERROR while opening file "<<fic_adjacence_clean<<endl;
      exit(EXIT_FAILURE);
   }
   vector<adjacence>::iterator it_adj;
   for(it_adj=adjacencies.begin();it_adj!=adjacencies.end();it_adj++){
      OffilAdj<<(*it_adj).gene1<<"\t"<<(*it_adj).gene2<<endl;
   }
   cout<<"Closing "<<fic_adjacence_clean<<" ... "<<flush;
   OffilAdj.close();
   cout<<"and clearing adjacencies vector ... "<<flush;
   adjacencies.clear();
   cout<<" DONE"<<endl<<endl;
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
}


int main(int argc, char* argv[]){
   time_t tbegin=time(NULL);              // get the current calendar time
   cout<<"\n\t#######################"<<endl;
   cout<<"\t###  Start Step1 !  ###"<<endl;
   cout<<"\t#######################"<<endl<<endl;
   
   //Read config file of DeCo/DeCo++
   lireFicConfig(fic_arbre,fic_gene,fic_especes,fic_adjacence,exp_name,directory,argc,argv);
   
   cout<<endl<<"File(s) used in Step1:"<<endl;
   cout<<"\t- Gene Trees file: "<<fic_arbre<<endl;
   cout<<"\t- Species Tree file: "<<fic_especes<<endl;
   cout<<"\t- Species-Gene file: "<<fic_gene<<endl<<endl;
   
   // Log file creation or opening
   if (file_log){
      //Get time
      time_t rawtime;
      struct tm datetime;
      char str_time[50];
      time(&rawtime);
      datetime = *localtime(&rawtime);
      strftime(str_time, 50, "%A %d %B %Y %H:%M:%S", &datetime);
      
      // Creation and/or opening log file
      if (is_readable(log_file)){
         Offile_log.open(log_file.c_str(), ios::out|ios::app);
         if (!Offile_log){
            cout<<"\nFI\tFrom Step1_reconciliation.cpp: ERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file open by Step1_reconciliation.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      else{
         Offile_log.open(log_file.c_str(), ios::out|ios::trunc);
         if (!Offile_log){
            cout<<"\nFI\tFrom Step1_reconciliation.cpp: ERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file created by Step1_reconciliation.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      Offile_log<<"\t#######################"<<endl;
      Offile_log<<"\t###  Start Step1 !  ###"<<endl;
      Offile_log<<"\t#######################"<<endl<<endl;
      Offile_log<<"\tFile(s) used in Step1:"<<endl;
      Offile_log<<"\t\t- Gene Trees file: "<<fic_arbre<<endl;
      Offile_log<<"\t\t- Species Tree file: "<<fic_especes<<endl;
      Offile_log<<"\t\t- Species-Gene file: "<<fic_gene<<endl<<endl;
   }
   
   // Variables
   //   TreeTemplate<Node> *S;
   //   S = new TreeTemplate<Node>;
   map<string,string> gene_species_EXT;
   vector<Tree*> Arbres;
   map<string,vector<int> > gene_GeneTreeID;
   
   
   //   //Store species tree in TreeTemplate<Node> * S
   //   StoreSpeciesTree(S);
   
   cout<<"Store Species tree in TreeTemplate<Node> * S..."<<flush;
   if(file_log)
      Offile_log<<"OS\tFrom Step1_reconciliation.cpp (main): Store Species tree in TreeTemplate<Node> * S..."<<flush;
   Newick * newickReaderEsp = new Newick(true,true);
   newickReaderEsp->enableExtendedBootstrapProperty(esp);
   TreeTemplate<Node> * S = newickReaderEsp->read(fic_especes);
   delete(newickReaderEsp);
   cout<<" DONE"<<endl<<endl;   
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;
   
   //Les noeuds internes sont lus comme des propriétés de branches
   //La f° ci-dessous permet de nommer les noeud internes et de faire passer la propriété au noeud
   if (file_log)
      Offile_log<<"OS\tFrom Step1_reconciliation.cpp (main:nommerNoeudsInternes): Affect species name to node in species tree ..."<<flush;
   nommerNoeudsInternes(S->getRootNode());
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
   
   //Création du fichier de sortie de correspondance des espèces
   if (!fic_SORTIE_EspeceDone){
      CreationOutputSpecies(S);
   }
   
   //Récupération & affectation d'info sur les branches des arbres de gènes avec l'arbre des espèces != Réconciliation!!!
   StoreGeneTrees(S,Arbres);
   
   //Associate Species_name with gene in map<string, string> gene_species_EXT
   AssociateGeneWithSpecies(gene_species_EXT);

   //Reconciliation of Gene trees with Species tree (Parsimony principle DL Model (Goodman, 1979)
   reconciliation(S,Arbres,gene_species_EXT);

/////////////////////////////////
/////   START Genes Clean   /////
/////////////////////////////////

   //Associe les gènes avec les numéros d'arbres réconciliés dans lesquels ils sont présents
   AssociateGeneWithReconciledTreeNb(Arbres,gene_GeneTreeID,gene_species_EXT);
   //   afficheMap(gene_GeneTreeID);
   
   //On libère la mémoire du vecteur d'arbres
   vector<Tree*>::iterator ita;
   for (ita=Arbres.begin();ita!=Arbres.end();ita++)
      delete(*ita);
   Arbres.clear();
   delete(S);
   
   //Épurage de la map gene_species_EXT pour ne garder que les gènes présents dans nos arbres
   int elim = 0;
   elim=GeneSpeciesClean(gene_GeneTreeID,gene_species_EXT);
   gene_GeneTreeID.clear();
   
///////////////////////////////
/////   END Genes Clean   /////
///////////////////////////////


/////////////////////////////////////////
/////   START Rewrite INPUT Files   /////
/////////////////////////////////////////

   if(elim>0){
	   //Réécriture du fichier d'association Gène-Espèce dans un NOUVEAU fichier et si on a éliminé des gènes
      RewriteSpeciesGeneFile(gene_species_EXT);
      //Store adjacencies in vector<adjacence> adjacencies & Rewrite adjacacencies file after removing of genes not present in gene trees.
      RewriteAdjFile(gene_species_EXT);
   }
   gene_species_EXT.clear();

///////////////////////////////////////
/////   END Rewrite INPUT Files   /////
///////////////////////////////////////


   time_t tend=time(NULL);                // get the current calendar time    
   // Compute execution time
   float texec=difftime(tend,tbegin);    // tend-tbegin (result in second)
   
   cout<<"Execution time of Step1 is: "<<texec<<"s"<<endl<<endl;
   cout<<"\t#####################"<<endl;
   cout<<"\t###  End Step1 !  ###"<<endl;
   cout<<"\t#####################"<<endl<<endl;
   
   if (file_log){
      Offile_log<<"Execution time of Step1 is: "<<texec<<"s"<<endl<<endl;
      Offile_log<<"\t#####################"<<endl;
      Offile_log<<"\t###  End Step1 !  ###"<<endl;
      Offile_log<<"\t#####################"<<endl<<endl;
   }
   Offile_log.close();

   return(0);
} 
