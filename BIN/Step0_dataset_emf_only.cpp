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

File: Step0_dataset_emf_only.cpp           Last modified on: 29/09/2015
Created by: Yoann Anselmetti               Created on: 29/10/2014
--------------------------------------------------------------------------
Specification:
Create Gene_Specie association and adjacences from Ensembl Database.
Using EMF files (files containing Gene trees in Ensembl Database).
=========================================================================*/
#include "Step0_dataset.h"

set<string> extantGenes(vector<Tree*> Arbres){
   vector<Tree *>::iterator it;
   set<string> genes_GT;
   string ID_gene;
   for(it=Arbres.begin();it!=Arbres.end();it++){
      vector<string> noms_feuilles;
      vector<string>::iterator it_feuilles;
      noms_feuilles=(*it)->getLeavesNames();
      for(it_feuilles=noms_feuilles.begin();it_feuilles!=noms_feuilles.end();it_feuilles++){
         ID_gene=(*it_feuilles).substr(0,(*it_feuilles).rfind(sep,(*it_feuilles).length()-1));
         if(ID_gene!=NomPer){
//            cout<<ID_gene<<endl;
            genes_GT.insert(ID_gene);
         }
      }
   }
   cout<<genes_GT.size()<<" extant genes in Gene trees."<<endl<<endl;
   return genes_GT;   
}

string getcwd_string(void){
   char buff[PATH_MAX];
   if(getcwd(buff,PATH_MAX)==0){
      cout<<"ERROR while quering current directory!!!"<<endl;
      exit(EXIT_FAILURE);
   }
   return string(buff);
}

int main(int argc, char* argv[]){
   time_t tbegin=time(NULL);              // get the current calendar time
   cout<<"\n\t#########################"<<endl;
   cout<<"\t###  Start Step0_emf  ###"<<endl;
   cout<<"\t#########################"<<endl<<endl;

   lireFicConfig(fic_arbre,fic_gene,fic_especes,fic_adjacence,exp_name,directory,argc,argv);

   cout<<"File(s) used in Step0:"<<endl;
   cout<<"\t- Selected species file: "<<file_spe<<endl;
   cout<<"\t- INPUT Species Tree file: "<<ST_input<<endl;
   cout<<"\t- OUTPUT Species Tree file: "<<fic_especes<<endl;
   cout<<"\t- Gene Trees file: "<<fic_arbre<<endl;
   cout<<"\t- Adjacencies file: "<<fic_adjacence<<endl;
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
            cout<<"\nFrom Step0_dataset_emf_only.cpp: ERROR while opening file "<<log_file<<endl<<endl;
            return(1);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file open by Step0_dataset_emf_only.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      else{
         Offile_log.open(log_file.c_str(), ios::out);
         if (!Offile_log){
            cout<<"\nFrom Step0_dataset_emf_only.cpp: ERROR while opening file "<<log_file<<endl<<endl;
            return(1);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file created by Step0_dataset_emf_only.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      Offile_log<<"\t#########################"<<endl;
      Offile_log<<"\t###  Start Step0_emf  ###"<<endl;
      Offile_log<<"\t#########################"<<endl<<endl;
      Offile_log<<"\tFile(s) used in Step0:"<<endl;
      Offile_log<<"\t\t- Selected species file: "<<file_spe<<endl;
      Offile_log<<"\t\t- INPUT Species Tree file: "<<ST_input<<endl;
      Offile_log<<"\t\t- OUTPUT Species Tree file: "<<fic_especes<<endl;
      Offile_log<<"\t\t- Gene Trees file: "<<fic_arbre<<endl;
      Offile_log<<"\t\t- Adjacencies file: "<<fic_adjacence<<endl;
      Offile_log<<"\t\t- Species-Gene file: "<<fic_gene<<endl<<endl;
   }


   // Variables
//   TreeTemplate<Node> S;
//   S = new TreeTemplate<Node>;
   map<string,string> gene_species_EXT;
   vector<Tree*> Arbres;


   //Récupération du chemin du répertoire BIN/ dans lequel sont situés les exécutables
   string dir_bin=getcwd_string();
//   string dir_bin=get_current_dir_name();
   const char * dirBIN=dir_bin.c_str();

   //Création du dossier EMBL/ pour stocker les fichiers EMBL
   cout<<"Creation of directory "<<dir_embl<<"..."<<flush;
   if(file_log)
      Offile_log<<"FI\tFrom Step0_dataset_emf_only.cpp: Creation of directory "<<dir_embl<<"..."<<flush;
   string commande="mkdir -p "+dir_embl;
   int ret=system(commande.c_str());
   if(ret<0){
      cout<<"\nFrom Step0_dataset_emf_only.cpp: ERROR while creating directory "<<dir_embl<<endl;
      if(file_log)
         Offile_log<<"FI\tFrom Step0_dataset_emf_only.cpp: ERROR while creating directory "<<dir_embl<<endl;
      return(1);
   }
   cout<<" DONE"<<endl<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;

   //Récupération du chemin jusqu'au dossier EMBL/
   int cd=chdir(dir_embl.c_str());
   if(cd<0){
      cout<<"ERROR while changing directory "<<dir_embl<<endl;
      return(1);
   }
   dir_embl=getcwd_string();
//   dir_embl=get_current_dir_name();

   //Replacement dans le répertoire BIN/
   cd=chdir(dirBIN);
   if(cd<0){
      cout<<"ERROR while changing directory "<<dirBIN<<endl;
      return(1);
   }

   //Création du répertoire EMF/ pour stocker fichier EMF
   cout<<"Creation of directory "<<dir_embl<<"... ";
   commande="mkdir -p "+dir_emf;
   ret=system(commande.c_str());
   if(ret<0)
   {
      cout<<"ERROR while creating directory "<<dir_emf<<endl;
      return(1);
   }
   cout<<"DONE"<<endl<<endl;
   //Récupération du chemin jusqu'au dossier EMF/
   cd=chdir(dir_emf.c_str());
   if(cd<0){
      cout<<"ERROR while changing directory "<<dir_emf<<endl;
      return(1);
   }
   dir_emf=getcwd_string();
//   dir_emf=get_current_dir_name();
   const char * dirEMF = dir_emf.c_str();

   //Replacement dans le répertoire BIN/
   cd=chdir(dirBIN);
   if(cd<0){
      cout<<"ERROR while changing directory "<<dirBIN<<endl;
      return(1);
   }


   //Création du vecteur d'espèces à prendre en entrée
   ifstream IffilGenome;
   IffilGenome.open(file_spe.c_str(), ios::in);
   cout<<"Reading file "<<file_spe<<endl;
   if (!IffilGenome){
      cout<<endl<<"ERROR while opening file "<<file_spe<<endl;
      return(1);
   }
   //Création de l'arbre étoile/buffer à partir duquel on va épurer l'arbre des espèces original
   ofstream OffilBufferTree;
   string ST_buffer="buffer_tree";
   OffilBufferTree.open(ST_buffer.c_str(), ios::out | ios::app);
   OffilBufferTree<<"(";
   cout<<"Species list: "<<endl;
   vector<string> choosen_species;
   int j=0;
   while (!IffilGenome.eof()){
      string genome;
      int nb_chr;
      IffilGenome>>genome;
      IffilGenome>>nb_chr;
      if(genome!=""){
         if (j==0)
            OffilBufferTree<<genome;
         else
            OffilBufferTree<<","<<genome;
      }
      j=1;
      choosen_species.push_back(genome);
      if(genome!="")
         cout<<"\t"<<genome<<"\t"<<nb_chr<<endl;
   }
   OffilBufferTree<<");";
   OffilBufferTree.close();
   IffilGenome.close();


   ///////////////////////////
   /// Formulaire de choix ///
   ///////////////////////////

   string data_ST;
   string nb_release;
   string data_trees;
   string gene_tree="";
   string annot;
   string species_gene_adj;
   string clean;
   string emf_file;

   cout<<endl<<"Do you want to generate your Species tree with list of choosen species and the referent Species tree (species_tree.nwk)? (Y/y (Yes) or N/n (No))"<<endl;
   cin>>data_ST;
   cout<<endl;
   if(data_ST=="N"||data_ST=="n"){
      cout<<"=> You have ever your Species tree."<<endl<<endl;
   }
   else if(data_ST=="Y"||data_ST=="y"){
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////   Création du fichier d'arbre d'espèces au format Newick à partir de la liste d'espèces et de l'arbre d'espèces de référence   ///////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      //Épurage de l'arbre des espèces avec l'arbre tampon
      cout<<"Creation Species tree file to the Newick format... "<<flush;
      string commande="./script/restrict2commonsDeLuxe.pl "+ST_input+" "+ST_buffer+"|head -1 >"+fic_especes;
      ret=system(commande.c_str());
      if(ret<0){
         cout<<"ERROR while refining Species tree "<<ST_input<<endl;
         exit(EXIT_FAILURE);
      }
      if(file_exists(ST_buffer.c_str())){
         ret=remove(ST_buffer.c_str());
         if(ret!=0){
            cout<<"ERROR while removing file "<<ST_buffer<<endl;
            exit(EXIT_FAILURE);
         }
      }
      cout<<"DONE"<<endl<<endl;
   }
   else{
      cout<<"ERROR you have to choose between Y or y (for Yes) and N or n (for No)"<<endl;
      return(1);
   }



   ///////////////////////////////////////////////////////////////////////////
   //////   Téléchargement de données de la Base de Données d'Ensembl   //////
   ///////////////////////////////////////////////////////////////////////////

   cout<<"Do you want to download EMF file (Gene trees) from Ensembl database?"<<endl;
   cin>>data_trees;
   cout<<endl;
   if(data_trees=="N"||data_trees=="n"){
      cout<<"=> You have ever EMF file. What is the name of your file (Necessary for the next step of the program and it has to be in the EMF/ directory:"<<dir_emf<<")"<<endl;
      cin>>emf_file;
      cout<<endl;
   }
   else if(data_trees=="Y"||data_trees=="y"){
      cout<<"Which Ensembl release number do you want?"<<endl;
      cin>>nb_release;
      cout<<endl;
      cout<<"Which Gene trees do you want to download? (prot/p (Protein trees) | rna/r (ncRNA trees) | all/a (Protein & ncRNA trees) | none/n (no Gene trees))"<<endl;
      cin>>gene_tree;
      cout<<endl;
      /////////////////////////////////////////////////////////////////////////
      ///////   Récupération des arbres de gènes via le FTP d'ensembl   ///////
      /////////////////////////////////////////////////////////////////////////

      //Récupération du fichier EMF de la Base de données d'Ensembl
      if(gene_tree!="none" && gene_tree!="n" && gene_tree!=""){
         //Suppression de l'ancien fichiers contenant les arbres de gènes
         if(file_exists(fic_arbre.c_str())){
            ret=remove(fic_arbre.c_str());
            if(ret!=0){
               cout<<"ERROR while removing file "<<fic_arbre<<endl;
               exit(EXIT_FAILURE);
            }
         }

         //Création du fichier d'arbres de gènes:
         cout<<"Opening file for Gene trees: "<<fic_arbre<<"..."<<flush;
         ofstream OffilGeneTree;
         OffilGeneTree.open(fic_arbre.c_str(), ios::out | ios::app);
         cout<<" DONE"<<endl<<endl;   
      
         //Placement dans le répertoire EMF/
         cd=chdir(dirEMF);
         if(cd<0){
            cout<<"ERROR while changing directory "<<dir_emf<<endl;
            return(1);
         }
         
         if(INPUT_FORMAT==1){
            // Si on récupère uniquement arbres protéiques
            if(gene_tree=="prot" || gene_tree=="p"){
               commande="rm *.emf; wget -np ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".protein.nhx.emf.gz -O proteins_trees.nhx.emf.gz; gunzip proteins_trees.nhx.emf.gz; chmod 755 proteins_trees.nhx.emf;";
               ret=system(commande.c_str());
               if(ret<0){
                  cout<<"ERROR while downloading file ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".protein.nhx.emf.gz"<<endl;
                  return(1);
               }
               emf_file="proteins_trees.nhx.emf";
               ifstream IffilEmf;
               IffilEmf.open(emf_file.c_str(), ios::in);
               if(!IffilEmf){
                  cout<<endl<<"ERROR while opening file "<<emf_file<<endl;
                  return(1);
               }
               cout<<"Creation of Gene trees file: "<<fic_arbre<<"... ";
               while(!IffilEmf.eof()){
                  string buffer;
                  IffilEmf>>buffer;
                  if(buffer=="DATA"){
                     IffilEmf>>buffer;
                     OffilGeneTree<<buffer<<endl;
                  }
               }
               OffilGeneTree.close();
               IffilEmf.close();
               cout<<"DONE"<<endl<<endl;
            }
            // Si on récupère uniquement arbres ncRNA
            else if(gene_tree=="rna" || gene_tree=="r"){
               commande="rm *.emf; wget -np ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".ncrna.nhx.emf.gz -O ncrna_trees.nhx.emf.gz; gunzip ncrna_trees.nhx.emf.gz;";
               ret=system(commande.c_str());
               if(ret<0){
                  cout<<"ERROR while downloading file ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".ncrna.nhx.emf.gz"<<endl;
                  return(1);
               }
               emf_file="ncrna_trees.nhx.emf";
               ifstream IffilEmf;
               IffilEmf.open(emf_file.c_str(), ios::in);
               if (!IffilEmf){
                  cout<<endl<<"ERROR while opening file "<<emf_file<<endl;
                  return(1);
               }
               cout<<"Creation of Gene trees file: "<<fic_arbre<<"... ";
               while (!IffilEmf.eof()){
                  string buffer;
                  IffilEmf>>buffer;
                  if(buffer=="DATA"){
                     IffilEmf>>buffer;
                     OffilGeneTree<<buffer<<endl;
                  }
               }
               OffilGeneTree.close();
               IffilEmf.close();
               cout<<"DONE"<<endl<<endl;
            }
            // Si on récupère arbres protéiques et arbres ncRNA
            else if(gene_tree=="all" || gene_tree=="a"){
               commande="rm *.emf; wget -np ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".protein.nhx.emf.gz -O proteins_trees.nhx.emf.gz; gunzip proteins_trees.nhx.emf.gz; wget -np ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".ncrna.nhx.emf.gz -O ncrna_trees.nhx.emf.gz; gunzip ncrna_trees.nhx.emf.gz;";
               ret=system(commande.c_str());
               if(ret<0){
                  cout<<"ERROR while downloading files ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".protein.nhx.emf.gz & ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".ncrna.nhx.emf.gz"<<endl;
                  return(1);
               }
               //Données sur les arbes de gènes protéiques
               string prot_file="proteins_trees.nhx.emf";
               ifstream IffilEmf;
               IffilEmf.open(prot_file.c_str(), ios::in);
               if(!IffilEmf){
                  cout<<endl<<"ERROR while opening file "<<prot_file<<endl;
                  return(1);
               }
               cout<<"Creation of Gene trees file: "<<fic_arbre<<"... ";
               while(!IffilEmf.eof()){
                  string buffer;
                  IffilEmf>>buffer;
                  if(buffer=="DATA"){
                     IffilEmf>>buffer;
                     OffilGeneTree<<buffer<<endl;
                  }
               }
               IffilEmf.close();

               //Données sur les arbes de gènes ncRNA
               string ncrna_file="ncrna_trees.nhx.emf";
               IffilEmf.open(ncrna_file.c_str(), ios::in);
               cout<<"File "<<ncrna_file<<" is open!"<<endl<<endl;
               if(!IffilEmf){
                  cout<<endl<<"ERROR while opening file "<<ncrna_file<<endl;
                  return(1);
               }
               while(!IffilEmf.eof()){
                  string buffer;
                  IffilEmf>>buffer;
                  if(buffer=="DATA"){
                     IffilEmf>>buffer;
                     OffilGeneTree<<buffer<<endl;
                  }
               }
               OffilGeneTree.close();
               IffilEmf.close();
               cout<<"DONE"<<endl<<endl;
               emf_file="all_trees.nhx.emf";
               commande="cat "+prot_file+" "+ncrna_file+">"+emf_file;
               ret=system(commande.c_str());
               if(ret<0){
                  cout<<"ERROR while concatening "<<prot_file<<" & "<<ncrna_file<<" to obtain "<<emf_file<<endl;
                  return(1);
               }
            }
         }
         if(INPUT_FORMAT==0){
            // Si on récupère uniquement arbres protéiques
            if(gene_tree=="prot" || gene_tree=="p"){
               commande="rm *.emf; wget -np ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".protein.nh.emf.gz -O proteins_trees.nh.emf.gz; gunzip proteins_trees.nh.emf.gz; chmod 755 proteins_trees.nh.emf;";
               ret=system(commande.c_str());
               if(ret<0){
                  cout<<"ERROR while downloading file ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".protein.nh.emf.gz"<<endl;
                  return(1);
               }
               emf_file="proteins_trees.nh.emf";
               ifstream IffilEmf;
               IffilEmf.open(emf_file.c_str(), ios::in);
               if(!IffilEmf){
                  cout<<endl<<"ERROR while opening file "<<emf_file<<endl;
                  return(1);
               }
               cout<<"Creation of Gene trees file: "<<fic_arbre<<"... ";
               while(!IffilEmf.eof()){
                  string buffer;
                  IffilEmf>>buffer;
                  if(buffer=="DATA"){
                     IffilEmf>>buffer;
                     OffilGeneTree<<buffer<<endl;
                  }
               }
               OffilGeneTree.close();
               IffilEmf.close();
               cout<<"DONE"<<endl<<endl;
            }
            // Si on récupère uniquement arbres ncRNA
            else if(gene_tree=="rna" || gene_tree=="r"){
               commande="rm *.emf; wget -np ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".ncrna.nh.emf.gz -O ncrna_trees.nh.emf.gz; gunzip ncrna_trees.nh.emf.gz;";
               ret=system(commande.c_str());
               if(ret<0){
                  cout<<"ERROR while downloading file ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".ncrna.nh.emf.gz"<<endl;
                  return(1);
               }
               emf_file="ncrna_trees.nh.emf";
               ifstream IffilEmf;
               IffilEmf.open(emf_file.c_str(), ios::in);
               if (!IffilEmf){
                  cout<<endl<<"ERROR while opening file "<<emf_file<<endl;
                  return(1);
               }
               cout<<"Creation of Gene trees file: "<<fic_arbre<<"... ";
               while (!IffilEmf.eof()){
                  string buffer;
                  IffilEmf>>buffer;
                  if(buffer=="DATA"){
                     IffilEmf>>buffer;
                     OffilGeneTree<<buffer<<endl;
                  }
               }
               OffilGeneTree.close();
               IffilEmf.close();
               cout<<"DONE"<<endl<<endl;
            }
            // Si on récupère arbres protéiques et arbres ncRNA
            else if(gene_tree=="all" || gene_tree=="a"){
               commande="rm *.emf; wget -np ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".protein.nh.emf.gz -O proteins_trees.nh.emf.gz; gunzip proteins_trees.nh.emf.gz; wget -np ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".ncrna.nh.emf.gz -O ncrna_trees.nh.emf.gz; gunzip ncrna_trees.nh.emf.gz;";
               ret=system(commande.c_str());
               if(ret<0){
                  cout<<"ERROR while downloading files ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".protein.nh.emf.gz & ftp://ftp.ensembl.org/pub/release-"+nb_release+"/emf/ensembl-compara/homologies/Compara."+nb_release+".ncrna.nh.emf.gz"<<endl;
                  return(1);
               }
               //Données sur les arbes de gènes protéiques
               string prot_file="proteins_trees.nh.emf";
               ifstream IffilEmf;
               IffilEmf.open(prot_file.c_str(), ios::in);
               if(!IffilEmf){
                  cout<<endl<<"ERROR while opening file "<<prot_file<<endl;
                  return(1);
               }
               cout<<"Creation of Gene trees file: "<<fic_arbre<<"... ";
               while(!IffilEmf.eof()){
                  string buffer;
                  IffilEmf>>buffer;
                  if(buffer=="DATA"){
                     IffilEmf>>buffer;
                     OffilGeneTree<<buffer<<endl;
                  }
               }
               IffilEmf.close();

               //Données sur les arbes de gènes ncRNA
               string ncrna_file="ncrna_trees.nh.emf";
               IffilEmf.open(ncrna_file.c_str(), ios::in);
               cout<<"File "<<ncrna_file<<" is open!"<<endl<<endl;
               if(!IffilEmf){
                  cout<<endl<<"ERROR while opening file "<<ncrna_file<<endl;
                  return(1);
               }
               while(!IffilEmf.eof()){
                  string buffer;
                  IffilEmf>>buffer;
                  if(buffer=="DATA"){
                     IffilEmf>>buffer;
                     OffilGeneTree<<buffer<<endl;
                  }
               }
               OffilGeneTree.close();
               IffilEmf.close();
               cout<<"DONE"<<endl<<endl;
               emf_file="all_trees.nh.emf";
               commande="cat "+prot_file+" "+ncrna_file+">"+emf_file;
               ret=system(commande.c_str());
               if(ret<0){
                  cout<<"ERROR while concatening "<<prot_file<<" & "<<ncrna_file<<" to obtain "<<emf_file<<endl;
                  return(1);
               }
            }
         }
      }
      else if(gene_tree=="none"||gene_tree=="n"||gene_tree==""){
         cout<<"No need to download Gene Trees from Ensembl Database"<<endl;
      }
      else{
         cout<<"ERROR you have to choose between prot / p (for protein trees) | rna / r (for ncRNA trees) | all / a (for protein & ncRNA trees) | none / n (for no gene trees)"<<endl;
         return(1);
      }
   }
   else{
      cout<<"ERROR you have to choose between Y or y (for Yes) and N or n (for No)"<<endl;
      return(1);
   }



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////   CRÉATION FICHIERS ANNOTATION GÉNOMIQUE & ASSOCIATION GÈNE-ESPÈCE ET ADJACENCES + RÉCONCILIATION      /////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   //Replacement dans le répertoire BIN/
   cd=chdir(dirBIN);
   if(cd<0){
      cout<<"ERROR while changing directory "<<dirBIN<<endl;
      return(1);
   }

   string emf_file_path=dir_emf+"/"+emf_file;
   string annotation_file=dir_embl+"/genomes_annotation";
   string buffer_annotation=dir_embl+"/buffer_annotation";
   cout<<"Do you want to create genome annotation for all species? (Y/y (Yes) or N/n (No))"<<endl;
   cin>>annot;
   cout<<endl;
   if(annot=="N"||annot=="n"){
      cout<<"=> You have ever your Genome annotation for all species."<<endl<<endl;
   }
   else if(annot=="Y"||annot=="y"){
      //////////////////////////////////////////////////////////////////////////////////////////////
      ///////   Création d'un fichier d'annotation génomique de toutes les espèces           ///////
      ///////   SANS nettoyage des gènes qui ne sont pas présents dans les arbres de gènes   ///////
      //////////////////////////////////////////////////////////////////////////////////////////////

      //Suppression genomes annotation file
      if(file_exists(annotation_file.c_str())){
         ret=remove(annotation_file.c_str());
         if(ret!=0){
            cout<<"ERROR while removing file "<<annotation_file<<endl;
            exit(EXIT_FAILURE);
         }
      }
      //Suppression buffer annotation
      if(file_exists(buffer_annotation.c_str())){
         ret=remove(buffer_annotation.c_str());
         if(ret!=0){
            cout<<"ERROR while removing file "<<buffer_annotation<<endl;
            exit(EXIT_FAILURE);
         }
      }

      cout<<"Creation of Genome annotation file from EMF file for all species from Ensembl database "<<buffer_annotation<<"... "<<flush;
      commande="grep SEQ "+emf_file_path+" | cut -c5- | sort -k1d,1d -k3d,3d -k4n,4n | awk -F\" \" 'BEGIN {OFS=\" \"} {print $1,$3,$4,$5,$7,$2}' | sed 's/ /\\t/g' | sed 's/^./\\u&/' > "+buffer_annotation;
      ret=system(commande.c_str());
      if(ret<0){
         cout<<"ERROR while creating buffer annotation file "<<buffer_annotation<<endl;
         return(1);
      }
      cout<<"DONE"<<endl<<endl;

      vector<string>::iterator it_spe;
      cout<<"Creation of Genome annotation file from EMF file for selected species "<<annotation_file<<"... DONE"<<endl<<endl;
      for(it_spe=choosen_species.begin();it_spe!=choosen_species.end();it_spe++){
         if((*it_spe)!=""){
            string spe=(*it_spe);
            cout<<spe<<endl;
            commande="grep "+spe+" "+buffer_annotation+">>"+annotation_file;
            ret=system(commande.c_str());
            if(ret<0){
               cout<<"ERROR while creating Genome annotation file "<<annotation_file<<endl;
               return(1);
            }
         }
         else
            break;
      }
   }
   else{
      cout<<"ERROR you have to choose between Y or y (for Yes) and N or n (for No)"<<endl;
      return(1);
   }
   choosen_species.clear();


   string buffer_annot=dir_embl+"/buffer_annot";
   cout<<endl<<"Do you want to generate species-gene association & adjacencies files from the genome annotation of all species? (Y/y (Yes) or N/n (No))"<<endl;
   cin>>species_gene_adj;
   cout<<endl;
   if(species_gene_adj=="N"||species_gene_adj=="n"){
      cout<<"=> You have ever your species-gene association & adjacencies files."<<endl<<endl;
   }
   else if(species_gene_adj=="Y"||species_gene_adj=="y"){
      /////////////////////////////////////////////////////////////////////////////////////////////////
      ///////   Création des fichiers d'entrée de DeCo (association espèce-gène & adjacences)   ///////
      /////////////////////////////////////////////////////////////////////////////////////////////////
      
      //Tri du fichier d'annotation génomique
      cout<<"Sorting Genome annotation file... "<<flush;
      commande="sort -k1d,1d -k2d,2d -k3n,3n "+annotation_file+" > "+buffer_annot;
      ret=system(commande.c_str());
      if(ret<0){
         cout<<"ERROR while sorting "<<annotation_file<<endl;
         return(1);
      }
      int rc = rename(buffer_annot.c_str(), annotation_file.c_str()); 
      if(rc){
         cout<<"Error while renaming "<<buffer_annot<<" in "<<annotation_file<<endl;
         exit(EXIT_FAILURE);
      }
      cout<<"DONE"<<endl<<endl;


      //Suppression ancien fichier d'adjacences:
      if(file_exists(fic_adjacence.c_str())){
         ret=remove(fic_adjacence.c_str());
         if(ret!=0){
            cout<<"ERROR while removing file "<<fic_adjacence<<endl;
            exit(EXIT_FAILURE);
         }
      }
      //Création du nouveau fichier d'adjacences:
      ofstream OffilAdj;
      OffilAdj.open(fic_adjacence.c_str(), ios::out | ios::app); 

      //Suppression ancien fichier d'association gène-espèce:
      if(file_exists(fic_gene.c_str())){
         ret=remove(fic_gene.c_str());
         if(ret!=0){
            cout<<"ERROR while removing file "<<fic_gene<<endl;
            exit(EXIT_FAILURE);
         }
      }
      //Création du nouveau fichier d'association gène-espèce:
      ofstream OffilGeneSpecie;
      OffilGeneSpecie.open(fic_gene.c_str(), ios::out | ios::app); 


      //Parcours de tous les fichiers d'annotations des génomes des espèces que l'on étudie pour créer le fichier Gène-Espèce d'entrée de DeCo
      cout<<"Creation of Gene-Species and Adjacencies files... "<<flush;
      ifstream Iffile (annotation_file.c_str(),ios::in);
      if(!Iffile){
            cout<<endl<<"ERROR while opening file "<<annotation_file<<endl<<endl;
            return(1);
      }
      else{
         string species="";
         string chr="";
         int startPos=0;
         int endPos=0;
   //      string direction="";
         string gene="";
         string rna_prot="";
         string bufferline="";
         string buffer="";
         string line;
         int bufferInt;
         bool overlap=false;
         while(getline(Iffile, line)){
            istringstream bufferline(line);
            bufferline>>buffer;
            if(!Iffile.eof()){
               if(buffer==species){
                  bufferline>>buffer;
                  if(buffer==chr){
                     OffilAdj<<rna_prot<<"\t";
                     bufferline>>bufferInt;
                     if(bufferInt<startPos){
                        cout<<"ERROR the file is not sorted in genome order position!!! "<<bufferInt<<"<"<<startPos<<endl;
   //                     return(1);
                        if(bufferInt<endPos){
   //                        cout<<"Overlapping of gene "<<gene<<" ("<<startPos<<"-"<<endPos<<") and gene ";
                           overlap=true;
                        }
                        startPos=bufferInt;
                        bufferline>>endPos;
                        bufferline>>gene;
                        if(overlap){
   //                        cout<<gene<<" ("<<startPos<<"-"<<endPos<<")"<<endl;
                           overlap=false;
                        }
                        bufferline>>rna_prot;
                        gene_species_EXT[rna_prot]=species;
                        OffilGeneSpecie<<species<<"\t"<<rna_prot<<endl;
                     }
                     else{
                        if(bufferInt<endPos){
   //                        cout<<"Overlapping of gene "<<gene<<" ("<<startPos<<"-"<<endPos<<") and gene ";
                           overlap=true;
                        }
                        startPos=bufferInt;
                        bufferline>>endPos;
                        bufferline>>gene;
                        if(overlap){
   //                        cout<<gene<<" ("<<startPos<<"-"<<endPos<<")"<<endl;
                           overlap=false;
                        }
                        bufferline>>rna_prot;
                        gene_species_EXT[rna_prot]=species;
                        OffilGeneSpecie<<species<<"\t"<<rna_prot<<endl;
                        OffilAdj<<rna_prot<<endl;
                     }
                  }
                  else{
                     chr=buffer;
                     bufferline>>startPos;
                     bufferline>>endPos;
                     bufferline>>gene;
                     bufferline>>rna_prot;
                     gene_species_EXT[rna_prot]=species;
                     OffilGeneSpecie<<species<<"\t"<<rna_prot<<endl;
                  }
               }
               else{
                  species=buffer;
                  bufferline>>chr;
                  bufferline>>startPos;
                  bufferline>>endPos;
                  bufferline>>gene;
                  bufferline>>rna_prot;
                  gene_species_EXT[rna_prot]=species;
                  OffilGeneSpecie<<species<<"\t"<<rna_prot<<endl;
               }
            }
         }
      }
      Iffile.close();
      OffilGeneSpecie.close();
      OffilAdj.close();
      cout<<"DONE"<<endl<<endl;
   }
   else{
      cout<<"ERROR you have to choose between Y or y (for Yes) and N or n (for No)"<<endl;
      return(1);
   }


   
   ////////////////////////////////////////
   ///////   STEP1_RECONCILIATION   ///////
   ////////////////////////////////////////

   // !!!!!   La réconciliation nécessite un  fichier d'association Gène-Espèce pour pouvoir sélectionner les arbres d'intérêts   !!!!!

   cout<<"*** START RECONCILIATION ***"<<endl<<endl;

   //Arbre des espèces
   cout<<"Storage Species tree in TreeTemplate<Node> * S..."<<flush;
   Newick * newickReaderEsp = new Newick(true,true);
   newickReaderEsp->enableExtendedBootstrapProperty(esp);
   TreeTemplate<Node> * S = newickReaderEsp->read(fic_especes);
   cout<<" DONE"<<endl<<endl;

   //Les noeuds internes sont lus comme des propriétés de branches
   //La f° ci-dessous permet de nommer les noeud internes et de faire passer la propriété au noeud
   nommerNoeudsInternes(S->getRootNode());
   
   //Création du fichier de sortie de correspondance des espèces
   if (!fic_SORTIE_EspeceDone){
      CreationOutputSpecies(S);
   }

   //Récupération & affectation d'info sur les branches des arbres de gènes avec l'arbre des espèces != Réconciliation!!!
   StoreGeneTrees(S,Arbres);

   //Réconciliation des arbres de gènes avec l'arbre des espèces
   reconciliation(S,Arbres,gene_species_EXT);
   delete(S);

   cout<<"*** END RECONCILIATION ***"<<endl<<endl;

   cout<<"Do you want to clean extant genes from genome species that are not present in Gene trees? (Avoid gap in adjacencies file) (Y/y (Yes) or N/n (No))"<<endl;
   cin>>clean;
   cout<<endl;
   if(clean=="N"||clean=="n"){
      cout<<"=> You don't want to clean extant Genes that are not present in Gene trees."<<endl<<endl;
   }
   else if(clean=="Y"||clean=="y"){
      //////////////////////////////////////////////////////////////////////////////////////////
      ///////   Création d'un fichier d'annotation génomique de toutes les espèces      ////////
      ///////   et épurage des gènes qui ne sont pas présents dans les arbres de gènes   ///////
      //////////////////////////////////////////////////////////////////////////////////////////

      //Récupération des gènes contenus dans les arbres de gènes dans vector<string> genes_GT
      set<string> genes_GT=extantGenes(Arbres);
      vector<Tree*>::iterator ita;
      for (ita=Arbres.begin();ita!=Arbres.end();ita++)
         delete(*ita);
      Arbres.clear();

      //Nettoyage du fichier d'annotation génomique des gènes qui ne sont pas contenus dans les arbres de gènes
      cout<<"Cleaning Genome annotation file..."<<flush;
      //Ouverture en écriture du fichier buffer_annot
      ofstream OffilBufAnnot;
      OffilBufAnnot.open(buffer_annot.c_str(), ios::out | ios::app);
      //Ouverture en lecture du fichier d'annotation génomique
      ifstream IffilGenomeAnnot(annotation_file.c_str(),ios::in);
      if(!IffilGenomeAnnot){
            cout<<endl<<"ERROR while opening file "<<annotation_file<<endl<<endl;
            return(1);
      }
      else{
         string gene="";
         string buffer="";
         string line;

         while(getline(IffilGenomeAnnot, line)){
            istringstream bufferline(line);
            bufferline>>buffer;
            bufferline>>buffer;
            bufferline>>buffer;
            bufferline>>buffer;
            bufferline>>buffer;
            bufferline>>gene;
            if(genes_GT.find(gene)!=genes_GT.end()){
               OffilBufAnnot<<line<<"\n";
            }
         }
         OffilBufAnnot.close();
         IffilGenomeAnnot.close();
      }

      cout<<"DONE"<<endl<<endl;

      int rc = rename(buffer_annot.c_str(), annotation_file.c_str()); 
      if(rc){
         cout<<"Error while renaming "<<buffer_annot<<" in "<<annotation_file<<endl;
         exit(EXIT_FAILURE);
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////
      ///////   Création des fichiers d'entrée de DeCo (association espèce-gène & adjacences)   ///////
      /////////////////////////////////////////////////////////////////////////////////////////////////

      //Suppression ancien fichier d'adjacences:
      if(file_exists(fic_adjacence.c_str())){
         ret=remove(fic_adjacence.c_str());
         if(ret!=0){
            cout<<"ERROR while removing file "<<fic_adjacence<<endl;
            exit(EXIT_FAILURE);
         }
      }
      //Création du nouveau fichier d'adjacences:
      ofstream OffilAdj;
      OffilAdj.open(fic_adjacence.c_str(), ios::out | ios::app); 

      //Suppression ancien fichier d'association gène-espèce:
      if(file_exists(fic_gene.c_str())){
         ret=remove(fic_gene.c_str());
         if(ret!=0){
            cout<<"ERROR while removing file "<<fic_gene<<endl;
            exit(EXIT_FAILURE);
         }
      }
      //Création du nouveau fichier d'association gène-espèce:
      ofstream OffilGeneSpecie;
      OffilGeneSpecie.open(fic_gene.c_str(), ios::out | ios::app); 

      //Parcours de tous les fichiers d'annotations des génomes des espèces que l'on étudie pour créer le fichier Gène-Espèce d'entrée de DeCo
      cout<<"Creation of Gene-Species and Adjacencies files... "<<flush;
      ifstream Iffile (annotation_file.c_str(),ios::in);         
      if(!Iffile){
            cout<<endl<<"ERROR while opening file "<<annotation_file<<endl<<endl;
            return(1);
      }
      else{
         string species="";
         string chr="";
         int startPos=0;
         int endPos=0;
   //      string direction="";
         string gene="";
         string rna_prot="";
         string bufferline="";
         string buffer="";
         string line;
         int bufferInt;
         bool overlap=false;
         while(getline(Iffile, line)){
            istringstream bufferline(line);
            bufferline>>buffer;
            if(!Iffile.eof()){            
               if(buffer==species){
                  bufferline>>buffer;
                  if(buffer==chr){
                     OffilAdj<<rna_prot<<"\t";
                     bufferline>>bufferInt;
                     if(bufferInt<startPos){
                        cout<<"ERROR the file is not sorted in genome order position!!! "<<bufferInt<<"<"<<startPos<<endl;
   //                     return(1);
                        if(bufferInt<endPos){
   //                        cout<<"Overlapping of gene "<<gene<<" ("<<startPos<<"-"<<endPos<<") and gene ";
                           overlap=true;
                        }
                        startPos=bufferInt;
                        bufferline>>endPos;
                        bufferline>>gene;
                        if(overlap){
   //                        cout<<gene<<" ("<<startPos<<"-"<<endPos<<")"<<endl;
                           overlap=false;
                        }
                        bufferline>>rna_prot;
                        gene_species_EXT[rna_prot]=species;
                        OffilGeneSpecie<<species<<"\t"<<rna_prot<<endl;
                     }
                     else{
                        if(bufferInt<endPos){
   //                        cout<<"Overlapping of gene "<<gene<<" ("<<startPos<<"-"<<endPos<<") and gene ";
                           overlap=true;
                        }
                        startPos=bufferInt;
                        bufferline>>endPos;
                        bufferline>>gene;
                        if(overlap){
   //                        cout<<gene<<" ("<<startPos<<"-"<<endPos<<")"<<endl;
                           overlap=false;
                        }
                        bufferline>>rna_prot;
                        gene_species_EXT[rna_prot]=species;
                        OffilGeneSpecie<<species<<"\t"<<rna_prot<<endl;
                        OffilAdj<<rna_prot<<endl;
                     }
                  }
                  else{
                     chr=buffer;
                     bufferline>>startPos;
                     bufferline>>endPos;
                     bufferline>>gene;
                     bufferline>>rna_prot;
                     gene_species_EXT[rna_prot]=species;
                     OffilGeneSpecie<<species<<"\t"<<rna_prot<<endl;
                  }
               }
               else{
                  species=buffer;
                  bufferline>>chr;
                  bufferline>>startPos;
                  bufferline>>endPos;
                  bufferline>>gene;
                  bufferline>>rna_prot;
                  gene_species_EXT[rna_prot]=species;
                  OffilGeneSpecie<<species<<"\t"<<rna_prot<<endl;
               }
            }
         }
      }
      Iffile.close();
      OffilGeneSpecie.close();
      OffilAdj.close();
      cout<<"DONE"<<endl<<endl;
   }
   else{
      cout<<"ERROR you have to choose between Y or y (for Yes) and N or n (for No)"<<endl;
      return(1);
   }
   gene_species_EXT.clear();

   time_t tend=time(NULL);                // get the current calendar time    
   // Compute execution time
   float texec=difftime(tend,tbegin);    // tend-tbegin (result in second)

   cout<<"Execution time of Step0_emf is: "<<texec<<"s"<<endl<<endl;
   cout<<"\t#######################"<<endl;
   cout<<"\t###  End Step0_emf  ###"<<endl;
   cout<<"\t#######################"<<endl<<endl;

   if (file_log){
      Offile_log<<"Execution time of Step0_emf is: "<<texec<<"s"<<endl<<endl;
      Offile_log<<"\t#######################"<<endl;
      Offile_log<<"\t###  End Step0_emf  ###"<<endl;
      Offile_log<<"\t#######################"<<endl<<endl;
   }
   Offile_log.close();

   return(0);
}
