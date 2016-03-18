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

File: Step3_DECO.cpp                              Last modified on: 22/05/2015
Created by: Yoann Anselmetti/Sèverine Bérard      Created on: 19/03/2014
--------------------------------------------------------------------------
Specification: 
File containing the DECO function dedicated to pairwise tree computations
=========================================================================*/
#include "Step3_DECO.h"

//Peut être utilisé sur les deux formats nhx et newick
//+ ajout de noeud à un seul fils pour les arbres d'adjacences !!
void afficheInfoNoeud(Node * n){
   convert<int> C;
   if (n->getNumberOfSons()==0)//(n->isLeaf())
      cout<<"Leaf with id "<<n->getId();
   else if (n->getNumberOfSons()==1)
      cout<<"Internal node with id "<<n->getId()<<"("<<n->getSon(0)->getId()<<")";
   else //c'est qu'il en a deux ! 
      cout<<"Internal node with id"<<n->getId()<<"("<<n->getSon(0)->getId()<<","<<n->getSon(1)->getId()<<")";
   cout<<"\t Name: ";
   if (n->hasName())
      cout<<n->getName();
   else
      cout<<"NONAME";
   cout<<"\t Specie: ";
   if (n->hasNodeProperty(esp)){
      BppString * prop = dynamic_cast<BppString*> (n->getNodeProperty(esp));
      cout<<C.from_s(prop->toSTL());
   }
   else
      cout<<"NOSPECIE";
   cout<<"\t Type: ";
   if (n->hasBranchProperty(typ)){
      BppString * type = dynamic_cast<BppString*> (n->getBranchProperty(typ));
      cout<<type->toSTL();
   }
   else
      cout<<"NOTYPE";
   cout<<"\t Father dist: ";
   if (n->hasDistanceToFather())
      cout<<n->getDistanceToFather();
   else
      cout<<"NODIST";
   cout<<endl;
}

//Peut être utilisé sur les deux formats nhx et newick ... à tester
void afficheInfoNoeuds(Node * n){
   if (n->getNumberOfSons()==0)//(n->isLeaf())
      afficheInfoNoeud(n);
   else if (n->getNumberOfSons()==1){
      afficheInfoNoeuds((n->getSons())[0]);
      afficheInfoNoeud(n);
   }
   else{
      afficheInfoNoeuds((n->getSons())[0]);
      afficheInfoNoeud(n);
      afficheInfoNoeuds((n->getSons())[1]);
   }
}



void ordonnePostfixe(Node * n, vector<Node *> & v){
   if (n->getNumberOfSons()==0) //(n->isLeaf())
      v.push_back(n);
   else{
      ordonnePostfixe(n->getSon(0),v);
      ordonnePostfixe(n->getSon(1),v);
      v.push_back(n);
   }
}

////Pour compter les feuilles en dessus du noeud n
//int compteFeuilles(Node * n){
//   if (n->getNumberOfSons()==0)
//      return 1;
//   else if (n->getNumberOfSons()==1)
//      return compteFeuilles(n->getSon(0));
//   else //2 fils 
//      return compteFeuilles(n->getSon(0))+compteFeuilles(n->getSon(1));
//}

////Écrit l'info sur noeud dans le fichier des adjacences ancestrales sous la forme :
////numEspece G1=num1|id G2=num2|id Nb_adj_filles_feuilles suivant des adj_filles_feuilles
////ATTENTION : que si le noeud est un noeud de spéciation !!
//void ecritAdjAncestrale(Node * n, int num1, int num2, TreeTemplate<Node> * A1, TreeTemplate<Node> * A2, ofstream & OfficAdj2){
//   convert<int> c;
//   //n° espece
//   BppString * Espece = dynamic_cast<BppString*> (n->getNodeProperty(esp));
//   OfficAdj2<<c.from_s(Espece->toSTL())<<"\t";

//   //Les deux gènes adjacents
//   string nom=n->getName();
//   string Id1=nom.substr(0,nom.find(sepAdj));
//   string Id2=nom.substr(nom.find(sepAdj)+1,nom.size());
//   OfficAdj2<<num1<<sep<<Id1<<"\t"<<num2<<sep<<Id2<<"\t";

//   //Leurs descendants dans l'arbre solution
//   int nb_adj_filles_feuilles=compteFeuilles(n);
//   OfficAdj2<<nb_adj_filles_feuilles;
//   vector<string> feuilles;
//   retourneFeuilles(n,feuilles);
//   for(int i=0;i<nb_adj_filles_feuilles;i++)
//      OfficAdj2<<"\t"<<feuilles[i];
//   OfficAdj2<<"\n";
//}

////Gère la récursion de l'écriture des info des noeuds de l'arbre dans le fic
//void ecritAdjAncestrales(Node * n, int num1, int num2, TreeTemplate<Node> * A1, TreeTemplate<Node> * A2, ofstream & OfficAdj2){
////On ne fait rien pour les feuilles
//   if (n->getNumberOfSons()!=0){
//      //Est-ce un noeud de spéciation ?
//      BppString * TYP = dynamic_cast<BppString*> (n->getBranchProperty(typ));
//      if(TYP->toSTL()==spe)
//         if (n->getNumberOfSons()==1){
//            ecritAdjAncestrale(n,num1,num2,A1,A2,OfficAdj2);
//            ecritAdjAncestrales(n->getSon(0),num1,num2,A1,A2,OfficAdj2);
//         }
//         else{ //2 fils
//            ecritAdjAncestrale(n,num1,num2,A1,A2,OfficAdj2);
//            ecritAdjAncestrales(n->getSon(0),num1,num2,A1,A2,OfficAdj2);
//            ecritAdjAncestrales(n->getSon(1),num1,num2,A1,A2,OfficAdj2);
//         }
//   }
//}

void ecritSORTIEDup(vector<adjacence> DUPLI, int num1, TreeTemplate<Node>* Arbre1, int num2, ofstream & OfficSORTIE_dup){
   vector<adjacence>::iterator it_adj;
   for(it_adj=DUPLI.begin();it_adj!=DUPLI.end();it_adj++){
      //n° espece (les deux gènes sont de la même espèce --> on travaille sur le 1er)
      convert<int> c;
      int Id1 = c.from_s((*it_adj).gene1); //c'est un numéro d'identifiant
      Node * n = Arbre1->getNode(Id1);
      BppString * Espece = dynamic_cast<BppString*> (n->getNodeProperty(esp));
      OfficSORTIE_dup<<c.from_s(Espece->toSTL())<<"\t";
      OfficSORTIE_dup<<c.to_s(num1)+sep+(*it_adj).gene1<<"\t";
      OfficSORTIE_dup<<c.to_s(num2)+sep+(*it_adj).gene2<<endl;
   }
}

//Écrit le fic SORTIE_adj (nouvelle version du fic des adjacences ancestrales) sous la forme :
//G1=num1|id G2=num2|id 
//ATTENTION : que si le noeud est un noeud de spéciation !!
void ecritSORTIEAdj(Node * n, int num1, int num2, ofstream & OfficAdj2){
   //On ne fait rien pour les feuilles
   if (n->getNumberOfSons()!=0){
      //Est-ce un noeud de spéciation ?
      BppString * TYP = dynamic_cast<BppString*> (n->getBranchProperty(typ));
      if(TYP->toSTL()==spe){
         convert<int> c;
         //n° espece
         BppString * Espece = dynamic_cast<BppString*> (n->getNodeProperty(esp));
         OfficAdj2<<c.from_s(Espece->toSTL())<<"\t";
         //Les deux gènes adjacents
         string nom=n->getName();
         string Id1=nom.substr(0,nom.find(sepAdj));
         string Id2=nom.substr(nom.find(sepAdj)+1,nom.size());
         OfficAdj2<<num1<<sep<<Id1<<"\t"<<num2<<sep<<Id2<<endl;
      }
   }
   if (n->getNumberOfSons()==1)
      ecritSORTIEAdj(n->getSon(0),num1,num2,OfficAdj2);
   else if (n->getNumberOfSons()==2){
      ecritSORTIEAdj(n->getSon(0),num1,num2,OfficAdj2);
      ecritSORTIEAdj(n->getSon(1),num1,num2,OfficAdj2);
   }
}

//Les deux arbres passés en paramètre n'ont pas forcément des numéros d'identifiants
//"complets" (ça peut être des sous-arbres !)
void DECO(TreeTemplate<Node>* S, TreeTemplate<Node>* Arbre1, TreeTemplate<Node>* Arbre2, vector<adjacence> * Adj_classe, ofstream& OfficAdj2, ofstream & OfficArbresAdj, ofstream & OfficSORTIE_dup, string nom_classe){
   if(affich_CalculCout||affich_ParcoursArriere){
      cout<<"\n\n*******************************\nOn entre dans la fonction DECO !!"<<endl;
      cout<<"On traite les adjacences "<<endl;
      vector<adjacence>::iterator it_adj;
      for(it_adj=Adj_classe->begin();it_adj!=Adj_classe->end();it_adj++)
         cout<<(*it_adj).gene1<<" - "<<(*it_adj).gene2<<endl;
      cout<<endl;
   }
   vector<Node *> noeudsA1,noeudsA2;
   vector<Node *>::iterator itv1,itv2;
   ordonnePostfixe(Arbre1->getRootNode(),noeudsA1);
   ordonnePostfixe(Arbre2->getRootNode(),noeudsA2);
   
   if(affich_CalculCout||affich_ParcoursArriere){
      cout<<"Affichage des noeuds ordonnés :"<<endl;
      BppString * espece= dynamic_cast <BppString *> (Arbre1->getRootNode()->getNodeProperty(esp));
      cout<<"\n\t Arbre 1 de racine d'id"<<Arbre1->getRootNode()->getId()<<", d'espèce "<<espece->toSTL()<<endl;
      afficheInfoNoeuds(Arbre1->getRootNode());   
      
      espece= dynamic_cast <BppString *> (Arbre2->getRootNode()->getNodeProperty(esp));
      cout<<"\n\t Arbre 2 de racine d'id"<<Arbre2->getRootNode()->getId()<<", d'espèce "<<espece->toSTL()<<endl;
      afficheInfoNoeuds(Arbre2->getRootNode());   
      
      cout<<"\n\n********************* DÉBUT CALCUL *********************"<<endl;
   }
   
   //Remise à zéro des map de coût qui sont des var. glob. pour ce nouveau calcul
   C0.clear();C1.clear();
   
   for(itv1=noeudsA1.begin();itv1!=noeudsA1.end();itv1++)
      for(itv2=noeudsA2.begin();itv2!=noeudsA2.end();itv2++)
         if (Espece(*itv1)==Espece(*itv2)){
            pair<Node *,Node *> p((*itv1),(*itv2));
            pair<Node *,Node *> p2((*itv2),(*itv1));
            C1[p]=calculeC1(*itv1,*itv2,Adj_classe);   //Calcul C1
            C1[p2]=C1[p];
            if(affich_CalculCout||affich_ParcoursArriere)
               cout<<"\t C1("<<(*itv1)->getId()<<","<<(*itv2)->getId()<<")="<<C1[p]<<endl;
            C0[p]=calculeC0(*itv1,*itv2,Adj_classe);   //Calcul C0
            C0[p2]=C0[p];
            if(affich_CalculCout||affich_ParcoursArriere)
               cout<<"\t C0("<<(*itv1)->getId()<<","<<(*itv2)->getId()<<")="<<C0[p]<<endl;
         }
   float c1=recupCoutC1(Arbre1->getRootNode(),Arbre2->getRootNode());
   float c0 =recupCoutC0(Arbre1->getRootNode(),Arbre2->getRootNode());
   
   if(affich_CalculCout||affich_ParcoursArriere){
      cout<<"\n**********************"<<endl;
      cout<<"C0 aux racines : "<<c0<<endl;
      cout<<"C1 aux racines (+Gain): "<<c1+Crea<<endl;
      cout<<"\n\n\n\n\n\n"<<endl;
   }

   if(affich_ParcoursArriere){
      cout<<"\n**********************"<<endl;
      cout<<"Parcours arrière :"<<endl;
   }
   id_arbres_adj=0;
   vector<Tree *> * ArbresDAdjacences= new vector<Tree *>();
   DUP.clear();
   
   //Pour avoir un nombre aléatoire : 
   float r=rand()/(float)RAND_MAX;
   bool avantageC1;
   if (r<=Adj_percentage) //en cas C0=C1 on choisit C1
      avantageC1=true;
   else
      avantageC1=false;

   if (c0>c1+Crea|| (c1+Crea==c0 && avantageC1)){
      //Création du premier arbre d'adjacence
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(Arbre1->getRootNode(),Arbre2->getRootNode(),Adj_classe,ArbresDAdjacences,true));
      ArbresDAdjacences->push_back(T);
   }
   else{
      parcoursArriereC0(Arbre1->getRootNode(),Arbre2->getRootNode(),Adj_classe,ArbresDAdjacences,true);
   }

   if(affich_CalculCout||affich_ParcoursArriere)
      cout<<"\nÉcriture des arbres d'adjacences et des adjacences ancestrales de cette classe dans les fichiers correspondants."<<endl;
   OfficArbresAdj<<"Adjacency tree(s) solution of class "<<nom_classe<<endl; //<<" ( c0(R1,R2)="<<c0<<" | c1(R1,R2)="<<c1+Crea<<" )"
   vector<Tree *>::iterator it;
   int i=1;
   convert<int> c;
   int num1; int num2;
   num1=c.from_s(nom_classe.substr(0,nom_classe.find(sep)));
   if (nom_classe.find(sep)==nom_classe.rfind(sep))
      num2=c.from_s(nom_classe.substr(nom_classe.find(sep)+1,nom_classe.size()-1));
   else   //cas où on a une classe du type 12|12|57
      num2=c.from_s(nom_classe.substr(nom_classe.find(sep)+1,nom_classe.rfind(sep)-1));

   //On range les gènes dupliqués ensemble dans le fichier SORTIE_dup
   ecritSORTIEDup(DUP,num1,Arbre1,num2,OfficSORTIE_dup);

   for(it=ArbresDAdjacences->begin();it!=ArbresDAdjacences->end();it++){
      Tree * Gp=*it;
      TreeTemplate<Node> *G = dynamic_cast <TreeTemplate<Node> *> (Gp);                  
      if(affich_ParcoursArriere){
         //Affichage écran
         cout<<"\n****\nArbre "<<i<<" :"<<endl;
         afficheInfoNoeuds(G->getRootNode());
         cout<<endl<<endl;
      }
      //Écriture dans le fichier des adjacences ancestrales
      OfficArbresAdj<<"\nTree "<<i<<":"<<endl;
      if(OUTPUT_FORMAT==0){
         Newick * newickWriter = new Newick(); 
         newickWriter->enableExtendedBootstrapProperty(esp);
         if (G->getNumberOfLeaves()!=0){
            affecteInfoSurBrancheNoeuds(G->getRootNode());
            newickWriter->write(*G,OfficArbresAdj);
         }
         else
            OfficArbresAdj<<G->getRootNode()->getName()<<endl;
         delete(newickWriter);
      }
      else if(OUTPUT_FORMAT==1){
         Nhx * nhxWriter = new Nhx(true);
         nhxWriter->registerProperty(Nhx::Property("Duplication", "D", true, 0)); 
         if (G->getNumberOfLeaves()!=0)
            nhxWriter->write(*G,OfficArbresAdj);
         else
            OfficArbresAdj<<G->getRootNode()->getName()<<endl;
         delete(nhxWriter);
      }
      else{
         cout<<"From DECO.cpp (DECO): ERROR wrong OUTPUT_FORMAT"<<endl;
         exit(EXIT_FAILURE);
      }
      i++;
      //Calcul et écriture dans le fichier des adjacences ancestrales
      ecritSORTIEAdj(G->getRootNode(),num1,num2,OfficAdj2);
   }
   if(affich_ParcoursArriere)
      cout<<"\n\n\n"<<endl;
   OfficArbresAdj<<"\n***************************************************************************"<<endl;

   //On libère la mémoire du vecteur d'arbres d'adjacences
   for(it=ArbresDAdjacences->begin();it!=ArbresDAdjacences->end();it++)
      delete(*it);
   ArbresDAdjacences->clear();
   delete(ArbresDAdjacences);
}

int main(int argc, char* argv[]){
   time_t tbegin=time(NULL);              // get the current calendar time
   cout<<"\n\t#########################"<<endl;
   cout<<"\t###  Start Step3 !!!  ###"<<endl;
   cout<<"\t#########################"<<endl<<endl;

   // Récupérer données fichier de config pour recéer gene_species_EXT et pour récupérer le préfixe pour rechercher le fichier de classes d'adjacences
   lireFicConfig(fic_arbre,fic_gene,fic_especes,fic_adjacence,exp_name,directory,argc,argv);
   //On regarde si un fichier d'adjacences nettoyées existe, et si oui, on l'utilise
   if (is_readable(fic_adjacence_clean))
      fic_adjacence=fic_adjacence_clean;
   //Idem pour le fichier de gènes
   if (is_readable(fic_gene_clean))
      fic_gene=fic_gene_clean;
   
   cout<<"\nFile(s) used in Step3:"<<endl;
   cout<<"\tReconciled trees file: "<<tree_reconciled_file<<endl;
   cout<<"\tAdjacencies classes file: "<<adj_class_file<<endl<<endl;


   if(file_log){
      //Get time
      time_t rawtime;
      struct tm datetime;
      char str_time[50];
      time(&rawtime);
      datetime = *localtime(&rawtime);
      strftime(str_time, 50, "%A %d %B %Y %H:%M:%S", &datetime);

      // Creation and/or opening log file
      if(is_readable(log_file)){
         Offile_log.open(log_file.c_str(), ios::out|ios::app);
         if(!Offile_log){
            cout<<"\tERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file created by Step3_DECO.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      else{
         Offile_log.open(log_file.c_str(), ios::out|ios::trunc);
         if(!Offile_log){
            cout<<"\tERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file created by Step3_DECO.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      Offile_log<<"\t#########################"<<endl;
      Offile_log<<"\t###  Start Step3 !!!  ###"<<endl;
      Offile_log<<"\t#########################"<<endl<<endl;
      Offile_log<<"\tFile(s) used in Step3:"<<endl;
      Offile_log<<"\t\t- Reconciled trees file: "<<tree_reconciled_file<<endl;
      Offile_log<<"\tAdjacencies classes file: "<<adj_class_file<<endl<<endl;
   }


   // Variables
//   TreeTemplate<Node> S;
//   S = new TreeTemplate<Node>;
   map<string,string> gene_species_EXT;
   vector<Tree*> Arbres;
   map<string,vector<int> > gene_GeneTreeID;
   map<string,vector<int> > classes_adjacences;
   vector<adjacence> adjacencies;
   map<int,vector<adjacence> > Adj_actuelles;


//   //Store species tree in TreeTemplate<Node> *S
//   StoreSpeciesTree(S);

   cout<<"Store Species tree in TreeTemplate<Node> * S..."<<flush;
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_classes.cpp (StoreSpeciesTree): Store Species tree in TreeTemplate<Node> * S..."<<flush;
   Newick * newickReaderEsp = new Newick(true,true);
   newickReaderEsp->enableExtendedBootstrapProperty(esp);
   TreeTemplate<Node> *S = newickReaderEsp->read(fic_especes);
   delete(newickReaderEsp);
   cout<<" DONE"<<endl<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;

   //Les noeuds internes sont lus comme des propriétés de branches
   //La f° ci-dessous permet de nommer les noeud internes et de faire passer la propriété au noeud
   if (file_log)
      Offile_log<<"OS\tFrom Step3_DECO.cpp (nommerNoeudsInternes): Affect species name to node in species tree ..."<<flush;
   nommerNoeudsInternes(S->getRootNode());
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;


/////////////////////////////////
/////   START Genes Clean   /////
/////////////////////////////////

   //Store arbres réconciliés dans un vecteur Arbres
   StoreReconciledTrees(Arbres);

   //Associate Species_name with gene in map<string, string> gene_species_EXT
   AssociateGeneWithSpecies(gene_species_EXT);

   //Associe les gènes avec les numéros d'arbres réconciliés dans lesquels ils sont présents
   AssociateGeneWithReconciledTreeNb(Arbres,gene_GeneTreeID,gene_species_EXT);
//   afficheMap(gene_GeneTreeID);

   //Nettoyage de la map gene_species_EXT pour ne garder que les gènes présents dans nos arbres
   int elim = 0;
   elim=GeneSpeciesClean(gene_GeneTreeID,gene_species_EXT);
   if(elim>0){
      cout<<"\nFrom Step3_DECO (main): WARNING while cleaning "<<elim<<" genes have been deleted"<<endl<<endl;
   }
   gene_GeneTreeID.clear();

///////////////////////////////
/////   END Genes Clean   /////
///////////////////////////////


   //Reconstruction de classes adjacences et à partir du fichier adj_class_file
   AdjClassVector(classes_adjacences);
//   afficheMap(classes_adjacences);

   //Construction du vecteur d'adjacences ordonnés (bonne numérotation) et de l'association Id espèce actuelle avec vector<adj>
   AdjVector_Esp_AdjExtant(S,gene_species_EXT,adjacencies,Adj_actuelles);
   //Affichage des adjacences actuelles / espèce
//   afficheMap(Adj_actuelles);
   Adj_actuelles.clear();


//////////////////////////////////
////   Opening OUTPUT files   ////
//////////////////////////////////

   // Ancestral adjacencies file
   string fic_adj_ances=prefixe+output_adj;
   ofstream OfficAdj2;

   // Adjacencies trees file
   string fic_arbresdAdjacences=prefixe+output_adj_trees_per_class;
   ofstream OfficArbresAdj;

   //fic sortie_dup des gènes dupliqués ensembles
   string fic_SORTIE_dup=prefixe+output_dup_gene_pairs;
   ofstream OfficSORTIE_dup;

//   //fic sortie sous arbres
//   string fic_sousArbresTraites=prefixe+output_treated_subT;
//   ofstream OfficSousArbres;

   if(!light_mode){
      // Ancestral adjacencies file
      OfficAdj2.open(fic_adj_ances.c_str(), ios::out|ios::trunc);
      if (!OfficAdj2){
         cout<<"\nFrom Step3_DECO.cpp: ERROR while opening file "<<fic_adj_ances<<endl;
         if (file_log)
            Offile_log<<"\nFI\tFrom Step3_DECO.cpp: ERROR while opening file "<<fic_adj_ances<<endl;
         exit(EXIT_FAILURE);
      }
      cout<<"File "<<fic_adj_ances<<" open"<<endl<<endl;
      if (file_log)
         Offile_log<<"FI\tFrom Step3_DECO.cpp: File "<<fic_adj_ances<<" open"<<endl;

      // Adjacencies trees file
      OfficArbresAdj.open(fic_arbresdAdjacences.c_str(), ios::out|ios::trunc);
      if (!OfficArbresAdj){
         cout<<"\nFrom Step3_DECO.cpp: ERROR while opening file "<<fic_arbresdAdjacences<<endl;
         if (file_log)
            Offile_log<<"\nFI\tFrom Step3_DECO.cpp: ERROR while opening file "<<fic_arbresdAdjacences<<endl;
         exit(EXIT_FAILURE);
      }
      cout<<"File "<<fic_arbresdAdjacences<<" open"<<endl<<endl;
      if (file_log)
         Offile_log<<"FI\tFrom Step3_DECO.cpp: File "<<fic_arbresdAdjacences<<" open"<<endl;

      //fic sortie_dup des gènes dupliqués ensembles
      OfficSORTIE_dup.open(fic_SORTIE_dup.c_str(), ios::out|ios::trunc);
      if (!OfficSORTIE_dup){
         cout<<"\nFrom Step3_DECO.cpp: ERROR while opening file "<<fic_SORTIE_dup<<endl;
         if (file_log)
            Offile_log<<"\nFI\tFrom Step3_DECO.cpp: ERROR while opening file "<<fic_SORTIE_dup<<endl;
         exit(EXIT_FAILURE);
      }
      cout<<"File "<<fic_SORTIE_dup<<" open"<<endl<<endl;
      if (file_log)
         Offile_log<<"FI\tFrom Step3_DECO.cpp: File "<<fic_SORTIE_dup<<" open"<<endl;

//      //fic sortie sous arbres
//      OfficSousArbres.open(fic_sousArbresTraites.c_str(), ios::out|ios::trunc);
//      if (!OfficSousArbres){
//         cout<<"\nFrom Step3_DECO.cpp: ERROR while opening file "<<fic_sousArbresTraites<<endl;
//         if (file_log)
//            Offile_log<<"\nFI\tFrom Step3_DECO.cpp: ERROR while opening file "<<fic_sousArbresTraites<<endl;
//         exit(EXIT_FAILURE);
//      }
//      cout<<"File "<<fic_sousArbresTraites<<" open"<<endl<<endl;
//      if (file_log)
//         Offile_log<<"FI\tFrom Step3_DECO.cpp: File "<<fic_sousArbresTraites<<" open"<<endl;
   }


////////////////////////////////////////////////
/////   Création des arbres d'adjacences   /////
////////////////////////////////////////////////

   cout<<"Creation of Adjacencies Trees from adjacencies classes ..."<<flush;
   if(file_log)
      Offile_log<<"OS\tFrom Step3_DECO.cpp: Creation of Adjacencies Trees from adjacencies classes ..."<<flush;

   // VARIABLES
   convert<int> c;
   TreeTemplate<Node> * Arbre1;//= new  TreeTemplate<Node>();
   TreeTemplate<Node> * Arbre2;//= new  TreeTemplate<Node>();
   int num_diff,num_ident;
   vector<adjacence> * Adj_classe= new vector<adjacence> ();
   num_diff=num_ident=0;
   //    1) Pour chaque classe
   vector<int>::iterator it_adj;
   map<string,vector<int> >::iterator it_cl;
   for(it_cl=classes_adjacences.begin();it_cl!=classes_adjacences.end();it_cl++){
      Adj_classe->clear();

      //         a) Récupérer arbres de gènes
      string nom_classe=(*it_cl).first;
      int num1=c.from_s(nom_classe.substr(0,nom_classe.find(sep)));
      int num2=c.from_s(nom_classe.substr(nom_classe.find(sep)+1,nom_classe.length()-1));
      if(affich_Classes)
         cout<<"On traite la classe "<<nom_classe<<endl;
      if(num1!=num2){
         if(affich_Classes){
            cout<<"\tNuméros différents."<<endl;
            cout<<"Extraction des identifiants des noeuds où couper :"<<endl;    
         }
         int coupe1,coupe2;
         coupe1=c.from_s(nom_classe.substr(nom_classe.rfind(sep)+1,nom_classe.find(sepAdj)));
         coupe2=c.from_s(nom_classe.substr(nom_classe.find(sepAdj)+1,nom_classe.length()-1));
         if(affich_Classes){
            cout<<"\tCoupe dans l'arbre 1 à l'id "<<coupe1<<endl;
            cout<<"\tCoupe dans l'arbre 2 à l'id "<<coupe2<<endl;
         }
         //Les manips ci-dessous sont-elles bien nécessaires ?
//         Arbre1= dynamic_cast <TreeTemplate<Node> *> (Arbres[num1]);
//         Arbre2= dynamic_cast <TreeTemplate<Node> *> (Arbres[num2]);
//         TreeTemplate<Node > A1=*Arbre1;
//         TreeTemplate<Node > A2=*Arbre2;
         TreeTemplate<Node > A1=*(Arbres[num1]);
         TreeTemplate<Node > A2=*(Arbres[num2]);
         Arbre1= new TreeTemplate<Node>();
         Arbre1=A1.cloneSubtree(coupe1);
         Arbre2= new  TreeTemplate<Node>();
         Arbre2=A2.cloneSubtree(coupe2);

         num_diff++;
      }
      else{
         if(affich_Classes)
            cout<<"\tMême numéros."<<endl;
         TreeTemplate<Node > ArbreUnique=*(Arbres[num1]);
         vector<int> Id_Genes;
         //parcours des adjacences de la classe pour récupérer les gènes puis leur identifiant dans l'arbre
         for(it_adj=(*it_cl).second.begin();it_adj!=(*it_cl).second.end();it_adj++){
            //gène 1
            string nom_gene =adjacencies[*it_adj].gene1+sep+gene_species_EXT[adjacencies[*it_adj].gene1];
            Id_Genes.push_back(ArbreUnique.getLeafId(nom_gene));
            //gène 2
            nom_gene =adjacencies[*it_adj].gene2+sep+gene_species_EXT[adjacencies[*it_adj].gene2];
            Id_Genes.push_back(ArbreUnique.getLeafId(nom_gene));
         }
         int id_LCA=TreeTools::getLastCommonAncestor(ArbreUnique,Id_Genes);
         Id_Genes=ArbreUnique.getSonsId(id_LCA);
         Arbre1= new  TreeTemplate<Node>();
         Arbre1=ArbreUnique.cloneSubtree(Id_Genes[0]);
         Arbre2= new  TreeTemplate<Node>();
         Arbre2=ArbreUnique.cloneSubtree(Id_Genes[1]); 
         num_ident++;
      }

       
      //         b) Récupérer les adjacences
      //OfficSousArbres<<"Adjacencies to be treated:"<<endl;
      for(it_adj=(*it_cl).second.begin();it_adj!=(*it_cl).second.end();it_adj++){
         //OfficSousArbres<<adjacencies[*it_adj].gene1<<" - "<<adjacencies[*it_adj].gene2<<endl;
         Adj_classe->push_back(adjacencies[*it_adj]);
      }

      //         c) Appliquer l'algo DéCo (produit des adjacences ancestrales qu'on accumule dans 
      //un fichier de même type que le fichier adjacences d'entrée mais où les noms d'espèces sont 
      //ceux des espèces ancestrales)
      if(INPUT_FORMAT==1){
         affecteInfoSurBrancheNoeuds(Arbre1->getRootNode());
         affecteInfoSurBrancheNoeuds(Arbre2->getRootNode());
      }
      // OfficSousArbres<<"\nSubtree 1 (tree "<<num1<<") :"<<endl;
      // if (Arbre1->getNumberOfLeaves()!=0)
      //    newickReaderEsp->write(*Arbre1,OfficSousArbres);
      // else
      //    OfficSousArbres<<Arbre1->getRootNode()->getName()<<endl;
      affecteInfoSurBrancheNoeuds(Arbre2->getRootNode());
      // OfficSousArbres<<"\nSubtree 2 (tree "<<num2<<") :"<<endl;
      // if (Arbre2->getNumberOfLeaves()!=0)
      //    newickReaderEsp->write(*Arbre2,OfficSousArbres);
      // else
      //    OfficSousArbres<<Arbre2->getRootNode()->getName()<<endl;
      // OfficSousArbres<<"\n***************************************************************************\n"<<endl;

      DECO(S,Arbre1,Arbre2,Adj_classe,OfficAdj2,OfficArbresAdj,OfficSORTIE_dup,nom_classe);
      delete(Arbre1);
      delete(Arbre2);
   }
   //On libère la mémoire du vecteur d'arbres
   vector<Tree*>::iterator ita;
   for (ita=Arbres.begin();ita!=Arbres.end();ita++)
      delete(*ita);
   Arbres.clear();
   gene_species_EXT.clear();
   classes_adjacences.clear();
   adjacencies.clear();
   if(!light_mode){
      OfficAdj2.close();
      OfficArbresAdj.close();
      OfficSORTIE_dup.close();
//      OfficSousArbres.close();
   }
   cout<<" DONE"<<endl<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;

   ////////////////////////////////////////////////////////////////////////////////////
   //STATS
   
   ////////////////////////////////////////////////////////////////////////////////////
   // Fichier de stats de DECO
   string DECO_stats_file=prefixe+output_DECO_stats;
   cout<<"Creation of the file for statistics: "<<DECO_stats_file<<" ..."<<flush;
   if(file_log)
      Offile_log<<"FI\tFrom Step3_DECOProba.cpp: Creation of the file for statistics: "<<DECO_stats_file<<" ..."<<flush;
   ofstream OffilStatsDECO;
   OffilStatsDECO.open(DECO_stats_file.c_str(), ios::out|ios::trunc);
   if (!OffilStatsDECO){
      cout<<"\nFrom Step3_DECOProba.cpp (main): ERROR while opening file "<<DECO_stats_file<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step3_DECOProba.cpp (main): ERROR while opening file "<<DECO_stats_file<<endl;
      exit(EXIT_FAILURE);
   }
   OffilStatsDECO<<"Statistics file of Step3-DECO"<<endl<<endl;
   OffilStatsDECO<<"Number of Adjacencies creation - nb_Crea= "<<nb_Crea<<endl;
   OffilStatsDECO<<"Number of Gene Duplication - nb_GDup= "<<nb_GDup<<endl;
   OffilStatsDECO<<"Number of Adjacencies Break - nb_Bk= "<<nb_Bk<<endl;
   OffilStatsDECO<<"Number of Adjacencies Duplication - nb_ADup= "<<nb_ADup;
   OffilStatsDECO.close();
   cout<<" DONE"<<endl<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;
   ////////////////////////////////////////////////////////////////////////////////////

   time_t tend=time(NULL);                // get the current calendar time
   // Compute execution time
   float texec=difftime(tend,tbegin);    // tend-tbegin (result in second)

   cout<<"Execution time of Step3 is: "<<texec<<"s"<<endl<<endl;
   cout<<"\t#######################"<<endl;
   cout<<"\t###  End Step3 !!!  ###"<<endl;
   cout<<"\t#######################"<<endl<<endl;

   if (file_log){
      Offile_log<<"Execution time of Step3 is: "<<texec<<"s"<<endl<<endl;
      Offile_log<<"\t#######################"<<endl;
      Offile_log<<"\t###  End Step3 !!!  ###"<<endl;
      Offile_log<<"\t#######################"<<endl<<endl;
   }
   Offile_log.close();

   return(0);
}
