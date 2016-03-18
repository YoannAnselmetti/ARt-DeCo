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

File: Step2_adjacencies_class.cpp               Last modified on: 22/05/2015
Created by: Yoann Anselmetti/Sèverine Bérard    Created on: 10/03/2014
--------------------------------------------------------------------------
Specification:
File use to create adjacencies classes wtih reconciled trees and extant adjacencies
=========================================================================*/
#include "Step2_adjacencies_class.h"



//renvoie vrai si le noeud d'id1 est ancêtre du noeud d'id2 dans A
bool isAncestor(int id1, int id2, TreeTemplate<Node> *A) 
{
   vector< int > Id_Genes;
   Id_Genes.push_back(id1);
   Id_Genes.push_back(id2);
   TreeTemplate<Node > A1=*A;
   int id_LCA=TreeTools::getLastCommonAncestor(A1,Id_Genes);
   return(id_LCA==id1);
}


//Renvoie le numéro d'identifiant du noeud d'espèce s
//au dessus du noeud n dans G
int trouveNoeudSpeSAuDessusDeN(int s, int id_n,TreeTemplate<Node> * G){
   Node * n = G->getNode(id_n);
   convert<int> C;
   BppString * ESPECE = dynamic_cast <BppString *> (n->getNodeProperty(esp));
   string espece = ESPECE->toSTL();
   while (n->hasFather()&&(C.from_s(espece)!=s)){
      n=n->getFather();
      ESPECE = dynamic_cast <BppString *> (n->getNodeProperty(esp));
      espece = ESPECE->toSTL();
   }
   //delete(n); //si on fait ça on détruit les données d'origine !!!
   //delete(ESPECE);
   return(n->getId());
}


void retourneIdSpeSAuDessousDeN(vector<int> & Ids, int s, Node * n,TreeTemplate<Node> * G){
   convert<int> C;
   BppString * ESPECE = dynamic_cast <BppString *> (n->getNodeProperty(esp));
   string espece = ESPECE->toSTL();
   if(C.from_s(espece)==s)
      Ids.push_back(n->getId());
   else if (n->getNumberOfSons()>1){
      retourneIdSpeSAuDessousDeN(Ids,s, n->getSon(0),G);
      retourneIdSpeSAuDessousDeN(Ids,s, n->getSon(1),G);
   }
   //delete(ESPECE);
}

//Adjacencies class Pre-construction
void AdjClassPreConstruction(TreeTemplate<Node> *S, vector<Tree*> &Arbres, map<string,vector<int> > &gene_GeneTreeID, map<string,string> &gene_species_EXT, map<string,vector<int> > &classes_adjacences, vector<adjacence> &adjacencies){
   
   ifstream IfficAdj(fic_adjacence.c_str(), ios::in);
   if (!IfficAdj){
      cout<<"\nFrom Step2_adjacencies_class.cpp (AdjClassPreConstruction): ERROR while opening file "<<fic_adjacence<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step2_adjacencies_class.cpp (AdjClassPreConstruction): ERROR while opening file "<<fic_adjacence<<endl;
      exit(EXIT_FAILURE);
   }
   cout<<"File "<<fic_adjacence<<" open"<<endl<<endl;
   if (file_log){
      Offile_log<<"FI\tFrom Step2_adjacencies_class.cpp: File "<<fic_adjacence<<" open"<<endl;
      Offile_log<<"OS\tFrom Step2_adjacencies_class.cpp: Store adjacencies in vector<adjacence> adjacences / Pre-construction of adjacencies classes in map<string,vector<int> > classes_adjacences ..."<<flush;
   }
   convert<int> c;
   string nom_classe;
   int j=0;   //j est le numéro des adjacences, on les numérote à 0
   adjacence adj;
   string gene1="";
   string gene2="";
   vector<int> ab1,ab2;
   vector<int>::iterator it_ab1,it_ab2;
   while(!IfficAdj.eof()){
      //on lit la ligne en entier;
      IfficAdj>>gene1;
      //Pour éviter un dernier tour de trop à cause d'un blanc en fin de fichier
      if(!IfficAdj.eof()){
         IfficAdj>>gene2;
         //cout<<"On traite l'adjacence "<< gene1<<" - "<<gene2;
         //on crée l'adjacence si les deux gènes sont dans gene_species_EXT (c-à-d dans nos arbres)
         if (gene_species_EXT.find(gene1)!=gene_species_EXT.end() && gene_species_EXT.find(gene2)!=gene_species_EXT.end()){
       //            cout<<", on la nomme "<<j<<endl;
            adj.gene1=gene1;
            adj.gene2=gene2;
            //on la stocke
            adjacencies.push_back(adj);

//INUTILE DE CRÉER LA MAP ADJ_ACTUELLES
//            //Storeage de l'adj actuelle dans son espèce dans Adj_actuelles
//            int esp_actuelle=S->getNode(gene_species_EXT[gene1])->getId();
//            if(Adj_actuelles.find(esp_actuelle)!=Adj_actuelles.end())
//               Adj_actuelles[esp_actuelle].push_back(adj);
//            else{
//               vector<adjacence> v;
//               v.push_back(adj);
//               Adj_actuelles[esp_actuelle]=v;
//            }

            //on s'intéresse aux classes :
            ab1=gene_GeneTreeID[gene1];
            ab2=gene_GeneTreeID[gene2];
       
            //Permet de parcourir les arbres dans lesquels les deux gènes de l'adjacence sont présents pour pouvoir les associer dans des classes d'adjacences
            //et améliorer la vitesse de calcul pour la reconstruction des arbres d'adjacences
            for(it_ab1=ab1.begin();it_ab1!=ab1.end();it_ab1++){
               for(it_ab2=ab2.begin();it_ab2!=ab2.end();it_ab2++){
                  //Attention tous les adjacences des classes entre même arbre 
                  //ne peuvent pas forcément être traitées ensemble !!
                  //Il faut qu'elles respectent même noeud ancêtre commun des extrémités !!!
                  //Nom de la classe x|x|idLCA_dans_x
                  //ATTENTION BIS idem pour les classes entre arbres de numéros différents
                  //on les subdivisera plus bas car il faut les trier en les comparant !!
        
                  if((*it_ab2)==(*it_ab1)){
                     TreeTemplate<Node> *A=dynamic_cast <TreeTemplate<Node> *> (Arbres[(*it_ab1)]);
                     TreeTemplate<Node> Arbre=*A;
           vector< int > Id_Genes;
                     string g1=gene1+sep+gene_species_EXT[gene1];
                     Id_Genes.push_back(Arbre.getLeafId(g1));
                     string g2=gene2+sep+gene_species_EXT[gene2];
                     Id_Genes.push_back(Arbre.getLeafId(g2));
                     int id_ancetre=TreeTools::getLastCommonAncestor(Arbre,Id_Genes);
                     nom_classe=c.to_s(*it_ab1)+sep+c.to_s(*it_ab2)+sep+c.to_s(id_ancetre);
           //delete(A);//provoque une erreur de segmentation MAIS OK car pointe vers Arbres[(*it_ab1)]
                  }
                  else 
                     if((*it_ab2)<(*it_ab1))
                        nom_classe=c.to_s(*it_ab2)+sep+c.to_s(*it_ab1);
                     else
                        nom_classe=c.to_s(*it_ab1)+sep+c.to_s(*it_ab2);
                  if (classes_adjacences.find(nom_classe)!=classes_adjacences.end()){
                     //cout<<nom_classe<<" est déjà dans classes_adjacences"<<endl;
                     //on vérifie qu'on a pas déjà ajouté cette adjacence
                     //cas où g1 et g2 appartiennent tous deux à ab1 et ab2
                     if (classes_adjacences[nom_classe].back()!=j)
                        classes_adjacences[nom_classe].push_back(j);
                  }
                  else{
                     //cout<<nom_classe<<" n'est PAS dans classes_adjacences"<<endl;
                     vector<int> v;
                     v.push_back(j);
                     classes_adjacences[nom_classe]=v;
                  }
               }
               j++;
            }
         }
         else{
       //            cout<<gene1<<" "<<gene2<<" --> on la jette."<<endl;
         }
      }
   }
   IfficAdj.close();
   // if(IfficAdj.fail())
   //    {
   //     cout<<"\nFrom Step2_adjacencies_class.cpp (AdjClassPreConstruction): ERROR while closing file "<<fic_adjacence<<endl;
   //     if (file_log)
   //        Offile_log<<"\nFI\tFrom Step2_adjacencies_class.cpp (AdjClassPreConstruction): ERROR while closing file "<<fic_adjacence<<endl;
   //     exit(EXIT_FAILURE);
   // }
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;
   
   cout<<"Adjacencies file: "<<fic_adjacence<<" read: "<<adjacencies.size()<<" adjacencies"<<endl<<endl;
   if(file_log)
      Offile_log<<"FI\tFrom Step2_adjacencies_class.cpp: Adjacencies file: "<<fic_adjacence<<" read: "<<adjacencies.size()<<" adjacencies"<<endl<<endl;
}


void AdjClassConstruction(TreeTemplate<Node> *S, vector<Tree*> &Arbres, map<string,string> &gene_species_EXT, map<string,vector<int> > &classes_adjacences, vector<adjacence> &adjacencies){
   cout<<"Adjacencies classes under construction..."<<flush;
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_class.cpp: Adjacencies classes under construction..."<<flush;
   
   convert<int> c;
   string nom_classe;
   int num1,num2;
   vector<int>::iterator it_adj;
   TreeTemplate<Node> * Arbre1;//= new  TreeTemplate<Node>();
   TreeTemplate<Node> * Arbre2;//= new  TreeTemplate<Node>();
   vector<string> classes_adjacences_A_SUPPR;
   vector<string>::iterator it_class;
   map<string,vector<int> > classes_adjacences_SUB;
   map<string, vector<int> >::iterator it_cl;
   map<string, vector<int> >::iterator it_cl_prec;
   string nom_gene1;
   string nom_gene2;
   cout<<"\nAvant traitement des classes :"<<endl;
   cout<<"classes_adjacences.size()="<<classes_adjacences.size()<<endl;
   int classe1=0;
   int classeGARD=0;
   int classeSUB=0;
   int classeSUPPR=0;
   int nb_tour=0;
   for(it_cl=classes_adjacences.begin();it_cl!=classes_adjacences.end();it_cl++){
      nb_tour++;
      nom_classe=(*it_cl).first;
      if(affich_Classes)
         cout<<"\nOn traite la classe "<<nom_classe<<" de taille "<<(*it_cl).second.size()<<endl;
      //Élimination des classes de taille 1
      if ((*it_cl).second.size()==1){
         classe1++;
         if(affich_Classes)
            cout<<" : on la jette (taille 1)."<<endl;
         if(it_cl!=classes_adjacences.begin()){
            it_cl_prec=it_cl;
            it_cl_prec--; 
         }
         else
            it_cl_prec=classes_adjacences.begin();
         classes_adjacences.erase(it_cl);
         it_cl=it_cl_prec;
      }
      else{
         //Subdivision éventuelle des classes entre deux arbres différents
         if(affich_Classes)
            cout<<" : on l'examine,";
         num1=c.from_s(nom_classe.substr(0,nom_classe.find(sep)));
         if (nom_classe.find(sep)==nom_classe.rfind(sep))
            num2=c.from_s(nom_classe.substr(nom_classe.find(sep)+1,nom_classe.size()-1));
         else
            //cas où on a une classe du type 12|12|57
            num2=c.from_s(nom_classe.substr(nom_classe.find(sep)+1,nom_classe.rfind(sep)-1));
         if(affich_Classes)
            cout<<" numéros extraits : "<<num1<<" et "<<num2<<endl;
         if(num1!=num2){
            if(affich_Classes)
               cout<<"\n******* On subdivise éventuellement la classe "<<nom_classe<<" de taille "<<(*it_cl).second.size()<<"."<<endl;
            //On va renommer cette classe quoi qu'il en soit
            classes_adjacences_A_SUPPR.push_back((*it_cl).first);
       classeSUPPR++;
            Arbre1= dynamic_cast <TreeTemplate<Node> *> (Arbres[num1]);
            Arbre2= dynamic_cast <TreeTemplate<Node> *> (Arbres[num2]);
            TreeTemplate<Node > A1=*Arbre1;
            TreeTemplate<Node > A2=*Arbre2;
       //On recherche les LCA des extremités des adj dans chaque arbre
            vector< int > Id_Genes1;
            vector< int > Id_Genes2;
            //parcours des adjacences
            for(unsigned int a=0;a!=(*it_cl).second.size();a++){
               //Le numéro de l'adjacance
               string nom_gene1,nom_gene2;
               int va=(*it_cl).second[a];
               //gène 1
               nom_gene1 =adjacencies[va].gene1+sep+gene_species_EXT[adjacencies[va].gene1];
               //gène 2
               nom_gene2 =adjacencies[va].gene2+sep+gene_species_EXT[adjacencies[va].gene2];
               if(affich_Classes)
                  cout<<"Adj n°"<<va<<" "<<nom_gene1<<" - "<<nom_gene2<<endl;
               try{
                  Id_Genes1.push_back(Arbre1->getLeafId(nom_gene1));
                  Id_Genes2.push_back(Arbre2->getLeafId(nom_gene2));
               }
               catch (Exception e){
                  Id_Genes1.push_back(Arbre1->getLeafId(nom_gene2));
                  Id_Genes2.push_back(Arbre2->getLeafId(nom_gene1));
               }
            }
            if(affich_Classes){
               cout<<"Id_Genes1 : ";
               for(unsigned int f=0;f<Id_Genes1.size();f++)
                  cout<<Id_Genes1[f]<<" ";
               cout<<endl;
               cout<<"Id_Genes2 : ";
               for(unsigned int f=0;f<Id_Genes2.size();f++)
                  cout<<Id_Genes2[f]<<" ";
               cout<<endl;
            }
            int id_LCA1=TreeTools::getLastCommonAncestor(A1,Id_Genes1);
            int id_LCA2=TreeTools::getLastCommonAncestor(A2,Id_Genes2);
            if(affich_Classes)
               cout<<"id_LCA1="<<id_LCA1<<" et id_LCA2="<<id_LCA2<<endl;
       
            BppString * EspLCA1 = dynamic_cast <BppString *> (Arbre1->getNode(id_LCA1)->getNodeProperty(esp));
            if(affich_Classes)
               cout<<"Espèce du LCA1="<<EspLCA1->toSTL()<<endl;
            BppString * EspLCA2 = dynamic_cast <BppString *> (Arbre2->getNode(id_LCA2)->getNodeProperty(esp));
            if(affich_Classes)
               cout<<"Espèce du LCA2="<<EspLCA2->toSTL()<<endl;
            //Si les deux ancêtres sont de même espèce alors pas de subdivision
            if(EspLCA1->toSTL()==EspLCA2->toSTL()){
               if(affich_Classes)
                  cout<<"***** Cas LCA1 et LCA2 de même espèce"<<endl;
               //On change juste le nom :
               string nom=c.to_s(num1)+sep+c.to_s(num2)+sep+c.to_s(id_LCA1)+sepAdj+c.to_s(id_LCA2);
               if(affich_Classes)
                  cout<<"Nouveau nom de la classe : "<<nom<<endl;
               //Une classe à ajouter : nom différent mais mêmes adjacences
               classes_adjacences_SUB[nom]=(*it_cl).second;
               classeSUB++;
            }
            else{
               if(affich_Classes)
                  cout<<"***** Cas LCA1 et LCA2 d'espèces différentes, on cherche les espèces des racines"<<endl;
               BppString * r1 = dynamic_cast <BppString *> (Arbre1->getRootNode()->getNodeProperty(esp));
               if(affich_Classes)
                  cout<<"Espèce de la racine de l'arbre 1 : r1="<<r1->toSTL()<<endl;
               BppString * r2 = dynamic_cast <BppString *> (Arbre2->getRootNode()->getNodeProperty(esp));
               if(affich_Classes)
                  cout<<"Espèce de la racine de l'arbre 2 : r2="<<r2->toSTL()<<endl;
               if(isAncestor(c.from_s(EspLCA1->toSTL()),c.from_s(EspLCA2->toSTL()),S)){
                  if(affich_Classes)
                     cout<<"\tLCA1 ancêtre de LCA2"<<endl;
                  if(isAncestor(c.from_s(r2->toSTL()),c.from_s(EspLCA1->toSTL()),S)){
                     if(affich_Classes)
                        cout<<"\tr2 ancêtre de LCA1"<<endl;
                     //il faut remonter dans arbre 2, à partir de LCA2 et couper au niveau de l'espèce de LCA1
                     int coupe2=trouveNoeudSpeSAuDessusDeN(c.from_s(EspLCA1->toSTL()), id_LCA2, Arbre2);
                     //On change le nom :
                     string nom=c.to_s(num1)+sep+c.to_s(num2)+sep+c.to_s(id_LCA1)+sepAdj+c.to_s(coupe2);
                     if(affich_Classes)
                        cout<<"Nouveau nom de la classe : "<<nom<<endl;
                     //Une classe à ajouter : nom différent mais mêmes adjacences
                     classes_adjacences_SUB[nom]=(*it_cl).second;
                     classeSUB++;
                  }
                  else{
                     if(affich_Classes)
                        cout<<"\tLCA1 > r2 => subdivision obligatoire"<<endl;//À PROUVER
                     //Il faut couper dans arbre 1 pour trouver les noeuds d'espèce LCA2
                     vector<int> Ids;
                     retourneIdSpeSAuDessousDeN(Ids,c.from_s(EspLCA2->toSTL()),Arbre1->getNode(id_LCA1),Arbre1);
                     if(Ids.size()==0)
                        cout<<"On n'a pas trouvé de noeud d'espèce LCA2 en dessous de LCA1, IMPOSSIBLE";
                     else if(Ids.size()==1)
                        cout<<"On a trouvé un seul noeud d'espèce LCA2 en dessous de LCA1, IMPOSSIBLE";
                     else{
                        if(affich_Classes)
                           cout<<"On doit subdiviser en "<<Ids.size()<<" classes."<<endl;
                        //Les nouveaux noms :
                        for(unsigned int k=0;k<Ids.size();k++){
                           string nom=c.to_s(num1)+sep+c.to_s(num2)+sep+c.to_s(Ids[k])+sepAdj+c.to_s(id_LCA2);
                           if(affich_Classes)
                              cout<<"Nouveau nom de la classe sub : "<<nom<<endl;
                           //Une classe à ajouter, on l'initialise on y rangera les adj plus tard
                           vector<int> a;
                           classes_adjacences_SUB[nom]=a;
                           classeSUB++;
                        }
                        //On range les adjacences dans les nouvelles classes
                        //Les identifiants des extremités des adj se trouvent dans l'ordre dans Id_Genes1 et Id_Genes2
                        for(unsigned int a=0;a!=(*it_cl).second.size();a++){
                           unsigned int cpt=0;
                           while(!isAncestor(Ids[cpt],Id_Genes1[a],Arbre1)&&cpt<Ids.size())
                              cpt++;
                           if(cpt==Ids.size())
                              cout<<"Problème, on n'a pas trouvé de classe pour l'adj "<<a<<endl;
                           else{
                              string nom=c.to_s(num1)+sep+c.to_s(num2)+sep+c.to_s(Ids[cpt])+sepAdj+c.to_s(id_LCA2);
                              classes_adjacences_SUB[nom].push_back((*it_cl).second[a]);
                           }
                        }
                     }
                  }
               }
               else{
                  if(affich_Classes)
                     cout<<"\tLCA2 ancêtre de LCA1"<<endl;
                  if(isAncestor(c.from_s(r1->toSTL()),c.from_s(EspLCA2->toSTL()),S)){
                     if(affich_Classes)
                        cout<<"\tr1 ancêtre de LCA2"<<endl;
                     //il faut remonter dans arbre 1, à partir de LCA1 et couper au niveau de l'espèce de LCA2
                     int coupe1=trouveNoeudSpeSAuDessusDeN(c.from_s(EspLCA2->toSTL()), id_LCA1, Arbre1);
                     //On change le nom :
                     string nom=c.to_s(num1)+sep+c.to_s(num2)+sep+c.to_s(coupe1)+sepAdj+c.to_s(id_LCA2);
                     if(affich_Classes)
                        cout<<"Nouveau nom de la classe : "<<nom<<endl;
                     //Une classe à ajouter : nom différent mais mêmes adjacences
                     classes_adjacences_SUB[nom]=(*it_cl).second;
                     classeSUB++;
                  }
                  else{
                     if(affich_Classes)
                        cout<<"\tLCA2 > r1 => subdivision obligatoire"<<endl;//À PROUVER
                     //Il faut couper dans arbre 2 pour trouver les noeuds d'espèce LCA1
                     vector<int> Ids;
                     retourneIdSpeSAuDessousDeN(Ids,c.from_s(EspLCA1->toSTL()),Arbre2->getNode(id_LCA2),Arbre2);
                     if(Ids.size()==0)
                        cout<<"On n'a pas trouvé de noeud d'espèce LCA1 en dessous de LCA2, IMPOSSIBLE";
                     else if(Ids.size()==1)
                        cout<<"On a trouvé un seul noeud d'espèce LCA1 en dessous de LCA2, IMPOSSIBLE";
                     else{
                        if(affich_Classes)
                           cout<<"On doit subdiviser en "<<Ids.size()<<" classes."<<endl;
                        //Les nouveaux noms :
                        for(unsigned int k=0;k<Ids.size();k++){
                           string nom=c.to_s(num1)+sep+c.to_s(num2)+sep+c.to_s(id_LCA1)+sepAdj+c.to_s(Ids[k]);
                           if(affich_Classes)
                              cout<<"Nouveau nom de la classe sub : "<<nom<<endl;
                           //Une classe à ajouter, on l'initialise on y rangera les adj plus tard
                           vector<int> a;
                           classes_adjacences_SUB[nom]=a;
                           classeSUB++;
                        }
                        //On range les adjacences dans les nouvelles classes
                        //Les identifiants des extremités des adj se trouvent dans l'ordre dans Id_Genes1 et Id_Genes2
                        for(unsigned int a=0;a!=(*it_cl).second.size();a++){
                           unsigned int cpt=0;
                           while(!isAncestor(Ids[cpt],Id_Genes2[a],Arbre2)&&cpt<Ids.size())
                              cpt++;
                           if(cpt==Ids.size())
                              cout<<"Problème, on n'a pas trouvé de classe pour l'adj "<<a<<endl;
                           else{
                              string nom=c.to_s(num1)+sep+c.to_s(num2)+sep+c.to_s(id_LCA1)+sepAdj+c.to_s(Ids[cpt]);
                              classes_adjacences_SUB[nom].push_back((*it_cl).second[a]);
                           }
                        }
                     }
                  }
               }
          // delete(r1);
          // delete(r2);          
            }
       // delete(EspLCA1);
       // delete(EspLCA2);
         }
         classeGARD++;
      }
      // delete(Arbre1); //idem, ne pas effacer car segfault sinon
      // delete(Arbre2);
   }
   
   if(file_log)
      Offile_log<<" DONE"<<endl;
   
   // cout<<"\n\n\nAprès traitement des classes ("<<nb_tour<<" tours de boucle) :"<<endl;
   // cout<<"classes_adjacences.size()"<<classes_adjacences.size()<<endl;
   // cout<<"classes_adjacences_SUB.size()"<<classes_adjacences_SUB.size()<<endl;
   // cout<<"classes_adjacences_A_SUPPR.size()"<<classes_adjacences_A_SUPPR.size()<<endl;
   // cout<<"classe1="<<classe1<<" classeGARD="<<classeGARD<<" classeSUB="<<classeSUB<<" classeSUPPR="<<classeSUPPR<<endl;
   // nb_tour=0;
   
   //Supprimer les anciennes
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_class.cpp: Remove old adjacencies classes in map<string,vector<int> > classes_adjacences..."<<flush;
   for(it_class=classes_adjacences_A_SUPPR.begin();it_class!=classes_adjacences_A_SUPPR.end();it_class++)
      {
    nb_tour++;
    classes_adjacences.erase(*it_class);
      }
   if(file_log)
      Offile_log<<" DONE"<<endl;
   
   // cout<<"\n\n\nAprès suppression des classes ("<<nb_tour<<" tours de boucle) :"<<endl;
   // cout<<"classes_adjacences.size()"<<classes_adjacences.size()<<endl;
   // nb_tour=0;

   //Ajouter les classes créées ou renommées de taille >=2 :
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_class.cpp: Add new adjacencies classes (created or renamed where size>=2) in map<string,vector<int> > classes_adjacences..."<<flush;
   for(it_cl=classes_adjacences_SUB.begin();it_cl!=classes_adjacences_SUB.end();it_cl++)
      {
    nb_tour++;
    if ((*it_cl).second.size()>=2)
       classes_adjacences[(*it_cl).first]=(*it_cl).second;
      }
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;
   
   // cout<<"\n\n\nAprès ajout des classes ("<<nb_tour<<" tours de boucle) :"<<endl;
   // cout<<"classes_adjacences.size()"<<classes_adjacences.size()<<endl;
   
   int nb_classes=classes_adjacences.size();
   cout<<"Number of classes of size >1 to be treated: "<<nb_classes<<endl;
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_class.cpp: Number of classes of size >1 to be treated: "<<nb_classes<<endl;
   int taille_tot= 0;
   for(it_cl=classes_adjacences.begin();it_cl!=classes_adjacences.end();it_cl++){
      int taille=(*it_cl).second.size();
      taille_tot+=taille;
   }
   cout<<"Classes average size = "<<taille_tot/float(nb_classes)<<endl;
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_class.cpp: Classes average size = "<<taille_tot/float(nb_classes)<<endl<<endl;
}


void WriteAdjClassFile(map<string,vector<int> > &classes_adjacences){
   string adj_class_file=prefixe+output_adj_class;
   cout<<"Creation of the file for adjacencies classes: "<<adj_class_file<<" ..."<<flush;
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_class.cpp (WriteAdjClassFile): Creation of the file for adjacencies classes: "<<adj_class_file<<" ..."<<flush;
   ofstream OffilAdjClass(adj_class_file.c_str(), ios::out|ios::trunc);
   if(!OffilAdjClass){      
      cout<<"\nFrom Step2_adjacencies_class.cpp (WriteAdjClassFile): ERROR while opening file "<<adj_class_file<<endl;
      if(file_log)
         Offile_log<<"\nFI\tFrom Step2_adjacencies_class.cpp: ERROR while opening file "<<adj_class_file<<endl;
      exit(EXIT_FAILURE);
   }
   map<string, vector<int> >::iterator it_cl;
   for(it_cl=classes_adjacences.begin();it_cl!=classes_adjacences.end();it_cl++){
      OffilAdjClass<<(*it_cl).first<<"\t"<<(*it_cl).second.size();
      for(unsigned int a=0;a!=(*it_cl).second.size();a++){
         int num_Adj=(*it_cl).second[a];
         OffilAdjClass<<"\t"<<num_Adj;
      }      
      OffilAdjClass<<endl;
   }
   OffilAdjClass.close();
   // if(OffilAdjClass.fail())
   //    {
   //     cout<<"\nFrom Step2_adjacencies_class.cpp (WriteAdjClassFile): ERROR while closing file "<<adj_class_file<<endl;
   //     if(file_log)
   //        Offile_log<<"\nFrom Step2_adjacencies_class.cpp (WriteAdjClassFile): ERROR while closing file "<<adj_class_file<<endl;
   //     exit(EXIT_FAILURE);
   //    }
   cout<<" DONE"<<endl<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;
}


void WriteAdjClassFileHumanRead(map<string,vector<int> > &classes_adjacences, vector<adjacence> adjacencies, map<string,string> gene_species_EXT){
   string adj_class_file_read=prefixe+output_adj_class_read;
   cout<<"Creation of the file for adjacencies classes (Human readable): "<<adj_class_file_read<<" ...";
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_class.cpp (WriteAdjClassFileHumanRead): Creation of the file for adjacencies classes (Human readable): "<<adj_class_file_read<<" ..."<<flush;

   ofstream OffilAdjClassRead(adj_class_file_read.c_str(), ios::out|ios::trunc);
   if(!OffilAdjClassRead){      
      cout<<"\nFrom Step2_adjacencies_class.cpp (WriteAdjClassFileHumanRead): ERROR while opening file "<<adj_class_file_read<<endl;
      if(file_log)
         Offile_log<<"\nFI\tFrom Step2_adjacencies_class.cpp (WriteAdjClassFileHumanRead): ERROR while opening file "<<adj_class_file_read<<endl;
      exit(EXIT_FAILURE);
   }
   
   map<string, vector<int> >::iterator it_cl;
   for(it_cl=classes_adjacences.begin();it_cl!=classes_adjacences.end();it_cl++){
      OffilAdjClassRead<<"Adjacencies list contained in adjacencies class "<<(*it_cl).first<<" of size "<<(*it_cl).second.size()<<" :"<<"\n\n";
      for(unsigned int a=0;a!=(*it_cl).second.size();a++){
         int va=(*it_cl).second[a];
         string adj=adjacencies[va].gene1+sep+gene_species_EXT[adjacencies[va].gene1]+sepAdj+adjacencies[va].gene2+sep+gene_species_EXT[adjacencies[va].gene2];
         OffilAdjClassRead<<"Adjacency "<<a<<" :"<<endl;
         OffilAdjClassRead<<adj<<"\n\n";
      }      
      OffilAdjClassRead<<"\n"<<"***********************************************************"<<"\n";
   }
   OffilAdjClassRead.close();
   // if(OffilAdjClassRead.fail())
   //    {
   //     cout<<"\nFrom Step2_adjacencies_class.cpp (WriteAdjClassFile): ERROR while closing file "<<adj_class_file_read<<endl;
   //     if(file_log)
   //        Offile_log<<"\nFrom Step2_adjacencies_class.cpp (WriteAdjClassFile): ERROR while closing file "<<adj_class_file_read<<endl;
   //     exit(EXIT_FAILURE);
   //    }
   cout<<" DONE"<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;
}


int main(int argc, char* argv[]){
   time_t tbegin=time(NULL);              // get the current calendar time
   cout<<"\n\t########################"<<endl;
   cout<<"\t###  Start Step2 !!  ###"<<endl;
   cout<<"\t########################"<<endl<<endl;
   
   lireFicConfig(fic_arbre,fic_gene,fic_especes,fic_adjacence,exp_name,directory,argc,argv);
   //On regarde si un fichier d'adjacences nettoyées existe, et si oui, on l'utilise
   if (is_readable(fic_adjacence_clean))
      fic_adjacence=fic_adjacence_clean;
   //Idem pour le fichier de gènes
   if (is_readable(fic_gene_clean))
      fic_gene=fic_gene_clean;
   
   cout<<endl<<"File(s) used in Step2:"<<endl;
   cout<<"\t- Reconciled trees file: "<<tree_reconciled_file<<endl;
   cout<<"\t- Adjacencies file: "<<fic_adjacence<<endl<<endl;
   
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
            cout<<endl<<"\nFrom Step2_adjacencies_class.cpp (main): ERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file open by Step2_adjacencies_class.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      else{
         Offile_log.open(log_file.c_str(), ios::out|ios::trunc);
         if(!Offile_log){
       cout<<endl<<"\nFrom Step2_adjacencies_class.cpp (main): ERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file created by Step2_adjacencies_class.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      Offile_log<<"\t########################"<<endl;
      Offile_log<<"\t###  Start Step2 !!  ###"<<endl;
      Offile_log<<"\t########################"<<endl<<endl;
      Offile_log<<"\tFile(s) used in Step2:"<<endl;
      Offile_log<<"\t\t- Reconciled trees file: "<<tree_reconciled_file<<endl;
      Offile_log<<"\t\t- Adjacencies file: "<<fic_adjacence<<endl<<endl;
   }
   
   // Variables
   //   TreeTemplate<Node> S;
   //   S = new TreeTemplate<Node>;
   map<string,string> gene_species_EXT;
   vector<Tree*> Arbres;
   map<string,vector<int> > gene_GeneTreeID;
   map<string,vector<int> > classes_adjacences;
   vector<adjacence> adjacencies;
   
   
   /////////////////////////////////
   /////   START Genes Clean   /////
   /////////////////////////////////
   

   ////////////////////////////////////////////////////////////////////////
   // Store species tree in TreeTemplate<Node> *S   
   cout<<"Store Species tree in TreeTemplate<Node> * S..."<<flush;
   if(file_log)
      Offile_log<<"OS\tFrom Step2_adjacencies_classes.cpp (main): Store Species tree in TreeTemplate<Node> * S..."<<flush;
   Newick * newickReaderEsp = new Newick(true,true);
   newickReaderEsp->enableExtendedBootstrapProperty(esp);
   TreeTemplate<Node> *S = newickReaderEsp->read(fic_especes);
   delete(newickReaderEsp);
   cout<<" DONE"<<endl<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;
   
   //Store arbres réconciliés dans un vecteur Arbres
   StoreReconciledTrees(Arbres);
   
   // A/ Correspondance entre adjacences et arbres de gène : création des classes 
   //      1) Parcours des arbres pour construire la correspondance g_x apparaît dans arbres y, z, t, ...
   
   //Associate Species_name with gene in map<string, string> gene_species_EXT
   AssociateGeneWithSpecies(gene_species_EXT);
   //   afficheMap(gene_species_EXT);
   
   //Associe les gènes avec les numéros d'arbres réconciliés dans lesquels ils sont présents
   AssociateGeneWithReconciledTreeNb(Arbres,gene_GeneTreeID,gene_species_EXT);
   //   afficheMap(gene_GeneTreeID);
   
   //Nettoyage de la map gene_species_EXT pour ne garder que les gènes présents dans nos arbres
   int elim = 0;
   elim=GeneSpeciesClean(gene_GeneTreeID,gene_species_EXT);
   if(elim>0){
      cout<<"\nFrom Step2_adjacencies_class (main): WARNING while cleaning "<<elim<<" genes have been deleted"<<endl<<endl;
   }

///////////////////////////////
/////   END Genes Clean   /////
///////////////////////////////



   //      2) Parcours des adjacences pour construire les classes
   //NB : si un ou les deux gènes de l'adjacence apparaît dans plusieurs arbres, alors l'adjacence apparaît dans plusieurs classes
   //Élimination des classes de taille 1
   //+ ATTENTION BIS subdivision des classes entre deux arbres différents pour que les extrémités des 
   //adjacences ait la même espèce ancestre dans les arbres
   //Nom de la classe x|x|n°espèce_ancestre_commun
   //+++ ATTENTION TER modif de la définition des classes, on ne subdivise plus les classes entre deux 
   //arbres, on regarde simplement qu'il y ait bien un noeud ancêtre de chaque extremité de 
   //même espèce dans chacun des arbres
   //==> l'élagage d'un ou des deux arbres suffit et c'est fait après dans le traitement de chaque classe
   //==> on commente tout le bazar que ne sert plus ...
   //++++ ATTENTION QUATRO après vérif on peut encore subdiviser ces classes :
   //soient r1 et r2 les espèces des racines des arbres 1 et 2
   //soient LCA1 et LCA2 les espèces des plus petits ancres communs des extrémités des adj dans arbre 1 et 2
   //si LCA1 = LCA2, on prend les sous arbres enracinés en ces noeuds
   //sinon, supposons LCA1 > LCA2
   //         soit LCA1<=r2, alors on remonte dans arbre2 jusqu'à un noeud d'espèce LC1 (tjs une seule classe)
   //         soit LCA1>r2, alors il faut couper dans arbre 1 au niveau des espèce de LCA1 (subdivision possible)
   //                             puis répartir les adj entre les classes subdivisées si nécessaire
   //Codage des noms de classes entre 2 arbres x et y (x!=y)
   //x|y|id1-id2 id1 et id2 sont les racines des sous-arbres des arbres 1 et 2 constituant la classe



   // Adjacencies class pre-construction
   AdjClassPreConstruction(S,Arbres,gene_GeneTreeID,gene_species_EXT,classes_adjacences,adjacencies);
   gene_GeneTreeID.clear();
//   afficheMap(classes_adjacences);

   //Adjacencies class construction
   AdjClassConstruction(S,Arbres,gene_species_EXT,classes_adjacences,adjacencies);
   //On libère la mémoire du vecteur d'arbres
   vector<Tree*>::iterator ita;
   for (ita=Arbres.begin();ita!=Arbres.end();ita++)
      delete(*ita);
   Arbres.clear();
   
   //Création du fichier de classes d'adjacences:
   WriteAdjClassFile(classes_adjacences);

   //   Version lisible   //
//   WriteAdjClassFileHumanRead(classes_adjacences,adjacencies,gene_species_EXT);
   classes_adjacences.clear();
   adjacencies.clear();


   time_t tend=time(NULL);                // get the current calendar time
   // Compute execution time
   float texec=difftime(tend,tbegin);    // tend-tbegin (result in second)
   
   cout<<"Execution time of Step2 is: "<<texec<<"s"<<endl<<endl;
   cout<<"\t######################"<<endl;
   cout<<"\t###  End Step2 !!  ###"<<endl;
   cout<<"\t######################"<<endl<<endl;

   if (file_log){
      Offile_log<<"Execution time of Step2 is: "<<texec<<"s"<<endl<<endl;
      Offile_log<<"\t######################"<<endl;
      Offile_log<<"\t###  End Step2 !!  ###"<<endl;
      Offile_log<<"\t######################"<<endl<<endl;
      Offile_log.close();
      if(Offile_log.fail()){
         cout<<"\nFrom Step2_adjacencies_class.cpp (main): ERROR while closing file "<<log_file<<endl<<endl;
         exit(EXIT_FAILURE);
      }
   }

   return(0);
}
