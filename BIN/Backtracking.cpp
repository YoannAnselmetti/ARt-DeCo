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

File: Backtracking.cpp               Last modified on: 20/11/2014
Created by: Sèverine Bérard            Created on: 10/01/2011
--------------------------------------------------------------------------
Specification: 
File containing the functions doing the backtracking part of the function
DECO in Step3_DECO.cpp
=========================================================================*/
#include "Backtracking.h"

Node * parcoursArriereCassure(Node * v1, Node * v2, bool v1dansA1){
   nb_Bk++;
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereCassure"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   convert<int> c;
   Node * n = new Node();
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   //Pour compter le nb de break par branche :
   if (espece_nb_break.find(c.from_s(Espece->toSTL()))!=espece_nb_break.end())
      espece_nb_break[c.from_s(Espece->toSTL())]++;
   else
      espece_nb_break[c.from_s(Espece->toSTL())]=1;
   
   n->setId(id_arbres_adj);
   id_arbres_adj++;
   n->setNodeProperty(esp,*Espece);
   n->setBranchProperty(typ,BppString(bk));
   n->setBranchProperty("D", BppString("?"));   
   if (v1dansA1)
      n->setName(NomCass+sep+c.to_s(v1->getId())+sepAdj+c.to_s(v2->getId()));
   else
      n->setName(NomCass+sep+c.to_s(v2->getId())+sepAdj+c.to_s(v1->getId()));
   return n;
}

/////////////////////////////////////////
//////      Cas pour Calcul C1      /////
/////////////////////////////////////////

// Cas 1 - E(v1)=Extant / E(v2)=Extant
Node * parcoursArriereC1ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe, vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC1ExtantExtant"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   Node * n = new Node();
   n->setId(id_arbres_adj);
   id_arbres_adj++;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   n->setNodeProperty(esp,*Espece);
   
   string nomV1 = v1->getName();
   string nomV2 = v2->getName();
   string nom_gene1=nomV1.substr(0,nomV1.rfind(sep,nomV1.length()-1));
   string nom_gene2=nomV2.substr(0,nomV2.rfind(sep,nomV2.length()-1));
   string espece=nomV2.substr(nomV2.rfind(sep,nomV2.length())+1,nomV2.length());
   //On cherche d'adj dans les deux sens g1-g2 ou g2-g1
   adjacence a12,a21;
   a12.gene1=nom_gene1;
   a12.gene2=nom_gene2;
   a21.gene1=nom_gene2;
   a21.gene2=nom_gene1;
   if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe)){
      n->setBranchProperty(typ,BppString(ga));
      n->setBranchProperty("D", BppString("?"));
      if (v1dansA1)
         n->setName(nom_gene1+sepAdj+nom_gene2+sep+espece);
      else
         n->setName(nom_gene2+sepAdj+nom_gene1+sep+espece);
   }
   else{
      n->setBranchProperty(typ,BppString("IMPOSSIBLE"));
      n->setBranchProperty("D", BppString("?"));
      n->setName("IMPOSSIBLE"+sep+espece);
   }
   return n;
}

// Cas 2 - E(v1)=GLos / E(v2)!=GLos
Node * parcoursArriereC1GLosAny(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC1GLosAny"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   convert<int> c;
   Node * n = new Node();
   n->setId(id_arbres_adj);
   id_arbres_adj++;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   n->setNodeProperty(esp,*Espece);
   n->setBranchProperty(typ,BppString(per));
   n->setBranchProperty("D", BppString("?"));
   if (v1dansA1)
      n->setName(NomPer+sep+c.to_s(v1->getId())+sepAdj+c.to_s(v2->getId()));
   else
      n->setName(NomPer+sep+c.to_s(v2->getId())+sepAdj+c.to_s(v1->getId()));
   return n;
}

// Cas 3 - E(v1)=GLos / E(v2)=GLos
Node * parcoursArriereC1GLosGLos(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC1GLosGLos"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   convert<int> c;
   Node * n = new Node();
   n->setId(id_arbres_adj);
   id_arbres_adj++;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   n->setNodeProperty(esp,*Espece);
   n->setBranchProperty(typ,BppString(Aper));
   n->setBranchProperty("D", BppString("?"));
   string nom;
   if(v1dansA1)
      nom=NomPerA+sep+c.to_s(v1->getId())+sepAdj+c.to_s(v2->getId());
   else
      nom=NomPerA+sep+c.to_s(v2->getId())+sepAdj+c.to_s(v1->getId());
   n->setName(nom);
   return n;
}

// Cas 4.1 - E(v1)=GDup / E(v2)=Spec
Node * parcoursArriereC1GDupExtantOrSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   nb_GDup++;
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC1GDupExtantOrSpec"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   convert<int> c;
   Node * n = new Node();
   n->setId(id_arbres_adj);
   id_arbres_adj++;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   n->setNodeProperty(esp,*Espece);
   n->setBranchProperty(typ,BppString(dupl));
   n->setBranchProperty("D", BppString("T"));
   if(v1dansA1)
      n->setName(c.to_s(v1->getId())+sepAdj+c.to_s(v2->getId()));
   else
      n->setName(c.to_s(v2->getId())+sepAdj+c.to_s(v1->getId()));
   float C = recupCoutC1(v1,v2);
   // (1)
   if(C==recupCoutC1(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)){    
      n->addSon(0,parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   // (2)
   else if(C==recupCoutC0(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)){    
      n->addSon(0,parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   // (3)
   else if(C==recupCoutC1(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)+Crea){    
      n->addSon(0,parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;

      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
   }
   // (4)
   else if(C==recupCoutC0(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)+Break){    
      n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2, v1dansA1));
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else 
      {
	 cout<<"From Backtracking.cpp (parcoursArriereC1GDupExtantOrSpec): ERROR none of the 4 costs correspond to the stored cost"<<endl;
	 exit(EXIT_FAILURE);
      }
   return n;
}

// Cas 5 - E(v1)=Spec / E(v2)=Spec
Node * parcoursArriereC1SpecSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC1SpecSpec"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   Node * fgV1=v1->getSon(0);
   Node * fdV1=v1->getSon(1);
   Node * fgV2;
   Node * fdV2;
   if (Espece(fgV1)==Espece(v2->getSon(0))){
      fgV2=v2->getSon(0);
      fdV2=v2->getSon(1);
   }
   else{
      fgV2=v2->getSon(1);
      fdV2=v2->getSon(0);
   }
   
   convert<int> c;
   Node * n = new Node();
   n->setId(id_arbres_adj);
   id_arbres_adj++;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   n->setNodeProperty(esp,*Espece);
   n->setBranchProperty(typ,BppString(spe));
   n->setBranchProperty("D", BppString("F"));
   if(v1dansA1)
      n->setName(c.to_s(v1->getId())+sepAdj+c.to_s(v2->getId()));
   else
      n->setName(c.to_s(v2->getId())+sepAdj+c.to_s(v1->getId()));
   float C = recupCoutC1(v1,v2);
   // (1)
   if(C==recupCoutC1(fgV1,fgV2) + recupCoutC1(fdV1,fdV2)){
      n->addSon(0,parcoursArriereC1(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      n->addSon(1,parcoursArriereC1(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1));
   }
   // (2)
   else if(C==recupCoutC1(fgV1,fgV2) + recupCoutC0(fdV1,fdV2) + Break){
      n->addSon(0,parcoursArriereC1(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      n->addSon(1,parcoursArriereCassure(fdV1,fdV2, v1dansA1));
      parcoursArriereC0(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   // (3)
   else if(C==recupCoutC0(fgV1,fgV2) + recupCoutC1(fdV1,fdV2) + Break){
      n->addSon(0,parcoursArriereC1(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      n->addSon(1,parcoursArriereCassure(fgV1,fgV2, v1dansA1));
      parcoursArriereC0(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   // (4)
   else if(C==recupCoutC0(fgV1,fgV2) + recupCoutC0(fdV1,fdV2) + 2*Break){
      n->addSon(0,parcoursArriereCassure(fgV1,fgV2, v1dansA1));
      n->addSon(1,parcoursArriereCassure(fdV1,fdV2, v1dansA1));
      parcoursArriereC0(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else 
      {
	 cout<<"From Backtracking.cpp (parcoursArriereC1SpecSpec): ERROR none of the 4 costs correspond to the stored cost"<<endl;
	 exit(EXIT_FAILURE);
      }
   return n;
}

// Cas 6 - E(v1)=GDup / E(v2)=GDup
Node * parcoursArriereC1GDupGDup(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC1GDupGDup"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   convert<int> c;
   Node * n = new Node();
   n->setId(id_arbres_adj);
   id_arbres_adj++;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   n->setNodeProperty(esp,*Espece);
   if(v1dansA1)
      n->setName(c.to_s(v1->getId())+sepAdj+c.to_s(v2->getId()));
   else
      n->setName(c.to_s(v2->getId())+sepAdj+c.to_s(v1->getId()));
   n->setBranchProperty("D", BppString("T"));
   
   float C = recupCoutC1(v1,v2);
   //Les cas "D1"
   if(C==recupCoutC1(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC0(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC1(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)+Crea){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
   }
   else if(C==recupCoutC0(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)+Break){
      nb_GDup++;
      cout<<"\t\t\tJe suis dans parcoursArriereC1GDupGDup D1, je passe dans le cas C0 + C0 + Break"<<endl;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2, v1dansA1));
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   //Les cas "D2"
   else if(C==recupCoutC1(v1,v2->getSon(0))+recupCoutC0(v1,v2->getSon(1))){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC0(v1,v2->getSon(0))+recupCoutC1(v1,v2->getSon(1))){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC1(v1,v2->getSon(0))+recupCoutC1(v1,v2->getSon(1))+Crea){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
   }
   else if(C==recupCoutC0(v1,v2->getSon(0))+recupCoutC0(v1,v2->getSon(1))+Break){
    nb_GDup++;
      cout<<"\t\t\tJe suis dans parcoursArriereC1GDupGDup D2, je passe dans le cas C0 + C0 + Break"<<endl;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereCassure(v1,v2->getSon(0), v1dansA1));
      parcoursArriereC0(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   //Les cas "D12"
   else{
      nb_ADup++;
      n->setBranchProperty(typ,BppString(Adupl));
      //Pour tracer les gènes qui se dupliquent ensemble :
      adjacence a;
      if(v1dansA1){
         a.gene1=c.to_s(v1->getId());
         a.gene2=c.to_s(v2->getId());
      }
      else{
         a.gene1=c.to_s(v2->getId());
         a.gene2=c.to_s(v1->getId());
      }
      DUP.push_back(a);
      //(1)
      if(C==recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      //(2)
      else if(C==recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Crea){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         //Création
         nb_Crea++;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]++;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=1;
            TreeTemplate<Node> * T = new TreeTemplate<Node>();
            T->setRootNode(parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
            ArbresDAdjacences->push_back(T);
         }
      //(3)
      else if(C==recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Crea){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         //Création
         nb_Crea++;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]++;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=1;
         TreeTemplate<Node> * T = new TreeTemplate<Node>();
         T->setRootNode(parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      //(15)
      else if(C==recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+2*Crea){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         //Création
         nb_Crea+=2;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]+=2;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=2;
         TreeTemplate<Node> * T1 = new TreeTemplate<Node>();
         T1->setRootNode(parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T1);
         //Création
         TreeTemplate<Node> * T2 = new TreeTemplate<Node>();
         T2->setRootNode(parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T2);
      }
      // (4)
      else if(C==recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         //Cassure
         n->addSon(1,parcoursArriereCassure(v1->getSon(1),v2->getSon(1), v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      // (13)
      else if(C==recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Break+Crea){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         //Cassure
         n->addSon(1,parcoursArriereCassure(v1->getSon(1),v2->getSon(1), v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         //Création
         nb_Crea++;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]++;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=1;
         TreeTemplate<Node> * T = new TreeTemplate<Node>();
         T->setRootNode(parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T);
      }
      //(14)
      else if(C==recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break+Crea){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         //Cassure
         n->addSon(1,parcoursArriereCassure(v1->getSon(1),v2->getSon(1), v1dansA1));
         //Création
         nb_Crea++;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]++;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=1;
         TreeTemplate<Node> * T = new TreeTemplate<Node>();
         T->setRootNode(parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      //(5) 
      else if(C==recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break){
         //Cassure    
         n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2->getSon(0), v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      // (11)   
      else if(C==recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Break+Crea){
         //Cassure    
         n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2->getSon(0), v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         //Création
         nb_Crea++;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]++;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=1;
         TreeTemplate<Node> * T = new TreeTemplate<Node>();
         T->setRootNode(parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T);
      }
      //(12)
      else if(C==recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break+Crea){
         //Cassure    
         n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2->getSon(0), v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         //Création
         nb_Crea++;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]++;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=1;
         TreeTemplate<Node> * T = new TreeTemplate<Node>();
         T->setRootNode(parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      // (6)
      else if(C==recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      }
      //(7)
      else if(C==recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Crea){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         //Création
         nb_Crea++;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]++;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=1;
         TreeTemplate<Node> * T = new TreeTemplate<Node>();
         T->setRootNode(parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T);
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      }
      //(8)
      else if(C==recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Crea){
         //Création
         nb_Crea++;
         //Pour compter le nb de gain par branche :
         if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
            espece_nb_gain[c.from_s(Espece->toSTL())]++;
         else
            espece_nb_gain[c.from_s(Espece->toSTL())]=1;
         TreeTemplate<Node> * T = new TreeTemplate<Node>();
         T->setRootNode(parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         ArbresDAdjacences->push_back(T);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      }
      //(10)
      else if(C==recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereCassure(v1->getSon(1),v2->getSon(0), v1dansA1));
      }
      // (9)
      else if(C==recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Break){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2->getSon(1), v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      }
      //(16) ajouté
      else if(C==recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+2*Break){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2->getSon(1), v1dansA1));
         n->addSon(1,parcoursArriereCassure(v1->getSon(1),v2->getSon(0), v1dansA1));
      }
      else 
         {
	    cout<<"From Backtracking.cpp (parcoursArriereC1GDupGDup): ERROR none of the 24 costs correspond to the stored cost"<<endl;
	    exit(EXIT_FAILURE);
	 }
   }
   return n;
}

/////////////////////////////////////////
//////      Cas pour Calcul C0      /////
/////////////////////////////////////////

// Cas 1 - E(v1)=Extant / E(v2)=Extant
void parcoursArriereC0ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe, vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0ExtantExtant"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   
   string nomV1 = v1->getName();
   string nomV2 = v2->getName();
   string nom_gene1=nomV1.substr(0,nomV1.rfind(sep,nomV1.length()-1));
   string nom_gene2=nomV2.substr(0,nomV2.rfind(sep,nomV2.length()-1)); 
   string espece=nomV2.substr(nomV2.rfind(sep,nomV2.length())+1,nomV2.length());
   //On cherche d'adj dans les deux sens g1-g2 ou g2-g1
   adjacence a12,a21;
   a12.gene1=nom_gene1;
   a12.gene2=nom_gene2;
   a21.gene1=nom_gene2;
   a21.gene2=nom_gene1;
   
   if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe)){
      //IMPOSSIBLE AVEC COÛT INFINI SAUF SI ARBRE RÉDUIT À UNE FEUILLE
      cout<<"From Backtracking.cpp (parcoursArriereC0ExtantExtant): IMPOSSIBLE CASE !!!"<<endl;
      exit(EXIT_FAILURE);
   }
   else{
      //On ne fait rien
   }
}

void parcoursArriereC0GLosAny(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   //On ne fait rien
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0GLosAny"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
}

void parcoursArriereC0GDupExtantOrSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   nb_GDup++;
   convert<int> c;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0GDupExtantOrSpec"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   float C = recupCoutC0(v1,v2);
   if(C==recupCoutC0(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)){
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC1(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)+Crea){
      nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC0(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)+Crea){
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
   }
   else if(recupCoutC1(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)+2*Crea){
      nb_Crea+=2;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]+=2;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=2;
      TreeTemplate<Node> * T1 = new TreeTemplate<Node>();
      T1->setRootNode(parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T1);
      TreeTemplate<Node> * T2 = new TreeTemplate<Node>();
      T2->setRootNode(parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T2);
   }
   else
      {
	 cout<<"From Backtracking.cpp (parcoursArriereC0GDupExtantOrSpec): ERROR none of the 4 costs correspond to the stored cost"<<endl;
	 exit(EXIT_FAILURE);
      }
}

void parcoursArriereC0GLosGLos(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   //On ne fait rien
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0GLosGLos"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
}

void parcoursArriereC0SpecSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0SpecSpec"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   Node * fgV1=v1->getSon(0);
   Node * fdV1=v1->getSon(1);
   Node * fgV2;
   Node * fdV2;
   convert<int> c;
   BppString * Especefg = dynamic_cast<BppString*> (fgV1->getNodeProperty(esp));
   BppString * Especefd = dynamic_cast<BppString*> (fdV1->getNodeProperty(esp));
   if (Espece(fgV1)==Espece(v2->getSon(0))){
      fgV2=v2->getSon(0);
      fdV2=v2->getSon(1);
   }
   else{
      fgV2=v2->getSon(1);
      fdV2=v2->getSon(0);
   }
   
   float C = recupCoutC0(v1,v2);
   if(C==recupCoutC0(fgV1,fgV2) + recupCoutC0(fdV1,fdV2)){
      parcoursArriereC0(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else  if(C==recupCoutC1(fgV1,fgV2) + recupCoutC0(fdV1,fdV2) + Crea){
      nb_Crea++;
      //Pour compter le nb de gain par branche fg :
      if (espece_nb_gain.find(c.from_s(Especefg->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Especefg->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Especefg->toSTL())]=1;
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
      parcoursArriereC0(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else  if(C==recupCoutC0(fgV1,fgV2) + recupCoutC1(fdV1,fdV2) + Crea){
      parcoursArriereC0(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1);
      nb_Crea++;
      //Pour compter le nb de gain par branche fd:
      if (espece_nb_gain.find(c.from_s(Especefd->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Especefd->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Especefd->toSTL())]=1;
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
   }
   else  if(C==recupCoutC1(fgV1,fgV2) + recupCoutC1(fdV1,fdV2) + 2*Crea){
      nb_Crea+=2;
      //Pour compter le nb de gain par branche fg :
      if (espece_nb_gain.find(c.from_s(Especefg->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Especefg->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Especefg->toSTL())]=1;
      //Pour compter le nb de gain par branche fd:
      if (espece_nb_gain.find(c.from_s(Especefd->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Especefd->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Especefd->toSTL())]=1;
      TreeTemplate<Node> * T1 = new TreeTemplate<Node>();
      T1->setRootNode(parcoursArriereC1(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T1);
      TreeTemplate<Node> * T2 = new TreeTemplate<Node>();
      T2->setRootNode(parcoursArriereC1(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T2);
   }
   else
      {
	 cout<<"From Backtracking.cpp (parcoursArriereC0SpecSpec): ERROR none of the 4 costs correspond to the stored cost"<<endl;
	 exit(EXIT_FAILURE);
      }
}

void parcoursArriereC0GDupGDup(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0GDupGDup"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<endl;
   float C = recupCoutC0(v1,v2);
   convert<int> c;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   
   //Les 2 cas originaux d'Éric
   if(C==recupCoutC0(v1->getSon(0),v2) + recupCoutC0(v1->getSon(1),v2)){
      nb_GDup++;
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC0(v1,v2->getSon(0)) + recupCoutC0(v1,v2->getSon(1))){
      nb_GDup++;
      parcoursArriereC0(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   //3 cas où v1 se duplique avant v2
   else if(C==recupCoutC0(v1->getSon(0),v2) + recupCoutC1(v1->getSon(1),v2)+Crea){
      nb_GDup++;nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
   }
   else if(C==recupCoutC1(v1->getSon(0),v2) + recupCoutC0(v1->getSon(1),v2)+Crea){
      nb_GDup++;nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC1(v1->getSon(0),v2) + recupCoutC1(v1->getSon(1),v2)+2*Crea){
      nb_GDup++;nb_Crea+=2;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]+=2;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=2;
      TreeTemplate<Node> * T1 = new TreeTemplate<Node>();
      T1->setRootNode(parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T1);
      TreeTemplate<Node> * T2 = new TreeTemplate<Node>();
      T2->setRootNode(parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T2);
   }
   //3 cas où v2 se duplique avant v1
   else if(C==recupCoutC0(v1,v2->getSon(0)) + recupCoutC1(v1,v2->getSon(1))+Crea){
      nb_GDup++;nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;
      parcoursArriereC0(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
   }
   else if(C==recupCoutC1(v1,v2->getSon(0)) + recupCoutC0(v1,v2->getSon(1))+Crea){
      nb_GDup++;nb_Crea++;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]++;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=1;
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T);
      parcoursArriereC0(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==recupCoutC1(v1,v2->getSon(0)) + recupCoutC1(v1,v2->getSon(1))+2*Crea){
      nb_GDup++;nb_Crea+=2;
      //Pour compter le nb de gain par branche :
      if (espece_nb_gain.find(c.from_s(Espece->toSTL()))!=espece_nb_gain.end())
         espece_nb_gain[c.from_s(Espece->toSTL())]+=2;
      else
         espece_nb_gain[c.from_s(Espece->toSTL())]=2;
      TreeTemplate<Node> * T1 = new TreeTemplate<Node>();
      T1->setRootNode(parcoursArriereC1(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T1);
      TreeTemplate<Node> * T2 = new TreeTemplate<Node>();
      T2->setRootNode(parcoursArriereC1(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
      ArbresDAdjacences->push_back(T2);
   }
   else
      {
	 cout<<"From Backtracking.cpp (parcoursArriereC1GDupGDup): ERROR none of the 8 costs correspond to the stored cost"<<endl;  
	 exit(EXIT_FAILURE);
      }
}

// parcoursArriereC1
Node * parcoursArriereC1(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   switch (EvtToInt(E(v1))){  
      case 1 :
         switch (EvtToInt(E(v2))){
            case 1 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Extant et v2 Extant"<<endl;
               return parcoursArriereC1ExtantExtant(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;
            case 4 :    
               if(affich_ParcoursArriere)
                  cout<<"v1 Extant et v2 GLos"<<endl;
               return parcoursArriereC1GLosAny(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);
               break;    
            case 2 :
               cout<<"From Backtracking.cpp (parcoursArriereC1): v1 Extant and v2 Spec IMPOSSIBLE"<<endl;
	       exit(EXIT_FAILURE);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Extant et v2 GDup"<<endl;
               return parcoursArriereC1GDupExtantOrSpec(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);    
               break;
            default: 
               cout<<"From Backtracking.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
	       exit(EXIT_FAILURE);
               break;  
         }   
         break;
      case 4 :
         switch (EvtToInt(E(v2))){
            case 1 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GLos et v2 Extant"<<endl;
               return parcoursArriereC1GLosAny(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;
            case 4 :    
               if(affich_ParcoursArriere)
                  cout<<"v1 GLos et v2 Glos"<<endl;
               return parcoursArriereC1GLosGLos(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;    
            case 2 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GLos et v2 Spec"<<endl;
               return parcoursArriereC1GLosAny(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GLos et v2 GDup"<<endl;
               return parcoursArriereC1GLosAny(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);    
               break;
            default: 
               cout<<"From Backtracking.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
	       exit(EXIT_FAILURE);
            break;  
         }
         break;
      case 2 :
         switch (EvtToInt(E(v2))){
            case 1 :
               cout<<"From Backtracking.cpp (parcoursArriereC1): v1 Spec and v2 Extant IMPOSSIBLE"<<endl;
	       exit(EXIT_FAILURE);
               break;
            case 4 :    
               if(affich_ParcoursArriere)
                  cout<<"v1 Spec et v2 Glos"<<endl;
               return parcoursArriereC1GLosAny(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);
               break;    
            case 2 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Spec et v2 Spec"<<endl;
               return parcoursArriereC1SpecSpec(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Spec et v2 GDup"<<endl;
               return parcoursArriereC1GDupExtantOrSpec(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);
               break;
            default: 
               cout<<"From Backtracking.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
	       exit(EXIT_FAILURE);
            break;  
         }
         break;
      case 3 :
         switch (EvtToInt(E(v2))){
            case 1 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GDup et v2 Extant"<<endl;
               return parcoursArriereC1GDupExtantOrSpec(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;
            case 4 :    
               if(affich_ParcoursArriere)
                  cout<<"v1 GDup et v2 Glos"<<endl;
               return parcoursArriereC1GLosAny(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);
               break;    
            case 2 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GDup et v2 Spec"<<endl;
               return parcoursArriereC1GDupExtantOrSpec(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GDup et v2 GDup"<<endl;
               return parcoursArriereC1GDupGDup(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);    
               break;
            default: 
               cout<<"From Backtracking.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
	       exit(EXIT_FAILURE);
            break;  
         }
         break;
      default: 
         cout<<"From Backtracking.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v1)) on node v1 with id "<<v1->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v1))<<endl;
	 exit(EXIT_FAILURE);
      break;  
   }
}

// parcoursArriereC0
void parcoursArriereC0(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   switch (EvtToInt(E(v1))){  
      case 1 :
         switch (EvtToInt(E(v2))){
            case 1 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Extant et v2 Extant"<<endl;
               parcoursArriereC0ExtantExtant(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;
            case 4 :    
               if(affich_ParcoursArriere)
                  cout<<"v1 Extant et v2 GLos"<<endl;
               parcoursArriereC0GLosAny(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);
               break;    
            case 2 :
               cout<<"From Backtracking.cpp (parcoursArriereC0): v1 Extant and v2 Spec IMPOSSIBLE"<<endl;
	       exit(EXIT_FAILURE);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Extant et v2 GDup"<<endl;
               parcoursArriereC0GDupExtantOrSpec(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);    
               break;
            default: 
               cout<<"From Backtracking.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
	       exit(EXIT_FAILURE);
            break;  
         }   
         break;
      case 4 :
         switch (EvtToInt(E(v2))){
            case 1 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GLos et v2 Extant"<<endl;
               parcoursArriereC0GLosAny(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;
            case 4 :    
               if(affich_ParcoursArriere)
                  cout<<"v1 GLos et v2 Glos"<<endl;
               parcoursArriereC0GLosGLos(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;    
            case 2 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GLos et v2 Spec"<<endl;
               parcoursArriereC0GLosAny(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GLos et v2 GDup"<<endl;
               parcoursArriereC0GLosAny(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);    
               break;
            default: 
               cout<<"From Backtracking.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
	       exit(EXIT_FAILURE);
            break;  
         }
         break;
      case 2 :
         switch (EvtToInt(E(v2))){
            case 1 :
               cout<<"From Backtracking.cpp (parcoursArriereC0): v1 Spec and v2 Extant IMPOSSIBLE"<<endl;
	       exit(EXIT_FAILURE);
               break;
            case 4 :    
               if(affich_ParcoursArriere)
                  cout<<"v1 Spec et v2 Glos"<<endl;
               parcoursArriereC0GLosAny(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);
               break;    
            case 2 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Spec et v2 Spec"<<endl;
               parcoursArriereC0SpecSpec(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Spec et v2 GDup"<<endl;
               parcoursArriereC0GDupExtantOrSpec(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);
               break;
            default: 
               cout<<"From Backtracking.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
	       exit(EXIT_FAILURE);
            break;  
         }
         break;
      case 3 :
         switch (EvtToInt(E(v2))){
            case 1 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GDup et v2 Extant"<<endl;
               parcoursArriereC0GDupExtantOrSpec(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;
            case 4 :    
               if(affich_ParcoursArriere)
                  cout<<"v1 GDup et v2 Glos"<<endl;
               parcoursArriereC0GLosAny(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);
               break;    
            case 2 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GDup et v2 Spec"<<endl;
               parcoursArriereC0GDupExtantOrSpec(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 GDup et v2 GDup"<<endl;
               parcoursArriereC0GDupGDup(v1,v2,Adj_classe,ArbresDAdjacences, v1dansA1);    
               break;
            default: 
               cout<<"From Backtracking.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
	       exit(EXIT_FAILURE);
            break;  
         }
         break;
      default: 
         cout<<"From Backtracking.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v1)) on node v1 with id "<<v1->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v1))<<endl;
	 exit(EXIT_FAILURE);
      break;  
   }
}
