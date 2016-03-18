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

File: BacktrackingProba.cpp                    Last modified on: 26/03/2015
Created by: Yoann Anselmetti/Sèverine Bérard   Created on: 10/01/2011
--------------------------------------------------------------------------
Specification: 
File containing the functions doing the backtracking part of the function
DECO in Step3_DECO.cpp
=========================================================================*/
#include "BacktrackingProba.h"

Node * parcoursArriereCassure(Node * v1, Node * v2, bool v1dansA1){
   nb_Bk++;
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereCassure"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
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
      cout <<"Je suis dans parcoursArriereC1ExtantExtant"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
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
   //On cherche l'adj dans les deux sens g1-g2 ou g2-g1
   adjacence a12,a21;
   a12.gene1=nom_gene1;
   a12.gene2=nom_gene2;
   a21.gene1=nom_gene2;
   a21.gene2=nom_gene1;

   if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe)){
//      cout<<"Existing adjacency!!!!"<<endl;
      n->setBranchProperty(typ,BppString(ga));
      n->setBranchProperty("D", BppString("?"));
      if (v1dansA1)
         n->setName(nom_gene1+sepAdj+nom_gene2+sep+espece);
      else
         n->setName(nom_gene2+sepAdj+nom_gene1+sep+espece);
   }
   else{
//      cout<<"There is 1 new adjacency!!!!"<<endl;
      convert<int> c;
      ofstream OffilAdj_file;
      if(!light_mode){
         string output_adj_file=prefixe+output_ext_adj_file;
         OffilAdj_file.open(output_adj_file.c_str(), ios::out | ios::app);
         if (!OffilAdj_file){
            cout<<"\nFrom BacktrackingProba.cpp: ERROR while opening file "<<output_adj_file<<endl;
            if (file_log)
               Offile_log<<"\nFI\tFrom BacktrackingProba.cpp: ERROR while opening file "<<output_adj_file<<endl;
            exit(EXIT_FAILURE);
         }
      }

      ofstream OffilNewAdj;
      string output_new_adj_file=prefixe+output_new_adj;
      OffilNewAdj.open(output_new_adj_file.c_str(), ios::out | ios::app);
      if (!OffilNewAdj){
         cout<<"\nFrom BacktrackingProba.cpp: ERROR while opening file "<<output_new_adj_file<<endl;
         if (file_log)
            Offile_log<<"\nFI\tFrom BacktrackingProba.cpp: ERROR while opening file "<<output_new_adj_file<<endl;
         exit(EXIT_FAILURE);
      }

      n->setBranchProperty(typ,BppString(ga));
      n->setBranchProperty("D", BppString("?"));
      //bool input=true;
      if (v1dansA1){
//         int adj_class_input=adjacency_class(nom_gene1, nom_gene2, input);
//         input=false;
//         int adj_class_output=adjacency_class(nom_gene1, nom_gene2, input);
         n->setName(nom_gene1+sepAdj+nom_gene2+sep+espece);
//         cout<<"New extant adjacency: "<<nom_gene1<<sepAdj<<nom_gene2<<sep<<espece<<endl;
         if(!light_mode){
            OffilAdj_file<<nom_gene1<<"\t"<<nom_gene2<<endl;
         }
//         OffilNewAdj<<nom_gene1<<"\t"<<nom_gene2<<"\t"<<c.to_s(adj_class_input)<<"\t"<<c.to_s(adj_class_output)<<endl;
         OffilNewAdj<<nom_gene1<<"\t"<<nom_gene2<<"\t"<<c.to_s(extant_gene_adj_nb[nom_gene1])<<"\t"<<c.to_s(extant_gene_adj_nb[nom_gene2])<<endl;
//         cout<<nom_gene1<<"\t"<<nom_gene2<<"\t"<<c.to_s(extant_gene_adj_nb[nom_gene1])<<"\t"<<c.to_s(extant_gene_adj_nb[nom_gene2])<<endl;
      }
      else{
//         int adj_class_input=adjacency_class(nom_gene2, nom_gene1, input);
//         input=false;
//         int adj_class_output=adjacency_class(nom_gene2, nom_gene1, input);
         n->setName(nom_gene2+sepAdj+nom_gene1+sep+espece);
//         cout<<"New extant adjacency: "<<nom_gene2<<sepAdj<<nom_gene1<<sep<<espece<<endl;
         if(!light_mode){
            OffilAdj_file<<nom_gene2<<"\t"<<nom_gene1<<endl;
         }
//         OffilNewAdj<<nom_gene2<<"\t"<<nom_gene1<<"\t"<<c.to_s(adj_class_input)<<"\t"<<c.to_s(adj_class_output)<<endl;
         OffilNewAdj<<nom_gene2<<"\t"<<nom_gene1<<"\t"<<c.to_s(extant_gene_adj_nb[nom_gene2])<<"\t"<<c.to_s(extant_gene_adj_nb[nom_gene1])<<endl;
//         cout<<nom_gene2<<"\t"<<nom_gene1<<"\t"<<c.to_s(extant_gene_adj_nb[nom_gene2])<<"\t"<<c.to_s(extant_gene_adj_nb[nom_gene1])<<endl;
      }
      if(!light_mode){
         OffilAdj_file.close();
      }
      OffilNewAdj.close();
   }
   return n;
}

// Cas 2 - E(v1)=GLos / E(v2)!=GLos
Node * parcoursArriereC1GLosAny(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC1GLosAny"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
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
      cout <<"Je suis dans parcoursArriereC1GLosGLos"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
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
      cout <<"Je suis dans parcoursArriereC1GDupExtantOrSpec"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
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
   int C = recupCasC1(v1,v2);
   // (1)
   if(C==141){    
      n->addSon(0,parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   // (2)
   else if(C==142){    
      n->addSon(0,parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   // (3)
   else if(C==143){    
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
   else if(C==144){    
      n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2, v1dansA1));
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else{
      cout<<"From BacktrackingProba.cpp (parcoursArriereC1GDupExtantOrSpec): ERROR none of the 4 cases correspond to the stored one : "<<C<<endl;
      exit(EXIT_FAILURE);}
   return n;
}

// Cas 5 - E(v1)=Spec / E(v2)=Spec
Node * parcoursArriereC1SpecSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC1SpecSpec"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;

   if (!v1dansA1)//dans les cas que l'on a calculés lors de la phase ascendante v1 et v2 sont inversés
      {
    Node * v3=v1;
    v1=v2;
    v2=v3;
    v1dansA1=true;
      }
   
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

   int C = recupCasC1(v1,v2);

 if(affich_ParcoursArriere)
    cout <<"C="<<C<<",id fgV1="<<fgV1->getId()<<",id fdV1="<<fdV1->getId()<<",id fgV2="<<fgV2->getId()<<",id fdV2="<<fdV2->getId()<<", v1dansA1="<<v1dansA1<<endl;

   // (1)
   if(C==151){
      n->addSon(0,parcoursArriereC1(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      n->addSon(1,parcoursArriereC1(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1));
   }
   // (2)
   else if(C==152){
      n->addSon(0,parcoursArriereC1(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      n->addSon(1,parcoursArriereCassure(fdV1,fdV2, v1dansA1));
      parcoursArriereC0(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   // (3)
   else if(C==153){
      n->addSon(0,parcoursArriereC1(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1));
      n->addSon(1,parcoursArriereCassure(fgV1,fgV2, v1dansA1));
      parcoursArriereC0(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   // (4)
   else if(C==154){
      n->addSon(0,parcoursArriereCassure(fgV1,fgV2, v1dansA1));
      n->addSon(1,parcoursArriereCassure(fdV1,fdV2, v1dansA1));
      parcoursArriereC0(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else 
      {
    cout<<"From BacktrackingProba.cpp (parcoursArriereC1SpecSpec): ERROR none of the 4 cases correspond to the stored one : "<<C<<" (v1dansA1 "<<v1dansA1<<")"<<endl;
    exit(EXIT_FAILURE);
      }
   return n;
}

// Cas 6 - E(v1)=GDup / E(v2)=GDup
Node * parcoursArriereC1GDupGDup(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere){
      cout <<"Je suis dans parcoursArriereC1GDupGDup"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
      cout<<"v1->getSon(0)="<<v1->getSon(0)->getId()<<" v1->getSon(1)="<<v1->getSon(1)->getId()<<endl;
      cout<<"v2->getSon(0)="<<v2->getSon(0)->getId()<<" v2->getSon(1)="<<v2->getSon(1)->getId()<<endl;
   }
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
   
   int C = recupCasC1(v1,v2);

   if (!v1dansA1)//dans les cas que l'on a calculés lors de la phase ascendante v1 et v2 sont inversés
      {
    Node * v3=v1;
    v1=v2;
    v2=v3;
    v1dansA1=true;
      }
   
   if(affich_ParcoursArriere)
      cout <<"C="<<C<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
   
   //Les cas "D1"
   if(C==161){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==162){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==163){
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
   else if(C==164){
      nb_GDup++;
      cout<<"\t\t\tJe suis dans parcoursArriereC1GDupGDup D1, je passe dans le cas C0 + C0 + Break"<<endl;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2, v1dansA1));
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   //Les cas "D2"
   else if(C==165){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==166){
      nb_GDup++;
      n->setBranchProperty(typ,BppString(dupl));
      n->addSon(0,parcoursArriereC1(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
      parcoursArriereC0(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==167){
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
   else if(C==168){
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
      if(C==1601){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      //(2)
      else if(C==1602){
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
      else if(C==1603){
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
      //(4)
      else if(C==1604){
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
      // (5)
      else if(C==1605){
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
         //Cassure
         n->addSon(1,parcoursArriereCassure(v1->getSon(1),v2->getSon(1), v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      // (6)
      else if(C==1606){
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
      //(7)
      else if(C==1607){
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
      //(8) 
      else if(C==1608){
         //Cassure    
         n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2->getSon(0), v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         parcoursArriereC0(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      }
      // (9)   
      else if(C==1609){
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
      //(10)
      else if(C==1610){
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
      // (11)
      else if(C==1611){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      }
      //(12)
      else if(C==1612){
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
      //(13)
      else if(C==1613){
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
      //(14)
      else if(C==1614){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereC1(v1->getSon(0),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1));
         n->addSon(1,parcoursArriereCassure(v1->getSon(1),v2->getSon(0), v1dansA1));
      }
      // (15)
      else if(C==1615){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2->getSon(1), v1dansA1));
         n->addSon(1,parcoursArriereC1(v1->getSon(1),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1));
      }
      //(16) ajouté
      else if(C==1616){
         parcoursArriereC0(v1->getSon(0),v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
         parcoursArriereC0(v1->getSon(1),v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
         n->addSon(0,parcoursArriereCassure(v1->getSon(0),v2->getSon(1), v1dansA1));
         n->addSon(1,parcoursArriereCassure(v1->getSon(1),v2->getSon(0), v1dansA1));
      }
      else 
         {
       cout<<"From BacktrackingProba.cpp (parcoursArriereC1GDupGDup): ERROR none of the 24 cases correspond to the stored one : "<<C<<endl;
       exit(EXIT_FAILURE);
    }
   }
   return n;
}

/////////////////////////////////////////
//////      Cas pour Calcul C0      /////
/////////////////////////////////////////

// Cas 1 - E(v1)=Extant / E(v2)=Extant
//   Modif   //
void parcoursArriereC0ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe, vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0ExtantExtant"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
   
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
      cout<<"From BacktrackingProba.cpp (parcoursArriereC0ExtantExtant): IMPOSSIBLE CASE !!!"<<endl;
      exit(EXIT_FAILURE);
   }
   else{
      //On ne fait rien
   }
}

void parcoursArriereC0GLosAny(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   //On ne fait rien
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0GLosAny"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
}

void parcoursArriereC0GDupExtantOrSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   nb_GDup++;
   convert<int> c;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0GDupExtantOrSpec"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
   int C = recupCasC0(v1,v2);
   if(C==41){
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==42){
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
   else if(C==43){
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
   else if(C==44){
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
    cout<<"From BacktrackingProba.cpp (parcoursArriereC0GDupExtantOrSpec) : ERROR none of the 4 cases correspond to the stored one : "<<C<<" (v1dansA1 "<<v1dansA1<<")"<<endl;
    exit(EXIT_FAILURE);
      }
}

void parcoursArriereC0GLosGLos(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   //On ne fait rien
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0GLosGLos"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
}

void parcoursArriereC0SpecSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0SpecSpec"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;

   if (!v1dansA1)//dans les cas que l'on a calculés lors de la phase ascendante v1 et v2 sont inversés
      {
    Node * v3=v1;
    v1=v2;
    v2=v3;
    v1dansA1=true;
      }
   
   
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
   
   int C = recupCasC0(v1,v2);
   
   if(affich_ParcoursArriere)
      cout <<"C="<<C<<",id fgV1="<<fgV1->getId()<<",id fdV1="<<fdV1->getId()<<",id fgV2="<<fgV2->getId()<<",id fdV2="<<fdV2->getId()<<", v1dansA1="<<v1dansA1<<endl;
   
   if(C==51){
      parcoursArriereC0(fgV1,fgV2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(fdV1,fdV2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else  if(C==52){
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
   else  if(C==53){
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
   else  if(C==54){
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
    cout<<"From BacktrackingProba.cpp (parcoursArriereC0SpecSpec): ERROR none of the 4 cases correspond to the stored one : "<<C<<" (v1dansA1 "<<v1dansA1<<")"<<endl;
    exit(EXIT_FAILURE);
      }
}

void parcoursArriereC0GDupGDup(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
   if(affich_ParcoursArriere)
      cout <<"Je suis dans parcoursArriereC0GDupGDup"<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
   int C = recupCasC0(v1,v2);
   convert<int> c;
   BppString * Espece = dynamic_cast<BppString*> (v1->getNodeProperty(esp));
   
   if (!v1dansA1)//dans les cas que l'on a calculés lors de la phase ascendante v1 et v2 sont inversés
      {
    Node * v3=v1;
    v1=v2;
    v2=v3;
    v1dansA1=true;
      }
   
   if(affich_ParcoursArriere)
      cout <<"C="<<C<<",id v1="<<v1->getId()<<",id v2="<<v2->getId()<<", v1dansA1="<<v1dansA1<<endl;
   
   //Les 2 cas originaux d'Éric
   if(C==61){
      nb_GDup++;
      parcoursArriereC0(v1->getSon(0),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1->getSon(1),v2,Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   else if(C==62){
      nb_GDup++;
      parcoursArriereC0(v1,v2->getSon(0),Adj_classe,ArbresDAdjacences, v1dansA1);
      parcoursArriereC0(v1,v2->getSon(1),Adj_classe,ArbresDAdjacences, v1dansA1);
   }
   //3 cas où v1 se duplique avant v2
      else if(C==63){
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
      else if(C==64){
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
      else if(C==65){
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
      else if(C==66){
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
      else if(C==67){
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
      else if(C==68){
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
    cout<<"From BacktrackingProba.cpp (parcoursArriereC0GDupGDup): ERROR none of the 8 cases correspond to the stored one : "<<C<<endl;      
    exit(EXIT_FAILURE);
      }
}

// parcoursArriereC1
Node * parcoursArriereC1(Node * v1, Node * v2, vector<adjacence> * Adj_classe, vector<Tree *> * ArbresDAdjacences, bool v1dansA1){
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
               cout<<"From BacktrackingProba.cpp (parcoursArriereC1): v1 Extant and v2 Spec IMPOSSIBLE"<<endl;
          exit(EXIT_FAILURE);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Extant et v2 GDup"<<endl;
               return parcoursArriereC1GDupExtantOrSpec(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);    
               break;
            default: 
               cout<<"From BacktrackingProba.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
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
               cout<<"From BacktrackingProba.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
          exit(EXIT_FAILURE);
          break;  
         }
         break;
      case 2 :
         switch (EvtToInt(E(v2))){
            case 1 :
               cout<<"From BacktrackingProba.cpp (parcoursArriereC1): v1 Spec and v2 Extant IMPOSSIBLE"<<endl;
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
               cout<<"From BacktrackingProba.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
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
               cout<<"From BacktrackingProba.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
          exit(EXIT_FAILURE);
          break;  
         }
         break;
   default: 
      cout<<"From BacktrackingProba.cpp (parcoursArriereC1): ERROR the call to EvtToInt(E(v1)) on node v1 with id "<<v1->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v1))<<endl;
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
               cout<<"From BacktrackingProba.cpp (parcoursArriereC0): v1 Extant and v2 Spec IMPOSSIBLE"<<endl;
          exit(EXIT_FAILURE);
               break;    
            case 3 :
               if(affich_ParcoursArriere)
                  cout<<"v1 Extant et v2 GDup"<<endl;
               parcoursArriereC0GDupExtantOrSpec(v2,v1,Adj_classe,ArbresDAdjacences, !v1dansA1);    
               break;
            default: 
               cout<<"From BacktrackingProba.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
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
               cout<<"From BacktrackingProba.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
          exit(EXIT_FAILURE);
            break;  
         }
         break;
      case 2 :
         switch (EvtToInt(E(v2))){
            case 1 :
               cout<<"From BacktrackingProba.cpp (parcoursArriereC0): v1 Spec and v2 Extant IMPOSSIBLE"<<endl;
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
               cout<<"From BacktrackingProba.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
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
               cout<<"From BacktrackingProba.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
          exit(EXIT_FAILURE);
            break;  
         }
         break;
      default: 
         cout<<"From BacktrackingProba.cpp (parcoursArriereC0): ERROR the call to EvtToInt(E(v1)) on node v1 with id "<<v1->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v1))<<endl;
    exit(EXIT_FAILURE);
      break;  
   }
}

// int adjacency_class(string nom_gene1, string nom_gene2, bool input){
//    if (extant_gene_adj_nb[nom_gene1]==0){
//       if (extant_gene_adj_nb[nom_gene2]==0){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 1;
//          }
//          else
//             return 1;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]==1){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 2;
//          }
//          else
//             return 2;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]==2){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 5;
//          }
//          else
//             return 5;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]>2){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 10;
//          }
//          else
//             return 10;
//       }
//       else{
//          cout<<"From BacktrackingProba.cpp (adjacency_class): ERROR the call to extant_gene_adj_nb[nom_gene2] on gene nom_gene2 "<<nom_gene2<<" doesn't return one of the know values, it returns:"<<extant_gene_adj_nb[nom_gene2]<<endl;
//          exit(EXIT_FAILURE);
//       }
//    }
//    else if (extant_gene_adj_nb[nom_gene1]==1){
//       if (extant_gene_adj_nb[nom_gene2]==0){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 3;
//          }
//          else
//             return 3;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]==1){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 4;
//          }
//          else
//             return 4;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]==2){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 7;
//          }
//          else
//             return 7;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]>2){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 12;
//          }
//          else
//             return 12;
//       }
//       else{
//          cout<<"From BacktrackingProba.cpp (adjacency_class): ERROR the call to extant_gene_adj_nb[nom_gene2] on gene nom_gene2 "<<nom_gene2<<" doesn't return one of the know values, it returns:"<<extant_gene_adj_nb[nom_gene2]<<endl;
//          exit(EXIT_FAILURE);
//       }
//    }
//    else if (extant_gene_adj_nb[nom_gene1]==2){
//       if (extant_gene_adj_nb[nom_gene2]==0){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 6;
//          }
//          else
//             return 6;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]==1){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 8;
//          }
//          else
//             return 8;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]==2){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 9;
//          }
//          else
//             return 9;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]>2){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 14;
//          }
//          else
//             return 14;
//       }
//       else{
//          cout<<"From BacktrackingProba.cpp (adjacency_class): ERROR the call to extant_gene_adj_nb[nom_gene2] on gene nom_gene2 "<<nom_gene2<<" doesn't return one of the know values, it returns:"<<extant_gene_adj_nb[nom_gene2]<<endl;
//          exit(EXIT_FAILURE);
//       }
//    }
//    else if (extant_gene_adj_nb[nom_gene1]>2){
//       if (extant_gene_adj_nb[nom_gene2]==0){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 11;
//          }
//          else
//             return 11;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]==1){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 13;
//          }
//          else
//             return 13;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]==2){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 15;
//          }
//          else
//             return 15;
//       }
//       else if (extant_gene_adj_nb[nom_gene2]>2){
//          if (input){
//             extant_gene_adj_nb[nom_gene1]+=1;
//             extant_gene_adj_nb[nom_gene2]+=1;
//             return 16;
//          }
//          else
//             return 16;
//       }
//       else{
//          cout<<"From BacktrackingProba.cpp (adjacency_class): ERROR the call to extant_gene_adj_nb[nom_gene2] on gene nom_gene2 "<<nom_gene2<<" doesn't return one of the know values, it returns:"<<extant_gene_adj_nb[nom_gene2]<<endl;
//          exit(EXIT_FAILURE);
//       }
//    }
//    else{
//       cout<<"From BacktrackingProba.cpp (adjacency_class): ERROR the call to extant_gene_adj_nb[nom_gene1] on gene nom_gene1 "<<nom_gene1<<" doesn't return one of the know values, it returns:"<<extant_gene_adj_nb[nom_gene1]<<endl;
//          exit(EXIT_FAILURE);
//    }
// }
