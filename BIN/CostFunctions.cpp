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

File: CostFunctions.cpp                     Last modified on: 12/01/2015
Created by: Sèverine Bérard               Created on: 10/01/2011
--------------------------------------------------------------------------
Specification: 
File containing the cost functions of the function DECO in Step3_DECO.cpp
=========================================================================*/
#include "CostFunctions.h"


//renvoie l'événement associé au noeud n
string E(Node * n){
   if (n->hasBranchProperty(typ)){
      BppString * Evt = dynamic_cast<BppString*> (n->getBranchProperty(typ));
      return Evt->toSTL();
   }
   else{
      cout<<"From CostFunctions.cpp (E): node with id "<<n->getId()<<" has no associated event !!"<<endl; 
      return "ND";
   }
}

//renvoie l'espèce associée au noeud n
int Espece(Node * n){
   convert<int> c;
   if (n->hasNodeProperty(esp)){
      BppString * Espece = dynamic_cast<BppString*> (n->getNodeProperty(esp));
      return c.from_s(Espece->toSTL());
   } 
   else{
      cout<<"From CostFunctions.cpp (Espece): node with id "<<n->getId()<<" has no associated specie !!"<<endl; 
      return -1;
   }
}

int EvtToInt(string e){
   if (e=="Extant")
      return 1;
   else if (e=="Spec")
      return 2;
   else if (e=="GDup")
      return 3;
   else if (e=="GLos")
      return 4;
   else if (e=="ADup")
      return 5;
   else if (e=="ALos")
      return 6;
   else if (e=="Crea")
      return 7;
   else if (e=="Break")
      return 8;
   else{
      cout<<"From CostFunctions.cpp (EvtToInt): event "<<e<<" is unknown"<<endl; 
      return -1;
   }
}

//
//min à plusieurs éléments
//
float min(float x, float y){
   return x > y ? y : x;
}
float min(float x, float y, float z){
   return x > y ? min(y, z) : min(x, z);
}
float min(float x, float y, float z, float a, float b){
   return min(min(x,y),min(z,a,b));
}
int min_pos(float x, float y){
   return x > y ? 2 : 1;
}

int min_pos(float x, float y, float z){
   return x < y ? (x < z ? 1 : min_pos(y, z)+1) : min_pos(y, z)+1;
}
int min_pos(float x, float y, float z, float a, float b)
{
   return min(x,y) < min(z,a,b) ? min_pos(x,y) : min_pos(z,a,b)+2;
}

bool appAdj(adjacence a, vector<adjacence> * Adjacences){ 
   bool trouve=false;
   vector<adjacence>::iterator it=Adjacences->begin();
   while (!trouve && it!=Adjacences->end()){
      if (a.gene1==(*it).gene1 && a.gene2==(*it).gene2)
         trouve=true;
      it++;
   }
   return trouve;
}

float recupCoutC1(Node * v1, Node * v2){
   pair<Node *,Node *> p(v1,v2);
   return C1[p];
}

float recupCoutC0(Node * v1, Node * v2){
   pair<Node *,Node *> p(v1,v2);
   return C0[p];
}

/////////////////////////////////////////
//////      Cas pour Calcul C1      /////
/////////////////////////////////////////

// Cas 1 - E(v1)=Extant / E(v2)=Extant
float C1ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe){
   if(affich_CalculCout)
      cout <<"Je suis dans C1ExtantExtant"<<endl;
   
   string nomV1 = v1->getName();
   string nomV2 = v2->getName();
   string nom_gene1=nomV1.substr(0,nomV1.rfind(sep,nomV1.length()-1));
   string nom_gene2=nomV2.substr(0,nomV2.rfind(sep,nomV2.length()-1));
   //On cherche l'adj dans les deux sens g1-g2 ou g2-g1
   adjacence a12,a21;
   a12.gene1=nom_gene1;
   a12.gene2=nom_gene2;
   a21.gene1=nom_gene2;
   a21.gene2=nom_gene1;

   if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe))
      return 0.0;
   else
      return INFINI;
}

float C1GLosAny(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1GLosAny"<<endl;
   return 0.0;
}

float C1GDupExtantOrSpec(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1GDupExtantOrSpec"<<endl;
   return min(min(recupCoutC1(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2),
        recupCoutC0(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2),
        recupCoutC1(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)+Crea),
        recupCoutC0(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)+Break);
}

float C1GLosGLos(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1GLosGLos"<<endl;
   return 0.0;
}

float C1SpecSpec(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1SpecSpec"<<endl;
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
   
   return min(min(recupCoutC1(fgV1,fgV2) + recupCoutC1(fdV1,fdV2),
        recupCoutC1(fgV1,fgV2) + recupCoutC0(fdV1,fdV2) + Break,
        recupCoutC0(fgV1,fgV2) + recupCoutC1(fdV1,fdV2) + Break),
         recupCoutC0(fgV1,fgV2) + recupCoutC0(fdV1,fdV2) + 2*Break);
}

float C1GDupGDup(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1GDupGDup"<<endl;
   float D1,D2,D12;
   
   //idem DupExtantOrSpec
   D1=min(min(recupCoutC1(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2),
         recupCoutC0(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2),
         recupCoutC1(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)+Crea),
         recupCoutC0(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)+Break);
   
   D2=min(min(recupCoutC1(v1,v2->getSon(0))+recupCoutC0(v1,v2->getSon(1)),
         recupCoutC0(v1,v2->getSon(0))+recupCoutC1(v1,v2->getSon(1)),
         recupCoutC1(v1,v2->getSon(0))+recupCoutC1(v1,v2->getSon(1))+Crea),
         recupCoutC0(v1,v2->getSon(0))+recupCoutC0(v1,v2->getSon(1))+Break);
   
   //un min5 de min3 => 15 lignes
   D12=min(min(
          //(1) (2) et (3)
          min(recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+
         recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0)),
         recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+
         recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Crea,
         recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+
         recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Crea),
          //(4) (5) et (6)
          min(recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+
         recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+2*Crea,
         recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+
         recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break,
         recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+
         recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Break+Crea),
          //(7) (8) et (9)   
          min(recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+
         recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break+Crea,
         recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+
         recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break,
         recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+
         recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Break+Crea),
          //(10) (11) et (12)
          min(recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+
         recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break+Crea,
         recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+
         recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0)),
         recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC1(v1->getSon(1),v2->getSon(1))+
         recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Crea),
          //(13) (14) et (15)
          min(recupCoutC1(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+
         recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Crea,
         recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+
         recupCoutC1(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+Break,
         recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+
         recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC1(v1->getSon(1),v2->getSon(0))+Break)),
      //(16) ajouté
      recupCoutC0(v1->getSon(0),v2->getSon(0))+recupCoutC0(v1->getSon(1),v2->getSon(1))+
      recupCoutC0(v1->getSon(0),v2->getSon(1))+recupCoutC0(v1->getSon(1),v2->getSon(0))+2*Break);
   
   return min(D1,D2,D12);
}

float C0ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe){
   if(affich_CalculCout)
      cout <<"Je suis dans C0ExtantExtant"<<endl;
   string nomV1 = v1->getName();
   string nomV2 = v2->getName();
   string nom_gene1=nomV1.substr(0,nomV1.rfind(sep,nomV1.length()-1));
   string nom_gene2=nomV2.substr(0,nomV2.rfind(sep,nomV2.length()-1));
   //On cherche d'adj dans les deux sens g1-g2 ou g2-g1
   adjacence a12,a21;
   a12.gene1=nom_gene1;
   a12.gene2=nom_gene2;
   a21.gene1=nom_gene2;
   a21.gene2=nom_gene1;
   if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe))
      return INFINI;
   else
      return 0.0;
}

float C0GLosAny(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0GLosAny"<<endl;
   return 0.0;
}

float C0GDupExtantOrSpec(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0GDupExtantOrSpec"<<endl;
   return min(recupCoutC0(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2),
         min(recupCoutC1(v1->getSon(0),v2)+recupCoutC0(v1->getSon(1),v2)+Crea,
        recupCoutC0(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)+Crea,
        recupCoutC1(v1->getSon(0),v2)+recupCoutC1(v1->getSon(1),v2)+2*Crea));
}

float C0GLosGLos(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0GLosGLos"<<endl;
   return 0.0;
}

float C0SpecSpec(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0SpecSpec"<<endl;
   Node * fgV1=v1->getSon(0);
   Node * fdV1=v1->getSon(1);
   Node * fgV2;
   Node * fdV2;
   if (Espece(fgV1)==Espece(v2->getSon(0)))
      {
    fgV2=v2->getSon(0);
    fdV2=v2->getSon(1);
      }
   else
      {
    fgV2=v2->getSon(1);
    fdV2=v2->getSon(0);
      }

   return min(min(recupCoutC0(fgV1,fgV2) + recupCoutC0(fdV1,fdV2),
        recupCoutC1(fgV1,fgV2) + recupCoutC0(fdV1,fdV2) + Crea,
        recupCoutC0(fgV1,fgV2) + recupCoutC1(fdV1,fdV2) + Crea),
         recupCoutC1(fgV1,fgV2) + recupCoutC1(fdV1,fdV2) + 2*Crea);
   
}

float C0GDupGDup(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0GDupGDup"<<endl;
   float D1;
   D1=min(//Les 2 cas originaux d'Éric
     min(recupCoutC0(v1->getSon(0),v2) + recupCoutC0(v1->getSon(1),v2),
         recupCoutC0(v1,v2->getSon(0)) + recupCoutC0(v1,v2->getSon(1))),
     //3 cas où v1 se duplique avant v2
     min(recupCoutC0(v1->getSon(0),v2) + recupCoutC1(v1->getSon(1),v2)+Crea,
         recupCoutC1(v1->getSon(0),v2) + recupCoutC0(v1->getSon(1),v2)+Crea,
         recupCoutC1(v1->getSon(0),v2) + recupCoutC1(v1->getSon(1),v2)+2*Crea),
     //3 cas où v2 se duplique avant v1
     min(recupCoutC0(v1,v2->getSon(0)) + recupCoutC1(v1,v2->getSon(1))+Crea,
         recupCoutC1(v1,v2->getSon(0)) + recupCoutC0(v1,v2->getSon(1))+Crea,
         recupCoutC1(v1,v2->getSon(0)) + recupCoutC1(v1,v2->getSon(1))+2*Crea));
   return D1;
}

// C1
//
//On suppose que les noeuds sont bien de la même espèce
float calculeC1(Node * v1, Node * v2, vector<adjacence> * Adj_classe)
{
   switch (EvtToInt(E(v1))){  
   case 1 :
      switch (EvtToInt(E(v2))){
      case 1 :
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 Extant"<<endl;
    return C1ExtantExtant(v1,v2,Adj_classe);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 GLos"<<endl;
    return C1GLosAny(v2,v1);
    break;    
      case 2 :
    cout<<"From CostFunctions.cpp (calculeC1): v1 Extant and v2 Spec IMPOSSIBLE"<<endl;
    exit(EXIT_FAILURE);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 GDup"<<endl;
    return C1GDupExtantOrSpec(v2,v1);    
    break;
      default: 
    cout<<"From CostFunctions.cpp (calculeC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }   
      break;
      
   case 4 :
      switch (EvtToInt(E(v2))){
      case 1 :
    if(affich_CalculCout)
       cout<<"v1 GLos et v2 Extant"<<endl;
    return C1GLosAny(v1,v2);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 GLos et v2 Glos"<<endl;
    return C1GLosGLos(v1,v2);
    break;    
      case 2 :
    if(affich_CalculCout)
       cout<<"v1 GLos et v2 Spec"<<endl;
    return C1GLosAny(v1,v2);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 GLos et v2 GDup"<<endl;
    return C1GLosAny(v1,v2);    
    break;
      default: 
    cout<<"From CostFunctions.cpp (calculeC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   case 2 :
      switch (EvtToInt(E(v2))){
      case 1 :
    cout<<"From CostFunctions.cpp (calculeC1): v1 Spec and v2 Extant IMPOSSIBLE"<<endl;
    exit(EXIT_FAILURE);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 Spec et v2 Glos"<<endl;
    return C1GLosAny(v2,v1);
    break;    
      case 2 :
    if(affich_CalculCout)
       cout<<"v1 Spec et v2 Spec"<<endl;
    return C1SpecSpec(v1,v2);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 Spec et v2 GDup"<<endl;
    return C1GDupExtantOrSpec(v2,v1);
    break;
      default: 
    cout<<"From CostFunctions.cpp (calculeC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   case 3 :
      switch (EvtToInt(E(v2))){
      case 1 :
    if(affich_CalculCout)
       cout<<"v1 GDup et v2 Extant"<<endl;
    return C1GDupExtantOrSpec(v1,v2);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 GDup et v2 Glos"<<endl;
    return C1GLosAny(v2,v1);
    break;    
      case 2 :
    if(affich_CalculCout)
       cout<<"v1 GDup et v2 Spec"<<endl;
    return C1GDupExtantOrSpec(v1,v2);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 GDup et v2 GDup"<<endl;
    return C1GDupGDup(v1,v2);    
    break;
      default: 
    cout<<"From CostFunctions.cpp (calculeC1): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   default: 
      cout<<"From CostFunctions.cpp (calculeC1): ERROR the call to EvtToInt(E(v1)) on node v1 with id "<<v1->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v1))<<endl;
      exit(EXIT_FAILURE);
      break;  
   }
}

// C0
//
//On suppose que les noeuds sont bien de la même espèce
float calculeC0(Node * v1, Node * v2, vector<adjacence> * Adj_classe)
{
   switch (EvtToInt(E(v1))){  
   case 1 :
      switch (EvtToInt(E(v2))){
      case 1 :
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 Extant"<<endl;
    return C0ExtantExtant(v1,v2,Adj_classe);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 GLos"<<endl;
    return C0GLosAny(v2,v1);
    break;    
      case 2 :
    cout<<"From CostFunctions.cpp (calculeC0): v1 Extant and v2 Spec IMPOSSIBLE"<<endl;
    exit(EXIT_FAILURE);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 GDup"<<endl;
    return C0GDupExtantOrSpec(v2,v1);    
    break;
      default: 
    cout<<"From CostFunctions.cpp (calculeC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }   
      break;
      
   case 4 :
      switch (EvtToInt(E(v2))){
      case 1 :
    if(affich_CalculCout)
       cout<<"v1 GLos et v2 Extant"<<endl;
    return C0GLosAny(v1,v2);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 GLos et v2 Glos"<<endl;
    return C0GLosGLos(v1,v2);
    break;    
      case 2 :
    if(affich_CalculCout)
       cout<<"v1 GLos et v2 Spec"<<endl;
    return C0GLosAny(v1,v2);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 GLos et v2 GDup"<<endl;
    return C0GLosAny(v1,v2);    
    break;
      default: 
    cout<<"From CostFunctions.cpp (calculeC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   case 2 :
      switch (EvtToInt(E(v2))){
      case 1 :
    cout<<"From CostFunctions.cpp (calculeC1): v1 Spec and v2 Extant IMPOSSIBLE"<<endl;
    exit(EXIT_FAILURE);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 Spec et v2 Glos"<<endl;
    return C0GLosAny(v2,v1);
    break;    
      case 2 :
    if(affich_CalculCout)
       cout<<"v1 Spec et v2 Spec"<<endl;
    return C0SpecSpec(v1,v2);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 Spec et v2 GDup"<<endl;
    return C0GDupExtantOrSpec(v2,v1);
    break;
      default: 
    cout<<"From CostFunctions.cpp (calculeC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   case 3 :
      switch (EvtToInt(E(v2))){
      case 1 :
    if(affich_CalculCout)
       cout<<"v1 GDup et v2 Extant"<<endl;
    return C0GDupExtantOrSpec(v1,v2);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 GDup et v2 Glos"<<endl;
    return C0GLosAny(v2,v1);
    break;    
      case 2 :
    if(affich_CalculCout)
       cout<<"v1 GDup et v2 Spec"<<endl;
    return C0GDupExtantOrSpec(v1,v2);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 GDup et v2 GDup"<<endl;
    return C0GDupGDup(v1,v2);    
    break;
      default: 
    cout<<"From CostFunctions.cpp (calculeC0): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   default: 
      cout<<"From CostFunctions.cpp (calculeC0): ERROR the call to EvtToInt(E(v1)) on node v1 with id "<<v1->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v1))<<endl;
      exit(EXIT_FAILURE);
      break;  
   }
}
