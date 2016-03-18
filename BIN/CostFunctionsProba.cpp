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

File: CostFunctionsProba.cpp                   Last modified on: 27/05/2015
Created by: Yoann Anselmetti/Sèverine Bérard   Created on: 10/01/2011
--------------------------------------------------------------------------
Specification: 
File containing the cost functions of the function DECO in Step3_DECO.cpp
=========================================================================*/
#include "CostFunctionsProba.h"

//renvoie l'événement associé au noeud n
string E(Node * n){
   if (n->hasBranchProperty(typ)){
      BppString * Evt = dynamic_cast<BppString*> (n->getBranchProperty(typ));
      return Evt->toSTL();
   }
   else{
      cout<<"From CostFunctionsProba.cpp (E): node with id "<<n->getId()<<" has no associated event !!"<<endl; 
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
      cout<<"From CostFunctionsProba.cpp (Espece): node with id "<<n->getId()<<" has no associated specie !!"<<endl; 
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
      cout<<"From CostFunctionsProba.cpp (EvtToInt): event "<<e<<" is unknown"<<endl; 
      return -1;
   }
}

//////////////////////////////////////////////////////////////////////////
//min à plusieurs éléments
//////////////////////////////////////////////////////////////////////////
//min2
pair<float,int> minPaire(float x, float y)
{
   pair<float,int> p;
   if (x<=y){
      p.first=x;
      p.second=1;
   }
   else{
      p.first=y;
      p.second=2;
   }
   // cout<<" x="<<x<<" y="<<y<<endl;
   // cout<<"minPaire 2 args retourne min="<<p.first<<" et pos="<<p.second<<endl;
   return p;
}

//min3
pair<float,int> minPaire(float x, float y, float z)
{
   pair<float,int> p;
   if (x<=y)
      if (x<=z){
    p.first=x;
    p.second=1;
      }
      else{
    p.first=z;
    p.second=3;
      }
   else
      if (y<=z){
    p.first=y;
    p.second=2;
      }
      else{
    p.first=z;
    p.second=3;
      }
   // cout<<" x="<<x<<" y="<<y<<" z="<<z<<endl;
   // cout<<"minPaire retourne 3 args min="<<p.first<<" et pos="<<p.second<<endl;
   return p;
}

//min4
pair<float,int> minPaire(float x, float y, float z, float a){
   pair<float,int> p1 = minPaire(x,y);
   pair<float,int> p2 = minPaire(z,a);
   pair<float,int> p;
   if (p1.first<=p2.first){
      p.first=p1.first;
      p.second=p1.second;
   }
   else{
      p.first=p2.first;
      p.second=2+p2.second;
   }
   // cout<<" x="<<x<<" y="<<y<<" z="<<z<<" a="<<a<<endl;
   // cout<<"minPaire 4 args retourne min="<<p.first<<" et pos="<<p.second<<endl;
   return p;
}

//min5 de 5 paires de min calculées sur 3 valeurs
pair<float,int> minPaire5x3(pair<float,int> x, pair<float,int> y, pair<float,int> z, pair<float,int> a, pair<float,int> b){
   pair<float,int> p1 = minPaire(x.first,y.first);
   pair<float,int> p2 = minPaire(z.first,a.first,b.first);
   pair<float,int> p;
   if (p1.first<=p2.first){
      p.first=p1.first;
      if (p1.second==1)
    p.second=x.second;
      else
    p.second=3+y.second;
   }
   else{
      p.first=p2.first;
      if (p2.second==1)
    p.second=6+z.second;
      else if (p2.second==2)
    p.second=9+a.second;
      else
    p.second=12+b.second;
   }
   // cout<<" x="<<x.first<<"-"<<x.second<<" y="<<y.first<<"-"<<y.second<<" z="<<z.first<<"-"<<z.second<<" a="<<a.first<<"-"<<a.second<<" b="<<a.first<<"-"<<b.second<<endl;
   // cout<<"minPaire5x3 retourne min="<<p.first<<" et pos="<<p.second<<endl;
   return p;
}

//min3 de 3 paires de min calculées sur 2 valeurs pour le 1re paire et 3 pour les 2 autres
pair<float,int> minPaire3x233(pair<float,int> x, pair<float,int> y, pair<float,int> z){
   pair<float,int> p1 = minPaire(x.first,y.first,z.first);
   pair<float,int> p;
   p.first=p1.first;
   if (p1.second==1)
      p.second=x.second;
   else if (p1.second==2)
      p.second=2+y.second;
   else
      p.second=5+z.second;

   // cout<<" x="<<x.first<<"-"<<x.second<<" y="<<y.first<<"-"<<y.second<<" z="<<z.first<<"-"<<z.second<<endl;
   // cout<<"minPaire3x233 retourne min="<<p.first<<" et pos="<<p.second<<endl;
   return p;
}


float min(float x, float y){
   return x > y ? y : x;
}
int min_pos(float x, float y){
   return x > y ? 2 : 1;
}

float min(float x, float y, float z){
   return x > y ? min(y, z) : min(x, z);
}
int min_pos(float x, float y, float z){
   return x < y ? (x < z ? 1 : min_pos(y, z)+1) : min_pos(y, z)+1;
}

float min(float x, float y, float z, float a, float b){
   return min(min(x,y),min(z,a,b));
}
int min_pos(float x, float y, float z, float a, float b)
{
   return min(x,y) < min(z,a,b) ? min_pos(x,y) : min_pos(z,a,b)+2;
}
//////////////////////////////////////////////////////////////////////////


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

float recupCoutC1P(Node * v1, Node * v2){
   pair<Node *,Node *> p(v1,v2);
   return C1P[p].first;
}

float recupCoutC0P(Node * v1, Node * v2){
   pair<Node *,Node *> p(v1,v2);
   return C0P[p].first;
}

int recupCasC1(Node * v1, Node * v2){
   pair<Node *,Node *> p(v1,v2);
   return C1P[p].second;
}

int recupCasC0(Node * v1, Node * v2){
   pair<Node *,Node *> p(v1,v2);
   return C0P[p].second;
}
float log_b(float base, float x){
   return log(x)/log(base);
}

/////////////////////////////////////////
//////      Cases to compute C1P      /////
/////////////////////////////////////////

// Case 1 - E(v1)=Extant / E(v2)=Extant
pair<float,int> C1ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe, TreeTemplate<Node> * S, map<string,string> &gene_species_EXT, map<string,vector<float> > &esp_nb_contig){
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

   if (proba_mode){   
      if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe)){
         pair<float,int> p (0.0,110);
         return p;
      }
      else{
         //Nombre voisins pour v1 et v2
         int nb_adj_gene1=extant_gene_adj_nb[nom_gene1];
         int nb_adj_gene2=extant_gene_adj_nb[nom_gene2];
         float rho=-1.0;
         //Détermination de rho
         if(nb_adj_gene1==2 || nb_adj_gene2==2)
            rho=0.0;
         else if(nb_adj_gene1==1 && nb_adj_gene2==1)
            rho=1.0;
         else if(nb_adj_gene1==0 && nb_adj_gene2==1)
            rho=2.0;
         else if(nb_adj_gene1==1 && nb_adj_gene2==0)
            rho=2.0;
         else if(nb_adj_gene1==0 && nb_adj_gene2==0)
            rho=4.0;
         else if(nb_adj_gene1>2 && nb_adj_gene2>2){
            cout<<"Gene "<<nom_gene1<<" has "<<nb_adj_gene1<<" neighboors!!!"<<endl;
            cout<<"Gene "<<nom_gene2<<" has "<<nb_adj_gene2<<" neighboors!!!"<<endl;
            rho=0.0;
         }
         else if(nb_adj_gene1>2){
            cout<<"Gene "<<nom_gene1<<" has "<<nb_adj_gene1<<" neighboors!!!"<<endl;
            rho=0.0;
         }
         else if(nb_adj_gene2>2){
            cout<<"Gene "<<nom_gene2<<" has "<<nb_adj_gene2<<" neighboors!!!"<<endl;
            rho=0.0;
         }
         if (rho==-1.0){
            cout<<"From CostFunctionsProba.cpp C1ExtantExtant : ERROR in rho determination"<<endl;
            exit(EXIT_FAILURE);
         }
         string species=gene_species_EXT[nom_gene1];
         float proba=rho*esp_nb_contig[species][3];
         float base_log=esp_nb_contig[species][2];

         pair<float,int> p (-log_b(base_log,proba),110);
         return p;
      }
   }
   else{
      if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe)){
         pair<float,int> p (0.0,110);
         return p;
      }
      else{
         pair<float,int> p (INFINI,110);
         return p;
      }
   }
}

// Case 2 - E(v1)=GLos / E(v2)!=GLos
pair<float,int>  C1GLosAny(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1GLosAny"<<endl;
   pair<float,int> p (0.0,120);
   return p;
}

// Case 3 - E(v1)=GLos / E(v2)=GLos
pair<float,int>  C1GLosGLos(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1GLosGLos"<<endl;
   pair<float,int> p (0.0,130);
   return p;
}

// Case 4 - E(v1)=GDup / E(v2)=Extant/Spec
pair<float,int>  C1GDupExtantOrSpec(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1GDupExtantOrSpec"<<endl;

   pair<float,int> a=minPaire(
     recupCoutC1P(v1->getSon(0),v2)+recupCoutC0P(v1->getSon(1),v2),//141
     recupCoutC0P(v1->getSon(0),v2)+recupCoutC1P(v1->getSon(1),v2),//142
     recupCoutC1P(v1->getSon(0),v2)+recupCoutC1P(v1->getSon(1),v2)+Crea,//143
     recupCoutC0P(v1->getSon(0),v2)+recupCoutC0P(v1->getSon(1),v2)+Break);//144

   a.second+=140;
   
   if(affich_CalculCout)
      cout <<"Je retourne la paire "<<a.first<< " "<<a.second<<endl;
   
   return a;
}

// Case 5 - E(v1)=Spec / E(v2)=Spec
pair<float,int>  C1SpecSpec(Node * v1, Node * v2){
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
   
   pair<float,int> a=minPaire(
     recupCoutC1P(fgV1,fgV2) + recupCoutC1P(fdV1,fdV2),//151
     recupCoutC1P(fgV1,fgV2) + recupCoutC0P(fdV1,fdV2) + Break,//152
     recupCoutC0P(fgV1,fgV2) + recupCoutC1P(fdV1,fdV2) + Break,//153
     recupCoutC0P(fgV1,fgV2) + recupCoutC0P(fdV1,fdV2) + 2*Break);//154
   
   a.second+=150;
   
   if(affich_CalculCout)
      cout <<"Je retourne la paire "<<a.first<< " "<<a.second<<endl;
   
   return a;
}

// Case 6 - E(v1)=GDup / E(v2)=GDup
pair<float,int>  C1GDupGDup(Node * v1, Node * v2){
   if(affich_CalculCout)
      cout <<"Je suis dans C1GDupGDup"<<endl;

   //idem DupExtantOrSpec
   pair<float,int> a1=minPaire(
     recupCoutC1P(v1->getSon(0),v2)+recupCoutC0P(v1->getSon(1),v2),//161
     recupCoutC0P(v1->getSon(0),v2)+recupCoutC1P(v1->getSon(1),v2),//162
     recupCoutC1P(v1->getSon(0),v2)+recupCoutC1P(v1->getSon(1),v2)+Crea,//163
     recupCoutC0P(v1->getSon(0),v2)+recupCoutC0P(v1->getSon(1),v2)+Break);//164

   pair<float,int> a2=minPaire(
     recupCoutC1P(v1,v2->getSon(0))+recupCoutC0P(v1,v2->getSon(1)),//165
     recupCoutC0P(v1,v2->getSon(0))+recupCoutC1P(v1,v2->getSon(1)),//166
     recupCoutC1P(v1,v2->getSon(0))+recupCoutC1P(v1,v2->getSon(1))+Crea,//167
     recupCoutC0P(v1,v2->getSon(0))+recupCoutC0P(v1,v2->getSon(1))+Break);//168
   
   pair<float,int> a12=minPaire5x3(
   minPaire(
      recupCoutC1P(v1->getSon(0),v2->getSon(0))+recupCoutC1P(v1->getSon(1),v2->getSon(1))+recupCoutC0P(v1->getSon(0),v2->getSon(1))+recupCoutC0P(v1->getSon(1),v2->getSon(0)),//1601
      recupCoutC1P(v1->getSon(0),v2->getSon(0))+recupCoutC1P(v1->getSon(1),v2->getSon(1))+recupCoutC0P(v1->getSon(0),v2->getSon(1))+recupCoutC1P(v1->getSon(1),v2->getSon(0))+Crea,//1602
      recupCoutC1P(v1->getSon(0),v2->getSon(0))+recupCoutC1P(v1->getSon(1),v2->getSon(1))+recupCoutC1P(v1->getSon(0),v2->getSon(1))+recupCoutC0P(v1->getSon(1),v2->getSon(0))+Crea),//1603
   //(4) (5) et (6)
   minPaire(
      recupCoutC1P(v1->getSon(0),v2->getSon(0))+recupCoutC1P(v1->getSon(1),v2->getSon(1))+recupCoutC1P(v1->getSon(0),v2->getSon(1))+recupCoutC1P(v1->getSon(1),v2->getSon(0))+2*Crea,//1604
      recupCoutC1P(v1->getSon(0),v2->getSon(0))+recupCoutC0P(v1->getSon(1),v2->getSon(1))+recupCoutC0P(v1->getSon(0),v2->getSon(1))+recupCoutC0P(v1->getSon(1),v2->getSon(0))+Break,//1605
      recupCoutC1P(v1->getSon(0),v2->getSon(0))+recupCoutC0P(v1->getSon(1),v2->getSon(1))+recupCoutC0P(v1->getSon(0),v2->getSon(1))+recupCoutC1P(v1->getSon(1),v2->getSon(0))+Break+Crea),//1606
   //(7) (8) et (9)   
   minPaire(
      recupCoutC1P(v1->getSon(0),v2->getSon(0))+recupCoutC0P(v1->getSon(1),v2->getSon(1))+recupCoutC1P(v1->getSon(0),v2->getSon(1))+recupCoutC0P(v1->getSon(1),v2->getSon(0))+Break+Crea,//1607
      recupCoutC0P(v1->getSon(0),v2->getSon(0))+recupCoutC1P(v1->getSon(1),v2->getSon(1))+recupCoutC0P(v1->getSon(0),v2->getSon(1))+recupCoutC0P(v1->getSon(1),v2->getSon(0))+Break,//1608
      recupCoutC0P(v1->getSon(0),v2->getSon(0))+recupCoutC1P(v1->getSon(1),v2->getSon(1))+recupCoutC0P(v1->getSon(0),v2->getSon(1))+recupCoutC1P(v1->getSon(1),v2->getSon(0))+Break+Crea),//1609
   //(10) (11) et (12)
   minPaire(
      recupCoutC0P(v1->getSon(0),v2->getSon(0))+recupCoutC1P(v1->getSon(1),v2->getSon(1))+recupCoutC1P(v1->getSon(0),v2->getSon(1))+recupCoutC0P(v1->getSon(1),v2->getSon(0))+Break+Crea,//1610
      recupCoutC0P(v1->getSon(0),v2->getSon(0))+recupCoutC0P(v1->getSon(1),v2->getSon(1))+recupCoutC1P(v1->getSon(0),v2->getSon(1))+recupCoutC1P(v1->getSon(1),v2->getSon(0)),//1611
      recupCoutC0P(v1->getSon(0),v2->getSon(0))+recupCoutC1P(v1->getSon(1),v2->getSon(1))+recupCoutC1P(v1->getSon(0),v2->getSon(1))+recupCoutC1P(v1->getSon(1),v2->getSon(0))+Crea),//1612
   //(13) (14) et (15)
   minPaire(
      recupCoutC1P(v1->getSon(0),v2->getSon(0))+recupCoutC0P(v1->getSon(1),v2->getSon(1))+recupCoutC1P(v1->getSon(0),v2->getSon(1))+recupCoutC1P(v1->getSon(1),v2->getSon(0))+Crea,//1613
      recupCoutC0P(v1->getSon(0),v2->getSon(0))+recupCoutC0P(v1->getSon(1),v2->getSon(1))+recupCoutC1P(v1->getSon(0),v2->getSon(1))+recupCoutC0P(v1->getSon(1),v2->getSon(0))+Break,//1614
      recupCoutC0P(v1->getSon(0),v2->getSon(0))+recupCoutC0P(v1->getSon(1),v2->getSon(1))+recupCoutC0P(v1->getSon(0),v2->getSon(1))+recupCoutC1P(v1->getSon(1),v2->getSon(0))+Break));//1615

   //(16) ajouté
   float cas16=recupCoutC0P(v1->getSon(0),v2->getSon(0))+recupCoutC0P(v1->getSon(1),v2->getSon(1))+recupCoutC0P(v1->getSon(0),v2->getSon(1))+recupCoutC0P(v1->getSon(1),v2->getSon(0))+2*Break;//1616

   pair<float,int> a=minPaire(a1.first,a2.first,a12.first,cas16);
   
   switch (a.second){
   case 1:
      a.second=160+a1.second;
      break;
      
   case 2:
      a.second=164+a2.second;
      break;
      
   case 3:
      a.second=1600+a12.second;
      break;
      
   case 4:
      a.second=1616;
      break;
      
   default:
      cout<<"\nFrom CostFunctionsProba.cpp (C1GDupGDup): ERROR in computing min case"<<endl;
      if (file_log)
         Offile_log<<"\nFrom CostFunctionsProba.cpp (C1GDupGDup): ERROR in computing min case"<<endl;
      exit(EXIT_FAILURE);
   }
   if(affich_CalculCout)
      cout <<"Je retourne la paire "<<a.first<< " "<<a.second<<endl;
   
   return a;
}

pair<float,int>  C0ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe, TreeTemplate<Node> * S, map<string,string> &gene_species_EXT, map<string,vector<float> > &esp_nb_contig){

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

   if(proba_mode){
      float proba;
      string species=gene_species_EXT[nom_gene1];
      if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe))
         proba=1;
      else{
         //Nombre voisins pour v1 et v2
         int nb_adj_gene1=extant_gene_adj_nb[nom_gene1];
         int nb_adj_gene2=extant_gene_adj_nb[nom_gene2];
         float rho=-1.0;
         //Détermination de rho
         if(nb_adj_gene1==2 || nb_adj_gene2==2)
            rho=0.0;
         else if(nb_adj_gene1==1 && nb_adj_gene2==1)
            rho=1.0;
         else if(nb_adj_gene1==0 && nb_adj_gene2==1)
            rho=2.0;
         else if(nb_adj_gene1==1 && nb_adj_gene2==0)
            rho=2.0;
         else if(nb_adj_gene1==0 && nb_adj_gene2==0)
            rho=4.0;
         else if(nb_adj_gene1>2 && nb_adj_gene2>2){
            cout<<"Gene "<<nom_gene1<<" has "<<nb_adj_gene1<<" neighboors!!!"<<endl;
            cout<<"Gene "<<nom_gene2<<" has "<<nb_adj_gene2<<" neighboors!!!"<<endl;
            rho=0.0;
         }
         else if(nb_adj_gene1>2){
            cout<<"Gene "<<nom_gene1<<" has "<<nb_adj_gene1<<" neighboors!!!"<<endl;
            rho=0.0;
         }
         else if(nb_adj_gene2>2){
            cout<<"Gene "<<nom_gene2<<" has "<<nb_adj_gene2<<" neighboors!!!"<<endl;
            rho=0.0;
         }
         if (rho==-1.0){
            cout<<"From CostFunctionsProba.cpp (C0ExtantExtant): ERROR in rho determination"<<endl;
            exit(EXIT_FAILURE);
         }

         proba=rho*esp_nb_contig[species][3];
      }
      float base_log=esp_nb_contig[species][2];
      pair<float,int> p (-log_b(base_log,(1.0-proba)),10);
      return p;
   }
   else{
      if (appAdj(a12,Adj_classe) || appAdj(a21,Adj_classe)){
         pair<float,int> p (INFINI,10);
         return p;
      }
      else{
         pair<float,int> p (0.0,10);
         return p;
      }
   }
}

pair<float,int>  C0GLosAny(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0GLosAny"<<endl;
   pair<float,int> p (0.0,20);
   return p;
}


pair<float,int>  C0GLosGLos(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0GLosGLos"<<endl;
   pair<float,int> p (0.0,30);
   return p;
}

pair<float,int>  C0GDupExtantOrSpec(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0GDupExtantOrSpec"<<endl;

   pair<float,int> a=minPaire(
     recupCoutC0P(v1->getSon(0),v2)+recupCoutC0P(v1->getSon(1),v2),//41
     recupCoutC1P(v1->getSon(0),v2)+recupCoutC0P(v1->getSon(1),v2)+Crea,//42
     recupCoutC0P(v1->getSon(0),v2)+recupCoutC1P(v1->getSon(1),v2)+Crea,//43
     recupCoutC1P(v1->getSon(0),v2)+recupCoutC1P(v1->getSon(1),v2)+2*Crea);//44

   a.second+=40;
   
   if(affich_CalculCout)
      cout <<"Je retourne la paire "<<a.first<< " "<<a.second<<endl;
   
   return a;
}

pair<float,int>  C0SpecSpec(Node * v1, Node * v2)
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

   if(affich_CalculCout)
      cout <<"fgV1="<<fgV1->getId()<<" fdV1="<<fdV1->getId()<<" fgV2="<<fgV2->getId()<<" fdV2="<<fdV2->getId()<<endl;

   pair<float,int> a=minPaire(
     recupCoutC0P(fgV1,fgV2) + recupCoutC0P(fdV1,fdV2),//51
     recupCoutC1P(fgV1,fgV2) + recupCoutC0P(fdV1,fdV2) + Crea,//52
     recupCoutC0P(fgV1,fgV2) + recupCoutC1P(fdV1,fdV2) + Crea,//53
     recupCoutC1P(fgV1,fgV2) + recupCoutC1P(fdV1,fdV2) + 2*Crea);//54

   a.second+=50;
   
   if(affich_CalculCout)
      cout <<"Je retourne la paire "<<a.first<< " "<<a.second<<endl;
   
   return a;
}

pair<float,int>  C0GDupGDup(Node * v1, Node * v2)
{
   if(affich_CalculCout)
      cout <<"Je suis dans C0GDupGDup"<<endl;

   pair<float,int> a=minPaire3x233(
     //Les 2 cas originaux d'Éric
     minPaire(recupCoutC0P(v1->getSon(0),v2) + recupCoutC0P(v1->getSon(1),v2),//61
         recupCoutC0P(v1,v2->getSon(0)) + recupCoutC0P(v1,v2->getSon(1))),//62
     //3 cas où v1 se duplique avant v2
     minPaire(recupCoutC0P(v1->getSon(0),v2) + recupCoutC1P(v1->getSon(1),v2)+Crea,//63
         recupCoutC1P(v1->getSon(0),v2) + recupCoutC0P(v1->getSon(1),v2)+Crea,//64
         recupCoutC1P(v1->getSon(0),v2) + recupCoutC1P(v1->getSon(1),v2)+2*Crea),//65
     //3 cas où v2 se duplique avant v1
     minPaire(recupCoutC0P(v1,v2->getSon(0)) + recupCoutC1P(v1,v2->getSon(1))+Crea,//66
         recupCoutC1P(v1,v2->getSon(0)) + recupCoutC0P(v1,v2->getSon(1))+Crea,//67
         recupCoutC1P(v1,v2->getSon(0)) + recupCoutC1P(v1,v2->getSon(1))+2*Crea));//68
   
   a.second+=60;
   
   if(affich_CalculCout)
      cout <<"Je retourne la paire "<<a.first<< " "<<a.second<<endl;
   
   return a;
}

// C1P
//
//On suppose que les noeuds sont bien de la même espèce
pair<float,int> calculeC1P(Node * v1, Node * v2, vector<adjacence> * Adj_classe, TreeTemplate<Node> * S, map<string,string> &gene_species_EXT, map<string,vector<float> > &esp_nb_contig){
   switch (EvtToInt(E(v1))){  
   case 1 :
      switch (EvtToInt(E(v2))){
      case 1 :
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 Extant"<<endl;
    return C1ExtantExtant(v1,v2,Adj_classe,S,gene_species_EXT,esp_nb_contig);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 GLos"<<endl;
    return C1GLosAny(v2,v1);
    break;    
      case 2 :
    cout<<"From CostFunctionsProba.cpp (calculeC1P): v1 Extant and v2 Spec IMPOSSIBLE"<<endl;
    exit(EXIT_FAILURE);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 GDup"<<endl;
    return C1GDupExtantOrSpec(v2,v1);    
    break;
      default: 
    cout<<"From CostFunctionsProba.cpp (calculeC1P): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
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
    cout<<"From CostFunctionsProba.cpp (calculeC1P): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   case 2 :
      switch (EvtToInt(E(v2))){
      case 1 :
    cout<<"From CostFunctionsProba.cpp (calculeC1P): v1 Spec and v2 Extant IMPOSSIBLE"<<endl;
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
    cout<<"From CostFunctionsProba.cpp (calculeC1P): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
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
    cout<<"From CostFunctionsProba.cpp (calculeC1P): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   default: 
      cout<<"From CostFunctionsProba.cpp (calculeC1P): ERROR the call to EvtToInt(E(v1)) on node v1 with id "<<v1->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v1))<<endl;
      exit(EXIT_FAILURE);
      break;  
   }
}

// C0P
//
//On suppose que les noeuds sont bien de la même espèce
pair<float,int> calculeC0P(Node * v1, Node * v2, vector<adjacence> * Adj_classe, TreeTemplate<Node> * S, map<string,string> &gene_species_EXT, map<string,vector<float> > &esp_nb_contig){
   switch (EvtToInt(E(v1))){  
   case 1 :
      switch (EvtToInt(E(v2))){
      case 1 :
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 Extant"<<endl;
    return C0ExtantExtant(v1,v2,Adj_classe,S,gene_species_EXT,esp_nb_contig);
    break;
      case 4 :    
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 GLos"<<endl;
    return C0GLosAny(v2,v1);
    break;    
      case 2 :
    cout<<"From CostFunctionsProba.cpp (calculeC0P): v1 Extant and v2 Spec IMPOSSIBLE"<<endl;
    exit(EXIT_FAILURE);
    break;    
      case 3 :
    if(affich_CalculCout)
       cout<<"v1 Extant et v2 GDup"<<endl;
    return C0GDupExtantOrSpec(v2,v1);    
    break;
      default: 
    cout<<"From CostFunctionsProba.cpp (calculeC0P): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
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
    cout<<"From CostFunctionsProba.cpp (calculeC0P): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   case 2 :
      switch (EvtToInt(E(v2))){
      case 1 :
    cout<<"From CostFunctionsProba.cpp (calculeC1P): v1 Spec and v2 Extant IMPOSSIBLE"<<endl;
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
    cout<<"From CostFunctionsProba.cpp (calculeC0P): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
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
    cout<<"From CostFunctionsProba.cpp (calculeC0P): ERROR the call to EvtToInt(E(v2)) on node v2 with id "<<v2->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v2))<<endl;
    exit(EXIT_FAILURE);
    break;  
      }
      break;
      
   default: 
      cout<<"From CostFunctionsProba.cpp (calculeC0P): ERROR the call to EvtToInt(E(v1)) on node v1 with id "<<v1->getId()<<" doesn't return one of the know events, it returns:"<<EvtToInt(E(v1))<<endl;
      exit(EXIT_FAILURE);
      break;  
   }
}
