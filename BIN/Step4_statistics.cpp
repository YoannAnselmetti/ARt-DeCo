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

File: Step4_statistics.cpp                     Last modified on: 16/06/2015
Created by: Yoann Anselmetti/Sèverine Bérard   Created on: 10/03/2014
--------------------------------------------------------------------------
Specification: 
File for statistics analysis
=========================================================================*/
#include "Step4_statistics.h"



//renvoie l'événement associé au noeud n
string E(Node * n){
   if (n->hasBranchProperty(typ)){
      BppString * Evt = dynamic_cast<BppString*> (n->getBranchProperty(typ));
      return Evt->toSTL();
   }
   else{
      cout<<"From Step4_statistics.cpp (E): node with id "<<n->getId()<<" has no associated event !!"<<endl; 
      return "ND";
   }
}

//Renvoie vrai si la compo connexe v est linéaire
//Càd si tous ses sommets sont de degrès 2 sauf 2 qui sont de degrès 1
bool isLinear(vector<string> v, map<string,int> deg){
   vector<string>::iterator it_v;
   bool linear = true;
   int nb_sommets_deg_1=0;
   for(it_v=v.begin();it_v!=v.end();it_v++){
      if (deg[*it_v]>2)
         linear=false;
      else if (deg[*it_v]==1)
         nb_sommets_deg_1++;
   }
   return linear&&(nb_sommets_deg_1==2);
}


vector<adjacence> composantesConnexes(vector<adjacence> adj_vector, map<int,int> &distrib_taille_compo_connexe, map<int, vector<string> > &numCC_listAdj, int i){
   //Pour calculer les degrés
   map<string,int> deg;
   vector<adjacence>::iterator it_v;
   map<string,int>::iterator it_m;
   map<string,int> corres_gene_numCompo;
   map<string,int>::iterator it_adj;
   map<string,int> corres_adj_gene_numCompo;
   string adj;
   bool compo_lineaire=true;

   for(it_v=adj_vector.begin();it_v!=adj_vector.end();it_v++){
      adj=(*it_v).gene1+"\t"+(*it_v).gene2;
      if (corres_gene_numCompo.find((*it_v).gene1)!=corres_gene_numCompo.end() && corres_gene_numCompo.find((*it_v).gene2)!=corres_gene_numCompo.end() && corres_gene_numCompo[(*it_v).gene2]!=corres_gene_numCompo[(*it_v).gene1]){
         //fusion des deux composantes connexes !
         int c1=corres_gene_numCompo[(*it_v).gene1];
         int c2=corres_gene_numCompo[(*it_v).gene2];
         for(it_m=corres_gene_numCompo.begin();it_m!=corres_gene_numCompo.end();it_m++)
            if ((*it_m).second==c1)
               (*it_m).second=c2;
         for(it_adj=corres_adj_gene_numCompo.begin();it_adj!=corres_adj_gene_numCompo.end();it_adj++)
            if ((*it_adj).second==c1)
               (*it_adj).second=c2;
         corres_adj_gene_numCompo[adj]=c2;
      }
      else if (corres_gene_numCompo.find((*it_v).gene1)!=corres_gene_numCompo.end() && corres_gene_numCompo.find((*it_v).gene2)!=corres_gene_numCompo.end()){
    //Coup de bol, elles ont déjà été classées dans la même compo connexe
    //Rien à faire !
      }    
      else if (corres_gene_numCompo.find((*it_v).gene1)!=corres_gene_numCompo.end()){
         corres_gene_numCompo[(*it_v).gene2]=corres_gene_numCompo[(*it_v).gene1];
         corres_adj_gene_numCompo[adj]=corres_gene_numCompo[(*it_v).gene1];
      }
      else if (corres_gene_numCompo.find((*it_v).gene2)!=corres_gene_numCompo.end()){
         corres_gene_numCompo[(*it_v).gene1]=corres_gene_numCompo[(*it_v).gene2];
         corres_adj_gene_numCompo[adj]=corres_gene_numCompo[(*it_v).gene2];
      }
      else{ //nouvelle compo connexe
         corres_adj_gene_numCompo[adj]=i;
         corres_gene_numCompo[(*it_v).gene1]=i;
         corres_gene_numCompo[(*it_v).gene2]=i;
         i++;
      }
      //Pour les degrés
      if (deg.find((*it_v).gene1)!=deg.end())
         deg[(*it_v).gene1]++;
      else
         deg[(*it_v).gene1]=1;
      if (deg.find((*it_v).gene2)!=deg.end())
         deg[(*it_v).gene2]++;
      else
         deg[(*it_v).gene2]=1;
   }
   
   // cout<<"Affichage des degrés :"<<endl;
   // afficheMap(deg);
   
   for(it_adj=corres_adj_gene_numCompo.begin();it_adj!=corres_adj_gene_numCompo.end();it_adj++){
      if(numCC_listAdj.find((*it_adj).second)!=numCC_listAdj.end())
         numCC_listAdj[(*it_adj).second].push_back((*it_adj).first);
      else{
         vector<string> v;
         v.push_back((*it_adj).first);
         numCC_listAdj[(*it_adj).second]=v;
      }
   }

   map<int,vector<string> > compoConnexe;
   map<int,vector<string> >::iterator it_c;
   for(it_m=corres_gene_numCompo.begin();it_m!=corres_gene_numCompo.end();it_m++){
      if(compoConnexe.find((*it_m).second)!=compoConnexe.end())
         compoConnexe[(*it_m).second].push_back((*it_m).first);
      else{
         vector<string> v;
         v.push_back((*it_m).first);
         compoConnexe[(*it_m).second]=v;
      }
   }
   int nb_compo_taille2=0;
   //Tous les gènes dans une composante de taille 2 vont être rangé dans la compo 0
   for(it_c=compoConnexe.begin();it_c!=compoConnexe.end();it_c++){
      int taille=((*it_c).second).size();
      if (taille==2){
         nb_compo_taille2++;
         for(int j=0;j<taille;j++)
            corres_gene_numCompo[((*it_c).second)[j]]=0;
      }
      if (distrib_taille_compo_connexe.find(taille)!=distrib_taille_compo_connexe.end())
         distrib_taille_compo_connexe[taille]++;
      else
         distrib_taille_compo_connexe[taille]=1;
      compo_lineaire=compo_lineaire&&isLinear((*it_c).second,deg);
   }
   
   cout<<"There are "<<compoConnexe.size()<<" connected components (of size >1)."<<endl;
   cout<<"\t"<<nb_compo_taille2<<" of size =2."<<endl;
   cout<<"\t"<<compoConnexe.size()-nb_compo_taille2<<" of size >2."<<endl;
   cout<<"Are all the components linear? "<<compo_lineaire<<endl;
   //On garde de adj_vector tous les gènes impliqués dans des composantes de taille > 2
   vector<adjacence> adj_vector_sup2;
   for(it_v=adj_vector.begin();it_v!=adj_vector.end();it_v++)
      if(corres_gene_numCompo[(*it_v).gene1]!=0)
         adj_vector_sup2.push_back(*it_v);   
   return adj_vector_sup2;
}


void File_DOT_for_Adj_Graph(map<int,vector<adjacence> > Adj_actuelles, map<int,vector<adjacence> > Adj_nouvelles, TreeTemplate<Node> * S){
   cout<<"OK"<<endl;
   convert<int> c;
   //Création dossier dans lequel seront stockés les graphes d'adjacence
   string dossier_graph=prefixe+output_graph+"/DOT";

   const char * DOT = dossier_graph.c_str();
   mkdir_rec(DOT);

//   string commande="mkdir -p "+dossier_graph;
//   int ret=system(commande.c_str());
//   if(ret<0){
//      cout<<"From Step4_statistics.cpp (File_DOT_for_Adj_Graph): ERROR while creating directory "<<dossier_graph<<endl;
//      if(file_log)
//    Offile_log<<"From Step4_statistics.cpp (File_DOT_for_Adj_Graph): ERROR while creating directory "<<dossier_graph<<endl;
//      exit(EXIT_FAILURE);
//   }
   
   map<int,vector<adjacence> >::iterator it_ori;
   //Parcours de la map Adj_actuelles pour créer le graphe des adjacences
   for(it_ori=Adj_actuelles.begin();it_ori!=Adj_actuelles.end();it_ori++){
      //Récupération nom_espèce
      string species=S->getNode((*it_ori).first)->getName();
      
      //Parcours du vector<adj> de l'esp X pour créer graphe d'Adj
      vector<adjacence> adj_ori=(*it_ori).second;
      vector<adjacence>::iterator it_adj_ori;
      map<string,int> gene_nodeID;
      map<string,int>::iterator it_name_id;
      vector<string> v_edge_dot, v_adj_dot;
      int id_edge=0;
      int id_node=0;
      for(it_adj_ori=adj_ori.begin();it_adj_ori!=adj_ori.end();it_adj_ori++){
         string common_edge_dot="";
         id_edge++;
         // Stockage gènes/noeuds et adjacences/branches dans map et vecteurs
         // 1er gène/noeud
         //Si déjà dans nos vecteurs
         if(gene_nodeID.find((*it_adj_ori).gene1)!=gene_nodeID.end()){
            common_edge_dot=c.to_s(gene_nodeID[(*it_adj_ori).gene1])+" -> ";
         }
         //Si pas encore dans nos vecteurs
         else{
            id_node++;
            gene_nodeID[(*it_adj_ori).gene1]=id_node;
            common_edge_dot=c.to_s(id_node)+" -> ";
         }
         // 1er gène/noeud
         //Si déjà dans nos vecteurs
         if(gene_nodeID.find((*it_adj_ori).gene2)!=gene_nodeID.end()){
            common_edge_dot+=c.to_s(gene_nodeID[(*it_adj_ori).gene2]);
         }
         //Si pas encore dans nos vecteurs
         else{
            id_node++;
            gene_nodeID[(*it_adj_ori).gene2]=id_node;
            common_edge_dot+=c.to_s(id_node);
         }
         v_edge_dot.push_back(common_edge_dot);
      }
      int first_new_node=id_node+1;
      try{
         vector<adjacence> adj_new=Adj_nouvelles[(*it_ori).first];
         vector<adjacence>::iterator it_adj_new;
         int nb_new_adj=0;
         for(it_adj_new=adj_new.begin();it_adj_new!=adj_new.end();it_adj_new++){
            id_edge++;
            nb_new_adj++;
            string new_edge_dot="";
            if(gene_nodeID.find((*it_adj_new).gene1)!=gene_nodeID.end()){
               new_edge_dot=c.to_s(gene_nodeID[(*it_adj_new).gene1])+" -> ";
            }
            else{
               id_node++;
               gene_nodeID[(*it_adj_new).gene1]=id_node;
               new_edge_dot=c.to_s(id_node)+" -> ";
            }
            if(gene_nodeID.find((*it_adj_new).gene2)!=gene_nodeID.end()){
               new_edge_dot+=c.to_s(gene_nodeID[(*it_adj_new).gene2]);
            }
            else{
               id_node++;
               gene_nodeID[(*it_adj_new).gene2]=id_node;
               new_edge_dot+=c.to_s(id_node);
            }
            v_adj_dot.push_back(new_edge_dot);
         }
         cout<<species<<" has "<<nb_new_adj<<" new adjacencies."<<endl;
         if(file_log)
            Offile_log<<"OS\t"<<species<<" has "<<nb_new_adj<<" new adjacencies."<<endl;
      }
      catch(exception e){
         //cout<<e.what();
         cout<<species<<" has not new adjacencies."<<endl;
      }
      
      ////////////////////////////////////////////////////////
      ///// Création du fichier de graphe au format DOT   ////
      ////////////////////////////////////////////////////////
      
      //Création du fichier DOT pour l'espèce "species"
      ofstream OffilDot;
      string fic_dot=dossier_graph+"/Adj_graph_"+species+".dot";
      OffilDot.open(fic_dot.c_str(), ios::out|ios::trunc);
      if (!OffilDot){
         cout<<"From Step4_statistics.cpp (File_DOT_for_Adj_Graph): ERROR while opening file "<<fic_dot<<endl;
    if (file_log)
       Offile_log<<"From Step4_statistics.cpp (File_DOT_for_Adj_Graph): ERROR while opening file "<<fic_dot<<endl;
    exit(EXIT_FAILURE);
      }
      
      //Auteur et commentaires
      OffilDot<<"// Author(s): \"Yoann Anselmetti\")"<<endl;
      OffilDot<<"// This file was generated by DeCo. Adjacency graph for species "<<species<<endl<<endl;
      //Date
      struct tm Today;
      time_t maintenant;
      time(&maintenant);
      Today = *localtime(&maintenant);
      OffilDot<<"// Date:";
      OffilDot<<c.to_s(Today.tm_year + 1900)<<"-"<<c.to_s(Today.tm_mon + 1)<<"-"<<c.to_s(Today.tm_mday)<<" "<<c.to_s(Today.tm_hour)<<":"<<c.to_s(Today.tm_min)<<":"<<c.to_s(Today.tm_sec)<<endl<<endl;

      //Entête fichier DOT
      OffilDot<<"digraph "+species+" { "<<endl;
      OffilDot<<"\tedge [color=black];"<<endl; 

      //Affectation des branches (adjacences)
      vector<string>::iterator it_edge_dot;
      for(it_edge_dot=v_edge_dot.begin();it_edge_dot!=v_edge_dot.end();it_edge_dot++){
         OffilDot<<"\t"<<(*it_edge_dot)<<";"<<endl;
      }
      // Affectation couleur rouge aux nouvelles branches (adjacences)
      vector<string>::iterator it_adj_dot;
      for(it_adj_dot=v_adj_dot.begin();it_adj_dot!=v_adj_dot.end();it_adj_dot++){
         OffilDot<<"\t"<<(*it_adj_dot)<<" [color=red, arrowsize=5];"<<endl;
      }
      //Affectation nom gène aux noeuds
      for(it_name_id=gene_nodeID.begin();it_name_id!=gene_nodeID.end();it_name_id++){
         if((*it_name_id).second>=first_new_node)
            OffilDot<<"\t"<<(*it_name_id).second<<" [color=red, style=filled, shape=box, label=\""<<(*it_name_id).first<<"\"];"<<endl;
         else            
            OffilDot<<"\t"<<(*it_name_id).second<<" [label=\""<<(*it_name_id).first<<"\"];"<<endl;
      }
      OffilDot<<"}";
      OffilDot.close();
   }
   cout<<"Creation of extant adjacencies graph for each species to DOT format ...DONE"<<endl<<endl;
}

void File_TULIP_for_Adj_Graph(map<int,vector<adjacence> > Adj_actuelles, map<int,vector<adjacence> > Adj_nouvelles, TreeTemplate<Node> * S){
   cout<<"OK"<<endl;
   convert<int> c;
   //Création dossier dans lequel seront stockés les graphes d'adjacences
   string dossier_graph=prefixe+output_graph+"/TULIP";
   
   const char * TULIP = dossier_graph.c_str();
   mkdir_rec(TULIP);

//   string commande="mkdir -p "+dossier_graph;
//   int ret=system(commande.c_str());
//   if(ret<0){
//      cout<<"From Step4_statistics (File_TULIP_for_Adj_Graph): ERROR while creating directory "<<dossier_graph<<endl;
//      if(file_log)
//    Offile_log<<"From Step4_statistics (File_TULIP_for_Adj_Graph): ERROR while creating directory "<<dossier_graph<<endl;
//      exit(EXIT_FAILURE);
//   }

   map<int,vector<adjacence> >::iterator it_ori;
   //Parcours de la map Adj_actuelles pour créer le graphe des adjacences
   for(it_ori=Adj_actuelles.begin();it_ori!=Adj_actuelles.end();it_ori++){
      //Récupération nom_espèce
      string species=S->getNode((*it_ori).first)->getName();
      
      //Parcours du vector<adj> de l'esp X pour créer graphe d'Adj
      vector<adjacence> adj_ori=(*it_ori).second;
      vector<adjacence>::iterator it_adj_ori;
      map<string,int> gene_nodeID;
      map<string,int>::iterator it_name_id;
      vector<string> v_edge_tlp, v_adj_tlp;
      int id_edge=0;
      int id_node=0;
      for(it_adj_ori=adj_ori.begin();it_adj_ori!=adj_ori.end();it_adj_ori++){
         string common_edge_tlp="";
         id_edge++;
         // Stockage gènes/noeuds et adjacences/branches dans map et vecteurs
         // 1er gène/noeud
         //Si déjà dans nos vecteurs
         if(gene_nodeID.find((*it_adj_ori).gene1)!=gene_nodeID.end()){
            common_edge_tlp="(edge "+c.to_s(id_edge)+" "+c.to_s(gene_nodeID[(*it_adj_ori).gene1])+" ";
         }
         //Si pas encore dans nos vecteurs
         else{
            id_node++;
            gene_nodeID[(*it_adj_ori).gene1]=id_node;
            common_edge_tlp="(edge "+c.to_s(id_edge)+" "+c.to_s(id_node)+" ";
         }
         // 1er gène/noeud
         //Si déjà dans nos vecteurs
         if(gene_nodeID.find((*it_adj_ori).gene2)!=gene_nodeID.end()){
            common_edge_tlp+=c.to_s(gene_nodeID[(*it_adj_ori).gene2])+")";
         }
         //Si pas encore dans nos vecteurs
         else{
            id_node++;
            gene_nodeID[(*it_adj_ori).gene2]=id_node;
            common_edge_tlp+=c.to_s(id_node)+")";
         }
         v_edge_tlp.push_back(common_edge_tlp);
      }
      int first_new_edge=id_edge+1;
      int first_new_node=id_node+1;
      try{
         vector<adjacence> adj_new=Adj_nouvelles[(*it_ori).first];
         vector<adjacence>::iterator it_adj_new;
         int nb_new_adj=0;
         for(it_adj_new=adj_new.begin();it_adj_new!=adj_new.end();it_adj_new++){
            id_edge++;
            nb_new_adj++;
            string new_edge_tlp="";
            if(gene_nodeID.find((*it_adj_new).gene1)!=gene_nodeID.end()){
               new_edge_tlp="(edge "+c.to_s(id_edge)+" "+c.to_s(gene_nodeID[(*it_adj_new).gene1])+" ";
            }
            else{
               id_node++;
               gene_nodeID[(*it_adj_new).gene1]=id_node;
               new_edge_tlp="(edge "+c.to_s(id_edge)+" "+c.to_s(id_node)+" ";
            }
            if(gene_nodeID.find((*it_adj_new).gene2)!=gene_nodeID.end()){
               new_edge_tlp+=c.to_s(gene_nodeID[(*it_adj_new).gene2])+")";
            }
            else{
               id_node++;
               gene_nodeID[(*it_adj_new).gene2]=id_node;
               new_edge_tlp+=c.to_s(id_node)+")";
            }
            v_adj_tlp.push_back(new_edge_tlp);
         }
         cout<<species<<" has "<<nb_new_adj<<" new adjacencies."<<endl;
         if(file_log)
            Offile_log<<"OS\t"<<species<<" has "<<nb_new_adj<<" new adjacencies."<<endl;
            
      }
      catch(exception e){
         //cout<<e.what();
         cout<<species<<" has not new adjacencies."<<endl;
      }
      
      ///////////////////////////////////////////////////////////
      /////   Création du fichier de graphe au format TULIP   ///// 
      ///////////////////////////////////////////////////////////
      
      //Création du fichier TULIP pour l'espèce "species"
      ofstream OffilTulip;
      string fic_tulip=dossier_graph+"/Adj_graph_"+species+".tlp";
      OffilTulip.open(fic_tulip.c_str(), ios::out|ios::trunc);
      if (!OffilTulip){
         cout<<"From Step4_statistics.cpp (File_TULIP_for_Adj_Graph): ERROR while opening file "<<fic_tulip<<endl;
         if (file_log)
            Offile_log<<"\nFI\tFrom Step4_statistics.cpp (File_TULIP_for_Adj_Graph): ERROR while opening file "<<fic_tulip<<endl;
         exit(EXIT_FAILURE);
      }

      //entête fichier TULIP
      OffilTulip<<"(tlp \"2.3\""<<endl;

      //Date
      struct tm Today;
      time_t maintenant;
      time(&maintenant);
      Today = *localtime(&maintenant);
      OffilTulip<<"(date \"";
      OffilTulip<<c.to_s(Today.tm_year + 1900)<<"-"<<c.to_s(Today.tm_mon + 1)<<"-"<<c.to_s(Today.tm_mday)<<" "<<c.to_s(Today.tm_hour)<<":"<<c.to_s(Today.tm_min)<<":"<<c.to_s(Today.tm_sec)<<"\")"<<endl;

      //Auteur et commentaires
      OffilTulip<<"(author \"Yoann Anselmetti\")"<<endl;
      OffilTulip<<"(comments \"This file was generated by DeCo. Adjacency graph for species "<<species<<"\")"<<endl;

      OffilTulip<<"(nb_nodes "<<id_node<<")"<<endl;
      OffilTulip<<";(nodes <node_id> <node_id> ...)"<<endl;
      OffilTulip<<"(nodes 1.."<<id_node<<")"<<endl;

      OffilTulip<<"(nb_edges "<<c.to_s(id_edge)<<")"<<endl;
      OffilTulip<<";(edge <edge_id> <source_id> <target_id>)"<<endl;
      vector<string>::iterator it_edge;
      for(it_edge=v_edge_tlp.begin();it_edge!=v_edge_tlp.end();it_edge++){
         OffilTulip<<(*it_edge)<<endl;
      }
      vector<string>::iterator it_adj;
      for(it_adj=v_adj_tlp.begin();it_adj!=v_adj_tlp.end();it_adj++){
         OffilTulip<<(*it_adj)<<endl;
      }

      //Propriété label pour la couleur des noeuds et des branches
      OffilTulip<<"(property  0 color \"viewColor\""<<endl;
      OffilTulip<<"(default \"(0,0,255,255)\" \"(0,0,255,255)\")"<<endl;
      int i,j;
      for(i=first_new_node;i<=id_node;i++){
         OffilTulip<<"(node "<<i<<" \"(255,0,0,255)\" \"(255,0,0,255)\")"<<endl;
      }
      for(j=first_new_edge;j<=id_edge;j++){
         OffilTulip<<"(edge "<<j<<" \"(255,0,0,255)\" \"(255,0,0,255)\")"<<endl;
      }
      OffilTulip<<")"<<endl<<endl;

      //La propriété label pour les noms des noeuds
      OffilTulip<<"(property  0 string \"viewLabel\""<<endl;
      OffilTulip<<"(default \"\" \"\" )"<<endl;
      for(it_name_id=gene_nodeID.begin();it_name_id!=gene_nodeID.end();it_name_id++)
         OffilTulip<<"(node "<<(*it_name_id).second<<" \""<<(*it_name_id).first<<"\")"<<endl;
      OffilTulip<<")"<<endl;
      OffilTulip<<")";
      OffilTulip.close();
   }
   cout<<"Creation of extant adjacencies graph for each species to TULIP format ...DONE"<<endl<<endl;
}

//Pour extraire le nom du gène du nom de la feuille de l'arbre
string nomGene(string nom){
   return nom.substr(0,nom.find(sep));
}

//Même fonction qu'au dessus mais on tronque le nom de la feuille pour ne garder que celui du gène
void retourneFeuillesNomsGenes(Node * n, vector<string> & feuilles){
   if (n->getNumberOfSons()==0){
      if (nomGene(n->getName())!=NomPer)
         feuilles.push_back(nomGene(n->getName()));
   }
   else if (n->getNumberOfSons()==1)
      retourneFeuillesNomsGenes(n->getSon(0),feuilles);
   else{ //2 fils
      retourneFeuillesNomsGenes(n->getSon(0),feuilles);
      retourneFeuillesNomsGenes(n->getSon(1),feuilles);
   }
}

//num espèce puis nom gène puis liste des descendants
void ecritGenes(Node * n, int num_arbre, ofstream & Offic_SORTIE_genes, map<string,int> &gene_species_ANC){
   convert<int> c;
   //Le numéro de l'espèce est son identifiant dans l'arbre des espèces
   BppString * ESPECE = dynamic_cast <BppString *> (n->getNodeProperty(esp));
   BppString * TYP = dynamic_cast<BppString*> (n->getBranchProperty(typ));
   
   //Nom si c'est une feuille ou code du gène si c'est un noeud de spéciation
   if(n->getNumberOfSons()==0){
      if (TYP->toSTL()!=per){
         Offic_SORTIE_genes<<c.from_s(ESPECE->toSTL())<<"\t";
         Offic_SORTIE_genes<<nomGene(n->getName());
         Offic_SORTIE_genes<<"\n";
      }
   }
   else{
      if(TYP->toSTL()==spe){
         stringstream ss;
         ss << num_arbre<<sep<<n->getId();
         string anc_gene = ss.str();
//         cout<<anc_gene<<endl;
         gene_species_ANC[anc_gene]=c.from_s(ESPECE->toSTL());
         Offic_SORTIE_genes<<c.from_s(ESPECE->toSTL())<<"\t";
         Offic_SORTIE_genes<<anc_gene;
         //Les descendants de ce gène
         vector<string> feuilles;
         retourneFeuillesNomsGenes(n,feuilles);
         for(unsigned int i=0;i<feuilles.size();i++)
            Offic_SORTIE_genes<<"\t"<<feuilles[i];
         Offic_SORTIE_genes<<"\n";
      }
   }
   
   if (n->getNumberOfSons()==1)
      ecritGenes(n->getSon(0),num_arbre,Offic_SORTIE_genes,gene_species_ANC);
   else if (n->getNumberOfSons()==2){
      ecritGenes(n->getSon(0),num_arbre,Offic_SORTIE_genes,gene_species_ANC);
      ecritGenes(n->getSon(1),num_arbre,Offic_SORTIE_genes,gene_species_ANC);
   }
}

//Pour compter le nombre de noeud ayant la propriété prop à la valeur val
int countNodesWithBranchProperty(Node * n, string prop, string val){
   int p=0;
   if (n->hasBranchProperty(prop)){
      BppString * VAL = dynamic_cast<BppString*> (n->getBranchProperty(prop)); 
      if(VAL->toSTL()==val)
         p=1;
   }
   if (n->getNumberOfSons()==0)
      return p;
   else if (n->getNumberOfSons()==1)
      return p+countNodesWithBranchProperty(n->getSon(0),prop,val);
   else //2 fils 
      return p
    +countNodesWithBranchProperty(n->getSon(0),prop,val)
    +countNodesWithBranchProperty(n->getSon(1),prop,val);
   
}

void lectureFicRelBin(string fic,map<int,vector<adjacence> > & miva){
   ifstream Iffic;
   string tampon;
   convert<int> c;
   Iffic.open(fic.c_str(), ios::in);
   if (!Iffic){
      cout<<"\nFrom Step4_statistics.cpp (lectureFicRelBin): ERROR while opening file "<<fic<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step4_statistics.cpp (lectureFicRelBin): ERROR while opening file "<<fic<<endl;
      exit(EXIT_FAILURE);
   } 
   while (!Iffic.eof()){
      //num espèce
      Iffic>>tampon;
      int e=c.from_s(tampon);
      //Pour éviter un dernier tour de boucle inutile à cause d'un blanc en fin de fichier
      if(!Iffic.eof()){
         //Les deux gènes adjacents
         adjacence Adj;
         Iffic>>tampon;
         Adj.gene1=tampon;
         Iffic>>tampon;
         Adj.gene2=tampon;
    
         //On entre ça dans la map 
         if (miva.find(e)!=miva.end())
            miva[e].push_back(Adj);
         else{
            vector<adjacence> v;
            v.push_back(Adj);
            miva[e]=v;
         }
      }  
   }
   Iffic.close();
}

int main(int argc, char* argv[]){
   cout<<"\n\t##########################"<<endl;
   cout<<"\t###  Start Step4 !!!!  ###"<<endl;
   cout<<"\t##########################"<<endl<<endl;
   time_t tbegin=time(NULL);              // get the current calendar time

   /////////////////////////////////////////////////////////////////////////////
   // Récupérer données fichier de config pour recréer gene_species_EXT
   // et pour récupérer le préfixe pour rechercher le fichier de
   // classes d'adjacences
   lireFicConfig(fic_arbre,fic_gene,fic_especes,fic_adjacence,exp_name,directory,argc,argv);
   if (is_readable(fic_adjacence_clean))
      fic_adjacence=fic_adjacence_clean;
   //Idem pour le fichier de gènes
   if (is_readable(fic_gene_clean))
      fic_gene=fic_gene_clean;
   
   if (file_log){
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
            cout<<endl<<"From Step4_statistics.cpp (main): ERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file open by Step4_statistics.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      else{
         Offile_log.open(log_file.c_str(), ios::out);
         if(!Offile_log){
            cout<<"\nERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file created by Step4_statistics.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      Offile_log<<"\t#########################"<<endl;
      Offile_log<<"\t###  Start Step4 !!!  ###"<<endl;
      Offile_log<<"\t#########################"<<endl<<endl;
      Offile_log<<"\tFile(s) used in Step3_proba:"<<endl;
      Offile_log<<"\t\t- Reconciled trees file: "<<tree_reconciled_file<<endl;
      //Offile_log<<"\t\t- Adjacencies file: "<<fic_adjacence<<endl;
      Offile_log<<"\t\t- Selected species file: "<<file_spe<<endl<<endl;
   }

   cout<<endl<<"File(s) used in Step4:"<<endl;
   cout<<"\t- Reconciled trees file: "<<tree_reconciled_file<<endl;
   //cout<<"\t- Adjacencies file: "<<fic_adjacence<<endl;
   cout<<"\t- Selected species file: "<<file_spe<<endl<<endl;


   // Variables
   //   TreeTemplate<Node> S;
   //   S = new TreeTemplate<Node>;
   map<string,string> gene_species_EXT;
   vector<Tree*> Arbres;
   map<string,vector<int> > gene_GeneTreeID;
   map<string,vector<int> > classes_adjacences;
   vector<adjacence> adjacencies;
   map<int,int> extant_species_gene_nb;
   map<int,vector<adjacence> > Adj_actuelles;
   
   ////////////////////////////////////////////////////////////////////////
   // Store species tree in TreeTemplate<Node> *S
   cout<<"Store Species tree in TreeTemplate<Node> * S..."<<flush;
   if(file_log)
      Offile_log<<"OS\tFrom Step4_statistics (main): Store Species tree in TreeTemplate<Node> * S..."<<flush;
   Newick * newickReaderEsp = new Newick(true,true);
   newickReaderEsp->enableExtendedBootstrapProperty(esp);
   TreeTemplate<Node> *S = newickReaderEsp->read(fic_especes);
   delete(newickReaderEsp);
   cout<<" DONE"<<endl<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;

   /////////////////////////////////
   /////   START Genes Clean   /////
   /////////////////////////////////


   ////////////////////////////////////////////////////////////////////////
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
      cout<<"\nFrom Step4_statistics (main): WARNING while cleaning "<<elim<<" genes have been deleted"<<endl<<endl;
   }
   gene_GeneTreeID.clear();

///////////////////////////////
/////   END Genes Clean   /////
///////////////////////////////


   //Reconstruction de classes adjacences et à partir du fichier adj_class_file
   AdjClassVector(classes_adjacences);
//   afficheMap(classes_adjacences);


   string output_adj_file=prefixe+output_ext_adj_file;
   if (is_readable(output_adj_file)){
      int j=0;
      string gene1="";
      string gene2="";
      adjacence adj;      
      ifstream IfficAdj (output_adj_file.c_str(), ios::in);
      if (!IfficAdj){
         cout<<endl<<"From General.cpp: ERROR while opening file "<<output_adj_file<<endl;
         exit(EXIT_FAILURE);
      }
      cout<<"File "<<output_adj_file<<" open"<<endl<<endl;
      cout<<"Creation of the sorted vector<adjacence> adjacences..."<<flush;
      while (!IfficAdj.eof()){
         //on lit la ligne en entier;
         IfficAdj>>gene1;      
         //Pour éviter un dernier tour de trop à cause d'un blanc en fin de fichier
         if(!IfficAdj.eof()){
            IfficAdj>>gene2;      
   //         cout<<"On traite l'adjacence "<< gene1<<" - "<<gene2;
            //on crée l'adjacence si les deux gènes sont dans gene_species_EXT (c-à-d dans nos arbres)
            if (gene_species_EXT.find(gene1)!=gene_species_EXT.end() && gene_species_EXT.find(gene2)!=gene_species_EXT.end()){
   //            cout<<", on la nomme "<<j<<endl;
               adj.gene1=gene1;
               adj.gene2=gene2;
               //on la stocke
               adjacencies.push_back(adj);
               //Stockage de l'adj actuelle dans son espèce dans Adj_actuelles
               int esp_actuelle=S->getNode(gene_species_EXT[gene1])->getId();
               if(Adj_actuelles.find(esp_actuelle)!=Adj_actuelles.end())
                  Adj_actuelles[esp_actuelle].push_back(adj);
               else{
                  vector<adjacence> v;
                  v.push_back(adj);
                  Adj_actuelles[esp_actuelle]=v;
               }
            }
         }
         j++;
      }
      IfficAdj.close();
      cout<<" DONE"<<endl<<endl;
   }
   else{
      //cout<<e.what();
      //Construction du vecteur d'adjacences ordonné et de l'association Id_espèce_actuelle avec vector<adjacences actuelles>
      AdjVector_Esp_AdjExtant(S,gene_species_EXT,adjacencies,Adj_actuelles);
      //   //Affichage des adjacences actuelles / espèce
      //   afficheMap(Adj_actuelles);
   }



// FONCTIONS SPÉCIFIQUES À STEP4_STATISTICS
   
   //Création des graphes d'adjacences actuelles pour chaque espèce
   string output_new_adj_file=prefixe+output_new_adj;
   if (is_readable(output_new_adj_file)){
      //Création de la map<int, vector<adjacence>> permettant de stocker les nouvelles adjacences avec l'identifiant espèce
      map<int,vector<adjacence> > Adj_nouvelles;
      string gene1="";
      string gene2="";
      adjacence adj;
      ifstream IffilNewAdj;
      IffilNewAdj.open(output_new_adj_file.c_str(), ios::in);
      if (!IffilNewAdj)
    {
       cout<<endl<<"From Step4_statistics.cpp: ERROR while opening file "<<output_new_adj_file<<endl;
       exit(EXIT_FAILURE);
    }
      cout<<"File "<<output_new_adj_file<<" open"<<endl<<endl;
      cout<<"Creation of the sorted map<int, vector<adjacence>> Adj_nouvelles"<<flush;
      while (!IffilNewAdj.eof()){
         //on lit la ligne en entier;
         IffilNewAdj>>gene1;      
         //Pour éviter un dernier tour de trop à cause d'un blanc en fin de fichier
         if(!IffilNewAdj.eof()){
            IffilNewAdj>>gene2;      
//            cout<<"On traite l'adjacence "<< gene1<<" - "<<gene2;
            //on crée l'adjacence si les deux gènes sont dans gene_species_EXT (c-à-d dans nos arbres)
            if (gene_species_EXT.find(gene1)!=gene_species_EXT.end() && gene_species_EXT.find(gene2)!=gene_species_EXT.end()){
//               cout<<", on la nomme "<<j<<endl;
               adj.gene1=gene1;
               adj.gene2=gene2;
               //on la stocke
               adjacencies.push_back(adj);
               //Stockage de l'adj actuelle dans son espèce dans Adj_actuelles
               int esp_actuelle=S->getNode(gene_species_EXT[gene1])->getId();
               if(Adj_nouvelles.find(esp_actuelle)!=Adj_nouvelles.end())
                  Adj_nouvelles[esp_actuelle].push_back(adj);
               else{
                  vector<adjacence> v;
                  v.push_back(adj);
                  Adj_nouvelles[esp_actuelle]=v;
               }
            }
         }
      }
      IffilNewAdj.close();
      cout<<" DONE"<<endl<<endl;


      File_TULIP_for_Adj_Graph(Adj_actuelles, Adj_nouvelles, S);
      File_DOT_for_Adj_Graph(Adj_actuelles, Adj_nouvelles, S);
   }


//Récupération EXT & ANC Spe
   //Récupération des espèces actuelles
   vector<string> leaves_name_S_tree;
   leaves_name_S_tree=S->getLeavesNames();

   //Récupération des espèces ancestrales
   vector<int> inner_node_S_tree;
   inner_node_S_tree=S->getInnerNodesId();

   //Récupération de toutes les espèces
   vector<int> node_S_tree;
   node_S_tree=S->getNodesId();


// Création du fichier de sortie OUTPUT_genes
   string output_genes_file=prefixe+output_genes;   
   ofstream Offic_SORTIE_genes;
   if(!fic_SORTIE_GeneDone){
      cout<<"Creation of the OUTPUT file for the matching of genes: "<<output_genes_file<<endl;
      Offic_SORTIE_genes.open(output_genes_file.c_str(), ios::out|ios::trunc);
      if (!Offic_SORTIE_genes){
         cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<output_genes_file<<endl;
         if (file_log)
            Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<output_genes_file<<endl;
         exit(EXIT_FAILURE);
      }
   }

   //On compte les noeuds de duplication, de spéciation et autres stats
   int nb_GDup_arbres=0;//comptés sur arbres réconciliés
   int nb_GDup_arbres_BIS=0; //devra être égal à nb_GDup_arbres, c'est juste qu'on n'utilisera pas la même fonction pour compter
   int nb_Spe_arbres=0;
   int nb_Extant_arbres=0;
   int nb_Per_arbres=0;
   map<int,int> especes_ancestrales_nb_gene;
   map<int,int> especes_nb_pertes;
   map<int,int> especes_nb_dup;
   map<string,int> gene_species_ANC;
   int i=0;
   vector<Tree *>::iterator it;
   for(it=Arbres.begin();it!=Arbres.end();it++){ 
      Tree * Gp=*it;
      TreeTemplate<Node> *G = dynamic_cast <TreeTemplate<Node> *> (Gp);
      nb_GDup_arbres+=countNodesWithBranchProperty(G->getRootNode(),typ,dupl);
      nb_Spe_arbres+=countNodesWithBranchPropertyBySpecies(G->getRootNode(),typ,spe,especes_ancestrales_nb_gene);
      nb_Extant_arbres+=countNodesWithBranchPropertyBySpecies(G->getRootNode(),typ,ga,extant_species_gene_nb);
      nb_GDup_arbres_BIS+=countNodesWithBranchPropertyBySpecies(G->getRootNode(),typ,dupl,especes_nb_dup);
      nb_Per_arbres+=countNodesWithBranchPropertyBySpecies(G->getRootNode(),typ,per,especes_nb_pertes);
      //Création du fichier de sortie de correspondance des gènes
      if(!fic_SORTIE_GeneDone)
         ecritGenes(G->getRootNode(),i,Offic_SORTIE_genes,gene_species_ANC);
      i++;
   }
   //On libère la mémoire du vecteur d'arbres
   vector<Tree*>::iterator ita;
   for (ita=Arbres.begin();ita!=Arbres.end();ita++)
      delete(*ita);
   Arbres.clear();

   if(!fic_SORTIE_GeneDone)
      Offic_SORTIE_genes.close();


// Recovery of stats from Step3_DECO.cpp for Stats Bloc ADJACENCIES TREES
   ifstream IffilStatsDECO (DECO_stats_file.c_str(), ios::in);
   if (IffilStatsDECO){
      string buffer;
      //int num_spe,nb_gene;
      while (!IffilStatsDECO.eof()){
         IffilStatsDECO>>buffer;
         if (buffer=="nb_Crea=")
            IffilStatsDECO>>nb_Crea;
         if (buffer=="nb_GDup=")
            IffilStatsDECO>>nb_GDup;
         if (buffer=="nb_Bk=")
            IffilStatsDECO>>nb_Bk;
         if (buffer=="nb_ADup=")
            IffilStatsDECO>>nb_ADup;
      }
   }
   else
      {
    cout<<endl<<"From Step4_statistics.cpp: ERROR while opening file "<<DECO_stats_file<<endl;
    exit(EXIT_FAILURE);
      }
   IffilStatsDECO.close();

////   Récupération des espèces dont on veut supprimer les données (TEST)
//   cout<<"List of species where data have to be deleted :"<<endl;
   vector<string> del_species;
//   ifstream IffilDelData (del_spe_data.c_str(),ios::in);
//   if(IffilDelData){
//      string buffer;
//      while (!IffilDelData.eof()){
//         IffilDelData>>buffer;
//         if(!IffilDelData.eof()){
//            del_species.push_back(buffer);
//            cout<<"\t"<<buffer<<endl;
//         }
//      }   
//   }
//   else
//      cout<<"\nFrom Step4_statistics.cpp: ERROR while opening file "<<del_spe_data<<endl;   
//   IffilDelData.close();

//EXTANT
   //Méthode permettant de récupérer les gènes actuels n'ayant pas d'adjacences et espèces actuelles contenant des gènes
   map<string,string>::iterator it_act;
   vector<adjacence>::iterator it_adj;
   vector<string>::iterator it_str;
   map<string,string> extant_gene_none_adj;
   for(it_act=gene_species_EXT.begin();it_act!=gene_species_EXT.end();it_act++){
      string gene_name=(*it_act).first;
      string specie_name=(*it_act).second;
   }
   adjacencies.clear();

//   for(it_act=extant_gene_none_adj.begin();it_act!=extant_gene_none_adj.end();it_act++){
//      if((*it_act).second=="Homo_sapiens")
//         cout<<(*it_act).first<<"\t";         
//   }
   extant_gene_none_adj.clear();


   //Méthode permettant de créer la liste des espèces actuels n'ayant pas de données sur les gènes et de supprimer les données des espèces dont on ne veut pas de données
   vector<string>::iterator it_str1;
   vector<string>::iterator it_str2;
   vector<string> extant_species_none_genes;
   for(it_str1=leaves_name_S_tree.begin();it_str1!=leaves_name_S_tree.end();it_str1++){
      int specie_id=S->getNode(*it_str1)->getId();

      for(it_str2=del_species.begin();it_str2!=del_species.end();it_str2++){
         if(*it_str2==*it_str1)
            extant_species_gene_nb.erase(specie_id);
      }

      if (extant_species_gene_nb.find(specie_id)==extant_species_gene_nb.end())
         extant_species_none_genes.push_back(*it_str1);
   }
   del_species.clear();

   //Méthode permettant de créer la liste des espèces ancestrales n'ayant pas de données sur les gènes
   vector<int>::iterator it_i;
   vector<int> ances_species_none_genes;
   for(it_i=inner_node_S_tree.begin();it_i!=inner_node_S_tree.end();it_i++){
      int specie_id=*it_i;

      if (especes_ancestrales_nb_gene.find(specie_id)==especes_ancestrales_nb_gene.end())
         ances_species_none_genes.push_back(specie_id);
   }
   inner_node_S_tree.clear();

//Pour statistiques sur adjacences actuelles
   int nb_ea_tot=0;
   int nb_ea_min=100000;
   int nb_ea_max=0;
   int nb_adj_par_geneE_tot=0;
   int nb_adj_par_geneE_min=100000;
   int nb_adj_par_geneE_max=0;
   int nb_adj_par_geneE_1=0;
   int nb_adj_par_geneE_2=0;
   int nb_adj_par_geneE_sup2=0;
   int nb_genes_e=0;
   //Pour stocker le nb de gènes par espèce actuelle
   map<string,int> extant_gene_adj_nb4;
   map<string,int>::iterator it_gana;
   map<int,int> nb_gene_extant_spe;
   //Stats pour le nb d'adj actuelless par espèce + nb gènes actuels
   map<int,vector<adjacence> >::iterator it_AA;
   for(it_AA=Adj_actuelles.begin();it_AA!=Adj_actuelles.end();it_AA++){
      vector<adjacence> v=(*it_AA).second;
      vector<adjacence>::iterator it_v;
      //nb d'adj actuelles par espèce      
      int n=v.size();
      nb_ea_tot+=n;
      if(nb_ea_min>n)
         nb_ea_min=n; 
      if(nb_ea_max<n)
         nb_ea_max=n;
      //nb gènes actuelles
      set<string> genes;
      for(it_v=v.begin();it_v!=v.end();it_v++){
         genes.insert((*it_v).gene1);
         genes.insert((*it_v).gene2);
         if(extant_gene_adj_nb4.find((*it_v).gene1)!=extant_gene_adj_nb4.end())
            extant_gene_adj_nb4[(*it_v).gene1]++;
         else
            extant_gene_adj_nb4[(*it_v).gene1]=1;
         if(extant_gene_adj_nb4.find((*it_v).gene2)!=extant_gene_adj_nb4.end())
            extant_gene_adj_nb4[(*it_v).gene2]++;
         else
            extant_gene_adj_nb4[(*it_v).gene2]=1;
      }
      nb_gene_extant_spe[(*it_AA).first]=genes.size();
      nb_genes_e+=genes.size();
   } 

   for(it_gana=extant_gene_adj_nb4.begin();it_gana!=extant_gene_adj_nb4.end();it_gana++){
      int n=(*it_gana).second;   
         nb_adj_par_geneE_tot+=n;
      if(nb_adj_par_geneE_min>n)
         nb_adj_par_geneE_min=n; 
      if(nb_adj_par_geneE_max<n)
         nb_adj_par_geneE_max=n;
      if(n==1)
         nb_adj_par_geneE_1++;
      if(n==2)
         nb_adj_par_geneE_2++;
      if(n>2){
         nb_adj_par_geneE_sup2++;
//         Affiche les gènes actuels qui ont plus de 2 adjacences
//         cout<<(*it_gana).first<<endl;
      }
   }

   map<int,int> distrib_nb_adj_extant_spec;
   map<string, int>::iterator p2;
   for(p2=extant_gene_adj_nb4.begin();p2!=extant_gene_adj_nb4.end();p2++){
      if (not distrib_nb_adj_extant_spec[p2->second])
         distrib_nb_adj_extant_spec[p2->second]=1;
      else
         distrib_nb_adj_extant_spec[p2->second]+=1;
   }

//Pour statistiques sur adjacences ancestrales
   //map contenant les adjacences ancestrales par espèces
   map<int,vector<adjacence> > Adj_ancestrales;
   int nb_aa_tot=0;
   int nb_aa_min=100000;
   int nb_aa_max=0;
   int nb_adj_par_geneA_tot=0;
   int nb_adj_par_geneA_min=100000;
   int nb_adj_par_geneA_max=0;
   int nb_adj_par_geneA_1=0;
   int nb_adj_par_geneA_2=0;
   int nb_adj_par_geneA_sup2=0;
   int nb_genes_a=0;

   string adj_ances_file=prefixe+output_adj;
   cout<<"\nReading the file of ancestral adjacencies: "<<adj_ances_file<<endl<<endl;
   lectureFicRelBin(adj_ances_file, Adj_ancestrales);

   //Pour stocker le nb de gènes par espèce ancestrale
   map<int,int> nb_gene_ances_spe;
   map<string,int> ances_gene_adj_nb;
   //Stats pour le nb d'adj ancestrales par espèce + nb gènes ancestraux
   for(it_AA=Adj_ancestrales.begin();it_AA!=Adj_ancestrales.end();it_AA++){
      vector<adjacence> v=(*it_AA).second;
      vector<adjacence>::iterator it_v;
      //nb d'adj ancestrales par espèce      
      int n=v.size();
      nb_aa_tot+=n;
      if(nb_aa_min>n)
         nb_aa_min=n; 
      if(nb_aa_max<n)
         nb_aa_max=n;
      //nb gènes ancestraux
      set<string> genes;
      for(it_v=v.begin();it_v!=v.end();it_v++){
         genes.insert((*it_v).gene1);
         genes.insert((*it_v).gene2);
         if(ances_gene_adj_nb.find((*it_v).gene1)!=ances_gene_adj_nb.end())
            ances_gene_adj_nb[(*it_v).gene1]++;
         else
            ances_gene_adj_nb[(*it_v).gene1]=1;
         if(ances_gene_adj_nb.find((*it_v).gene2)!=ances_gene_adj_nb.end())
            ances_gene_adj_nb[(*it_v).gene2]++;
         else
            ances_gene_adj_nb[(*it_v).gene2]=1;
      }
      nb_gene_ances_spe[(*it_AA).first]=genes.size();
      nb_genes_a+=genes.size();
   }

   for(it_gana=ances_gene_adj_nb.begin();it_gana!=ances_gene_adj_nb.end();it_gana++){
      int n=(*it_gana).second;   
         nb_adj_par_geneA_tot+=n;
      if(nb_adj_par_geneA_min>n)
         nb_adj_par_geneA_min=n; 
      if(nb_adj_par_geneA_max<n)
         nb_adj_par_geneA_max=n;
      if(n==1)
         nb_adj_par_geneA_1++;
      if(n==2)
         nb_adj_par_geneA_2++;
      if(n>2){
         nb_adj_par_geneA_sup2++;
//         Affiche les gènes ancestraux qui ont plus de 2 adjacences
//         cout<<(*it_gana).first<<endl;
      }
   }

   map<int,int> distrib_nb_adj_ances_spec;
   for(p2=ances_gene_adj_nb.begin();p2!=ances_gene_adj_nb.end();p2++){
      if (not distrib_nb_adj_ances_spec[p2->second])
         distrib_nb_adj_ances_spec[p2->second]=1;
      else
         distrib_nb_adj_ances_spec[p2->second]+=1;
   }

//Pour statistiques sur Arbres d'adjacences
   string adj_trees_file=prefixe+output_adj_trees_per_class;
   ifstream Iffile_OUTPUT_adjacencies_trees(adj_trees_file.c_str(), ios::in);

   vector<Tree *> Adjacencies_trees;
   int nb_adj_trees_per_class_min=100000;
   int nb_adj_trees_per_class_max=0;
   int tot_nb_adj_trees_per_class=0;
   map<string,int> adj_trees_per_class;
   map<int,int> distrib_nb_adj_tree_class;
   if (Iffile_OUTPUT_adjacencies_trees){
      string tampon, adj_class, buffer;
      while (!Iffile_OUTPUT_adjacencies_trees.eof()){
         Iffile_OUTPUT_adjacencies_trees>>tampon;
         if (tampon=="class"){
            Iffile_OUTPUT_adjacencies_trees>>adj_class;
            adj_trees_per_class[adj_class]=0;
         }
         if (tampon=="Tree"){   
            adj_trees_per_class[adj_class]+=1;

            Iffile_OUTPUT_adjacencies_trees>>tampon;
            Iffile_OUTPUT_adjacencies_trees>>tampon;
            istringstream buffer(tampon);
            if(OUTPUT_FORMAT==0){
               Newick * newickReader = new Newick(true,true);
               newickReader->enableExtendedBootstrapProperty(esp);
               newickReader->read(buffer,Adjacencies_trees);
               delete(newickReader);
            }
            else if(OUTPUT_FORMAT==1){
               Nhx * nhxReader = new Nhx(true);
               nhxReader->read(buffer,Adjacencies_trees);
               delete(nhxReader);
            }
            else{
               cout<<"From Step4_statistics.cpp: ERROR wrong OUTPUT_FORMAT"<<endl;
               exit(EXIT_FAILURE);
            }
//            if(adj_class=="103|3755|64-50")
//               cout<<adj_class<<" : "<<Adjacencies_trees.back()->getNumberOfNodes()<<" "<<endl;
//            if(adj_class=="1042|2355|58-64")
//               cout<<adj_class<<" : "<<Adjacencies_trees.back()->getNumberOfNodes()<<" "<<endl;
//            if(adj_class=="1050|1050|51")
//               cout<<adj_class<<" : "<<Adjacencies_trees.back()->getNumberOfNodes()<<" "<<endl;
//            if(adj_class=="1057|3898|111-77")
//               cout<<adj_class<<" : "<<Adjacencies_trees.back()->getNumberOfNodes()<<" "<<endl;
//            if(adj_class=="1045|2043|50-111")
//               cout<<adj_class<<" : "<<Adjacencies_trees.back()->getNumberOfNodes()<<" "<<endl;
         }
      }
      Iffile_OUTPUT_adjacencies_trees.close();
   }
   else{
      cout<<endl<<"From Step4_statistics.cpp: ERROR while opening file "<<adj_trees_file<<endl;
      exit(EXIT_FAILURE);
   }
   
   map<string, int>::iterator p;
   for(p=adj_trees_per_class.begin();p!=adj_trees_per_class.end();p++){
      int class_size=(*p).second;
      tot_nb_adj_trees_per_class+=class_size;
      if(nb_adj_trees_per_class_min>class_size)
         nb_adj_trees_per_class_min=class_size; 
      if(nb_adj_trees_per_class_max<class_size)
         nb_adj_trees_per_class_max=class_size;

      if (not distrib_nb_adj_tree_class[p->second])
         distrib_nb_adj_tree_class[p->second]=1;
      else
         distrib_nb_adj_tree_class[p->second]+=1;
   }


   ///////////////////////////////////////////////////////////////////////////
   // Création du fichier de sortie OUTPUT_stats_human_readable
   string OUTPUT_human_stats=prefixe+output_human_stats;   
   ofstream Offile_OUTPUT_stats_human;
   cout<<"Creation of the OUTPUT file for the matching of genes: "<<OUTPUT_human_stats<<endl;
   Offile_OUTPUT_stats_human.open(OUTPUT_human_stats.c_str(), ios::out|ios::trunc);
   if(!Offile_OUTPUT_stats_human)
      {
    cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<OUTPUT_human_stats<<endl;
    if(file_log)
Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<OUTPUT_human_stats<<endl;
    exit(EXIT_FAILURE);
      }
   ///////////////////////////////////////////////////////////////////////////
   // Création du fichier de sortie OUTPUT_stats_machine
   string OUTPUT_machine_stats=prefixe+output_machine_stats;   
   ofstream Offile_OUTPUT_stats_machine;
   cout<<"Creation of the OUTPUT file for the matching of genes: "<<OUTPUT_machine_stats<<endl<<endl;
   Offile_OUTPUT_stats_machine.open(OUTPUT_machine_stats.c_str(), ios::out|ios::trunc);
   if(!Offile_OUTPUT_stats_machine)
      {
    cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<OUTPUT_machine_stats<<endl;
    if(file_log)
       Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<OUTPUT_machine_stats<<endl;
    exit(EXIT_FAILURE);
      }
   

   //   DÉBUT SORTIE DES DONNÉES STATISTIQUES   //
   cout<<"###################################################"<<endl;
   cout<<"######  Computation done -- STATISTIC OUTPUT  #####"<<endl;
   cout<<"###################################################"<<endl<<endl;
   
   //BLOC DE STATS SUR LES ESPÈCES
   cout<<"###################################"<<endl;
   cout<<"#####  Statistics on SPECIES  #####"<<endl;
   cout<<"###################################"<<endl<<endl;
   
   Offile_OUTPUT_stats_human<<"#####  Statistics on SPECIES  #####"<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"@SPECIES"<<endl<<endl;
   
   //Variables pour espèces actuelles
   int nb_extant_spe_tot=leaves_name_S_tree.size();
   int nb_extant_spe_genes=extant_species_gene_nb.size();
   int nb_extant_spe_none_genes=nb_extant_spe_tot-nb_extant_spe_genes;
   //Variables pour espèces ancestrales
   int nb_ances_spe_tot=node_S_tree.size()-leaves_name_S_tree.size();
   leaves_name_S_tree.clear();
   int nb_ances_spe_genes=especes_ancestrales_nb_gene.size();
   int nb_ances_spe_none_genes=nb_ances_spe_tot-nb_ances_spe_genes;
   //Variables pour espèces totales
   map<int,int>::iterator it_m;
   int tot_nb_species=node_S_tree.size();
   node_S_tree.clear();
   int tot_nb_species_genes=nb_ances_spe_genes+nb_extant_spe_genes;
   int tot_nb_species_none_genes=nb_ances_spe_none_genes+nb_extant_spe_none_genes;
   
   //OUTPUT screen du bloc ESPÈCES
   cout<<"Number of species :\t\t\t"<<tot_nb_species<<"\t(EXT "<<nb_extant_spe_tot<<" | ANC "<<nb_ances_spe_tot<<")"<<endl;
   cout<<"Number of species + Gene data :\t\t"<<tot_nb_species_genes<<"\t(EXT "<<nb_extant_spe_genes<<" | ANC "<<nb_ances_spe_genes<<")"<<endl;
   cout<<"Number of species - Gene data :\t\t"<<tot_nb_species_none_genes<<"\t(Ratio "<<float(tot_nb_species_none_genes*100.00)/float(tot_nb_species)<<"%)"<<endl;
   cout<<"Number of EXT species - Gene data :\t"<<nb_extant_spe_none_genes<<"\t(Ratio "<<float(nb_extant_spe_none_genes*100.00)/float(nb_extant_spe_tot)<<"%)"<<endl;
   cout<<"Number of ANC species - Gene data :\t"<<nb_ances_spe_none_genes<<"\t(Ratio "<<float(nb_ances_spe_none_genes*100.00)/float(nb_ances_spe_tot)<<"%)"<<endl<<endl;
   //ÉCRITURE dans fichier stats sortie HUMAN READABLE
   Offile_OUTPUT_stats_human<<"Number of species :\t\t\t"<<tot_nb_species<<"\t(EXT "<<nb_extant_spe_tot<<" | ANC "<<nb_ances_spe_tot<<")"<<endl;
   Offile_OUTPUT_stats_human<<"Number of species + Gene data :\t\t"<<tot_nb_species_genes<<"\t(EXT "<<nb_extant_spe_genes<<" | ANC "<<nb_ances_spe_genes<<")"<<endl;
   Offile_OUTPUT_stats_human<<"Number of species - Gene data :\t\t"<<tot_nb_species_none_genes<<"\t(Ratio "<<float(tot_nb_species_none_genes*100.00)/float(tot_nb_species)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of EXT species - Gene data :\t"<<nb_extant_spe_none_genes<<"\t(Ratio "<<float(nb_extant_spe_none_genes*100.00)/float(nb_extant_spe_tot)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of ANC species - Gene data :\t"<<nb_ances_spe_none_genes<<"\t(Ratio "<<float(nb_ances_spe_none_genes*100.00)/float(nb_ances_spe_tot)<<"%)"<<endl<<endl;
   
   Offile_OUTPUT_stats_machine<<"TOT\t"<<tot_nb_species<<"\t"<<nb_extant_spe_tot<<"\t"<<nb_ances_spe_tot<<endl;
   Offile_OUTPUT_stats_machine<<"+GENE\t"<<tot_nb_species_genes<<"\t"<<nb_extant_spe_genes<<"\t"<<nb_ances_spe_genes<<endl;
   Offile_OUTPUT_stats_machine<<"-GENE\t"<<tot_nb_species_none_genes<<"\t"<<nb_extant_spe_none_genes<<"\t"<<nb_ances_spe_none_genes<<endl<<endl;
//   cout<<"Species with Gene data :"<<endl;
//   cout<<"-----------------------------------"<<endl;
//   //Extant species

//   for(it_m=extant_species_gene_nb.begin();it_m!=extant_species_gene_nb.end();it_m++)
//   {
//      int id=(*it_m).first;
//      string name=S->getNode(id)->getName();
//      cout<<"EXT Specie #"<<id<<"\t("<<name<<")"<<endl;
//   }
//   //Ancestral species
//   for(it_m=especes_ancestrales_nb_gene.begin();it_m!=especes_ancestrales_nb_gene.end();it_m++)
//   {
//      cout<<"ANC Specie #"<<(*it_m).first<<endl;
//   }
//   cout<<endl;
   cout<<"Species without Gene data :"<<endl;
   cout<<"-----------------------------------"<<endl;

   Offile_OUTPUT_stats_human<<"Species without Gene data :"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   //Extant species
   for(it_str=extant_species_none_genes.begin();it_str!=extant_species_none_genes.end();it_str++){
      int id=S->getNode(*it_str)->getId();
      cout<<"EXT Specie #"<<id<<"\t("<<(*it_str)<<")"<<endl;
      Offile_OUTPUT_stats_human<<"EXT Specie #"<<id<<"\t("<<(*it_str)<<")"<<endl;
      Offile_OUTPUT_stats_machine<<"EXT\t#"<<id<<"\t"<<(*it_str)<<endl;
   }
   //Ancestral species
   for(it_i=ances_species_none_genes.begin();it_i!=ances_species_none_genes.end();it_i++){
      cout<<"ANC Specie #"<<(*it_i)<<endl;
      Offile_OUTPUT_stats_human<<"ANC Specie #"<<(*it_i)<<endl;
      Offile_OUTPUT_stats_machine<<"ANC\t#"<<(*it_i)<<endl;
   }
   cout<<endl<<endl<<endl;
   Offile_OUTPUT_stats_human<<endl<<endl<<endl;
   Offile_OUTPUT_stats_machine<<endl<<endl<<endl;

//BLOC DE STATS SUR LES GÈNES
   cout<<"#################################"<<endl;
   cout<<"#####  Statistics on GENES  #####"<<endl;
   cout<<"#################################"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"#####  Statistics on GENES  #####"<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"@GENES"<<endl<<endl;
//Variables pour gènes actuels
   int nb_extant_gene_min=100000;
   int nb_extant_gene_max=0;
   //   nb_extant_genes == nb_extant_gene
   int nb_extant_gene=0;
   //int nb_extant_genes=gene_species_EXT.size();
   int nb_extant_gene_adj=extant_gene_adj_nb4.size();

   for(it_m=extant_species_gene_nb.begin();it_m!=extant_species_gene_nb.end();it_m++){
      int id=(*it_m).first;
      int nb=extant_species_gene_nb[id];
      nb_extant_gene+=nb;
      if(nb_extant_gene_min>nb)
         nb_extant_gene_min=nb; 
      if(nb_extant_gene_max<nb)
         nb_extant_gene_max=nb; 
   }
   float aver_nb_gene_ext_spe=nb_extant_gene/float(nb_extant_spe_genes);
   int nb_extant_gene_none_adj=nb_extant_gene-nb_extant_gene_adj;
//Variables pour gènes ancestraux
   int nb_ances_gene_min=100000;
   int nb_ances_gene_max=0;
   int nb_ances_gene=0;
   int nb_ances_gene_adj=ances_gene_adj_nb.size();
   //Nombre de gènes/espèces sur l'ensemble des gènes ancestraux (map especes_ancestrales_nb_gene)
   for(it_m=especes_ancestrales_nb_gene.begin();it_m!=especes_ancestrales_nb_gene.end();it_m++){
      int id=(*it_m).first;
      int nb=especes_ancestrales_nb_gene[id];
      nb_ances_gene+=nb;
      if(nb_ances_gene_min>nb)
         nb_ances_gene_min=nb; 
      if(nb_ances_gene_max<nb)
         nb_ances_gene_max=nb;
   }
   float aver_nb_gene_ances_spe=nb_ances_gene/float(nb_ances_spe_genes);
   int nb_ances_gene_none_adj=nb_ances_gene-nb_ances_gene_adj;
//Variables pour gènes totaux
   int min_genes=100000;
   int max_genes=0;
   int tot_nb_genes=nb_extant_gene+nb_ances_gene;
   int nb_gene_adj=nb_extant_gene_adj+nb_ances_gene_adj;
   int nb_gene_none_adj=tot_nb_genes-nb_gene_adj;
   if (nb_extant_gene_min>nb_ances_gene_min)
      min_genes=nb_ances_gene_min; 
   else
      min_genes=nb_extant_gene_min;
   if (nb_extant_gene_max>nb_ances_gene_max)
      max_genes=nb_extant_gene_max;
   else
      max_genes=nb_ances_gene_max;
   float aver_nb_gene_spe=(tot_nb_genes)/float(tot_nb_species_genes);

//OUTPUT screen du bloc GÈNES
   cout<<"Number of genes : "<<tot_nb_genes<<" (EXT "<<nb_extant_gene/*<<" (extant_species_gene_nb) | "<<nb_extant_genes<<" (=gene_species_EXT.size())"*/<<" | ANC "<<nb_ances_gene<<")"<<endl<<endl;
   Offile_OUTPUT_stats_human<<"Number of genes : "<<tot_nb_genes<<" (EXT "<<nb_extant_gene/*<<" (extant_species_gene_nb) | "<<nb_extant_genes<<" (=gene_species_EXT.size())"*/<<" | ANC "<<nb_ances_gene<<")"<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"TOT\t"<<tot_nb_genes<<"\t"<<nb_extant_gene<<"\t"<<nb_ances_gene<<endl<<endl;
//   cout<<"Number of genes by specie :"<<endl;
//   cout<<"-----------------------------------"<<endl;
//   for(it_m=extant_species_gene_nb.begin();it_m!=extant_species_gene_nb.end();it_m++)
//   {
//      int id=(*it_m).first;
//      string name=S->getNode(id)->getName();
//      cout<<"EXT Specie #"<<id<<" : "<<extant_species_gene_nb[id]<<" genes\t("<<name<<")"<<endl;
//   }
//   for(it_m=especes_ancestrales_nb_gene.begin();it_m!=especes_ancestrales_nb_gene.end();it_m++)
//   {
//      int id=(*it_m).first;
//      cout<<"ANC Specie #"<<id<<" : "<<especes_ancestrales_nb_gene[id]<<" genes"<<endl;
//   }
//   cout<<endl;
   cout<<"Average number of genes per specie : "<<aver_nb_gene_spe<<" (EXT "<<aver_nb_gene_ext_spe<<" | ANC "<<aver_nb_gene_ances_spe<<")"<<endl;   
   cout<<"Minimum number : "<<min_genes<<" (EXT "<<nb_extant_gene_min<<" | ANC "<<nb_ances_gene_min<<")"<<endl;
   cout<<"Maximum number : "<<max_genes<<" (EXT "<<nb_extant_gene_max<<" | ANC "<<nb_ances_gene_max<<")"<<endl<<endl;
   cout<<"Number of genes belonging to adjacencies :\t"<<nb_gene_adj<<" ("<<nb_gene_none_adj<<" genes which have no adjacency | Ratio "<<float(nb_gene_none_adj*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   cout<<"Number of EXT genes belonging to adjacencies :\t"<<nb_extant_gene_adj<<" ("<<nb_extant_gene_none_adj<<" genes which have no adjacency | Ratio "<<float(nb_extant_gene_none_adj*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   cout<<"Number of ANC genes belonging to adjacencies :\t"<<nb_ances_gene_adj<<" ("<<nb_ances_gene_none_adj<<" genes which have no adjacency | Ratio "<<float(nb_ances_gene_none_adj*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Average number of genes per specie : "<<aver_nb_gene_spe<<" (EXT "<<aver_nb_gene_ext_spe<<" | ANC "<<aver_nb_gene_ances_spe<<")"<<endl;
   Offile_OUTPUT_stats_human<<"Minimum number : "<<min_genes<<" (EXT "<<nb_extant_gene_min<<" | ANC "<<nb_ances_gene_min<<")"<<endl;
   Offile_OUTPUT_stats_human<<"Maximum number : "<<max_genes<<" (EXT "<<nb_extant_gene_max<<" | ANC "<<nb_ances_gene_max<<")"<<endl<<endl;
   Offile_OUTPUT_stats_human<<"Number of genes belonging to adjacencies :\t"<<nb_gene_adj<<" ("<<nb_gene_none_adj<<" genes which have no adjacency | Ratio "<<float(nb_gene_none_adj*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of EXT genes belonging to adjacencies :\t"<<nb_extant_gene_adj<<" ("<<nb_extant_gene_none_adj<<" genes which have no adjacency | Ratio "<<float(nb_extant_gene_none_adj*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of ANC genes belonging to adjacencies :\t"<<nb_ances_gene_adj<<" ("<<nb_ances_gene_none_adj<<" genes which have no adjacency | Ratio "<<float(nb_ances_gene_none_adj*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl<<endl<<endl;

   Offile_OUTPUT_stats_machine<<"AV\t"<<aver_nb_gene_spe<<"\t"<<aver_nb_gene_ext_spe<<"\t"<<aver_nb_gene_ances_spe<<endl;
   Offile_OUTPUT_stats_machine<<"MIN\t"<<min_genes<<"\t"<<nb_extant_gene_min<<"\t"<<nb_ances_gene_min<<endl;
   Offile_OUTPUT_stats_machine<<"MAX\t"<<max_genes<<"\t"<<nb_extant_gene_max<<"\t"<<nb_ances_gene_max<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"+ADJ\t"<<nb_gene_adj<<"\t"<<nb_extant_gene_adj<<"\t"<<nb_ances_gene_adj<<endl;
   Offile_OUTPUT_stats_machine<<"-ADJ\t"<<nb_gene_none_adj<<"\t"<<nb_extant_gene_none_adj<<"\t"<<nb_ances_gene_none_adj<<endl<<endl<<endl<<endl;

//Liste des gènes sans adjacences (EXT/ANC) (Voir ligne 351) => Utile???


//BLOC DE STATS SUR LES ADJACENCES
   cout<<"#######################################"<<endl;
   cout<<"#####  Statistics on Adjacencies  #####"<<endl;
   cout<<"#######################################"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"#####  Statistics on Adjacencies  #####"<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"@ADJACENCIES"<<endl<<endl;

//Variables pour adjacences actuelles (Voir plus haut l.413-482)
   float aver_nb_adj_extant_spe=nb_ea_tot/float(nb_extant_spe_genes);
   float aver_nb_adj_extant_gene=nb_adj_par_geneE_tot/float(nb_extant_gene);
//Variables pour adjacences ancestrales (Voir plus haut l.482-556)
   float aver_nb_adj_ances_spe=nb_aa_tot/float(nb_ances_spe_genes);
   float aver_nb_adj_ances_gene=nb_adj_par_geneA_tot/float(nb_ances_gene);
//Variables pour adjacences totales
   int tot_nb_adj=nb_ea_tot+nb_aa_tot;
   float aver_nb_adj_spe=(tot_nb_adj)/float(tot_nb_species_genes);
   int tot_nb_adj_par_gene=nb_adj_par_geneE_tot+nb_adj_par_geneA_tot;
   float aver_nb_adj_gene=(tot_nb_adj_par_gene)/float(tot_nb_genes);

   int min_adj,max_adj;
   if (nb_ea_min>nb_aa_min)
      min_adj=nb_aa_min; 
   else
      min_adj=nb_ea_min;
   if (nb_ea_max>nb_aa_max)
      max_adj=nb_ea_max;
   else
      max_adj=nb_aa_max;
   int min_adj_par_gene,max_adj_par_gene;
   if(nb_ances_gene_none_adj>0)
      nb_adj_par_geneA_min=0;
   if(nb_extant_gene_none_adj>0)
      nb_adj_par_geneE_min=0;
   if (nb_adj_par_geneE_min>nb_adj_par_geneA_min)
      min_adj_par_gene=nb_adj_par_geneA_min; 
   else
      min_adj_par_gene=nb_adj_par_geneE_min;
   if (nb_adj_par_geneE_max>nb_adj_par_geneA_max)
      max_adj_par_gene=nb_adj_par_geneE_max;
   else
      max_adj_par_gene=nb_adj_par_geneA_max;
   int nb_adj_par_gene_1=nb_adj_par_geneA_1+nb_adj_par_geneE_1;
   int nb_adj_par_gene_2=nb_adj_par_geneA_2+nb_adj_par_geneE_2;
   int nb_adj_par_gene_sup2=nb_adj_par_geneA_sup2+nb_adj_par_geneE_sup2;

//OUTPUT screen du bloc ADJACENCES
   cout<<"Number of adjacencies : "<<tot_nb_adj<<" (EXT "<<nb_ea_tot<<" | ANC "<<nb_aa_tot<<")"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Number of adjacencies : "<<tot_nb_adj<<" (EXT "<<nb_ea_tot<<" | ANC "<<nb_aa_tot<<")"<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"TOT\t"<<tot_nb_adj<<"\t"<<nb_ea_tot<<"\t"<<nb_aa_tot<<endl<<endl;

   cout<<"Number of genes and adjacencies by specie"<<endl;
   cout<<"-----------------------------------"<<endl;

   Offile_OUTPUT_stats_human<<"Number of genes and adjacencies by specie"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;

   for(it_m=extant_species_gene_nb.begin();it_m!=extant_species_gene_nb.end();it_m++)
      try{
         int id=(*it_m).first;
         string name=S->getNode(id)->getName();
         cout<<"EXT Specie #"<<id<<" : "<<extant_species_gene_nb[id]<<" genes - "<<Adj_actuelles[id].size()<<" adjacencies\t("<<name<<")"<<endl;
         Offile_OUTPUT_stats_human<<"EXT Specie #"<<id<<" : "<<extant_species_gene_nb[id]<<" genes - "<<Adj_actuelles[id].size()<<" adjacencies\t("<<name<<")"<<endl;
         Offile_OUTPUT_stats_machine<<"EXT\t#"<<id<<"\t"<<extant_species_gene_nb[id]<<"\t"<<Adj_actuelles[id].size()<<"\t"<<name<<endl;
      }
      catch(exception e){
         //cout<<e.what();
         cout<<"EXT Specie #"<<(*it_m).first<<" : "<<extant_species_gene_nb[(*it_m).first]<<" genes - "<<Adj_actuelles[(*it_m).first].size()<<" adjacencies"<<endl;
         Offile_OUTPUT_stats_human<<"EXT Specie #"<<(*it_m).first<<" : "<<extant_species_gene_nb[(*it_m).first]<<" genes - "<<Adj_actuelles[(*it_m).first].size()<<" adjacencies"<<endl;
         Offile_OUTPUT_stats_machine<<"EXT\t#"<<(*it_m).first<<"\t"<<extant_species_gene_nb[(*it_m).first]<<"\t"<<Adj_actuelles[(*it_m).first].size()<<endl;
      }
   for(it_m=especes_ancestrales_nb_gene.begin();it_m!=especes_ancestrales_nb_gene.end();it_m++)
      try{
         int id=(*it_m).first;
         string name=S->getNode(id)->getName();
         cout<<"ANC Specie #"<<id<<" : "<<(*it_m).second<<" genes - "<<Adj_ancestrales[id].size()<<" adjacencies\t("<<name<<")"<<endl;
         Offile_OUTPUT_stats_human<<"ANC Specie #"<<id<<" : "<<(*it_m).second<<" genes - "<<Adj_ancestrales[id].size()<<" adjacencies\t("<<name<<")"<<endl;
         Offile_OUTPUT_stats_machine<<"ANC\t#"<<id<<"\t"<<(*it_m).second<<"\t"<<Adj_ancestrales[id].size()<<"\t"<<name<<endl;
      }
      catch(exception e){
         //cout<<e.what();
         cout<<"ANC Specie #"<<(*it_m).first<<" : "<<(*it_m).second<<" genes - "<<Adj_ancestrales[(*it_m).first].size()<<" adjacencies"<<endl;
         Offile_OUTPUT_stats_human<<"ANC Specie #"<<(*it_m).first<<" : "<<(*it_m).second<<" genes - "<<Adj_ancestrales[(*it_m).first].size()<<" adjacencies"<<endl;
         Offile_OUTPUT_stats_machine<<"ANC\t#"<<(*it_m).first<<"\t"<<(*it_m).second<<"\t"<<Adj_ancestrales[(*it_m).first].size()<<endl;
      }
   cout<<endl;
   cout<<"Number of adjacencies by specie"<<endl;
   cout<<"-----------------------------------"<<endl;
   cout<<"Average number of adjacencies per specie : "<<aver_nb_adj_spe<<" (EXT "<<aver_nb_adj_extant_spe<<" | ANC "<<aver_nb_adj_ances_spe<<")"<<endl;
   cout<<"Minimum number : "<<min_adj<<" (EXT "<<nb_ea_min<<" | ANC "<<nb_aa_min<<")"<<endl; 
   cout<<"Maximum number : "<<max_adj<<" (EXT "<<nb_ea_max<<" | ANC "<<nb_aa_max<<")"<<endl<<endl;


   cout<<"Number of adjacencies by gene"<<endl;
   cout<<"-----------------------------------"<<endl;
   cout<<"Average number of adjacencies per gene : "<<aver_nb_adj_gene<<" (EXT "<<aver_nb_adj_extant_gene<<" | ANC "<<aver_nb_adj_ances_gene<<")"<<endl;
   cout<<"Minimum number : "<<min_adj_par_gene<<" (EXT "<<nb_adj_par_geneE_min<<" | ANC "<<nb_adj_par_geneA_min<<")"<<endl; 
   cout<<"Maximum number : "<<max_adj_par_gene<<" (EXT "<<nb_adj_par_geneE_max<<" | ANC "<<nb_adj_par_geneA_max<<")"<<endl<<endl;

   cout<<"Number of genes which have no adjacency :\t"<<nb_gene_none_adj<<"\t(Ratio "<<float(nb_gene_none_adj*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   cout<<"Number of EXT genes which have no adjacency :\t"<<nb_extant_gene_none_adj<<"\t(Ratio "<<float(nb_extant_gene_none_adj*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   cout<<"Number of ANC genes which have no adjacency :\t"<<nb_ances_gene_none_adj<<"\t(Ratio "<<float(nb_ances_gene_none_adj*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl;

   cout<<"Number of genes with stricly 1 adjacencie :\t"<<nb_adj_par_gene_1<<"\t(Ratio "<<float(nb_adj_par_gene_1*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   cout<<"Number of EXT genes with stricly 1 adjacencie :\t"<<nb_adj_par_geneE_1<<"\t(Ratio "<<float(nb_adj_par_geneE_1*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   cout<<"Number of ANC genes with stricly 1 adjacencie :\t"<<nb_adj_par_geneA_1<<"\t(Ratio "<<float(nb_adj_par_geneA_1*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl;

   cout<<"Number of genes with stricly 2 adjacencies :\t\t"<<nb_adj_par_gene_2<<"\t(Ratio "<<float(nb_adj_par_gene_2*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   cout<<"Number of EXT genes with stricly 2 adjacencies :\t"<<nb_adj_par_geneE_2<<"\t(Ratio "<<float(nb_adj_par_geneE_2*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   cout<<"Number of ANC genes with stricly 2 adjacencies :\t"<<nb_adj_par_geneA_2<<"\t(Ratio "<<float(nb_adj_par_geneA_2*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl;

   cout<<"Number of genes with stricly more than 2 adjacencies :\t\t"<<nb_adj_par_gene_sup2<<"\t(Ratio "<<float(nb_adj_par_gene_sup2*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   cout<<"Number of EXT genes with stricly more than 2 adjacencies :\t"<<nb_adj_par_geneE_sup2<<"\t(Ratio "<<float(nb_adj_par_geneE_sup2*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   cout<<"Number of ANC genes with stricly more than 2 adjacencies :\t"<<nb_adj_par_geneA_sup2<<"\t(Ratio "<<float(nb_adj_par_geneA_sup2*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl;


   Offile_OUTPUT_stats_human<<endl;
   Offile_OUTPUT_stats_human<<"Number of adjacencies by specie"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   Offile_OUTPUT_stats_human<<"Average number of adjacencies per specie : "<<aver_nb_adj_spe<<" (EXT "<<aver_nb_adj_extant_spe<<" | ANC "<<aver_nb_adj_ances_spe<<")"<<endl;
   Offile_OUTPUT_stats_human<<"Minimum number : "<<min_adj<<" (EXT "<<nb_ea_min<<" | ANC "<<nb_aa_min<<")"<<endl; 
   Offile_OUTPUT_stats_human<<"Maximum number : "<<max_adj<<" (EXT "<<nb_ea_max<<" | ANC "<<nb_aa_max<<")"<<endl<<endl;


   Offile_OUTPUT_stats_human<<"Number of adjacencies by gene"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   Offile_OUTPUT_stats_human<<"Average number of adjacencies per gene : "<<aver_nb_adj_gene<<" (EXT "<<aver_nb_adj_extant_gene<<" | ANC "<<aver_nb_adj_ances_gene<<")"<<endl;
   Offile_OUTPUT_stats_human<<"Minimum number : "<<min_adj_par_gene<<" (EXT "<<nb_adj_par_geneE_min<<" | ANC "<<nb_adj_par_geneA_min<<")"<<endl; 
   Offile_OUTPUT_stats_human<<"Maximum number : "<<max_adj_par_gene<<" (EXT "<<nb_adj_par_geneE_max<<" | ANC "<<nb_adj_par_geneA_max<<")"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Number of genes which have no adjacency :\t"<<nb_gene_none_adj<<"\t(Ratio "<<float(nb_gene_none_adj*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of EXT genes which have no adjacency :\t"<<nb_extant_gene_none_adj<<"\t(Ratio "<<float(nb_extant_gene_none_adj*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of ANC genes which have no adjacency :\t"<<nb_ances_gene_none_adj<<"\t(Ratio "<<float(nb_ances_gene_none_adj*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Number of genes with stricly 1 adjacencie :\t"<<nb_adj_par_gene_1<<"\t(Ratio "<<float(nb_adj_par_gene_1*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of EXT genes with stricly 1 adjacencie :\t"<<nb_adj_par_geneE_1<<"\t(Ratio "<<float(nb_adj_par_geneE_1*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of ANC genes with stricly 1 adjacencie :\t"<<nb_adj_par_geneA_1<<"\t(Ratio "<<float(nb_adj_par_geneA_1*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Number of genes with stricly 2 adjacencies :\t\t"<<nb_adj_par_gene_2<<"\t(Ratio "<<float(nb_adj_par_gene_2*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of EXT genes with stricly 2 adjacencies :\t"<<nb_adj_par_geneE_2<<"\t(Ratio "<<float(nb_adj_par_geneE_2*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of ANC genes with stricly 2 adjacencies :\t"<<nb_adj_par_geneA_2<<"\t(Ratio "<<float(nb_adj_par_geneA_2*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Number of genes with stricly more than 2 adjacencies :\t\t"<<nb_adj_par_gene_sup2<<"\t(Ratio "<<float(nb_adj_par_gene_sup2*100.00)/float(tot_nb_genes)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of EXT genes with stricly more than 2 adjacencies :\t"<<nb_adj_par_geneE_sup2<<"\t(Ratio "<<float(nb_adj_par_geneE_sup2*100.00)/float(nb_extant_gene)<<"%)"<<endl;
   Offile_OUTPUT_stats_human<<"Number of ANC genes with stricly more than 2 adjacencies :\t"<<nb_adj_par_geneA_sup2<<"\t(Ratio "<<float(nb_adj_par_geneA_sup2*100.00)/float(nb_ances_gene)<<"%)"<<endl<<endl;

   Offile_OUTPUT_stats_machine<<endl;
   Offile_OUTPUT_stats_machine<<"AV\t"<<aver_nb_adj_spe<<"\t"<<aver_nb_adj_extant_spe<<"\t"<<aver_nb_adj_ances_spe<<endl;
   Offile_OUTPUT_stats_machine<<"MIN\t"<<min_adj<<"\t"<<nb_ea_min<<"\t"<<nb_aa_min<<endl; 
   Offile_OUTPUT_stats_machine<<"MAX\t"<<max_adj<<"\t"<<nb_ea_max<<"\t"<<nb_aa_max<<endl<<endl;

   Offile_OUTPUT_stats_machine<<"AV\t"<<aver_nb_adj_gene<<"\t"<<aver_nb_adj_extant_gene<<"\t"<<aver_nb_adj_ances_gene<<endl;
   Offile_OUTPUT_stats_machine<<"MIN\t"<<min_adj_par_gene<<"\t"<<nb_adj_par_geneE_min<<"\t"<<nb_adj_par_geneA_min<<endl; 
   Offile_OUTPUT_stats_machine<<"MAX\t"<<max_adj_par_gene<<"\t"<<nb_adj_par_geneE_max<<"\t"<<nb_adj_par_geneA_max<<endl<<endl;

   Offile_OUTPUT_stats_machine<<"-ADJ\t"<<nb_gene_none_adj<<"\t"<<nb_extant_gene_none_adj<<"\t"<<nb_ances_gene_none_adj<<endl;
   Offile_OUTPUT_stats_machine<<"ADJ=1\t"<<nb_adj_par_gene_1<<"\t"<<nb_adj_par_geneE_1<<"\t"<<nb_adj_par_geneA_1<<endl;
   Offile_OUTPUT_stats_machine<<"ADJ=2\t"<<nb_adj_par_gene_2<<"\t"<<nb_adj_par_geneE_2<<"\t"<<nb_adj_par_geneA_2<<endl;
   Offile_OUTPUT_stats_machine<<"ADJ>2\t"<<nb_adj_par_gene_sup2<<"\t"<<nb_adj_par_geneE_sup2<<"\t"<<nb_adj_par_geneA_sup2<<endl<<endl;

   cout<<"Distribution Number of genes sharing the same number of adjacencies"<<endl;
   cout<<"-----------------------------------"<<endl;
   Offile_OUTPUT_stats_human<<"Distribution Number of genes sharing the same number of adjacencies"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;

   cout<<"EXT 0 Adj - "<<nb_extant_gene_none_adj<<" genes"<<endl;
   Offile_OUTPUT_stats_human<<"EXT 0 Adj - "<<nb_extant_gene_none_adj<<" genes"<<endl;
   Offile_OUTPUT_stats_machine<<"EXT\t0\t"<<nb_extant_gene_none_adj<<endl;
   for(it_m=distrib_nb_adj_extant_spec.begin();it_m!=distrib_nb_adj_extant_spec.end();it_m++){
      cout<<"EXT "<<(*it_m).first<<" Adj - "<<(*it_m).second<<" genes"<<endl;
      Offile_OUTPUT_stats_human<<"EXT "<<(*it_m).first<<" Adj - "<<(*it_m).second<<" genes"<<endl;
      Offile_OUTPUT_stats_machine<<"EXT\t"<<(*it_m).first<<"\t"<<(*it_m).second<<endl;
   }
   cout<<endl<<"ANC "<<"0 Adj - "<<nb_ances_gene_none_adj<<" genes"<<endl;
   Offile_OUTPUT_stats_human<<endl<<"ANC "<<"0 Adj - "<<nb_ances_gene_none_adj<<" genes"<<endl;
   Offile_OUTPUT_stats_machine<<"ANC\t0\t"<<nb_ances_gene_none_adj<<endl;
   for(it_m=distrib_nb_adj_ances_spec.begin();it_m!=distrib_nb_adj_ances_spec.end();it_m++){
      cout<<"ANC "<<(*it_m).first<<" Adj - "<<(*it_m).second<<" genes"<<endl;
      Offile_OUTPUT_stats_human<<"ANC "<<(*it_m).first<<" Adj - "<<(*it_m).second<<" genes"<<endl;
      Offile_OUTPUT_stats_machine<<"ANC\t"<<(*it_m).first<<"\t"<<(*it_m).second<<endl;
   }
   cout<<endl;
   cout<<"Distribution Number of genes sharing the same number of adjacencies by specie"<<endl;
   cout<<"-----------------------------------"<<endl;

   Offile_OUTPUT_stats_human<<endl;
   Offile_OUTPUT_stats_human<<"Distribution Number of genes sharing the same number of adjacencies by specie"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
//Stats on gene neighborhood for extant species
   map<int,int> adj_gene_extant;
   for(it_AA=Adj_actuelles.begin();it_AA!=Adj_actuelles.end();it_AA++){
      vector<adjacence> v=(*it_AA).second;
      vector<adjacence>::iterator it_v;
      extant_gene_adj_nb4.clear();
      adj_gene_extant.clear();
      for(it_v=v.begin();it_v!=v.end();it_v++){
         if(extant_gene_adj_nb4.find((*it_v).gene1)!=extant_gene_adj_nb4.end())
            extant_gene_adj_nb4[(*it_v).gene1]++;
         else
            extant_gene_adj_nb4[(*it_v).gene1]=1;
         if(extant_gene_adj_nb4.find((*it_v).gene2)!=extant_gene_adj_nb4.end())
            extant_gene_adj_nb4[(*it_v).gene2]++;            
         else
            extant_gene_adj_nb4[(*it_v).gene2]=1;
      }
      for(it_gana=extant_gene_adj_nb4.begin();it_gana!=extant_gene_adj_nb4.end();it_gana++){
         if (not adj_gene_extant[it_gana->second])
            adj_gene_extant[it_gana->second]=1;
         else
            adj_gene_extant[it_gana->second]+=1;
      }
      string specie_name=S->getNode((*it_AA).first)->getName();
//      cout<<"EXT Specie #"<<(*it_AA).first<<" ("<<specie_name<<") :";
      cout<<"EXT Specie #"<<(*it_AA).first<<":";

//      Offile_OUTPUT_stats_human<<"EXT Specie #"<<(*it_AA).first<<" ("<<specie_name<<") :";
      Offile_OUTPUT_stats_human<<"EXT Specie #"<<(*it_AA).first<<":";
//      Offile_OUTPUT_stats_machine<<endl<<"EXT\t#"<<(*it_AA).first<<"\t"<<specie_name;
      Offile_OUTPUT_stats_machine<<endl<<"EXT\t#"<<(*it_AA).first;

      int nb_ext_gene_none_adj_spe=(extant_species_gene_nb[(*it_AA).first])-(nb_gene_extant_spe[(*it_AA).first]);
      cout<<"\t"<<nb_ext_gene_none_adj_spe;
      Offile_OUTPUT_stats_human<<"\t"<<nb_ext_gene_none_adj_spe;
      Offile_OUTPUT_stats_machine<<"\t"<<nb_ext_gene_none_adj_spe;
      i=0;   
      for(it_m=adj_gene_extant.begin();it_m!=adj_gene_extant.end();it_m++){
         i++;
         if(i==(*it_m).first){
//            cout<<"\t"<<((*it_m).second)<<" - "<<((*it_m).first);
            cout<<"\t"<<((*it_m).second);
//            Offile_OUTPUT_stats_human<<"\t"<<((*it_m).second)<<" - "<<((*it_m).first);
            Offile_OUTPUT_stats_human<<"\t"<<((*it_m).second);
//            Offile_OUTPUT_stats_machine<<"\t"<<((*it_m).second)<<"-"<<((*it_m).first);
            Offile_OUTPUT_stats_machine<<"\t"<<((*it_m).second);
         }
         else{
            while(i<(*it_m).first){
//               cout<<"\t0 - "<<i;
               cout<<"\t0";
//               Offile_OUTPUT_stats_human<<"\t0 - "<<i;
               Offile_OUTPUT_stats_human<<"\t0";
//               Offile_OUTPUT_stats_machine<<"\t0-"<<i;
               Offile_OUTPUT_stats_machine<<"\t0";
               i++;
            }
//            cout<<"\t"<<((*it_m).second)<<" - "<<((*it_m).first);
            cout<<"\t"<<((*it_m).second);
//            Offile_OUTPUT_stats_human<<"\t"<<((*it_m).second)<<" - "<<((*it_m).first);
            Offile_OUTPUT_stats_human<<"\t"<<((*it_m).second);
//            Offile_OUTPUT_stats_machine<<"\t"<<((*it_m).second)<<"-"<<((*it_m).first);
            Offile_OUTPUT_stats_machine<<"\t"<<((*it_m).second);
         }
      }
      cout<<"\t("<<specie_name<<")"<<endl;
      Offile_OUTPUT_stats_human<<"\t("<<specie_name<<")"<<endl;
   }
   Offile_OUTPUT_stats_machine<<endl;

//Stats on gene neighborhood for ancestral species
   map<int,int> adj_gene_ances;
   for(it_AA=Adj_ancestrales.begin();it_AA!=Adj_ancestrales.end();it_AA++){
      vector<adjacence> v=(*it_AA).second;
      vector<adjacence>::iterator it_v;
      ances_gene_adj_nb.clear();
      adj_gene_ances.clear();
      for(it_v=v.begin();it_v!=v.end();it_v++){
         if(ances_gene_adj_nb.find((*it_v).gene1)!=ances_gene_adj_nb.end())
            ances_gene_adj_nb[(*it_v).gene1]++;
         else
            ances_gene_adj_nb[(*it_v).gene1]=1;
         if(ances_gene_adj_nb.find((*it_v).gene2)!=ances_gene_adj_nb.end())
            ances_gene_adj_nb[(*it_v).gene2]++;            
         else
            ances_gene_adj_nb[(*it_v).gene2]=1;
      }
      for(it_gana=ances_gene_adj_nb.begin();it_gana!=ances_gene_adj_nb.end();it_gana++){
         if (not adj_gene_ances[it_gana->second])
            adj_gene_ances[it_gana->second]=1;
         else
            adj_gene_ances[it_gana->second]+=1;
      }
      cout<<"ANC Specie #"<<(*it_AA).first<<":";
      Offile_OUTPUT_stats_human<<"ANC Specie #"<<(*it_AA).first<<":";
      Offile_OUTPUT_stats_machine<<"ANC\t#"<<(*it_AA).first<<":";
      int nb_anc_gene_none_adj_spe=(especes_ancestrales_nb_gene[(*it_AA).first])-(nb_gene_ances_spe[(*it_AA).first]);
      cout<<"\t"<<nb_anc_gene_none_adj_spe;
      Offile_OUTPUT_stats_human<<"\t"<<nb_anc_gene_none_adj_spe;
      Offile_OUTPUT_stats_machine<<"\t"<<nb_anc_gene_none_adj_spe;
      i=0;
      for(it_m=adj_gene_ances.begin();it_m!=adj_gene_ances.end();it_m++){
         i++;
         if(i==(*it_m).first){
//            cout<<"\t"<<((*it_m).second)<<" - "<<((*it_m).first);
            cout<<"\t"<<((*it_m).second);
//            Offile_OUTPUT_stats_human<<"\t"<<((*it_m).second)<<" - "<<((*it_m).first);
            Offile_OUTPUT_stats_human<<"\t"<<((*it_m).second);
//            Offile_OUTPUT_stats_machine<<"\t"<<((*it_m).second)<<"-"<<((*it_m).first);
            Offile_OUTPUT_stats_machine<<"\t"<<((*it_m).second);
         }
         else{
            while(i<(*it_m).first){
//               cout<<"\t0 - "<<i;
               cout<<"\t0";
//               Offile_OUTPUT_stats_human<<"\t0 - "<<i;
               Offile_OUTPUT_stats_human<<"\t0";
//               Offile_OUTPUT_stats_machine<<"\t0-"<<i;
               Offile_OUTPUT_stats_machine<<"\t0";
               i++;
            }
//            cout<<"\t"<<((*it_m).second)<<" - "<<((*it_m).first);
            cout<<"\t"<<((*it_m).second);
//            Offile_OUTPUT_stats_human<<"\t"<<((*it_m).second)<<" - "<<((*it_m).first);
            Offile_OUTPUT_stats_human<<"\t"<<((*it_m).second);
//            Offile_OUTPUT_stats_machine<<"\t"<<((*it_m).second)<<"-"<<((*it_m).first);
            Offile_OUTPUT_stats_machine<<"\t"<<((*it_m).second);
         }
      }
      cout<<endl;
      Offile_OUTPUT_stats_human<<endl;
      Offile_OUTPUT_stats_machine<<endl;
   }
   cout<<endl<<endl<<endl;
   Offile_OUTPUT_stats_human<<endl<<endl<<endl;
      Offile_OUTPUT_stats_machine<<endl<<endl<<endl;

//BLOC DE STATS SUR LES CLASSES D'ADJACENCES
   cout<<"###############################################"<<endl;
   cout<<"#####  Statistics on Adjacencies Classes  #####"<<endl;
   cout<<"###############################################"<<endl<<endl;


   Offile_OUTPUT_stats_human<<"#####  Statistics on Adjacencies Classes  #####"<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"@ADJACENCIES_CLASSES"<<endl<<endl;
//Variables pour classes d'adjacences
   int nb_classes=classes_adjacences.size();
   int tot_nb_adj_per_class= 0;
   int nb_adj_class_min=1000000000;
   int nb_adj_class_max=0;
   map<string,vector<int> >::iterator it_cl;
   for(it_cl=classes_adjacences.begin();it_cl!=classes_adjacences.end();it_cl++){
      int class_size=(*it_cl).second.size();      
      tot_nb_adj_per_class+=class_size;
      if(nb_adj_class_min>class_size)
         nb_adj_class_min=class_size; 
      if(nb_adj_class_max<class_size)
         nb_adj_class_max=class_size;
   }

   cout<<"Number of adjacencies classes : "<<nb_classes<<endl;
   cout<<"Average adjacencies class size : "<<tot_nb_adj_per_class/float(nb_classes)<<endl;
   cout<<"Minimum number : "<<nb_adj_class_min<<endl; 
   cout<<"Maximum number : "<<nb_adj_class_max<<endl<<endl;

   cout<<"Distribution adjacencies classes size"<<endl;
   cout<<"-----------------------------------"<<endl;

   Offile_OUTPUT_stats_human<<"Number of adjacencies classes : "<<nb_classes<<endl;
   Offile_OUTPUT_stats_human<<"Average adjacencies class size : "<<tot_nb_adj_per_class/float(nb_classes)<<endl;
   Offile_OUTPUT_stats_human<<"Minimum number : "<<nb_adj_class_min<<endl; 
   Offile_OUTPUT_stats_human<<"Maximum number : "<<nb_adj_class_max<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Distribution adjacencies classes size"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   map<int,int> distrib_nb_adj_per_class;
   for(it_cl=classes_adjacences.begin();it_cl!=classes_adjacences.end();it_cl++){
      if (not distrib_nb_adj_per_class[it_cl->second.size()])
         distrib_nb_adj_per_class[it_cl->second.size()]=1;
      else
         distrib_nb_adj_per_class[it_cl->second.size()]+=1;
   }
   classes_adjacences.clear();
   for(it_m=distrib_nb_adj_per_class.begin();it_m!=distrib_nb_adj_per_class.end();it_m++){
      cout<<(*it_m).first<<" Adj - "<<(*it_m).second<<" Class(es)"<<endl;
      Offile_OUTPUT_stats_human<<(*it_m).first<<" Adj - "<<(*it_m).second<<" Class(es)"<<endl;
      Offile_OUTPUT_stats_machine<<(*it_m).first<<"\t"<<(*it_m).second<<endl;
   }
   cout<<endl<<endl<<endl;
   Offile_OUTPUT_stats_human<<endl<<endl<<endl;
   Offile_OUTPUT_stats_machine<<endl<<endl<<endl;

//BLOC DE STATS SUR LES ARBRES D'ADJACENCES
   cout<<"#############################################"<<endl;
   cout<<"#####  Statistics on Adjacencies Trees  #####"<<endl;
   cout<<"#############################################"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"#####  Statistics on Adjacencies Trees  #####"<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"@ADJACENCIES_TREES"<<endl<<endl;

//Variables pour abres d'adjacences (Voir plus haut l.558-639)
   int nb_adj_trees=Adjacencies_trees.size();
   cout<<"Number of adjacencies trees : "<<nb_adj_trees/*<<" (=Adjacencies_trees.size()) | "<<tot_nb_adj_trees_per_class<<" (tot_nb_adj_trees_per_class)"*/<<endl<<endl;
   cout<<"Adjacencies classes size (Size=number of adjacencies)"<<endl;
   cout<<"-----------------------------------"<<endl;
   cout<<"Average adjacencies classes size : "<<nb_adj_trees/float(nb_classes)<<endl;
   cout<<"Minimum number : "<<nb_adj_trees_per_class_min<<endl;
   cout<<"Maximum number : "<<nb_adj_trees_per_class_max<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Number of adjacencies trees : "<<nb_adj_trees/*<<" (=Adjacencies_trees.size()) | "<<tot_nb_adj_trees_per_class<<" (tot_nb_adj_trees_per_class)"*/<<endl<<endl;
   Offile_OUTPUT_stats_human<<"Adjacencies classes size (Size=number of adjacencies)"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   Offile_OUTPUT_stats_human<<"Average adjacencies classes size : "<<nb_adj_trees/float(nb_classes)<<endl;
   Offile_OUTPUT_stats_human<<"Minimum number : "<<nb_adj_trees_per_class_min<<endl;
   Offile_OUTPUT_stats_human<<"Maximum number : "<<nb_adj_trees_per_class_max<<endl<<endl;

   Offile_OUTPUT_stats_machine<<"TOT\t"<<nb_adj_trees/*<<" (=Adjacencies_trees.size()) | "<<tot_nb_adj_trees_per_class<<" (tot_nb_adj_trees_per_class)"*/<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"Adjacencies_classes_size"<<endl;
   Offile_OUTPUT_stats_machine<<"AV\t"<<nb_adj_trees/float(nb_classes)<<endl;
   Offile_OUTPUT_stats_machine<<"MIN\t"<<nb_adj_trees_per_class_min<<endl;
   Offile_OUTPUT_stats_machine<<"MAX\t"<<nb_adj_trees_per_class_max<<endl<<endl;

   cout<<"Distribution of adjacencies classes size"<<endl;
   cout<<"-----------------------------------"<<endl;
   Offile_OUTPUT_stats_human<<"Distribution of adjacencies classes size"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   for(it_m=distrib_nb_adj_tree_class.begin();it_m!=distrib_nb_adj_tree_class.end();it_m++){
      cout<<(*it_m).first<< " Adj tree(s) - "<<(*it_m).second<<" Adj classe(s)"<<endl;
      Offile_OUTPUT_stats_human<<(*it_m).first<< " Adj tree(s) - "<<(*it_m).second<<" Adj classe(s)"<<endl;
      Offile_OUTPUT_stats_machine<<(*it_m).first<< "\t"<<(*it_m).second<<endl;
   }
   cout<<endl;
   Offile_OUTPUT_stats_human<<endl;
   Offile_OUTPUT_stats_machine<<endl;

   cout<<"Adjacencies trees size (Size=number of adjacencies)"<<endl;
   cout<<"-----------------------------------"<<endl;
   Offile_OUTPUT_stats_human<<"Adjacencies trees size (Size=number of adjacencies)"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   vector<Tree *>::iterator it_ATree;
   int adj_tree_size_min=100000;
   int adj_tree_size_max=0;
   int tot_nb_adj_adj_tree=0;
   map<int,int> distrib_nb_adj_per_adj_tree;
   int nb_adj_tree_leaves_min=100000;
   int nb_adj_tree_leaves_max=0;
   int tot_nb_adj_tree_leaves=0;
   map<int,int> distrib_nb_leaves_per_adj_tree;
   for(it_ATree=Adjacencies_trees.begin();it_ATree!=Adjacencies_trees.end();it_ATree++){
      int nb_adj_tree_leaves=(*it_ATree)->getNumberOfLeaves();
      int adj_tree_size=(*it_ATree)->getNumberOfNodes();
      Tree * Gp=*it_ATree;
      TreeTemplate<Node> *G = dynamic_cast <TreeTemplate<Node> *> (Gp);

      //Pour compenser les erreurs de calcul de la library Bpp lorsqu'il y a un seul noeud (1 adj seule) => Bpp calcul l'arbre comme un arbre de taille 2
      //Lorsqu'il y a un GDup à la racine, ArbreAdj : (X-Y,(.....)[&&NHX:D=T:Ev=GDup:S=n°S:ND=0]; => Mauvais calcul de la library Bio++ qui compte un noeud supplémentaire
      if(adj_tree_size==2){
         nb_adj_tree_leaves=1;
         adj_tree_size=1;
      }
      if (E(G->getRootNode())=="GDup"){
         nb_adj_tree_leaves--;
         adj_tree_size--;
      }

      //Stats on adjacencies tree size   
      tot_nb_adj_adj_tree+=adj_tree_size;
      if(adj_tree_size_min>adj_tree_size)
         adj_tree_size_min=adj_tree_size; 
      if(adj_tree_size_max<adj_tree_size)
         adj_tree_size_max=adj_tree_size;
      if (not distrib_nb_adj_per_adj_tree[adj_tree_size])
         distrib_nb_adj_per_adj_tree[adj_tree_size]=1;
      else
         distrib_nb_adj_per_adj_tree[adj_tree_size]+=1;

      //Stats on number of adjacencies tree leaves
      tot_nb_adj_tree_leaves+=nb_adj_tree_leaves;
      if(nb_adj_tree_leaves_min>nb_adj_tree_leaves)
         nb_adj_tree_leaves_min=nb_adj_tree_leaves; 
      if(nb_adj_tree_leaves_max<nb_adj_tree_leaves)
         nb_adj_tree_leaves_max=nb_adj_tree_leaves;
      if (not distrib_nb_leaves_per_adj_tree[nb_adj_tree_leaves])
         distrib_nb_leaves_per_adj_tree[nb_adj_tree_leaves]=1;
      else
         distrib_nb_leaves_per_adj_tree[nb_adj_tree_leaves]+=1;
   }
   cout<<"Average adjacencies trees size : "<<tot_nb_adj_adj_tree/float(nb_adj_trees)<<endl;
   cout<<"Minimum number : "<<adj_tree_size_min<<endl;
   cout<<"Maximum number : "<<adj_tree_size_max<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Average adjacencies trees size : "<<tot_nb_adj_adj_tree/float(nb_adj_trees)<<endl;
   Offile_OUTPUT_stats_human<<"Minimum number : "<<adj_tree_size_min<<endl;
   Offile_OUTPUT_stats_human<<"Maximum number : "<<adj_tree_size_max<<endl<<endl;

   Offile_OUTPUT_stats_machine<<"Adjacencies_trees_size"<<endl;
   Offile_OUTPUT_stats_machine<<"AV\t"<<tot_nb_adj_adj_tree/float(nb_adj_trees)<<endl;
   Offile_OUTPUT_stats_machine<<"MIN\t"<<adj_tree_size_min<<endl;
   Offile_OUTPUT_stats_machine<<"MAX\t"<<adj_tree_size_max<<endl<<endl;

   cout<<"Distribution of adjacencies trees size"<<endl;
   cout<<"-----------------------------------"<<endl;

   Offile_OUTPUT_stats_human<<"Distribution of adjacencies trees size"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   map<int,int> distrib_inter_adj_per_adj_tree;
   int j=0;
   for(it_m=distrib_nb_adj_per_adj_tree.begin();it_m!=distrib_nb_adj_per_adj_tree.end();it_m++){
      int nb_adj_per_adj_tree=(*it_m).first;
      if(nb_adj_per_adj_tree<=5)
         distrib_inter_adj_per_adj_tree[nb_adj_per_adj_tree]=(*it_m).second;
      if(nb_adj_per_adj_tree>5 && nb_adj_per_adj_tree<=10)
         distrib_inter_adj_per_adj_tree[7]+=(*it_m).second;

      if(nb_adj_per_adj_tree>10){
         int begin=(j*10+1);
         int end=(j+1)*10;
         if (nb_adj_per_adj_tree>end){
            j+=1;
            begin=(j*10+1);
            end=(j+1)*10;
            distrib_inter_adj_per_adj_tree[begin]+=(*it_m).second;
         }
         else
            distrib_inter_adj_per_adj_tree[begin]+=(*it_m).second;
      }
   }
   j=0;

//   Si on veut afficher sous forme d'intervalle:
//   for(it_m=distrib_inter_adj_per_adj_tree.begin();it_m!=distrib_inter_adj_per_adj_tree.end();it_m++)
//   {
//      int nb_adj=(*it_m).first;
//      j+=(*it_m).second;
//      if(nb_adj<=5)
//      {
//         cout<<(*it_m).first<<" Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_human<<(*it_m).first<<" Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_machine<<(*it_m).first<<"\t"<<(*it_m).second<<endl;
//      }
//      if(nb_adj==7)
//      {
//         cout<<"[6-10] Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         cout<<"[1-10] Adj - "<<j<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_human<<"[6-10] Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_human<<"[1-10] Adj - "<<j<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_machine<<"6-10\t"<<(*it_m).second<<endl;
//         Offile_OUTPUT_stats_machine<<"1-10\t"<<j<<endl;
//      }
//      if(nb_adj>10)
//      {
//         cout<<"["<<nb_adj<<"-"<<nb_adj+9<<"] Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_human<<"["<<nb_adj<<"-"<<nb_adj+9<<"] Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_machine<<nb_adj<<"-"<<nb_adj+9<<"\t"<<(*it_m).second<<endl;
//      }         
//   }

   for(it_m=distrib_nb_adj_per_adj_tree.begin();it_m!=distrib_nb_adj_per_adj_tree.end();it_m++){
      cout<<(*it_m).first<<" Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
      Offile_OUTPUT_stats_human<<(*it_m).first<<" Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
      Offile_OUTPUT_stats_machine<<(*it_m).first<<"\t"<<(*it_m).second<<endl;
   }
   cout<<endl;
   Offile_OUTPUT_stats_human<<endl;
   Offile_OUTPUT_stats_machine<<endl;

   cout<<"Number of extant adjacencies per adjacencies tree"<<endl;
   cout<<"-----------------------------------"<<endl;
   Offile_OUTPUT_stats_human<<"Number of extant adjacencies per adjacencies tree"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;

   cout<<"Average number of adjacencies trees leaves : "<<tot_nb_adj_tree_leaves/float(nb_adj_trees)<<endl;
   cout<<"Minimum number : "<<nb_adj_tree_leaves_min<<endl;
   cout<<"Maximum number : "<<nb_adj_tree_leaves_max<<endl<<endl;

   Offile_OUTPUT_stats_human<<"Average number of adjacencies trees leaves : "<<tot_nb_adj_tree_leaves/float(nb_adj_trees)<<endl;
   Offile_OUTPUT_stats_human<<"Minimum number : "<<nb_adj_tree_leaves_min<<endl;
   Offile_OUTPUT_stats_human<<"Maximum number : "<<nb_adj_tree_leaves_max<<endl<<endl;

   Offile_OUTPUT_stats_machine<<"Extant_adjacencies_trees_size"<<endl;
   Offile_OUTPUT_stats_machine<<"AV\t"<<tot_nb_adj_tree_leaves/float(nb_adj_trees)<<endl;
   Offile_OUTPUT_stats_machine<<"MIN\t"<<nb_adj_tree_leaves_min<<endl;
   Offile_OUTPUT_stats_machine<<"MAX\t"<<nb_adj_tree_leaves_max<<endl<<endl;

   cout<<"Distribution of number of extant adjacencies per adjacencies trees"<<endl;
   cout<<"-----------------------------------"<<endl;
   Offile_OUTPUT_stats_human<<"Distribution of number of extant adjacencies per adjacencies trees"<<endl;
   Offile_OUTPUT_stats_human<<"-----------------------------------"<<endl;
   map<int,int> distrib_inter_ext_adj_per_adj_tree;
   j=0;
   for(it_m=distrib_nb_leaves_per_adj_tree.begin();it_m!=distrib_nb_leaves_per_adj_tree.end();it_m++){
      int nb_ext_adj_per_adj_tree=(*it_m).first;
      if(nb_ext_adj_per_adj_tree<=5)
         distrib_inter_ext_adj_per_adj_tree[nb_ext_adj_per_adj_tree]=(*it_m).second;
      if(nb_ext_adj_per_adj_tree>5 && nb_ext_adj_per_adj_tree<=10)
         distrib_inter_ext_adj_per_adj_tree[7]+=(*it_m).second;

      if(nb_ext_adj_per_adj_tree>10){
         int begin=(j*10+1);
         int end=(j+1)*10;
         if (nb_ext_adj_per_adj_tree>end){
            j+=1;
            begin=(j*10+1);
            end=(j+1)*10;
            distrib_inter_ext_adj_per_adj_tree[begin]+=(*it_m).second;
         }
         else
            distrib_inter_ext_adj_per_adj_tree[begin]+=(*it_m).second;
      }
   }
   j=0;

//   Si on veut afficher sous forme d'intervalle:
//   for(it_m=distrib_inter_ext_adj_per_adj_tree.begin();it_m!=distrib_inter_ext_adj_per_adj_tree.end();it_m++)
//   {
//      int nb_adj=(*it_m).first;
//      j+=(*it_m).second;
//      if(nb_adj<=5)
//      {
//         cout<<(*it_m).first<<" EXT Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_human<<(*it_m).first<<" EXT Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_machine<<(*it_m).first<<"\t"<<(*it_m).second<<endl;
//      }
//      if(nb_adj==7)
//      {
//         cout<<"[6-10] EXT Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         cout<<"[1-10] EXT Adj - "<<j<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_human<<"[6-10] EXT Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_human<<"[1-10] EXT Adj - "<<j<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_machine<<"6-10\t"<<(*it_m).second<<endl;
//         Offile_OUTPUT_stats_machine<<"1-10\t"<<j<<endl;
//      }
//      if(nb_adj>10)
//      {
//         cout<<"["<<nb_adj<<"-"<<nb_adj+9<<"] EXT Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_human<<"["<<nb_adj<<"-"<<nb_adj+9<<"] EXT Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
//         Offile_OUTPUT_stats_machine<<nb_adj<<"-"<<nb_adj+9<<"\t"<<(*it_m).second<<endl;   
//      }   
//   }

   for(it_m=distrib_nb_leaves_per_adj_tree.begin();it_m!=distrib_nb_leaves_per_adj_tree.end();it_m++){
      cout<<(*it_m).first<<" EXT Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
      Offile_OUTPUT_stats_human<<(*it_m).first<<" EXT Adj - "<<(*it_m).second<<" Adj tree(s)"<<endl;
      Offile_OUTPUT_stats_machine<<(*it_m).first<<"\t"<<(*it_m).second<<endl;
   }
   cout<<endl<<endl<<endl;
   Offile_OUTPUT_stats_human<<endl<<endl<<endl;
   Offile_OUTPUT_stats_machine<<endl<<endl<<endl;

   cout<<"#######################################"<<endl;
   cout<<"#####  Statistics on Events node  #####"<<endl;
   cout<<"#######################################"<<endl<<endl;

   Offile_OUTPUT_stats_human<<"#####  Statistics on Events node  #####"<<endl<<endl;
   Offile_OUTPUT_stats_machine<<"@EVOLUTION_EVENTS"<<endl<<endl;

   cout<<"Creation number : "<<nb_Crea<<" ; (nb_arbres_adj_tot - nb_classes= "<<nb_adj_trees - nb_classes<<")"<<endl;
   cout<<"Break number : "<<nb_Bk<<endl;
   cout<<"Gene Duplication number : "<<nb_GDup<<" (compté par DeCo) ; "<<nb_GDup_arbres<<" (ds les arbres réconciliés) "<<nb_GDup_arbres_BIS<<" (pour les compter par espèces)"<<endl;
   cout<<"Adjacencie Duplication number : "<<nb_ADup<<endl;

   Offile_OUTPUT_stats_human<<"Creation number : "<<nb_Crea<<" ; (nb_arbres_adj_tot - nb_classes= "<<nb_adj_trees - nb_classes<<")"<<endl;
   Offile_OUTPUT_stats_human<<"Break number : "<<nb_Bk<<endl;
   Offile_OUTPUT_stats_human<<"Gene Duplication number : "<<nb_GDup<<" (compté par DeCo) ; "<<nb_GDup_arbres<<" (ds les arbres réconciliés) "<<nb_GDup_arbres_BIS<<" (pour les compter par espèces)"<<endl;
   Offile_OUTPUT_stats_human<<"Adjacencie Duplication number : "<<nb_ADup<<endl;

   Offile_OUTPUT_stats_machine<<"CREA\t"<<nb_Crea<<" ; (nb_arbres_adj_tot - nb_classes= "<<nb_adj_trees - nb_classes<<")"<<endl;
   Offile_OUTPUT_stats_machine<<"BREAK\t"<<nb_Bk<<endl;
   Offile_OUTPUT_stats_machine<<"GDUP\t"<<nb_GDup<<" (compté par DeCo) ; "<<nb_GDup_arbres<<" (ds les arbres réconciliés) "<<nb_GDup_arbres_BIS<<" (pour les compter par espèces)"<<endl;
   Offile_OUTPUT_stats_machine<<"ADUP\t"<<nb_ADup<<endl;





//////////////////////////////
/////   COMPO CONNEXES   /////
//////////////////////////////

////   Pour les pertes ne pas prendre en compte les espèces qui n'ont pas de données sur les gènes!!!!!
//   cout<<"\nEn ce qui concerne les pertes ("<<nb_Per_arbres<<" au total) :"<<endl;
//   cout<<"-------------------------------------------"<<endl;
//   for(it_m=especes_nb_pertes.begin();it_m!=especes_nb_pertes.end();it_m++)
//      try{
//         string nom=S->getNode((*it_m).first)->getName();
//         cout<<nom<<" : "<<especes_nb_pertes[(*it_m).first]<<" pertes dans arbres réconciliés."<<endl;
//      }
//      catch(exception e){
//         //cout<<e.what();
//         cout<<"Espèce n°"<<(*it_m).first<<" : "<<especes_nb_pertes[(*it_m).first]<<" pertes dans arbres réconciliés."<<endl;
//      }
//   especes_nb_pertes.clear();

//   map<string, int>::iterator it_msi;
//   for (it_msi=ances_gene_adj_nb.begin();it_msi!=ances_gene_adj_nb.end();it_msi++){
//      cout<<(*it_msi).first<< " : ";
//      cout<<(*it_msi).second<<" adjacences"<<endl;
//   }

//   map<int,vector<adjacence> > Dup_ancestrales;
//   string fic_SORTIE_dup=prefixe+output_dup_gene_pairs;
//   cout<<"Lecture du fichier de SORTIE des gènes dupliqués ensemble "<<fic_SORTIE_dup<<endl;
//   lectureFicRelBin(fic_SORTIE_dup, Dup_ancestrales);

////   //Normalisation par le fils (nb_gene) et autres ... :
//   cout<<"Ordre des colonnes : ";
//   cout<<"1) n° espèce ";
//   cout<<"2) nom de l'espèce ";
//   cout<<"3) Nb_dup ";
//   cout<<"4) Nb_gene ";
//   cout<<"5) Nb_dup_adj ";
//   cout<<"6) Nb Gains ";
//   TreeTemplate<Node> * S6 = new TreeTemplate<Node>(*S);
//   cout<<"7) Nb Cassure ";
//   TreeTemplate<Node> * S7 = new TreeTemplate<Node>(*S);
//   cout<<"8) Nb_dup_adj/Nb_gene ";
//   TreeTemplate<Node> * S8 = new TreeTemplate<Node>(*S);
//   cout<<"9) Nb_dup_adj/Nb_dup ";
//   TreeTemplate<Node> * S9 = new TreeTemplate<Node>(*S);
//   cout<<"10) (Nb_dup/(Nb_dup - Nb_dup_adj)) - 1";
//   TreeTemplate<Node> * S10 = new TreeTemplate<Node>(*S);
//   cout<<endl;
//   int nb_nodes=S->getNumberOfNodes();
//   for(i=0;i<nb_nodes;i++){
//      //1)
//      cout<<i<<" ";
//      Node * n=S->getNode(i);
//      //2)
//      try{
//         cout<<n->getName()<<" ";
//      }
//      catch(exception e){
//         cout<<"noname"<<" ";
//      }
//      float perDup8;
//      float perDup9;
//      float perDup10;
//      //3)
//      cout<<especes_nb_dup[i]<<" ";
//      if (n->isLeaf()){ //espèce actuelle
//         //4)
//         cout<<extant_species_gene_nb[i]<<" ";
//         if (Dup_ancestrales.find(i)!=Dup_ancestrales.end()){
//            //5)
//            cout<<Dup_ancestrales[i].size()<<" ";

//            perDup8=1.0*Dup_ancestrales[i].size()/extant_species_gene_nb[i];
//            perDup9=1.0*Dup_ancestrales[i].size()/especes_nb_dup[i];
//            perDup10=(1.0*especes_nb_dup[i]/(especes_nb_dup[i]-Dup_ancestrales[i].size()))-1.0;
//         }
//         else{
//            //5)
//            cout<<"0 ";

//            perDup8=0;
//            perDup9=0;
//            perDup10=0;
//         }
//      }
//      else{ //espèce ancestrale
//         //4)
//         cout<<especes_ancestrales_nb_gene[i]<<" ";
//         if (Dup_ancestrales.find(i)!=Dup_ancestrales.end()){
//            //5)
//            cout<<Dup_ancestrales[i].size()<<" ";

//            perDup8=1.0*Dup_ancestrales[i].size()/especes_ancestrales_nb_gene[i];
//            perDup9=1.0*Dup_ancestrales[i].size()/especes_nb_dup[i];
//            perDup10=(1.0*especes_nb_dup[i]/(especes_nb_dup[i]-Dup_ancestrales[i].size()))-1.0;
//         }
//         else{
//            //5)
//            cout<<"0 ";

//            perDup8=0;
//            perDup9=0;
//            perDup10=0;
//         }
//      }
//      //6)
//      cout<<espece_nb_gain[i]<<" ";
//      S6->getNode(i)->setDistanceToFather(espece_nb_gain[i]);
//      //7)
//      cout<<espece_nb_break[i]<<" ";
//      S7->getNode(i)->setDistanceToFather(espece_nb_break[i]);
//      //8)
//      cout<<perDup8<<" ";
//      S8->getNode(i)->setDistanceToFather(perDup8);
//      //9)
//      cout<<perDup9<<" ";
//      S9->getNode(i)->setDistanceToFather(perDup9);
//      //10)
//      cout<<perDup10<<endl;
//      S10->getNode(i)->setDistanceToFather(perDup10);
//   }
//   especes_nb_dup.clear();
//   especes_ancestrales_nb_gene.clear();
//   extant_species_gene_nb.clear();


////   cout<<"\nÉcriture des fichiers d'arbre d'espèce avec nb dupli ou autres infos sur branches"<<endl;
////   //NEWICK
//   string fic_especes6=fic_especes+"_NbGain";
//   string fic_especes7=fic_especes+"_NbBreak";
//   string fic_especes8=fic_especes+"_NbDupAdjSurNbGene";
//   string fic_especes9=fic_especes+"_NbDupAdjSurNbDupGene";
//   string fic_especes10=fic_especes+"_NbDupGeneSurNbDupGeneMoinsNbDupAdjLeToutMoinsUn";
//   Newick * newickWriter = new Newick(); 
//   newickWriter->enableExtendedBootstrapProperty(esp);
//   newickWriter->write(*S6,fic_especes6.c_str());
//   newickWriter->write(*S7,fic_especes7.c_str());
//   newickWriter->write(*S8,fic_especes8.c_str());
//   newickWriter->write(*S9,fic_especes9.c_str());
//   newickWriter->write(*S10,fic_especes10.c_str());
//   delete(newickWriter);

////   //NHX   //nhxWriter->registerProperty(Nhx::Property("Nb_Adj_Dup", NAD, true, 2)); 
//   string fic_especes6nhx=fic_especes6+".nhx";
//   string fic_especes7nhx=fic_especes7+".nhx";
//   string fic_especes8nhx=fic_especes8+".nhx";
//   string fic_especes9nhx=fic_especes9+".nhx";
//   string fic_especes10nhx=fic_especes10+".nhx";
//   Nhx * nhxWriter = new Nhx(true);
//   nhxWriter->write(*S6,fic_especes6nhx.c_str());
//   nhxWriter->write(*S7,fic_especes7nhx.c_str());
//   nhxWriter->write(*S8,fic_especes8nhx.c_str());
//   nhxWriter->write(*S9,fic_especes9nhx.c_str());
//   nhxWriter->write(*S10,fic_especes10nhx.c_str());
//   delete(nhxWriter);

//   convert<int> c;
//   map<int,int> distrib_CC_Dup_anc;
//   map<int,vector<adjacence> > cc_list_adj_Dup_anc;
//   for(it_AA=Dup_ancestrales.begin();it_AA!=Dup_ancestrales.end();it_AA++){
//      vector<adjacence> v=(*it_AA).second;
//      string nom;
//      try{
//         nom=S->getNode((*it_AA).first)->getName();
//      }
//      catch(exception e){
//         //cout<<e.what();
//         nom="Espece "+c.to_s((*it_AA).first);
//      }
//      cout<<endl<<nom<<" : ";
//      (*it_AA).second = composantesConnexes(v,distrib_CC_Dup_anc,cc_list_adj_Dup_anc);
//   }
//   distrib_CC_Dup_anc.clear();
//   cc_list_adj_Dup_anc.clear();
//   cout<<endl;
//   cout<<"Affichage distribution taille des composantes connexes :"<<endl;
//   afficheMap(distrib_CC_Dup_anc);



////////////////////////////////////////////////////////////////////
//////////   CONNECTED COMPONENTS FOR ART-DECO and DECO   //////////
////////////////////////////////////////////////////////////////////

   cout<<endl<<endl<<"##################################################"<<endl;
   cout<<"### START construction of Connected Components ###"<<endl;
   cout<<"##################################################"<<endl<<endl;
   cout<<"Get list of ancestral genes without adjacency... "<<flush;
// Remove genes present in ancestral adjacencies from gene_species_ANC to keep only alone ANC gene
   for(it_AA=Adj_ancestrales.begin();it_AA!=Adj_ancestrales.end();it_AA++){
      vector<adjacence> v=(*it_AA).second;
      for(it_adj=v.begin();it_adj!=v.end();it_adj++){
         if (gene_species_ANC.find((*it_adj).gene1)!=gene_species_ANC.end()){
            gene_species_ANC.erase((*it_adj).gene1);
         }
         if (gene_species_ANC.find((*it_adj).gene2)!=gene_species_ANC.end()){
            gene_species_ANC.erase((*it_adj).gene2); 
         }
      }
   }
   cout<<"DONE"<<endl;

//  Compo connexe pour les adjacences ancestrales inferrées par DeCo/ARt-DeCo   
   map<int,int> distrib_CC_adj_ANC;
   map<int,vector<string> > cc_list_adj_ANC;
   string output_CC_Adj_ANC="";
   if (is_readable(output_new_adj_file))
      output_CC_Adj_ANC=prefixe+"_Compo_connexe/CC_Adj_ANC_ARt-DeCo/";
   else
      output_CC_Adj_ANC=prefixe+"_Compo_connexe/CC_Adj_ANC_DeCo/";
   map<string,int>::iterator it_strint;
   cout<<"Creation of connected components files for each species in directory "<<output_CC_Adj_ANC<<" for ancestral adjacencies infered..."<<flush;
   for(it_AA=Adj_ancestrales.begin();it_AA!=Adj_ancestrales.end();it_AA++){
      vector<adjacence> v=(*it_AA).second;
      string ID_species=convertInt((*it_AA).first);
      string species_dir_CC="";
      if (is_readable(output_new_adj_file))
         species_dir_CC=output_CC_Adj_ANC+"CC_Adj_ANC_ARt-DeCo_"+ID_species;
      else
         species_dir_CC=output_CC_Adj_ANC+"CC_Adj_ANC_DeCo_"+ID_species;
      mkdir_rec(species_dir_CC.c_str());
      i=1;
//      Add alone genes in Connected Components files for species (name_species)
      for(it_strint=gene_species_ANC.begin();it_strint!=gene_species_ANC.end();it_strint++){
         if(convertInt((*it_strint).second)==ID_species){
		      ofstream Offile_OUTPUT_CC;
		      string CC_file="";
		      if (is_readable(output_new_adj_file))
		         CC_file=species_dir_CC+"/CC_Adj_ANC_ARt-DeCo_"+ID_species+"_"+convertInt(i);
		      else
		         CC_file=species_dir_CC+"/CC_Adj_ANC_DeCo_"+ID_species+"_"+convertInt(i);
		      Offile_OUTPUT_CC.open(CC_file.c_str(), ios::out|ios::trunc);
		      if(!Offile_OUTPUT_CC){
		         cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
		         if(file_log)
		            Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
		         exit(EXIT_FAILURE);
		      }
	         Offile_OUTPUT_CC<<(*it_strint).first<<endl;
            i++;
         }
      }
      cout<<endl<<"For ANC species "<<ID_species<<":"<<endl;
//      Creation of connected components for gene with adjacencies
      (*it_AA).second = composantesConnexes(v,distrib_CC_adj_ANC,cc_list_adj_ANC,i);
      map<int,vector<string> >::iterator it_aa;
//      int size_CC=cc_list_adj_ANC.size();
//      cout<<size_CC<<endl;
      for(it_aa=cc_list_adj_ANC.begin();it_aa!=cc_list_adj_ANC.end();it_aa++){
         ofstream Offile_OUTPUT_CC;
         string CC_file="";
         if (is_readable(output_new_adj_file))
            CC_file=species_dir_CC+"/CC_Adj_ANC_ARt-DeCo_"+ID_species+"_"+convertInt((*it_aa).first);
         else
            CC_file=species_dir_CC+"/CC_Adj_ANC_DeCo_"+ID_species+"_"+convertInt((*it_aa).first);
         Offile_OUTPUT_CC.open(CC_file.c_str(), ios::out|ios::trunc);
         if(!Offile_OUTPUT_CC){
            cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
            if(file_log)
               Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
            exit(EXIT_FAILURE);
         }
         vector<string>::iterator it_vs;
         for(it_vs=((*it_aa).second).begin();it_vs!=((*it_aa).second).end();it_vs++){
            Offile_OUTPUT_CC<<(*it_vs)<<endl;
         }
      }
      cc_list_adj_ANC.clear();
      distrib_CC_adj_ANC.clear();
   }
   gene_species_ANC.clear();
   Adj_ancestrales.clear();
   cout<<" DONE"<<endl<<endl;



//   Creation map Adj_actuelles_ORI to store extant adjacencies from INPUT adjacencies file
   map<int,vector<adjacence> > Adj_actuelles_ORI;
   string gene1="";
   string gene2="";
   adjacence adj;      
   ifstream IfficAdj (fic_adjacence.c_str(), ios::in);
   if (!IfficAdj){
      cout<<endl<<"From Step4_statistics.cpp: ERROR while opening file "<<fic_adjacence<<endl;
      exit(EXIT_FAILURE);
   }
   cout<<"Creation of map Adj_actuelles_ORI to associate species_ID to vector of adajacencies..."<<flush;
   while (!IfficAdj.eof()){
      //on lit la ligne en entier;
      IfficAdj>>gene1;      
      //Pour éviter un dernier tour de trop à cause d'un blanc en fin de fichier
      if(!IfficAdj.eof()){
         IfficAdj>>gene2;      
//         cout<<"On traite l'adjacence "<< gene1<<" - "<<gene2;
         //on crée l'adjacence si les deux gènes sont dans gene_species_EXT (c-à-d dans nos arbres)
         if (gene_species_EXT.find(gene1)!=gene_species_EXT.end() && gene_species_EXT.find(gene2)!=gene_species_EXT.end()){
//            cout<<", on la nomme "<<j<<endl;
            adj.gene1=gene1;
            adj.gene2=gene2;
            //on la stocke
            adjacencies.push_back(adj);
            //Stockage de l'adj actuelle dans son espèce dans Adj_actuelles_ORI
            int esp_actuelle=S->getNode(gene_species_EXT[gene1])->getId();
            if(Adj_actuelles_ORI.find(esp_actuelle)!=Adj_actuelles_ORI.end())
               Adj_actuelles_ORI[esp_actuelle].push_back(adj);
            else{
               vector<adjacence> v;
               v.push_back(adj);
               Adj_actuelles_ORI[esp_actuelle]=v;
            }
         }
      }
      j++;
   }
   IfficAdj.close();
   adjacencies.clear();
   cout<<" DONE"<<endl<<endl;
//   afficheMap(Adj_actuelles_ORI);

//      Remove genes present in ORI adjacencies from gene_species_EXT map to keep only alone gene
      for(it_AA=Adj_actuelles_ORI.begin();it_AA!=Adj_actuelles_ORI.end();it_AA++){
         vector<adjacence> v=(*it_AA).second;
         for(it_adj=v.begin();it_adj!=v.end();it_adj++){
            if (gene_species_EXT.find((*it_adj).gene1)!=gene_species_EXT.end()){
               gene_species_EXT.erase((*it_adj).gene1);
            }
            if (gene_species_EXT.find((*it_adj).gene2)!=gene_species_EXT.end()){
               gene_species_EXT.erase((*it_adj).gene2); 
            }
         }
      }

//  Compo connexe pour les adjacences actuelles présentes dans les fichiers d'adjacences en entrée de DeCo/ARt-DeCo
   map<int,int> distrib_CC_adj_Ext_ORI;
   map<int,vector<string> > cc_list_adj_Ext_ORI;
   string output_CC_Adj_ori=prefixe+"_Compo_connexe/CC_Adj_EXT_ORI/";
   cout<<"Creation of connected components files for each species in directory "<<output_CC_Adj_ori<<" for input adjacencies file..."<<flush;
   for(it_AA=Adj_actuelles_ORI.begin();it_AA!=Adj_actuelles_ORI.end();it_AA++){
      vector<adjacence> v=(*it_AA).second;
      string name_species=S->getNode((*it_AA).first)->getName();
      string species_dir_CC=output_CC_Adj_ori+"CC_Adj_EXT_ORI_"+name_species;
      mkdir_rec(species_dir_CC.c_str());
      i=1;
//      Add alone genes in Connected Components files for species (name_species)
      map<string,string>::iterator it_ss;
      for(it_ss=gene_species_EXT.begin();it_ss!=gene_species_EXT.end();it_ss++){
         if((*it_ss).second==name_species){
		      ofstream Offile_OUTPUT_CC;
		      string CC_file=species_dir_CC+"/CC_Adj_EXT_ORI_"+name_species+"_"+convertInt(i);
		      Offile_OUTPUT_CC.open(CC_file.c_str(), ios::out|ios::trunc);
		      if(!Offile_OUTPUT_CC){
		         cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
		         if(file_log)
		            Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
		         exit(EXIT_FAILURE);
		      }
	         Offile_OUTPUT_CC<<(*it_ss).first<<endl;
            i++;
         }
      }
      cout<<endl<<"For EXT species "<<name_species<<":"<<endl;
//      Creation of connected components for gene with adjacencies
      (*it_AA).second = composantesConnexes(v,distrib_CC_adj_Ext_ORI,cc_list_adj_Ext_ORI,i);
      map<int,vector<string> >::iterator it_aa;
//      int size_CC=cc_list_adj_Ext_ORI.size();
//      cout<<size_CC<<endl;
      for(it_aa=cc_list_adj_Ext_ORI.begin();it_aa!=cc_list_adj_Ext_ORI.end();it_aa++){
         ofstream Offile_OUTPUT_CC;
         string CC_file=species_dir_CC+"/CC_Adj_EXT_ORI_"+name_species+"_"+convertInt((*it_aa).first);
         Offile_OUTPUT_CC.open(CC_file.c_str(), ios::out|ios::trunc);
         if(!Offile_OUTPUT_CC){
            cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
            if(file_log)
               Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
            exit(EXIT_FAILURE);
         }
         vector<string>::iterator it_vs;
         for(it_vs=((*it_aa).second).begin();it_vs!=((*it_aa).second).end();it_vs++){
            Offile_OUTPUT_CC<<(*it_vs)<<endl;
         }
      }
      cc_list_adj_Ext_ORI.clear();
      distrib_CC_adj_Ext_ORI.clear();
   }
   Adj_actuelles_ORI.clear();
   cout<<" DONE"<<endl<<endl;


   if (is_readable(output_new_adj_file)){
//      Remove genes present in new adjacencies from gene_species_EXT map to keep only alone gene
      for(it_AA=Adj_actuelles.begin();it_AA!=Adj_actuelles.end();it_AA++){
         vector<adjacence> v=(*it_AA).second;
         for(it_adj=v.begin();it_adj!=v.end();it_adj++){
            if (gene_species_EXT.find((*it_adj).gene1)!=gene_species_EXT.end()){
               gene_species_EXT.erase((*it_adj).gene1);
            }
            if (gene_species_EXT.find((*it_adj).gene2)!=gene_species_EXT.end()){
               gene_species_EXT.erase((*it_adj).gene2); 
            }
         }
      }

//      Compo connexe pour les adjacences actuelles présentes dans les fichiers d'adjacences en entrée de DeCo/ARt-DeCo + les nouvelles adjacences inferrées par ARt-DeCo
      map<int,int> distrib_CC_adj_Ext;
      map<int,vector<string> > cc_list_adj_Ext;
      string output_CC_ARt_DeCo=prefixe+"_Compo_connexe/CC_Adj_EXT_ARt-DeCo/";
      cout<<"Creation of connected components files for each species in directory "<<output_CC_ARt_DeCo<<" for output adjacencies file of ARt-DeCo:"<<endl;
      for(it_AA=Adj_actuelles.begin();it_AA!=Adj_actuelles.end();it_AA++){
         vector<adjacence> v=(*it_AA).second;
         string name_species=S->getNode((*it_AA).first)->getName();
         string species_dir_CC=output_CC_ARt_DeCo+"CC_Adj_EXT_ARt-DeCo_"+name_species;
         mkdir_rec(species_dir_CC.c_str());
		   i=1;
//         Add alone genes in Connected Components files for species (name_species)
         map<string,string>::iterator it_ss;
		   for(it_ss=gene_species_EXT.begin();it_ss!=gene_species_EXT.end();it_ss++){
		      if((*it_ss).second==name_species){
				   ofstream Offile_OUTPUT_CC;
				   string CC_file=species_dir_CC+"/CC_Adj_EXT_ARt-DeCo_"+name_species+"_"+convertInt(i);
				   Offile_OUTPUT_CC.open(CC_file.c_str(), ios::out|ios::trunc);
				   if(!Offile_OUTPUT_CC){
				      cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
				      if(file_log)
				         Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
				      exit(EXIT_FAILURE);
				   }
			      Offile_OUTPUT_CC<<(*it_ss).first<<endl;
		         i++;
		      }
		   }
         cout<<endl<<"For EXT species "<<name_species<<":"<<endl;
//         Creation of connected components for gene with adjacencies
         (*it_AA).second = composantesConnexes(v,distrib_CC_adj_Ext,cc_list_adj_Ext,i);
         map<int,vector<string> >::iterator it_aa;
//         int size_CC=cc_list_adj_Ext.size();
//         cout<<size_CC<<endl;
         for(it_aa=cc_list_adj_Ext.begin();it_aa!=cc_list_adj_Ext.end();it_aa++){
            ofstream Offile_OUTPUT_CC;
            string CC_file=species_dir_CC+"/CC_Adj_EXT_ARt-DeCo_"+name_species+"_"+convertInt((*it_aa).first);
            Offile_OUTPUT_CC.open(CC_file.c_str(), ios::out|ios::trunc);
            if(!Offile_OUTPUT_CC){
               cout<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
               if(file_log)
                  Offile_log<<"From Step4_statistics.cpp (main): ERROR while opening file "<<CC_file<<endl;
               exit(EXIT_FAILURE);
            }
            vector<string>::iterator it_vs;
            for(it_vs=((*it_aa).second).begin();it_vs!=((*it_aa).second).end();it_vs++){
               Offile_OUTPUT_CC<<(*it_vs)<<endl;
            }
         }
         cc_list_adj_Ext.clear();
         distrib_CC_adj_Ext.clear();
      }
   }
   Adj_actuelles.clear();
   gene_species_EXT.clear();
   delete(S);
   cout<<" DONE"<<endl<<endl;


   cout<<"##################################################"<<endl;
   cout<<"### END construction of Connected Components ###"<<endl;
   cout<<"##################################################"<<endl<<endl;

   // FIN
   //

   Offile_OUTPUT_stats_human.close();
   Offile_OUTPUT_stats_machine.close();

   time_t tend=time(NULL);                // get the current calendar time    
   // Compute execution time
   float texec=difftime(tend,tbegin);    // tend-tbegin (result in second)

   cout<<"Execution time of Step4 is: "<<texec<<"s"<<endl<<endl;
      cout<<"\t########################"<<endl;
      cout<<"\t###  End Step4 !!!!  ###"<<endl;
      cout<<"\t########################"<<endl<<endl;

   if (file_log){
      Offile_log<<"Execution time of Step4 is: "<<texec<<"s"<<endl<<endl;
      Offile_log<<"\t########################"<<endl;
      Offile_log<<"\t###  End Step4 !!!!  ###"<<endl;
      Offile_log<<"\t########################"<<endl<<endl;
   }
   Offile_log.close();

   return(0);
}
