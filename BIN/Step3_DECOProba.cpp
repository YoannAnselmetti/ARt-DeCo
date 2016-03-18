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

File: Step3_DECOProba.cpp                      Last modified on: 20/08/2015
Created by: Yoann Anselmetti/Sèverine Bérard   Created on: 06/11/2014
--------------------------------------------------------------------------
Specification: 
File containing the DECO function dedicated to pairwise tree computations.
Using P(v1~v2) => Function in CostFunctionsProba.cpp
=========================================================================*/
#include "Step3_DECOProba.h"



void Nb_chr_species(map<string,int> &esp_nb_chr){
   //Ouverture du fichier d'espèces
   ifstream IffilGenome(file_spe.c_str(), ios::in);
   cout<<"Store number of chromosome by species in map<string,int> esp_nb_chr from file "<<file_spe<<"..."<<flush;
   if (file_log)
      Offile_log<<"FI\tFrom Step3_DECOProba.cpp (Nb_chr_species): Store number of chromosome by species in map<string,int> esp_nb_chr from file "<<file_spe<<"..."<<flush;
   if (!IffilGenome){
      cout<<endl<<"ERROR while opening file "<<file_spe<<endl;
      if (file_log)
         Offile_log<<"FI\tFrom Step3_DECOProba.cpp (Nb_chr_species): ERROR while opening file "<<file_spe<<endl<<endl;
      exit(EXIT_FAILURE);
   }
   while (!IffilGenome.eof()){
      string genome;
      int nb_chr;
      IffilGenome>>genome;
      IffilGenome>>nb_chr;
      if(genome!=""){
         esp_nb_chr[genome]=nb_chr;
//         cout<<"genome: "<<genome<<" | chr: "<<nb_chr<<endl;
      }
   }
   IffilGenome.close();
   cout<<" DONE"<<endl<<endl;
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
}


//float CalculPv1_v2(float p, float n){
//   //Nombre de chr et contig dans l'espèce de v1 et v2
//   float Proba,i,j;
//   float sum=0.0;


//   
//   
//   //Calcul de la Somme (sum) de x=2 à x=(n-p+1)
//   for(float x=2;x<=n-p+1;x++){
//      //cout<<"Pour x="<<x<<" :"<<endl;
//      //Calcul du 1er Produit (prodA) de i=1 à i=(x-1)
//      float prodA=1.0;
//      for(i=1.0;i<=(x-1.0);i++){
//         prodA=i/(n-i)*prodA;         
//         // cout<<"\t"<<"Pour i="<<i<<" :"<<endl;
//         // cout<<"\t"<<"\t"<<"i/(n-i)="<<i/(n-i)<<endl;
//      }
//      //Calcul du 2nd Produit (prodB) de j=(n-p+2) à j=n
//      float prodB=1.0;
//      for(j=(n-p+2.0);j<=n;j++){
//         prodB=(j-x)/j*prodB;
//         // cout<<"\t"<<"Pour j="<<j<<" :"<<endl;
//         // cout<<"\t"<<"\t"<<"(j-x)/j="<<(j-x)/j<<endl;
//      }

//      //SEV : ATTENTION DIVISION PAR 0 si nb chromo p=1, x va jusqu'à n dans la boucle
//      //=> VERIFIER LA FORMULE PAPIER POUR CE CAS 
//      sum+=prodA*prodB/((n-x)*(n-p+1));
//      
//      // cout<<"\t"<<"prodA="<<prodA<<endl;
//      // cout<<"\t"<<"prodB="<<prodB<<endl;
//      // cout<<"\t"<<"(n-x)*(n-p+1)="<<(n-x)*(n-p+1)<<endl;
//      // cout<<"\t"<<"prodA*prodB/(n-x)*(n-p+1)="<<prodA*prodB/((n-x)*(n-p+1))<<endl<<endl;
//   }

//   // cout<<"-------------------------------"<<endl<<endl;
//   
//   // cout<<"sum="<<sum<<endl<<endl;
//   
//   //Calcul P(v1~v2)
//   Proba=(p*(p-1)/2)*sum;
//   
//   // cout<<"-------------------------------"<<endl;
//   // cout<<"-------------------------------"<<endl<<endl;

//   return Proba;
//}

void Association_species_chr_contig_nb(TreeTemplate<Node> * S, vector<Tree*> &Arbres, map<string,int> &esp_nb_chr, map<string,vector<float> > &esp_nb_contig, map<int,int> &extant_species_gene_nb, map<int,vector<adjacence> > Adj_actuelles){
   vector<float> chr_contig_base_log;

   //Création map<nom_esp,nb_chr> esp_nb_chr
   Nb_chr_species(esp_nb_chr);
//   afficheMap(esp_nb_chr);
   vector<Tree *>::iterator it;
   //Création de la map<int,int> extant_species_gene_nb <Id_esp,nb_gene_actuel>
   int y=0;
   for(it=Arbres.begin();it!=Arbres.end();it++){ 
      Tree * Gp=*it;
      TreeTemplate<Node> *G = dynamic_cast <TreeTemplate<Node> *> (Gp);
      y+=countNodesWithBranchPropertyBySpecies(G->getRootNode(),typ,ga,extant_species_gene_nb);
   }

   string species_stats_file=prefixe+output_species_stats;
   ofstream Offile_species_stats;
   if(!light_mode){
      //Déterminer nombre de gènes avec 0 et 1 voisin dans l'espèce dans laquelle se situe v1 et v2
      //map<int,vector<adjacence>> Adj_actuelles <Id_esp,vector<adjacence>>
      string species_stats_file=prefixe+output_species_stats;
      Offile_species_stats.open(species_stats_file.c_str(), ios::out|ios::trunc);
      if (!Offile_species_stats){
         cout<<"\nFrom Step3_DECOProba.cpp (Association_species_chr_contig_nb()): ERROR while opening file "<<species_stats_file<<endl;
         if (file_log)
            Offile_log<<"\nFI\tFrom Step3_DECOProba.cpp (Association_species_chr_contig_nb()): ERROR while opening file "<<species_stats_file<<endl;
         exit(EXIT_FAILURE);
      }
      Offile_species_stats<<"Species_name\t#G_tot\t#G_0A\t#G_1A\t#G_2A\t#G_>2A\t#A_tot\t#Chr\t#Contig\tP(Adj)\tc0(0,0)\tc0(1,0/0,1)\tc0(1,1)\tc1(0,0)\tc1(1,0/0,1)\tc1(1,1)\tBase_log"<<endl;
   }
   map<string,int> extant_gene_adj_nb3; //le 3 pour différencier des variables de meme nom NON! sinon perdu pour calcul des couts !!!
   map<string,int>::iterator it3;
   map<int,vector<adjacence> >::iterator it_AA;
   for(it_AA=Adj_actuelles.begin();it_AA!=Adj_actuelles.end();it_AA++){
      vector<adjacence> v=(*it_AA).second;
      vector<adjacence>::iterator it_v;
      // int gene_nb_0_adj=0;
      // int gene_nb_1_adj=0;
      // int gene_nb_2_adj=0;
      int nb_gene_with_adj_esp=0;
      string species=S->getNode((*it_AA).first)->getName();
      extant_gene_adj_nb3.clear();//remettre à zéro la map à chaque nouvelle espèce pour la variable locale !!! 
      for(it_v=v.begin();it_v!=v.end();it_v++){
//         cout<<(*it_v).gene1<<"-"<<(*it_v).gene2<<endl;
         if(extant_gene_adj_nb3.find((*it_v).gene1)!=extant_gene_adj_nb3.end()){
            extant_gene_adj_nb3[(*it_v).gene1]++;
            extant_gene_adj_nb[(*it_v).gene1]++;
            //gene_nb_2_adj++;
         }
         else{
            extant_gene_adj_nb3[(*it_v).gene1]=1;
            extant_gene_adj_nb[(*it_v).gene1]=1;
            //nb_gene_with_adj_esp++;
         }
         if(extant_gene_adj_nb3.find((*it_v).gene2)!=extant_gene_adj_nb3.end()){
            extant_gene_adj_nb3[(*it_v).gene2]++;
            extant_gene_adj_nb[(*it_v).gene2]++;
            //gene_nb_2_adj++;
         }
         else{
            extant_gene_adj_nb3[(*it_v).gene2]=1;
            extant_gene_adj_nb[(*it_v).gene2]=1;
            //nb_gene_with_adj_esp++;
         }
      }
      nb_gene_with_adj_esp=extant_gene_adj_nb3.size();
      cout<<"\n**************** STATS for "<<species<<" ********************"<<endl;
      // cout<<"**** gene_nb_2_adj="<<gene_nb_2_adj<<endl;
      // cout<<"**** nb_gene_with_dj_esp="<<nb_gene_with_adj_esp<<endl;
      // cout<<"**** Total gene number="<<extant_species_gene_nb[(*it_AA).first]<<endl;
      cout<<"**** Gene number with at least one adjacency="<<nb_gene_with_adj_esp<<endl;
      if (file_log){
		   Offile_log<<"\n**************** STATS for "<<species<<" ********************"<<endl;
		   // Offile_log<<"**** gene_nb_2_adj="<<gene_nb_2_adj<<endl;
		   // Offile_log<<"**** nb_gene_with_adj_esp="<<nb_gene_with_adj_esp<<endl;
		   // Offile_log<<"**** Total gene number="<<extant_species_gene_nb[(*it_AA).first]<<endl;
		   Offile_log<<"**** Gene number with at least one adjacency="<<nb_gene_with_adj_esp<<endl;
      }
      
      int nb0=0,nb1=0,nb2=0,nbplus=0;
      for(it3=extant_gene_adj_nb3.begin(); it3!=extant_gene_adj_nb3.end();it3++)
         if ((*it3).second==1)
            nb1++;
         else if ((*it3).second==2)
            nb2++;
         else
            nbplus++;
      cout<<"**** nb1="<<nb1<<" nb2="<<nb2<<" nbplus="<<nbplus<<endl;
      if (file_log)
	      Offile_log<<"**** nb1="<<nb1<<" nb2="<<nb2<<" nbplus="<<nbplus<<endl;
         
      if(nbplus>0){
         cout<<"FI\tFrom Step3_DECOProba.cpp (Association_species_chr_contig_nb): ERROR "<<nbplus<<" genes with more than 2 adjacencies"<<endl;
         if (file_log)
            Offile_log<<"FI\tFrom Step3_DECOProba.cpp (Association_species_chr_contig_nb): ERROR "<<nbplus<<" genes with more than 2 adjacencies"<<endl;
         exit(EXIT_FAILURE);
      }
      // if (species=="Homo_sapiens")
      //     for(it3=extant_gene_adj_nb3.begin(); it3!=extant_gene_adj_nb3.end();it3++)
      //        cout<<(*it3).first<<endl;

      //Nombre de gènes avec 0 voisin dans l'esp de v1 et v2 = nb de gènes de l'espèce présents dans les arbres - ceux avec des adj
      nb0=(extant_species_gene_nb[(*it_AA).first])-nb_gene_with_adj_esp;

      //Nombre de gènes avec 1 voisin dans l'esp de v1 et v2
      //gene_nb_1_adj=nb_gene_with_adj_esp-gene_nb_2_adj;
      if (nb1%2!=0){
         cout<<"\tERROR the number of genes having only 1 adjacency isn't peer (Number of contig/chromosome extremities isn't peer)!!!"<<endl;
         if (file_log)
            Offile_log<<"\nOS\tFrom Step3_DECOProba.cpp (Association_species_chr_contig_nb): ERROR the number of genes having only 1 adjacency isn't peer (Number of contig/chromosome extremities isn't peer)!!!"<<endl;
         exit(EXIT_FAILURE);
      }
      
      //Nombre de contigs dans l'esp de v1 et v2
      cout<<"*************************** Calcul nb contig"<<endl;
      cout<<"gene_nb_0_adj : "<<nb0<<endl;
      cout<<"gene_nb_1_adj : "<<nb1<<endl;
      int contig=nb0+(nb1/2);
      cout<<"contig : "<<contig<<endl;
      if (file_log){
		   Offile_log<<"*************************** Calcul nb contig"<<endl;
		   Offile_log<<"gene_nb_0_adj : "<<nb0<<endl;
		   Offile_log<<"gene_nb_1_adj : "<<nb1<<endl;
		   Offile_log<<"contig : "<<contig<<endl;

      }
      
      //Nombre de chr dans l'esp de v1 et v2
      int chr=esp_nb_chr[species];

      if(chr>contig){
         cout<<"\tERROR nb_chr: "<<chr<<" > "<<contig<<" (=nb_contig in species "<<species<<")"<<endl;
         if (file_log)
            Offile_log<<"\nOS\tFrom Step3_DECOProba.cpp (Association_species_chr_contig_nb): ERROR nb_chr: "<<chr<<" > "<<contig<<" (=nb_contig in species "<<species<<")"<<endl;
         exit(EXIT_FAILURE);
      }
      
      //Compute P(v1~v2) for rho(v1~v2)=1
      float Proba=float(contig-chr)/float(2*contig*(contig-1));
//      float Proba=CalculPv1_v2(chr,contig);

      //Compute ln(base log) to compute c0P(v1~v2) & c1P(v1~v2)      
      float base_log;
      if (Proba==0)
         base_log=2.0;
      else
         base_log=ceil(pow((1.0-Proba)/Proba,1.0/Break))*coeff_b;   //e^(x*ln(a))=a^(x)
      
      float c1P_11=-log(Proba)/log(base_log);
      float c1P_01=-log(2*Proba)/log(base_log);
      float c1P_00=-log(4*Proba)/log(base_log);
      float c0P_11=-log(1.0-Proba)/log(base_log);
      float c0P_01=-log(1.0-(2*Proba))/log(base_log);
      float c0P_00=-log(1.0-(4*Proba))/log(base_log);

      
      //Remplissage de la map<string,vector<float> > esp_nb_contig
      chr_contig_base_log.clear();
      chr_contig_base_log.push_back(chr);
      chr_contig_base_log.push_back(contig);
      chr_contig_base_log.push_back(base_log);
      chr_contig_base_log.push_back(Proba);
      
      esp_nb_contig[species]=chr_contig_base_log;

      //Affichage nom_esp : nb_chr / nb_contig / base_log   / P(v1~v2)
      cout<<species<<": "<<endl;
      cout<<"\tnb_chr:"<<esp_nb_contig[species][0]<<endl;
      cout<<"\tnb_contig:"<<esp_nb_contig[species][1]<<endl;
      cout<<"\tbase_log:"<<base_log<<endl;
      cout<<"\tP(v1~v2):"<<esp_nb_contig[species][3]<<endl;
      cout<<"\tc0P:"<<c0P_11<<endl;
      cout<<"\tc1P:"<<c1P_11<<endl;
      if(file_log){
		   Offile_log<<species<<": "<<endl;
		   Offile_log<<"\tnb_chr:"<<esp_nb_contig[species][0]<<endl;
		   Offile_log<<"\tnb_contig:"<<esp_nb_contig[species][1]<<endl;
		   Offile_log<<"\tbase_log:"<<base_log<<endl;
		   Offile_log<<"\tP(v1~v2):"<<esp_nb_contig[species][3]<<endl;
		   Offile_log<<"\tc0P:"<<c0P_11<<endl;
		   Offile_log<<"\tc1P:"<<c1P_11<<endl;
      }


      if(!light_mode)
         Offile_species_stats<<species<<"\t"<<extant_species_gene_nb[(*it_AA).first]<<"\t"<<nb0<<"\t"<<nb1<<"\t"<<nb2<<"\t"<<nbplus<<"\t"<<v.size()<<"\t"<<esp_nb_contig[species][0]<<"\t"<<esp_nb_contig[species][1]<<"\t"<<esp_nb_contig[species][3]<<"\t"<<c0P_00<<"\t"<<c0P_01<<"\t"<<c0P_11<<"\t"<<c1P_00<<"\t"<<c1P_01<<"\t"<<c1P_11<<"\t"<<base_log<<endl;
   }
   if(!light_mode)
      Offile_species_stats.close();
}

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
//      return(1);
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
void DECO(TreeTemplate<Node>* S, TreeTemplate<Node>* Arbre1, TreeTemplate<Node>* Arbre2, vector<adjacence> * Adj_classe, ofstream& OfficAdj2, ofstream & OfficArbresAdj, ofstream & OfficSORTIE_dup, string nom_classe, map<string,string> &gene_species_EXT, map<string,vector<float> > &esp_nb_contig){
   if(file_log)
      Offile_log<<"\nOn traite la classe "<<nom_classe<<endl;
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
   C0P.clear();C1P.clear();

   if(file_log)
      Offile_log<<"\t Calcul de la matrice de couts ... "<<flush;
   for(itv1=noeudsA1.begin();itv1!=noeudsA1.end();itv1++)
      for(itv2=noeudsA2.begin();itv2!=noeudsA2.end();itv2++)
         if (Espece(*itv1)==Espece(*itv2)){
            pair<Node *,Node *> p((*itv1),(*itv2));
            pair<Node *,Node *> p2((*itv2),(*itv1));
            C1P[p]=calculeC1P(*itv1,*itv2,Adj_classe,S,gene_species_EXT,esp_nb_contig);   //Calcul C1P
            C1P[p2]=C1P[p];
            if(affich_CalculCout||affich_ParcoursArriere){
               cout<<"\t C1P("<<(*itv1)->getId()<<","<<(*itv2)->getId()<<").first="<<C1P[p].first<<endl;
               cout<<"\t C1P("<<(*itv1)->getId()<<","<<(*itv2)->getId()<<").second="<<C1P[p].second<<endl;
            }
            C0P[p]=calculeC0P(*itv1,*itv2,Adj_classe,S,gene_species_EXT,esp_nb_contig);   //Calcul C0P
            C0P[p2]=C0P[p];
            if(affich_CalculCout||affich_ParcoursArriere){
               cout<<"\t C0P("<<(*itv1)->getId()<<","<<(*itv2)->getId()<<").first="<<C0P[p].first<<endl;
               cout<<"\t C0P("<<(*itv1)->getId()<<","<<(*itv2)->getId()<<").second="<<C0P[p].second<<endl;
            }
         }
   if(file_log)
      Offile_log<<"OK"<<endl;
   float c1P=recupCoutC1P(Arbre1->getRootNode(),Arbre2->getRootNode());
   float c0P=recupCoutC0P(Arbre1->getRootNode(),Arbre2->getRootNode());
       
   if(affich_CalculCout||affich_ParcoursArriere){
      cout<<"\n**********************"<<endl;
      cout<<"C0P aux racines : "<<c0P<<endl;
      cout<<"C1P aux racines (+Gain): "<<c1P+Crea<<endl;
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
   bool avantageC1=false;
   if (r<=Adj_percentage) //en cas C0P=C1P on choisit C1P
       avantageC1=true;
    else
       avantageC1=false;

   if(file_log)
      Offile_log<<"\t Parcours arriere ... "<<flush;
   if (c0P>c1P+Crea || (c1P+Crea==c0P && avantageC1)){
      //if (c0P>c1P && c0P-c1P > 0.00001){ //Pour éviter d'entrer ici en cas d'égalité entre les flottants
      //Création du premier arbre d'adjacence
      TreeTemplate<Node> * T = new TreeTemplate<Node>();
      T->setRootNode(parcoursArriereC1(Arbre1->getRootNode(),Arbre2->getRootNode(),Adj_classe,ArbresDAdjacences,true));
      ArbresDAdjacences->push_back(T);
   }
   else{
      parcoursArriereC0(Arbre1->getRootNode(),Arbre2->getRootNode(),Adj_classe,ArbresDAdjacences,true);
   }
   if(file_log)
      Offile_log<<"OK ..."<<flush;

   vector<Tree *>::iterator it;
   if(!light_mode){
      if(affich_CalculCout||affich_ParcoursArriere)
         cout<<"\nÉcriture des arbres d'adjacences et des adjacences ancestrales de cette classe dans les fichiers correspondants."<<endl;
      OfficArbresAdj<<"Adjacency tree(s) solution of class "<<nom_classe<<endl;   //<<" ( c0(R1,R2)="<<c0P<<" | c1(R1,R2)="<<c1P+Crea<<" )"

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
   }

   //On libère la mémoire du vecteur d'arbres d'adjacences
   for(it=ArbresDAdjacences->begin();it!=ArbresDAdjacences->end();it++)
      delete(*it);
   ArbresDAdjacences->clear();
   delete(ArbresDAdjacences);

   if(file_log)
      Offile_log<<"\t COMPLET"<<endl;
}

int main(int argc, char* argv[]){
   cout<<"\n\t###############################"<<endl;
   cout<<"\t###  Start Step3_proba !!!  ###"<<endl;
   cout<<"\t###############################"<<endl<<endl;
   time_t tbegin=time(NULL);              // get the current calendar time

   // Récupérer données fichier de config pour recéer gene_species_EXT et pour récupérer le préfixe pour rechercher le fichier de classes d'adjacences
   lireFicConfig(fic_arbre,fic_gene,fic_especes,fic_adjacence,exp_name,directory,argc,argv);
   if (is_readable(fic_adjacence_clean))
      fic_adjacence=fic_adjacence_clean;
   //Idem pour le fichier de gènes
   if (is_readable(fic_gene_clean))
      fic_gene=fic_gene_clean;
   
   cout<<"File(s) used in Step3_proba:"<<endl;
   cout<<"\tReconciled trees file: "<<tree_reconciled_file<< " or "<<fic_arbre<<endl;
   cout<<"\tAdjacencies classes file: "<<adj_class_file<<endl;
   cout<<"\tSelected species file: "<<file_spe<<endl<<endl;

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
            cout<<"\nERROR while opening file "<<log_file<<endl<<endl;
            exit(EXIT_FAILURE);
         }
         else{
            Offile_log<<"####################################################################################"<<endl;
            Offile_log<<"Log file open by Step3_DECOProba.cpp at "<<str_time<<endl;
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
            Offile_log<<"Log file created by Step3_DECOProba.cpp at "<<str_time<<endl;
            Offile_log<<"####################################################################################"<<endl<<endl;
         }
      }
      Offile_log<<"\t###############################"<<endl;
      Offile_log<<"\t###  Start Step3_proba !!!  ###"<<endl;
      Offile_log<<"\t###############################"<<endl<<endl;
      Offile_log<<"\tFile(s) used in Step3_proba:"<<endl;
      Offile_log<<"\t\t- Reconciled trees file: "<<tree_reconciled_file<< " or "<<fic_arbre<<endl;
      Offile_log<<"\tAdjacencies classes file: "<<adj_class_file<<endl;
      Offile_log<<"\t\t- Selected species file: "<<file_spe<<endl<<endl;
   }


   // Variables
//   TreeTemplate<Node> S;
//   S = new TreeTemplate<Node>;
   map<string,string> gene_species_EXT;
   vector<Tree*> Arbres;
   map<string,vector<int> > gene_GeneTreeID;
   map<string,vector<int> > classes_adjacences;
   vector<adjacence> adjacencies;
   map<string,int> esp_nb_chr;
   map<string,vector<float> > esp_nb_contig;
   map<int,int> extant_species_gene_nb;
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
      cout<<"\nFrom Step3_DECOProba (main): WARNING while cleaning "<<elim<<" genes have been deleted"<<endl<<endl;
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

   //Création map<nom_esp,<nb_chr,nb_contig>>
   Association_species_chr_contig_nb(S,Arbres,esp_nb_chr,esp_nb_contig,extant_species_gene_nb,Adj_actuelles);
   Adj_actuelles.clear();
   esp_nb_chr.clear();
   extant_species_gene_nb.clear();

   ////////////////////////////////////////////////////////////////////////////////////
   //Suprression fichier sortie des nouvelles adjacences et nouveau fichier des adjacences
   string output_new_adj_file=prefixe+output_new_adj;
   ofstream NewAdjFile (output_new_adj_file.c_str(), ios::out|ios::trunc);
   if(!NewAdjFile){
      cout<<"\nFrom Step3_DECOProba.cpp (main): ERROR while removing file "<<output_new_adj_file<<endl<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step3_DECOProba.cpp (main): ERROR while removing file "<<output_new_adj_file<<endl<<endl;
      exit(EXIT_FAILURE);
   }
   else{
      cout<<endl<<"File "<<output_new_adj_file<<" is created or crushed"<<endl<<endl;
   }
   NewAdjFile.close();

   ////////////////////////////////////////////////////////////////////////////////////
   //Copy file fic_adjacence in output_adj_file
   if(!light_mode){
      string output_adj_file=prefixe+output_ext_adj_file;
      if (file_log)
         Offile_log<<"\nFI\tFrom Step3_DECOProba.cpp: Copy file "<<fic_adjacence<<" in "<<output_adj_file<<"..."<<flush;
      ifstream  src(fic_adjacence.c_str(), ios::binary);
      ofstream  dst(output_adj_file.c_str(), ios::binary|ios::trunc);
      if(!src||!dst)
         {
       cout<<"From Step3_DECOProba.cpp (main): ERROR copy file fic_adjacence in output_adj_file"<<endl;
       if (file_log)
          Offile_log<<"From Step3_DECOProba.cpp (main): ERROR copy file fic_adjacence in output_adj_file"<<endl;
       exit(EXIT_FAILURE);
         }
      dst << src.rdbuf();
      if (file_log)
         Offile_log<<" DONE"<<endl<<endl;
      src.close();
      dst.close();
   }

  
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
         cout<<"\tOn traite la classe "<<nom_classe<<endl;
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
         TreeTemplate<Node > A1=*(Arbres[num1]);
         TreeTemplate<Node > A2=*(Arbres[num2]);
         Arbre1= new TreeTemplate<Node>();
         Arbre1=A1.cloneSubtree(coupe1);
         Arbre2= new TreeTemplate<Node>();
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
         Arbre1= new TreeTemplate<Node>();
         Arbre1=ArbreUnique.cloneSubtree(Id_Genes[0]);
         Arbre2= new TreeTemplate<Node>();
         Arbre2=ArbreUnique.cloneSubtree(Id_Genes[1]); 

         num_ident++;
      }
      
      
      //         b) Récupérer les adjacences
      for(it_adj=(*it_cl).second.begin();it_adj!=(*it_cl).second.end();it_adj++){
         Adj_classe->push_back(adjacencies[*it_adj]);
      }
      
      //         c) Appliquer l'algo DéCo (produit des adjacences ancestrales 
      //qu'on accumule dans un fichier de même type que le fichier adjacences 
      //d'entrée mais où les noms d'espèces sont ceux des espèces ancestrales)
      if(INPUT_FORMAT==1){
         affecteInfoSurBrancheNoeuds(Arbre1->getRootNode());
         affecteInfoSurBrancheNoeuds(Arbre2->getRootNode());
      }
      affecteInfoSurBrancheNoeuds(Arbre2->getRootNode());
     
      /////////////////////////////////////////////////////////////////
      ///////////////  APPEL DE L'ALGO DECO ///////////////////////////
      /////////////////////////////////////////////////////////////////
      DECO(S,Arbre1,Arbre2,Adj_classe,OfficAdj2,OfficArbresAdj,OfficSORTIE_dup,nom_classe,gene_species_EXT,esp_nb_contig);
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
   esp_nb_contig.clear();
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
   
   if(!light_mode){
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
   }
   

   time_t tend=time(NULL);                // get the current calendar time
   // Compute execution time
   float texec=difftime(tend,tbegin);    // tend-tbegin (result in second)

   cout<<"Execution time of Step3_proba is: "<<texec<<"s"<<endl<<endl;
   cout<<"\t#############################"<<endl;
   cout<<"\t###  End Step3_proba !!!  ###"<<endl;
   cout<<"\t#############################"<<endl<<endl;

   if (file_log){
      Offile_log<<"Execution time of Step3_proba is: "<<texec<<"s"<<endl<<endl;
      Offile_log<<"\t#############################"<<endl;
      Offile_log<<"\t###  End Step3_proba !!!  ###"<<endl;
      Offile_log<<"\t#############################"<<endl<<endl;
   }
   Offile_log.close();

   return(0);
}
