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

File: General.cpp                              Last modified on: 20/08/2015
Created by: Yoann Anselmetti/Sèverine Bérard   Created on: 28/02/2014
--------------------------------------------------------------------------
Specification:
General file of the project to store global variables and methods
=========================================================================*/
#include "General.h"

// Parameters to open Config File   //
string fic_arbre;
string fic_gene;
string fic_especes;
string fic_adjacence;
string exp_name;
string directory;
string prefixe;
string ST_input;
string file_spe;
string dir_embl;
string dir_emf;
string del_spe_data;
string log_file;
ofstream Offile_log;

//Parameters from Config File
char sep;
char sepAdj='-';
int INPUT_FORMAT; //For gene trees (Species tree on Newick format)
int OUTPUT_FORMAT;//0 for Newick ; 1 for NHX 
float Adj_percentage;
float Crea;
float Break;
float coeff_b;
//float Spec;
//float GDup;
//float GLos;
//float ADup;
//float ALos;
bool ArbresRec;
string fic_gene_clean;
string fic_adjacence_clean;

// OUTPUT and intermediate files   //
string tree_reconciled_file;
string adj_class_file;
string DECO_stats_file;

//   BackTracking + Step3 + Step4   //
int nb_GDup=0;//comptés, comme les autres, dans Backtracking
int nb_ADup=0;
int nb_Crea=0;
int nb_Bk=0;


////////////////
//  BOOLEAN   //
////////////////

// Boolean for programm mode
bool proba_mode=true;
bool light_mode=false;   // If light_mode = true: Put bool fic_SORTIE_EspeceDone = true, bool fic_SORTIE_EspeceDone = true & bool file_log = false

// Boolean for output screen
bool affich_CalculCout=false;
bool affich_ParcoursArriere=false;
bool affich_Classes=false;

// Boolean for files creation
bool fic_SORTIE_EspeceDone=false;
bool fic_SORTIE_GeneDone=false;
bool file_log=true;


////////////////////////
///   MAP & VECTOR   ///
////////////////////////

// Step3_proba, Step_4 & Backtracking(Proba)
map<int,int> espece_nb_gain; //ces deux maps pour compter les événements sur les branches
map<int,int> espece_nb_break; // de l'arbre des espèces   => Pas utilisé Dans Backtracking(Proba)
map<string,int> extant_gene_adj_nb;   //string: nom gène | int: nb_adj du gène
vector<adjacence> DUP; //Pour recueillir les gènes dupliqués ensemble (sert pour des calculs intermédiaires)
int id_arbres_adj;

// Extern declaration for CostFunctions.cpp
map<pair<Node *,Node *>, float > C0;
map<pair<Node *,Node *>, float > C1;
map<pair<Node *,Node *>, pair<float,int> > C0P;
map<pair<Node *,Node *>, pair<float,int> > C1P;
map<pair<string,string>, float> P;   //=> Pour le modèle P(v1~v2) de DeCo




////   METHODS   ////

// BASIC METHODS:

void mkdir_rec(const char *dir){
   char tmp[256];
   char *p = NULL;
   size_t len;

   snprintf(tmp, sizeof(tmp),"%s",dir);
   len = strlen(tmp);
   if(tmp[len - 1] == '/')
      tmp[len - 1] = 0;
   for(p = tmp + 1; *p; p++)
      if(*p == '/') {
         *p = 0;
         mkdir(tmp, S_IRWXU);
         *p = '/';
      }
   mkdir(tmp, S_IRWXU);
}


string convertInt(int number){
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

string convertFloat(float number){
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

bool file_exists(const string& name){
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

//Pour comparer deux flottants, remplace le == qui ne fonctionne pa à cause des imprecisions
bool areEqual(float a,float b){
   float err=0.00001;
   //cout<<"a="<<a<<" b="<<b<<"abs(a - b)="<<abs(a - b)<<"  max(abs(a), abs(b))="<< max(abs(a), abs(b))<<" calcul="<< err * max(abs(a), abs(b)) * numeric_limits<float>::epsilon()<<endl;
   return abs(a - b) <=  err;
   //return a==b;
}

// Bool to know if a file exists
bool is_readable(string& file){
   ifstream fichier(file.c_str(), ios::in);
   bool ouverture_reussie=!fichier.fail();
   fichier.close();
   return ouverture_reussie;
}

//METHODS to visualize Map (afficheMap)

//Pour afficher une map <string, vector<int>>
void afficheMap(map<string, vector<int> > m){
   map<string, vector<int> >::iterator it;
   for (it=m.begin();it!=m.end();it++){
      cout<<(*it).first<< " :"<<flush;
      vector<int> v=(*it).second;
      vector<int>::iterator it2;
      for (it2=v.begin();it2!=v.end();it2++)
          cout<<" "<<*it2<<flush;
      cout<<endl;
   }
   cout<<endl;
}

//Pour afficher une map <int, vector<adjacence>>
void afficheMap(map<int, vector<adjacence> > m){
   map<int, vector<adjacence> >::iterator it;
   for (it=m.begin();it!=m.end();it++){
      cout<<(*it).first<< " :"<<flush;
      vector<adjacence> v=(*it).second;
      vector<adjacence>::iterator it2;
      for (it2=v.begin();it2!=v.end();it2++)
         cout<<" "<<(*it2).gene1<<sepAdj<<(*it2).gene2<<flush;
      cout<<endl;
   }
   cout<<endl;
}

void afficheMap(map<int, vector<string> > m){
map<int, vector<string> >::iterator it;
   for (it=m.begin();it!=m.end();it++){
      cout<<(*it).first<< " :"<<flush;
      vector<string> v=(*it).second;
      vector<string>::iterator it2;
      for (it2=v.begin();it2!=v.end();it2++)
         cout<<" "<<(*it2)<<flush;
      cout<<endl;
   }
   cout<<endl;
}

//Pour afficher une map <string, string>
void afficheMap(map<string, string> m){
   map<string, string>::iterator it;
   for (it=m.begin();it!=m.end();it++){
      cout<<(*it).first<< " : "<<flush;
      cout<<(*it).second<<endl;
   }
}

//Pour afficher une map <string, int>
void afficheMap(map<string, int> m){
   map<string, int>::iterator it;
   for (it=m.begin();it!=m.end();it++){
      cout<<(*it).first<< " : "<<flush;
      cout<<(*it).second<<endl;
   }
}

//Pour afficher une map <int, int>
void afficheMap(map<int, int> m){
   map<int, int>::iterator it;
   for (it=m.begin();it!=m.end();it++){
      cout<<(*it).first<< " : "<<flush;
      cout<<(*it).second<<endl;
   }
}




//   METHODS to read DeCo config file //

//Pour lire les noms de fichiers dans le fichier de config
int lireNomFic(ifstream& IfficConf, string& fic, string nomfic, char** argv){
   string tampon;
   convert<int> c;
   static const size_t npos = -1;
   IfficConf>>tampon;
   if (tampon==nomfic){
      IfficConf>>fic;
      if (fic.find("arg")!=npos){
         int a=c.from_s(fic.substr(3,fic.length()-1));
         fic=argv[a];
      }
      return(1);
   }
   else{
      cout<<"\nFrom General.cpp (lireNomFic): ERROR in the configuration file at "<<nomfic<<endl;
//      if(file_log)
//         Offile_log<<"\nFI\tFrom General.cpp (lireNomFic): ERROR in the configuration file at "<<nomfic<<endl;
      return(0);
   }
}

//Pour lire les coûts dans le fichier de config
int lireCout(ifstream& IfficConf, float& Cout, string nomcout){
   string tampon;
   convert<float> c;

   IfficConf>>tampon;
   if (tampon==nomcout){
      IfficConf>>tampon;
      Cout=c.from_s(tampon);
      return(1);
   }
   else{
      cout<<"\nFrom General.cpp (lireCout): ERROR in the configuration file at "<<nomcout<<endl;
//      if(file_log)
//         Offile_log<<"\nFI\tFrom General.cpp (lireCout): ERROR in the configuration file at "<<nomcout<<endl;
      return(0);
   }
}

void lireFicConfig(string& fic_arbre, string& fic_gene, string& fic_especes, string& fic_adjacence, string& exp_name, string& directory,int argc, char** argv){
   string ficConfig;
   if (argc==1)
      ficConfig="DeCo.conf";
   else if (argc==2)
      ficConfig=argv[1];
   else{
      cout<<"\nFrom General.cpp (lireFicConfig): Misuse of de "<<argv[0]<<", try:"<<endl;
      cout<<"\t"<<argv[0]<<endl;
      cout<<"(the default configuration file DeCo.conf is used)"<<endl;
      cout<<"or"<<endl;
      cout<<"\t"<<argv[0]<<" ConfigFile"<<endl;
      cout<<"(for another configuration file)"<<endl;
      exit(EXIT_FAILURE);
   }
   ifstream IfficConf (ficConfig.c_str(), ios::in);
   if (!IfficConf){
      cout<<"\nFrom General.cpp (lireFicConfig): ERROR while opening configuration file "<<ficConfig<<endl;
//      if(file_log)
//         Offile_log<<"\nFI\tFrom General.cpp (lireFicConfig): ERROR while opening configuration file "<<ficConfig<<endl;
      exit(EXIT_FAILURE);
   } 
   cout<<"Reading configuration file "<<ficConfig<<endl<<endl;
   
   int bp=0;
   //Les noms des fichiers
   bp+=lireNomFic(IfficConf,fic_arbre,"trees_file",argv);
   bp+=lireNomFic(IfficConf,fic_gene,"genes_file",argv);
   bp+=lireNomFic(IfficConf,fic_especes,"species_file",argv);
   bp+=lireNomFic(IfficConf,fic_adjacence,"adjacencies_file",argv);
   bp+=lireNomFic(IfficConf,exp_name,"exp_name",argv);
   bp+=lireNomFic(IfficConf,directory,"directory",argv);
   if (bp!=6){
      cout<<"From General.cpp (lireFicConfig):  ERROR in the configuration file missing line(s) in configuration file "<<ficConfig<<endl;
      exit(EXIT_FAILURE);
   }
   
   //Attribution des chemins pour les fichiers de sortie
   prefixe=directory+"/"+exp_name;
   log_file=prefixe+"_log_file";
   tree_reconciled_file=prefixe+output_rec;
   adj_class_file=prefixe+output_adj_class;
   DECO_stats_file=prefixe+output_DECO_stats;
   
   //Pour les fichiers d'entrée après nettoyage
   fic_gene_clean=prefixe+"_"+fic_gene.substr(fic_gene.rfind('/')+1,fic_gene.length()-1)+"_CLEAN";
   fic_adjacence_clean=prefixe+"_"+fic_adjacence.substr(fic_adjacence.rfind('/')+1,fic_adjacence.length()-1)+"_CLEAN";

   convert<int> c;
   //arbres déjà réconciliés ou pas
   string tampon;
   IfficConf>>tampon;
   if (tampon=="ReconcilDone"){
      IfficConf>>tampon;
      if (tampon=="true")
         ArbresRec=true;
      else if (tampon=="false")
         ArbresRec=false;
      else
         cout<<"From General.cpp (lireFicConfig):  ERROR in the configuration file at ReconcilDone (Wrong value: should be 'true' or 'false'"<<endl;
   }
   else{
      cout<<"From General.cpp (lireFicConfig):  ERROR in the configuration file at ReconcilDone"<<endl;
      exit(EXIT_FAILURE);
   }
   
   //les formats d'entrée/sortie
   IfficConf>>tampon;
   if (tampon=="INPUT_FORMAT"){
      IfficConf>>tampon;
      INPUT_FORMAT=c.from_s(tampon);
   }
   else{
      cout<<"From General.cpp (lireFicConfig):  ERROR in the configuration file at INPUT_FORMAT"<<endl;
      exit(EXIT_FAILURE);
   }
   IfficConf>>tampon;
   if (tampon=="OUTPUT_FORMAT"){
      IfficConf>>tampon;
      OUTPUT_FORMAT=c.from_s(tampon);
   }
   else{
      cout<<"From General.cpp (lireFicConfig):  ERROR in the configuration file at OUTPUT_FORMAT"<<endl;
      exit(EXIT_FAILURE);
   }
   
   //le séparateur
   IfficConf>>tampon;
   if (tampon=="sep"){
      IfficConf>>tampon;
      sep=tampon[0];
      }
   else{
      cout<<"From General.cpp (lireFicConfig):  ERROR in the configuration file at sep"<<endl;
      exit(EXIT_FAILURE);
   }

   // Création dossier pour les fichiers de sortie
   DIR *dir=opendir(directory.c_str());
   if (!dir){
      string commande="mkdir -p "+directory;
      int ret=system(commande.c_str());
      if(ret<0){
         cout<<"ERROR while creating directory "<<directory<<endl;
         exit(EXIT_FAILURE);
      }
      cout<<"OUTPUT directory: "<<directory<<endl;   
   }
   free(dir);
   
   //Les coûts
   bp=0;
   bp+=lireCout(IfficConf,Adj_percentage,"Adj_percentage");
   bp+=lireCout(IfficConf,Crea,"Gain");
   bp+=lireCout(IfficConf,Break,"Break");

   cout<<"Config file content :"<<endl;
   cout<<"Gene Trees file: "<<fic_arbre<<endl;
   cout<<"Species Tree file: "<<fic_especes<<endl;
   cout<<"Species-gene file: "<<fic_gene<<endl;
   cout<<"Prefixe: "<<prefixe<<endl;
   cout<<"Trees are they reconciled? "<<ArbresRec<<"   (0->No / 1->Yes)"<<endl;
   cout<<"Separator: "<<sep<<endl;
   cout<<"INPUT format: "<<INPUT_FORMAT<<"      (0->Newick Format / 1->NHX format)"<<endl;
   cout<<"OUTPUT format: "<<OUTPUT_FORMAT<<"   (0->Newick Format / 1->NHX format)"<<endl;
   cout<<"Adj_percentage: "<<Adj_percentage<<endl;
   cout<<"Gain cost: "<<Crea<<endl;
   cout<<"Break cost: "<<Break<<endl<<endl;

   if (bp!=3)
      exit(EXIT_FAILURE);

   bp=0;
   bp+=lireNomFic(IfficConf,file_spe,"genomes_file",argv);
   bp+=lireNomFic(IfficConf,ST_input,"species_tree_input",argv);
   bp+=lireNomFic(IfficConf,dir_embl,"embl_dir",argv);
   bp+=lireNomFic(IfficConf,dir_emf,"emf_dir",argv);
   bp+=lireCout(IfficConf,coeff_b,"coeff_b");
   IfficConf.close();

   
   if (bp==0){
      cout<<"No params to generate your data from the Ensembl database (Normal if you don't use Step0)"<<endl;
   }
   else if (bp==1 || bp==4 || bp==5){
      // OK nothing to do
   }
   else
      exit(EXIT_FAILURE);
}


// Step1 -> Step4 METHODS

void CreationOutputSpecies(TreeTemplate<Node> *S){
   string output_species_file=prefixe+output_spe;
   cout<<"Creation of OUTPUT file for association between ancestral species ID and extant species name: "<<output_species_file<<" ..."<<flush;
   if (file_log)
      Offile_log<<"FI\tFrom Step1 (CreationOutputSpecies): Creation of OUTPUT file for association between ancestral species ID and extant species name: "<<output_species_file<<" ..."<<flush;
   ofstream Offic_SORTIE_esp(output_species_file.c_str(), ios::out|ios::trunc);
   if (!Offic_SORTIE_esp){
      cout<<"\nFrom Step1 (CreationOutputSpecies):ERROR while opening file "<<output_species_file<<endl<<endl;
      if (file_log)
    Offile_log<<"\nFI\tFrom Step1 (CreationOutputSpecies):ERROR while opening file "<<output_species_file<<endl<<endl;
      exit(EXIT_FAILURE);
   }
   ecritEspeces(S->getRootNode(),Offic_SORTIE_esp);
   Offic_SORTIE_esp.close();
   cout<<" DONE"<<endl<<endl;
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
}

//Associate Gene with Species name
void AssociateGeneWithSpecies(map<string,string> &gene_species_EXT){
   string espece="";
   string gene="";
   cout<<"Creation of the map<string,string> gene_species_EXT from species-gene association file "<<fic_gene<<" ..."<<flush;
   if (file_log)
      Offile_log<<"FI\tFrom General.cpp (AssociateGeneWithSpecies): Creation of the map<string,string> gene_species_EXT from species-gene association file "<<fic_gene<<" ..."<<flush;
   ifstream IfficGene (fic_gene.c_str(), ios::in);
   if (!IfficGene) {
      cout<<"\nERROR while opening file "<<fic_gene<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom General.cpp (AssociateGeneWithSpecies): ERROR while opening file "<<adj_class_file<<endl;
      exit(EXIT_FAILURE);
   }
   while (!IfficGene.eof()){
      IfficGene>>espece;
      IfficGene>>gene;
      gene_species_EXT[gene]=espece;
   }   
   IfficGene.close();
   cout<<" DONE"<<endl<<endl;
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
}

// Step0(_emf) + Step1 + Step3(_proba)
//Rappatrie les noms des espèces aux noeuds + mets les distance à 1
//Utilisé sur l'arbre des espèces uniquement
void nommerNoeudsInternes(Node * n){
   convert<int> C;
   n->setNodeProperty(esp, BppString(C.to_s(n->getId())));
   if (n->isLeaf()){
      n->setDistanceToFather(1);
   }
   else{
      if (n->hasFather())
         n->setDistanceToFather(1);
      if (n->hasBranchProperty(esp)){
         BppString * Nom = dynamic_cast<BppString*> (n->getBranchProperty(esp));
         n->setName(Nom->toSTL());
      }
      else{
//         cout<<"\nFrom General.cpp (nommerNoeudsInternes): node "<<n->getId()<<" has no name !!!"<<endl;
//         if (file_log)
//            Offile_log<<"\nOS\tFrom General.cpp (nommerNoeudsInternes): node "<<n->getId()<<" has no name !!!"<<endl;
//         exit(EXIT_FAILURE);
      }
      nommerNoeudsInternes((n->getSons())[0]);
      nommerNoeudsInternes((n->getSons())[1]);
   }
}

//À n'utiliser que si arbre au format Newick a priori
void recupInfoSurBrancheNoeud(Node * n){
   convert<int> C;
   if (n->hasBranchProperty(esp)){
      BppString * PROP = dynamic_cast<BppString*> (n->getBranchProperty(esp)); 
      string prop=PROP->toSTL();

      //int Id =C.from_s(prop.substr(2, prop.find(sep)));

      prop=prop.substr(prop.find(sep)+1,prop.size());

      string name=prop.substr(0, prop.find(sep));
      n->setName(name);

      prop=prop.substr(prop.find(sep)+1,prop.size());

      string Espece=prop.substr(1, prop.find(sep)-1);
      n->setNodeProperty(esp,BppString(Espece));

      prop=prop.substr(prop.find(sep)+1,prop.size());

      string type=prop.substr(0, prop.find(sep));
      n->setBranchProperty(typ,BppString(type));
      //BppString * TYP = dynamic_cast<BppString*> (n->getBranchProperty(typ));

      prop=prop.substr(prop.find(sep)+1,prop.size());

      string D = prop.substr(1, prop.find(sep)-1);
      if (D!="NODIST"){
         int d=C.from_s(D);
         n->setDistanceToFather(d);
      }
   }
   else{
      cout<<"From General.cpp (recupInfoSurBrancheNoeud): WARNING no property on branches of node "<<n->getId()<<" !!!"<<endl;
      if (file_log)
         Offile_log<<"OS\tFrom General.cpp (recupInfoSurBrancheNoeud): WARNING no property on branches of node "<<n->getId()<<" !!!"<<endl;
      exit(EXIT_FAILURE);
   }
}

//À n'utiliser que si arbre au format Newick a priori
void recupInfoSurBrancheNoeuds(Node * n, TreeTemplate<Node> * S){
   if (n->isLeaf()){
      //Traitement spécial car seul son nom et son Id ont été écrit par le writer newick
      convert<int> C;
      string nom=n->getName();
      string gene=nom.substr(0,nom.find(sep));
      if (gene==NomPer)
         n->setBranchProperty(typ,BppString(per));
      else
         n->setBranchProperty(typ,BppString(ga));
      string nomEsp = nom.substr(nom.rfind(sep,nom.length()-1)+1, nom.length()-1);
      string Id = C.to_s((S->getNode(nomEsp))->getId());
      n->setNodeProperty(esp, BppString(Id));
   }
   else{
      recupInfoSurBrancheNoeuds((n->getSons())[0],S);
      recupInfoSurBrancheNoeud(n);
      recupInfoSurBrancheNoeuds((n->getSons())[1],S);
   }
}

//À n'utiliser que si arbre au format Newick a priori
void affecteInfoSurBrancheNoeud(Node * n){
   convert<int> C;
   string prop;

   prop+="Id";
   prop+=C.to_s(n->getId());
   prop+=sep;

   if (n->hasName())
      prop+=n->getName();
   else
      prop+="NONAME";
   prop+=sep;

   prop+="E";
   if (n->hasNodeProperty(esp)){
      BppString * espece = dynamic_cast<BppString*> (n->getNodeProperty(esp));
      prop+=espece->toSTL();
   }
   else
      prop+="NOSPECIE";
   prop+=sep;

   if (n->hasBranchProperty(typ)){
      BppString * type = dynamic_cast<BppString*> (n->getBranchProperty(typ));
      prop+=type->toSTL();
   }
   else
      prop+="NOTYPE";
   prop+=sep;

   prop+="D";
   if (n->hasDistanceToFather())
      prop+=C.to_s(n->getDistanceToFather());
   else
      prop+="NODIST";
   n->setBranchProperty(esp,BppString(prop));
}

//À n'utiliser que si arbre au format Newick a priori
//cas 1 fils pour les arbres d'adjacences
void affecteInfoSurBrancheNoeuds(Node * n){
   if (n->getNumberOfSons()==0)//isLeaf())
      affecteInfoSurBrancheNoeud(n);
   else if (n->getNumberOfSons()==1){
      affecteInfoSurBrancheNoeuds((n->getSons())[0]);
      affecteInfoSurBrancheNoeud(n);
   }
   else{
      affecteInfoSurBrancheNoeuds((n->getSons())[0]);
      affecteInfoSurBrancheNoeud(n);
      affecteInfoSurBrancheNoeuds((n->getSons())[1]);
   }
}

//Récupération & affectation d'info sur les branches des arbres de gènes avec l'arbre des espèces et stockage arbres de gènes dans un vector<Tree *> Arbres
void StoreGeneTrees(TreeTemplate<Node> * S, vector<Tree *> &Arbres){
   vector<Tree *>::iterator it;
   cout<<"Store Gene trees in vector<Tree *> Arbres..."<<flush;
   if (file_log){
      Offile_log<<"OS\tFrom General.cpp (StoreGeneTrees): Store Gene trees in vector<Tree *> Arbres..."<<flush;
   }
   if(INPUT_FORMAT==0){
      Newick * newickReader = new Newick(true,true);
      newickReader->enableExtendedBootstrapProperty(esp);
      newickReader->read(fic_arbre.c_str(),Arbres);
      delete(newickReader);
      if (ArbresRec){
         for(it=Arbres.begin();it!=Arbres.end();it++){
            Tree * Gp=*it;
            TreeTemplate<Node> *G = dynamic_cast <TreeTemplate<Node> *> (Gp);
            recupInfoSurBrancheNoeuds(G->getRootNode(),S);
            //Pour remettre à jour d'identifiant
            affecteInfoSurBrancheNoeuds(G->getRootNode());
         }
      }
   }
   else if(INPUT_FORMAT==1){
      Nhx * nhxReader = new Nhx(true);
      nhxReader->registerProperty(Nhx::Property("Duplication", "D", true, 0));
      nhxReader->read(fic_arbre.c_str(),Arbres);
      delete(nhxReader);
   }
   else{
      cout<<"\tFrom General.cpp (StoreGeneTrees): ERROR wrong INPUT_FORMAT in config file: "<<INPUT_FORMAT<<" (Should be 0 or 1)"<<endl;
      if (file_log){
         Offile_log<<"\tFrom General.cpp (StoreGeneTrees): ERROR wrong INPUT_FORMAT in config file: "<<INPUT_FORMAT<<" (Should be 0 or 1)"<<endl;
      }
      exit(EXIT_FAILURE);
   }
   cout<<" DONE"<<endl<<endl;
   if (file_log){
      Offile_log<<" DONE"<<endl<<endl;
   }
}


//   METHODS Step0(_emf) + Step1   //

void retourneFeuilles(Node * n, vector<string> & feuilles){
   if (n->getNumberOfSons()==0)
      feuilles.push_back(n->getName());
   else if (n->getNumberOfSons()==1)
      retourneFeuilles(n->getSon(0),feuilles);
   else{ //2 fils
      retourneFeuilles(n->getSon(0),feuilles);
      retourneFeuilles(n->getSon(1),feuilles);
   }
}

void ecritEspeces(Node * n, ofstream & Offic_SORTIE_esp){
   //Le numéro de l'espèce 
   Offic_SORTIE_esp<<n->getId();
   //Les descendants de cette espèce
   vector<string> feuilles;
   retourneFeuilles(n,feuilles);
   for(unsigned int i=0;i<feuilles.size();i++)
      Offic_SORTIE_esp<<"\t"<<feuilles[i];
   
   Offic_SORTIE_esp<<"\n";

   if (n->getNumberOfSons()==1)
      ecritEspeces(n->getSon(0),Offic_SORTIE_esp);
   else if (n->getNumberOfSons()==2){
      ecritEspeces(n->getSon(0),Offic_SORTIE_esp);
      ecritEspeces(n->getSon(1),Offic_SORTIE_esp);
   }
}

//Pour ajouter les noms des espèces aux feuilles des arbres de gènes si nécessaire
void reetiquette(Node * n, map<string,string> &gene_species_EXT){
   static const size_t npos = -1;
   if (n->getNumberOfSons()==0){
      if (n->getName().find(sep)==npos){
         string species_name="";
         if (gene_species_EXT.find(n->getName())!=gene_species_EXT.end()){
            species_name=gene_species_EXT[n->getName()];
         }
         string nom=n->getName()+sep+species_name;
         n->setName(nom);
      }
   }
   else{
      reetiquette(n->getSon(0),gene_species_EXT);
      reetiquette(n->getSon(1),gene_species_EXT);
   }
}

//VÉRIFIER TOUS LES GETNODE AU CAS OÙ UNE ESPÈCE DE L'ARBRE DE GÈNE NE
//SOIT PAS DANS L'ARBRES DES ESPÈCES --> ON ÉLAGUE L'ARBRE DES GÈNES
Node * preAffecteEspece(Node * n, TreeTemplate<Node> * G, TreeTemplate<Node> * S){
   //n est un noeud de G
   convert<int> C;
   if (n->getNumberOfSons()==0){
      string nomGene = n->getName();
      string nomEsp = nomGene.substr(nomGene.rfind(sep,nomGene.length()-1)+1, nomGene.length()-1);
      string Id = ""; 
      bool feuille=false;
      try{
         Id = C.to_s( S->getNode(nomEsp)->getId());
         if(S->getNode(nomEsp)->isLeaf()){
            feuille=true;
            n->setNodeProperty(esp, BppString("CONNUE"));
         }
      }
      catch (Exception e) {
         //perror(e.what()); 
      }
      if (Id=="" || !feuille){
         n->setNodeProperty(esp, BppString("INCONNUE"));
      }
      return n;
   }
   else{
      Node * fg = n->getSon(0);
      Node * fd = n->getSon(1);
      fg=preAffecteEspece(fg,G,S);
      fd=preAffecteEspece(fd,G,S); 
    
   //Après ces pre affectation n peut avoir perdu ses fils !!! (forcément les deux)
      if (!n->isLeaf()){
         BppString * EspeceFG = dynamic_cast<BppString*> (fg->getNodeProperty(esp));
         BppString * EspeceFD = dynamic_cast<BppString*> (fd->getNodeProperty(esp));

         if((EspeceFG->toSTL()=="INCONNUE") && (EspeceFD->toSTL()=="INCONNUE")){
            //cout<<"Les deux fils de "<<n->getId()<<" sont d'espèces INCONNUES"<<endl;
            n->setNodeProperty(esp, BppString("INCONNUE"));
            n->removeSons();
         }
         else if (EspeceFG->toSTL()=="INCONNUE"){
            //cout<<"Le fils gauche de "<<n->getId()<<" est d'espèce INCONNUE"<<endl;
            if (n->hasFather()){
               if(n->getFather()->getSon(0)->getId()==n->getId())
                  n->getFather()->setSon(0,fd);
               else
                  n->getFather()->setSon(1,fd);
               n=fd;
            }   
            else
               G->setRootNode(fd);
            //cout<<"Le fils gauche de "<<n->getId()<<" est d'espèce INCONNUE ----> TRAITÉ"<<endl;
         }
         else if (EspeceFD->toSTL()=="INCONNUE"){
         //cout<<"Le fils droit de "<<n->getId()<<" est d'espèce INCONNUE"<<endl;
               if (n->hasFather()){
                  if(n->getFather()->getSon(0)->getId()==n->getId())
                     n->getFather()->setSon(0,fg);
                  else
                     n->getFather()->setSon(1,fg);
                     n=fg;
               }
               else
                  G->setRootNode(fg);
               //cout<<"Le fils droit de "<<n->getId()<<" est d'espèce INCONNUE ----> TRAITÉ"<<endl;
         }
         else{
            //cout<<"Les deux fils de "<<n->getId()<<" sont d'espèces CONNUES"<<endl;
            n->setNodeProperty(esp, BppString("CONNUE"));
         }
      }
   }
   return n;
}

//Fonctionne pour des arbres binaires dont les feuilles ont été
//étiquetées. On affecte les Id des noeuds de l'arbre des espèces (au
//cas où on ne connaisse pas les noms)
void affecteEspece(Node * n, TreeTemplate<Node> * S){
   //n est un noeud de G
   convert<int> C;
   if (n->isLeaf()){
      string nomGene = n->getName();
      string nomEsp = nomGene.substr(nomGene.rfind(sep,nomGene.length()-1)+1, nomGene.length()-1);
      string Id = C.to_s( S->getNode(nomEsp)->getId());
      BppString * prop = new BppString(Id);
      n->setNodeProperty(esp, *prop);
      delete(prop);
   }
   else{
      Node * fg = (n->getSons())[0];
      Node * fd = (n->getSons())[1];

      affecteEspece(fg,S); 
      affecteEspece(fd,S);  
      BppString * EspeceFG = dynamic_cast<BppString*> (fg->getNodeProperty(esp));
      BppString * EspeceFD = dynamic_cast<BppString*> (fd->getNodeProperty(esp));

      vector< int > IdFils;
      IdFils.push_back(C.from_s(EspeceFG->toSTL()));
      IdFils.push_back(C.from_s(EspeceFD->toSTL()));
      int IdLCA = TreeTools::getLastCommonAncestor(* S, IdFils);
      n->setNodeProperty(esp, BppString(C.to_s(IdLCA)));
   }
}

//Apparemment il existe une fonction bio++ qui fait ça tout seul
void renumeroteId(Node * n, int & id){
   if (n->isLeaf()){
      n->setId(id);
      id++;
   }
   else{
      renumeroteId(n->getSon(0),id);
      renumeroteId(n->getSon(1),id);
      n->setId(id);
      id++;
   }
}

//Devrait fonctionner pour les deux formats
void annoteType(Node * n, TreeTemplate<Node> * S){
   convert<int> c;
   if (n->isLeaf()){
      n->setBranchProperty(typ, BppString(ga));
      n->setBranchProperty("D", BppString("?"));
   }
   else if (n->getNumberOfSons()!=2){
      cout<<"From General.cpp (annoteType): NON BINARY gene tree !!!"<<endl;
      if (file_log)
         Offile_log<<"OS\tFrom General.cpp (annoteType): NON BINARY gene tree !!!"<<endl;
      exit(EXIT_FAILURE);
   }
   else{
      Node * fg = (n->getSons())[0];
      Node * fd = (n->getSons())[1];
      annoteType(fg,S); 
      annoteType(fd,S);  
      BppString * EspeceN = dynamic_cast<BppString*> (n->getNodeProperty(esp));
      BppString * EspeceFG = dynamic_cast<BppString*> (fg->getNodeProperty(esp));
      BppString * EspeceFD = dynamic_cast<BppString*> (fd->getNodeProperty(esp));
    
      if (c.from_s(EspeceN->toSTL()) != c.from_s(EspeceFD->toSTL()) && c.from_s(EspeceN->toSTL()) != c.from_s(EspeceFG->toSTL())){
         n->setBranchProperty(typ, BppString(spe));
         n->setBranchProperty("D", BppString("F"));
      }
      else{
         n->setBranchProperty(typ, BppString(dupl));
         n->setBranchProperty("D", BppString("T"));
      }
   }
}

//Renvoie un noeud dont les fils sont f_G et p_G (une perte) 
//d'espèces correspondant aux fils du noeud x_S
Node * ajouteNoeudPerte(Node * x_S,Node * f_G,  int &id,TreeTemplate<Node> * S){
   convert<int> C;
   Node * p_G = new Node(); //futur noeud perte à ajouter
   p_G->setDistanceToFather(1);
   p_G->setId(id);id++;
   p_G->setBranchProperty(typ,BppString(per));
   p_G->setBranchProperty("D", BppString("?"));
   
   BppString * EspeceFils1 = dynamic_cast<BppString*> (x_S->getSon(1)->getNodeProperty(esp));
   BppString * EspeceF_G = dynamic_cast <BppString*> (f_G->getNodeProperty(esp));
   BppString * EspeceFils0 = dynamic_cast<BppString*> (x_S->getSon(0)->getNodeProperty(esp));
   if (EspeceFils0->toSTL()==EspeceF_G->toSTL()){
      p_G->setNodeProperty(esp,*EspeceFils1); 
      string NomPerte = NomPer;
      NomPerte+=sep;
      if (S->getNode(C.from_s(EspeceFils1->toSTL()))->hasName())
         NomPerte+=S->getNode(C.from_s(EspeceFils1->toSTL()))->getName();
      else
         NomPerte+=C.to_s(S->getNode(C.from_s(EspeceFils1->toSTL()))->getId());
      p_G->setName(NomPerte);
   }
   else{
      p_G->setNodeProperty(esp,*EspeceFils0);    
      string NomPerte = NomPer;
      NomPerte+=sep;
      if (S->getNode(C.from_s(EspeceFils0->toSTL()))->hasName())
         NomPerte+=S->getNode(C.from_s(EspeceFils0->toSTL()))->getName();
      else
         NomPerte+=C.to_s(S->getNode(C.from_s(EspeceFils0->toSTL()))->getId());
      p_G->setName(NomPerte);
   }     
   
   Node * t_G = new Node(); //futur noeud père de f_G et de la perte
   t_G->setDistanceToFather(1);
   
   t_G->setId(id);id++;
   //j'ajouterai l'id de T_G après avoir raccordé le tout sinon G->getNextId() me donne le même res que précédemment
   t_G->setBranchProperty(typ,BppString(spe));
   t_G->setBranchProperty("D", BppString("F"));
   t_G->setNodeProperty(esp,*(x_S->getNodeProperty(esp))); 
   
   t_G->addSon(f_G);
   t_G->addSon(p_G);
   
   return t_G;
}

void ajoutePerte(Node * n, int &id, TreeTemplate<Node> * S){
   convert<int> C;
   if (n->isLeaf()){
      //On ne fait rien
   }
   else if (n->getNumberOfSons()!=2){
      cout<<"From General.cpp (ajoutePerte): NON BINARY gene tree !!!"<<endl;
      if (file_log)
         Offile_log<<"OS\tFrom General.cpp (ajoutePerte): NON BINARY gene tree !!!"<<endl;
      exit(EXIT_FAILURE);
   }
   else{
      Node * fg = n->getSon(0);
      Node * fd = n->getSon(1);
      ajoutePerte(fg,id,S);//Normalement id est modifiée à l'issue de cet appel
      ajoutePerte(fd,id,S);

      BppString * EspeceN = dynamic_cast<BppString*> (n->getNodeProperty(esp));
      Node * NoeudEspN = S->getNode(C.from_s(EspeceN->toSTL()));
      BppString * typeN = dynamic_cast<BppString*> (n->getBranchProperty(typ));

      //Il faut traiter 2 arêtes : N-FG et N-FD  
      BppString * EspeceFG = dynamic_cast<BppString*> (fg->getNodeProperty(esp));
      Node * NoeudEspFG = S->getNode(C.from_s(EspeceFG->toSTL()));
      BppString * EspeceFD = dynamic_cast<BppString*> (fd->getNodeProperty(esp));
      Node * NoeudEspFD = S->getNode(C.from_s(EspeceFD->toSTL()));

      //VÉRIF
      vector< int > IdFils;
      IdFils.push_back(S->getNode(C.from_s(EspeceFG->toSTL()))->getId());
      IdFils.push_back(S->getNode(C.from_s(EspeceFD->toSTL()))->getId());
      IdFils.push_back(S->getNode(C.from_s(EspeceN->toSTL()))->getId());
      if (S->getNode(C.from_s(EspeceN->toSTL()))->getId()!=TreeTools::getLastCommonAncestor(* S, IdFils)){
         cout<<"From General.cpp (ajoutePerte): ERROR "<<n->getId()<<" is not an ancestor of its sons in the species tree"<<endl;
         if (file_log)
            Offile_log<<"OS\tFrom General.cpp (ajoutePerte): ERROR "<<n->getId()<<" is not an ancestor of its sons in the species tree"<<endl;
         exit(EXIT_FAILURE);
      }
    
      //
      //N-FG 
      //    
      if(EspeceFG!=EspeceN){
         //toutes les distances sont à 1 après appel de nommerNoeudsInternes
         int nbAretes=TreeTools::getDistanceBetweenAnyTwoNodes(* S, NoeudEspN->getId(),NoeudEspFG->getId());

         //Ajout d'un noeud supplémentaire en cas de duplication
         if(typeN->toSTL()==dupl)
            nbAretes++;

         Node * x_S= NoeudEspFG->getFather();//x_S parcours la branche dans l'arbre des espèces
         for(int i=0; i<nbAretes-1;i++){
            //cout<<"PERTES DÉTECTÉES !!!"<<endl;
            fg=ajouteNoeudPerte(x_S,fg,id,S); //fg modifié suite à cet appel fg=t_G;
            x_S=x_S->getFather();
         }
         //"on raccorde"
         n->setSon(0,fg); //Car t_G n'est pas connu en dehors
         //des boucles alors que fg oui
      }

      //
      //N-FD 
      //
      if(EspeceFD!=EspeceN){
         //toutes les distances sont à 1 après appel de nommerNoeudsInternes
         int nbAretes=TreeTools::getDistanceBetweenAnyTwoNodes(* S, NoeudEspN->getId(),NoeudEspFD->getId());

         //Ajout d'un noeud supplémentaire en cas de duplication
         if(typeN->toSTL()==dupl){
            nbAretes++;
         }
         Node * x_S= NoeudEspFD->getFather();//x_S parcours la branche dans l'arbre des espèces
         for(int i=0; i<nbAretes-1;i++){
            //cout<<"PERTES DÉTECTÉES !!!"<<endl;
            fd=ajouteNoeudPerte(x_S,fd,id,S); //fd modifié suite à cet appel fd=t_G;
            x_S=x_S->getFather();
         }
         //"on raccorde"
         n->setSon(1,fd); //Car t_G n'est pas connu en dehors
         //des boucle alors que fd oui
      }
   }
}

void reconcil(TreeTemplate<Node > * S,TreeTemplate<Node >* G){
   //
   // 1/ À chaque noeud de l'arbre de gène, on attribue un noeud de l'abre des espèces. 
   //
   preAffecteEspece(G->getRootNode(), G, S);   
   if (G->getRootNode()->getNumberOfSons()>0){
      //Suite à l'éventuel élagage de G on renumérote tout pour éviter
      //d'avoir des trous dans la numérotation des identifiants
      int id=0;
      renumeroteId(G->getRootNode(),id);
      affecteEspece(G->getRootNode(),S);
      //
      // 2/ On annote les duplications et spéciations 
      //
      annoteType(G->getRootNode(),S);
      //
      // 3/ On ajoute les pertes 
      //
      id=G->getNumberOfNodes();
      ajoutePerte(G->getRootNode(),id,S); 
   }
}

void reconciliation(TreeTemplate<Node > *S, vector<Tree *> &Arbres, map<string,string> &gene_species_EXT){
   int nb_arbres=Arbres.size();
   if (nb_arbres==0){
      cout<<"From General.cpp (reconciliation): No trees in vector<Tree *> Arbres... Stopping the program"<<endl;
      if (file_log)
         Offile_log<<"OS\tFrom General.cpp (reconciliation): No trees in vector<Tree *> Arbres... Stopping the program"<<endl;
      exit(EXIT_FAILURE);
   }
   cout<<"The "<<nb_arbres<<" trees of file "<<fic_arbre<<" are loaded"<<endl;
   cout<<"Reconciliation of Gene trees with Species tree... "<<flush;
   if (file_log){
      Offile_log<<"OS\tFrom General.cpp (reconciliation): The "<<nb_arbres<<" trees of file "<<fic_arbre<<" are loaded"<<endl;
      Offile_log<<"OS\tFrom General.cpp (reconciliation): Reconciliation of Gene trees with Species tree... "<<flush;
   }
   int i=1;
   if (!ArbresRec){
      vector<Tree *>::iterator it;
      for(it=Arbres.begin();it!=Arbres.end();it++){
//         if (file_log)
//            Offile_log<<"OS\tFrom General.cpp (reconciliation): Arbres réconciliés "<<i<<endl;
         Tree * Gp=*it;
         TreeTemplate<Node> *G = dynamic_cast <TreeTemplate<Node> *> (Gp);
         G->setBranchLengths(1);
         if (!TreeTemplateTools::isMultifurcating(*(G->getRootNode()))){
            // if (file_log)
            //    Offile_log<<"OS\tFrom General.cpp (reconciliation): Before reetiquette"<<endl;
            reetiquette(G->getRootNode(),gene_species_EXT);
            // if (file_log)
            //    Offile_log<<"OS\tFrom General.cpp (reconciliation): After reetiquette & Before reconcil"<<endl;
            reconcil(S,G);
            //Renumérote les id des noeuds pour qu'il n'y ait
            //pas de différence pour un même jeu de données
            //AVEC ou SANS étape de réconciliation :
            // if (file_log)
            //    Offile_log<<"OS\tFrom General.cpp (reconciliation): After reconcil & Before G->resetNodesId()"<<endl;
            G->resetNodesId();
            // if (file_log)
            //    Offile_log<<"OS\tFrom General.cpp (reconciliation): After G->resetNodesId() & Before affecteInfoSurBracheNoeuds(G->getRootNode())"<<endl;
            affecteInfoSurBrancheNoeuds(G->getRootNode());
         }
         if (TreeTemplateTools::isMultifurcating(*(G->getRootNode())) || G->getRootNode()->getNumberOfSons()==0){//aucune espèce dans notre arbre d'espèces
//            if (file_log)
//               Offile_log<<"OS\tFrom General.cpp (reconciliation): Tree "<<i<<" should be erased"<<endl;
            it=Arbres.erase(it);
            it--;
            //if(TreeTemplateTools::isMultifurcating(*(G->getRootNode())))
            // cout<<"\t Tree "<<i<<" treated, NON BINARY !!! => ERASED."<<endl;
            //else
            // cout<<"\t Tree "<<i<<" treated and ERASED."<<endl;
         }
         //else{
         //   cout<<"\t Tree "<<i<<" treated."<<endl;
         //}
         i++;
      }
      nb_arbres=Arbres.size();
      if (nb_arbres==0){
         cout<<"\nFrom General.cpp (reconciliation): ERROR No trees kept... stopping the program"<<endl;
         if (file_log)
            Offile_log<<"\nOS\tFrom General.cpp (reconciliation): ERROR No trees kept... stopping the program"<<endl;
         exit(EXIT_FAILURE);
      }
      cout<<" DONE: "<<nb_arbres<<" trees kept."<<endl<<endl;
      cout<<"Write reconciled trees in file: "<<tree_reconciled_file<<" ..."<<flush;
      if (file_log){
         Offile_log<<" DONE: "<<nb_arbres<<" trees kept."<<endl;
         Offile_log<<"OS\tFrom General.cpp (reconciliation): Write reconciled trees in file: "<<tree_reconciled_file<<" ..."<<flush;
      }
      if(OUTPUT_FORMAT==0){
         Newick * newickWriter = new Newick(true,true); 
         newickWriter->enableExtendedBootstrapProperty(esp);
         newickWriter->write(Arbres,tree_reconciled_file.c_str());
         delete(newickWriter);
      }
      else if(OUTPUT_FORMAT==1){
         if (file_log)
            Offile_log<<"OS\tFrom General.cpp (reconciliation): OUTPUT_FORMAT==1"<<endl;
         Nhx * nhxWriter = new Nhx(true);
         if (file_log)
            Offile_log<<"OS\tFrom General.cpp (reconciliation): After Nhx * nhxWriter = new Nhx(true); Before nhxWriter->registerProperty(Nhx::Property(\"Duplication\", \"D\", true, 0));)"<<endl;
         nhxWriter->registerProperty(Nhx::Property("Duplication", "D", true, 0)); 
         if (file_log)
            Offile_log<<"OS\tFrom General.cpp (reconciliation): After nhxWriter->registerProperty(Nhx::Property(\"Duplication\", \"D\", true, 0)); & Before nhxWriter->write(Arbres,tree_reconciled_file.c_str());"<<endl;
         nhxWriter->write(Arbres,tree_reconciled_file.c_str());
         if (file_log)
            Offile_log<<"OS\tFrom General.cpp (reconciliation): After nhxWriter->write(Arbres,tree_reconciled_file.c_str());"<<endl;
         delete(nhxWriter);
      }
      else{
         cout<<"From General.cpp (reconciliation): ERROR wrong OUTPUT_FORMAT"<<endl;
         if (file_log)
            Offile_log<<"OS\tFrom General.cpp (reconciliation): ERROR wrong OUTPUT_FORMAT"<<endl;
         exit(EXIT_FAILURE);
      }
      cout<<" DONE"<<endl<<endl;
      if (file_log)
         Offile_log<<" DONE"<<endl<<endl;
   }
   else{
      cout<<"Trees are already reconciled (boolean ReconcilDone is true in the configuration file)."<<endl;
      if (file_log)
         Offile_log<<"OS\tFrom General.cpp (reconciliation): Trees are already reconciled (boolean ReconcilDone is true in the configuration file)."<<endl;
   }
}

// Store Species tree in TreeTemplate<Node> * S (Step0 -> Step4)
void StoreSpeciesTree(TreeTemplate<Node> * S){
//   S = new TreeTemplate<Node>;
   cout<<"Store Species tree in TreeTemplate<Node> * S..."<<flush;
   if(file_log)
      Offile_log<<"OS\tFrom General.cpp (StoreSpeciesTree): Store Species tree in TreeTemplate<Node> * S..."<<flush;
   Newick * newickReaderEsp = new Newick(true,true);
   newickReaderEsp->enableExtendedBootstrapProperty(esp);
   S = newickReaderEsp->read(fic_especes);
   delete(newickReaderEsp);
   cout<<" DONE"<<endl<<endl;
   if(file_log)
      Offile_log<<" DONE"<<endl<<endl;
}

void StoreReconciledTrees(vector<Tree *> &Arbres){
   cout<<"Store reconciled trees in vector<Tree *> Arbres..."<<flush;
   if (file_log){
      Offile_log<<"OS\tFrom General.cpp (StoreReconciledTrees): Store reconciled trees in vector<Tree> Arbres..."<<flush;
   }
   
   string file_to_load;
   
   if (is_readable(tree_reconciled_file))
      file_to_load=tree_reconciled_file;
   else if (ArbresRec && is_readable(fic_arbre))
      file_to_load=fic_arbre;   
   else{
      cout<<"\nFrom General.cpp (StoreReconciledTrees): ERROR reconciled trees file "<<tree_reconciled_file<<" is not readable and bool ArbresRec is false"<<endl;
      if (file_log)
         Offile_log<<"\nFrom General.cpp (StoreReconciledTrees): ERROR reconciled trees file "<<tree_reconciled_file<<" is not readable and bool ArbresRec is false"<<endl;
      exit(EXIT_FAILURE);
   }    
   
   cout<<" from file "<<file_to_load<<flush;
   if (file_log){
      Offile_log<<" from file "<<file_to_load<<flush;
   }

   if(OUTPUT_FORMAT==0){
      Newick * newickReader = new Newick(true,true); 
      newickReader->enableExtendedBootstrapProperty(esp);
      newickReader->read(file_to_load.c_str(),Arbres);
      delete(newickReader);
   }
   else if(OUTPUT_FORMAT==1){
      Nhx * nhxReader = new Nhx(true);
      nhxReader->registerProperty(Nhx::Property("Duplication", "D", true, 0));
      nhxReader->read(file_to_load.c_str(),Arbres);
      delete(nhxReader);
   }
   else{
      cout<<"\nFrom General.cpp (StoreReconciledTrees): ERROR wrong OUTPUT_FORMAT"<<endl;
      if (file_log)
    Offile_log<<"\nOS\tFrom General.cpp (StoreReconciledTrees): ERROR wrong OUTPUT_FORMAT"<<endl;
      exit(EXIT_FAILURE);
   }
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
   cout<<" DONE"<<endl<<endl;
}

void AssociateGeneWithReconciledTreeNb(vector<Tree *> &Arbres, map<string,vector<int> > &gene_GeneTreeID, map<string,string> &gene_species_EXT){
   cout<<"Associate Genes with reconciled trees Id..."<<flush;
   if (file_log){
      Offile_log<<"OS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): Associate Genes with reconciled trees Id..."<<flush;
   }
   int i=0;
   string nom_gene="";
   set<string> GENES;
   vector<string> noms_feuilles;
   vector<string>::iterator it_feuilles;
   vector<Tree *>::iterator it;
   for(it=Arbres.begin();it!=Arbres.end();it++){
      //      cout<<"\nCourse tree "<<i<<"\n\tBefore getLeavesName"<<endl;
      //      if (file_log)
      //         Offile_log<<"\nOS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): Course tree "<<i<<"\nBefore getLeavesName"<<endl;
      noms_feuilles=(*it)->getLeavesNames();
//      cout<<"\tAfter getLeavesName"<<endl;
//      if (file_log)
//         Offile_log<<"\nOS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): After getLeavesName"<<endl;
      for(it_feuilles=noms_feuilles.begin();it_feuilles!=noms_feuilles.end();it_feuilles++){
         nom_gene=(*it_feuilles).substr(0,(*it_feuilles).rfind(sep,(*it_feuilles).length()-1));
//         cout<<"Analyse gene: "<<nom_gene<<endl;
//         if (file_log)
//            Offile_log<<"\nOS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): Analyse gene: "<<nom_gene<<endl;
         if(nom_gene!=NomPer){
            //Je relève les noms de feuilles :
//            cout<<"Before GENES insert"<<endl;
//            if (file_log)
//               Offile_log<<"\nOS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): Before GENES insert"<<endl;
            GENES.insert(nom_gene);
//            cout<<"After GENES insert"<<endl;
//            if (file_log)
//               Offile_log<<"\nOS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): After GENES insert"<<endl;
            //vérif que ces noms de feuilles est bien dans gene_species_EXT
            if (gene_species_EXT.find(nom_gene)!=gene_species_EXT.end()){
//               cout<<nom_gene<<" is WELL in gene_species_EXT"<<endl;
//               if (file_log)
//                  Offile_log<<"OS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): "<<nom_gene<<" is WELL in gene_species_EXT"<<endl;
               if (gene_GeneTreeID.find(nom_gene)!=gene_GeneTreeID.end()){
                  cout<<"\nWARNING "<<nom_gene<<" appears in multiple gene trees"<<endl;
                  if (file_log)
                     Offile_log<<"\nOS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): WARNING "<<nom_gene<<" appears in multiple gene trees"<<endl;
                  gene_GeneTreeID[nom_gene].push_back(i);
               }
               else{
//                  cout<<nom_gene<<" is NOT in gene_GeneTreeID"<<endl;
//                  if (file_log)
//                     Offile_log<<"OS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): "<<nom_gene<<" is NOT in gene_GeneTreeID"<<endl;
                  vector<int> v;
                  v.push_back(i);
                  gene_GeneTreeID[nom_gene]=v;
               }
            }
            else{
               cout<<"\nERROR "<<nom_gene<<" is NOT a known gene (=> Present in gene trees but not in Species-Gene association file)"<<endl;
               if (file_log)
                  Offile_log<<"\nOS\tFrom General.cpp (AssociateGeneWithReconciledTreeNb): WARNING "<<nom_gene<<" is NOT a known gene (=> Present in gene trees but not in Species-Gene association file)"<<endl;
               exit(EXIT_FAILURE);
            }
         }
      }
      i++;
   }
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
   cout<<" DONE"<<endl<<endl;
}

int GeneSpeciesClean(map<string,vector<int> > &gene_GeneTreeID, map<string,string> &gene_species_EXT){
   cout<<"Cleaning map gene_species_EXT of size "<<gene_species_EXT.size()<<" with map gene_GeneTreeID of size "<<gene_GeneTreeID.size()<<endl;
   cout<<"Remove genes that are not in Gene trees..."<<flush;
   if (file_log){
      Offile_log<<"OS\tFrom General.cpp (GeneSpeciesClean): Remove genes that are not in Gene trees..."<<flush;
   }
   map<string,string>::iterator it_ge;
   map<string,string>::iterator it_ge_prec;
   int elim=0;
   int gard=0;
   string gene="";
   map<string,string> GARDES;

   for(it_ge=gene_species_EXT.begin();it_ge!=gene_species_EXT.end();it_ge++){
      if (gene_GeneTreeID.find((*it_ge).first)==gene_GeneTreeID.end()){
         elim++;
      }
      else{
         GARDES[(*it_ge).first]=(*it_ge).second;
         gard++;
      }
   }
   gene_species_EXT.clear();
   gene_species_EXT=GARDES;
   //   cout<<"GARDES.size()="<<GARDES.size()<<" et maintenant gene_species_EXT.size()="<<gene_species_EXT.size()<<endl;
   GARDES.clear();
   cout<<" DONE"<<endl;
   cout<<elim<<" genes have been eliminated from the gene_species_EXT file cause not present in the reconciled trees"<<endl;
   cout<<gard<<" genes have been kept"<<endl<<endl;
   if (file_log){
      Offile_log<<" DONE"<<endl;
      Offile_log<<"OS\tFrom General.cpp (GeneSpeciesClean): "<<elim<<" genes have been eliminated from the gene_species_EXT file cause not present in the reconciled trees"<<endl;
      Offile_log<<"OS\tFrom General.cpp (GeneSpeciesClean): "<<gard<<" genes have been kept"<<endl<<endl;
   }
   return elim;
}



// Step3(_emf) -> Step4 METHODS

//Pour compter le nombre de noeud ayant la propriété prop à la valeur val et en stockant les résultats dans la map especes_nb_genes
int countNodesWithBranchPropertyBySpecies(Node * n, string prop, string val, map<int,int>& especes_nb_genes){
   convert<int> c;
   int p=0;
   if (n->hasBranchProperty(prop)){
      BppString * VAL = dynamic_cast<BppString*> (n->getBranchProperty(prop)); 
      if(VAL->toSTL()==val){
         p=1;
         //On rempli la map especes_nb_genes si le noeud est associé à une espèce
         if (n->hasNodeProperty(esp)){
            BppString * ESP = dynamic_cast<BppString*> (n->getNodeProperty(esp));
            int num_esp=c.from_s(ESP->toSTL()); 
            if(especes_nb_genes.find(num_esp)!=especes_nb_genes.end())
               especes_nb_genes[num_esp]++;
            else
               especes_nb_genes[num_esp]=1;
         }
      }
   }
   if (n->getNumberOfSons()==0)
      return p;
   else if (n->getNumberOfSons()==1)
      return p+countNodesWithBranchPropertyBySpecies(n->getSon(0),prop,val,especes_nb_genes);
   else //2 fils 
      return p
   +countNodesWithBranchPropertyBySpecies(n->getSon(0),prop,val,especes_nb_genes)
   +countNodesWithBranchPropertyBySpecies(n->getSon(1),prop,val,especes_nb_genes);  
}

void AdjClassVector(map<string,vector<int> > &classes_adjacences){
   cout<<"Creation of the map<string,vector<int> > classes_adjacences from adjacencies classes file "<<adj_class_file<<" ..."<<flush;
   if (file_log)
      Offile_log<<"FI\tFrom General.cpp (AdjClassVector): Creation of the map<string,vector<int> > classes_adjacences from adjacencies classes file "<<adj_class_file<<" ..."<<flush;
   ifstream IffilAdjClass (adj_class_file.c_str(), ios::in);
   if (!IffilAdjClass){
      cout<<"\nFrom General.cpp (AdjClassVector): ERROR while opening file "<<adj_class_file<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom General.cpp (AdjClassVector): ERROR while opening file "<<adj_class_file<<endl;
      exit(EXIT_FAILURE);
   }
   string line,buffer,name_class;
   int size_class,adj;
   while (getline(IffilAdjClass, line)){
      istringstream buffer(line);
      buffer>>name_class;
      buffer>>size_class;
//      cout<<name_class<<" "<<size_class;
      while (!buffer.eof()){
         buffer>>adj;
         classes_adjacences[name_class].push_back(adj);
//         cout<<" "<<adj;            
      }         
//      cout<<endl;
   }
   IffilAdjClass.close();
   cout<<" DONE"<<endl<<endl;
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
}

//Création de la map <int,vector<adjacence>> Adj_actuelles <Id_esp,adjacences de l'esp> !!! Ne contient pas les gènes qui n'ont pas d'adj !!!
void AdjVector_Esp_AdjExtant(TreeTemplate<Node> * S, map<string,string> &gene_species_EXT, vector<adjacence> &adjacencies, map<int,vector<adjacence> > &Adj_actuelles){
   cout<<"Creation of the sorted vector<adjacence> adjacences from adjacencies file "<<fic_adjacence<<" ..."<<flush;
   if (file_log)
      Offile_log<<"FI\tFrom General.cpp (AdjVector_Esp_AdjExtant): Creation of the sorted vector<adjacence> adjacences from adjacencies file "<<fic_adjacence<<" ..."<<flush;
   int j=0;
   string gene1="";
   string gene2="";
   string ori1="";
   string ori2="";
//   float dist=0.0;
//   float score=0.0;
//   string adj_tag="";
   adjacence adj;
   ifstream IfficAdj (fic_adjacence.c_str(), ios::in);
   if (!IfficAdj){
      cout<<"\nFrom General.cpp (AdjVector_Esp_AdjExtant): ERROR while opening file "<<fic_adjacence<<endl;
      if (file_log)
         Offile_log<<"\nFI\tFrom General.cpp (AdjVector_Esp_AdjExtant): ERROR while opening file "<<fic_adjacence<<endl;
      exit(EXIT_FAILURE);
   }
   while (!IfficAdj.eof()){
      //on lit la ligne en entier;
      IfficAdj>>gene1;      
      //Pour éviter un dernier tour de trop à cause d'un blanc en fin de fichier
      if(!IfficAdj.eof()){
         IfficAdj>>gene2;
//         if(IfficAdj>>ori1){
//            IfficAdj>>ori2;
//            IfficAdj>>dist;
//            IfficAdj>>score;
//            IfficAdj>>adj_tag;
//         }
//         cout<<gene1+" "+gene2+" "+ori1+" "+ori2+" "+convertFloat(dist)+" "+convertFloat(score)+" "+adj_tag<<endl;
//         cout<<"On traite l'adjacence "<< gene1<<" - "<<gene2;
         //on crée l'adjacence si les deux gènes sont dans gene_species_EXT (c-à-d dans nos arbres)
         if (gene_species_EXT.find(gene1)!=gene_species_EXT.end() && gene_species_EXT.find(gene2)!=gene_species_EXT.end()){
//            cout<<", on la nomme "<<j<<endl;
            adj.gene1=gene1;
            adj.gene2=gene2;
            //on la stocke
            adjacencies.push_back(adj);
            //Store extant adjacencies with associate species in Adj_actuelles
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
   if (file_log)
      Offile_log<<" DONE"<<endl<<endl;
}
