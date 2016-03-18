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

File: General.h                                Last modified on: 22/05/2015
Created by: Yoann Anselmetti/Sèverine Bérard   Created on: 28/02/2014
--------------------------------------------------------------------------
Specification:
Main file of the project
=========================================================================*/
#ifndef _GENERAL_H_
#define _GENERAL_H_

#define LOG_PRINT(...) log_print(__FILE__, __LINE__, __VA_ARGS__ )

#include <Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Phyl/Tree.h>
#include <Phyl/TreeTemplate.h>
#include <Phyl/TreeTemplateTools.h>
#include <Phyl/NodeTemplate.h>
#include <Phyl/Node.h>

using namespace bpp;

// The STL library:
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <time.h>
#include <cstdlib>//Pour le rand
//#include <stdlib.h>
#include <ctime>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <cmath> // pour std::abs 
#include <limits> // pour std::numeric_limits 

using namespace std;

// Use in BackT / CostFct / General / Step1 / Step3
template <class T> class convert{  
   public:
      T from_s(string s1){
         T res;
         istringstream ftmp(s1);
         ftmp >> res;
         return (res);
      }   
      string to_s(const T& f1){
         ostringstream ftmp;
         ftmp << f1;
         return (ftmp.str());
      }
};

// Use in BackT / CostFct / General / Step1 / Step3
struct adjacence{
   string gene1;
   string gene2;
};

//   CONSTANTES   //
const string esp = "S";  //nom de la propriété des noeuds contenant leur espèce
const string typ = "Ev";  //nom de la propriété des noeuds contenant leur type
const string NAD = "L";  //nom de la propriété des branches contenant le nb d'adj dupl.
//Les noms des événements :
const string dupl = "GDup"; 
const string Adupl = "ADup"; 
const string spe = "Spec";   
const string bk = "Break"; 
const string ga = "Extant";  
const string Aper = "ALos";   
const string per = "GLos";   
const string NomPer = "Loss";  
const string NomPerA = "ADJLoss";  
const string NomCass = "Breakage";  
const float INFINI=100000000.0;      // Use only in CostFunctions.cpp
//Noms des fichiers de sortie
const string output_genes="_OUTPUT_genes";
const string output_spe="_OUTPUT_species";
const string output_rec="_reconcil";
const string output_adj_class="_OUTPUT_adjacencies_classes";
const string output_adj_class_read="_OUTPUT_adjacencies_classes_vis";
const string output_adj_trees_per_class="_OUTPUT_adjacencies_trees";
const string output_adj="_OUTPUT_adjacencies";
const string output_dup_gene_pairs="_OUTPUT_duplicated_gene_pairs";
const string output_treated_subT="_OUTPUT_treated_subtrees";
const string output_reconcil_stats="_OUTPUT_stats_reconcil";
const string output_DECO_stats="_OUTPUT_stats_DECO";
const string output_human_stats="_OUTPUT_stats_human_readable";
const string output_machine_stats="_OUTPUT_stats_machine";
const string output_new_adj="_OUTPUT_new_adjacencies";
const string output_ext_adj_file="_OUTPUT_extant_adjacencies_file";
const string output_graph="_OUTPUT_adjacencies_graph";
const string output_species_stats="_OUTPUT_species_stats";
const string output_AdjClass_with_GT="_OUTPUT_adjacencies_classes_with_gene_trees";


//   VARIABLES GLOBALES   //

// Parameters to open config file   //
extern string fic_arbre;
extern string fic_gene;
extern string fic_especes;
extern string fic_adjacence;
extern string exp_name;
extern string directory;
extern string prefixe;
extern string ST_input;
extern string file_spe;
extern string dir_embl;
extern string dir_emf;
extern string del_spe_data;
extern string log_file;
extern ofstream Offile_log;

//Parameters from config file
extern char sep;
extern char sepAdj;
extern int INPUT_FORMAT; //valable pour les arbres de gènes (l'arbre des espèces étant en newick)
extern int OUTPUT_FORMAT;//0 pour newick ; 1 pour NHX 
extern float Adj_percentage;
extern float Crea;
extern float Break;
extern float coeff_b;
//extern float Spec;
//extern float GDup;
//extern float GLos;
//extern float ADup;
//extern float ALos;
extern bool ArbresRec;


//Parameters after cleaning genes and adjacencies
extern string fic_gene_clean;
extern string fic_adjacence_clean;

//   Fichiers intermédiaires et de sortie   //
extern string tree_reconciled_file;
extern string adj_class_file;
extern string DECO_stats_file;

//For statistics (Step4) => initialiser dans Step3_DECO(Proba)
extern int nb_GDup;
extern int nb_ADup;
extern int nb_Crea;
extern int nb_Bk;

////////////////
//  BOOLEAN   //
////////////////

// Boolean for programm mode
extern bool proba_mode;
extern bool light_mode;

// Boolean for output screen
extern bool affich_CalculCout;
extern bool affich_ParcoursArriere;
extern bool affich_Classes;

// Boolean for files creation
extern bool fic_SORTIE_EspeceDone;
extern bool fic_SORTIE_GeneDone;
extern bool file_log;

////////////////////////
///   MAP & VECTOR   ///
////////////////////////

// Step3_proba, Step_4 & Backtracking(Proba)
extern map<int,int> espece_nb_gain; //ces deux maps pour compter les événements sur les branches
extern map<int,int> espece_nb_break; // de l'arbre des espèces => Pas utilisé Dans Backtracking(Proba)
extern map<string,int> extant_gene_adj_nb;   //string: nom gène | int: nb_adj du gène
extern vector<adjacence> DUP; //Pour recueillir les gènes dupliqués ensemble (sert pour des calculs intermédiaires)
extern int id_arbres_adj;

// Pour la déclaration en externe dans CostFunctions.cpp
extern map<pair<Node *,Node *>, float > C0;
extern map<pair<Node *,Node *>, float > C1;
extern map<pair<Node *,Node *>, pair<float,int> > C0P;
extern map<pair<Node *,Node *>, pair<float,int> > C1P;
extern map<pair<string,string>, float> P;


/////////////////////
////   METHODS   ////
/////////////////////

void mkdir_rec(const char *dir);
string convertInt(int number);
bool file_exists(const string& name);
bool areEqual(float a,float b);
bool is_readable(string& file);

//Pour affichage de map
void afficheMap(map<string, vector<int> > m);
void afficheMap(map<int, vector<adjacence> > m);
void afficheMap(map<int, vector<string> > m);
void afficheMap(map<string, string> m);
void afficheMap(map<string, int> m);
void afficheMap(map<int, int> m);

//Pour lecture du fichier de config
int lireNomFic(ifstream& IfficConf, string& fic, string nomfic, char** argv);
int lireCout(ifstream& IfficConf, float& Cout, string nomcout);
void lireFicConfig(string& fic_arbre, string& fic_gene, string& fic_especes, string& fic_adjacence, string& exp_name, string& directory,int argc, char** argv);

//   METHODS Step0 -> Step4   //
void CreationOutputSpecies(TreeTemplate<Node> *S);
void AssociateGeneWithSpecies(map<string,string> &gene_species_EXT);
void nommerNoeudsInternes(Node * n);
void recupInfoSurBrancheNoeud(Node * n);
void recupInfoSurBrancheNoeuds(Node * n, TreeTemplate<Node> * S);
void affecteInfoSurBrancheNoeud(Node * n);
void affecteInfoSurBrancheNoeuds(Node * n); 
//   METHODS Step0(_emf) + Step1   //
void StoreGeneTrees(TreeTemplate<Node> *S, vector<Tree *> &Arbres);
void retourneFeuilles(Node * n, vector<string> & feuilles);
void ecritEspeces(Node * n, ofstream & Offic_SORTIE_esp);
void reetiquette(Node * n, map<string,string> &gene_species_EXT);
Node * preAffecteEspece(Node * n, TreeTemplate<Node> * G, TreeTemplate<Node> * S);
void affecteEspece(Node * n, TreeTemplate<Node> * S);
void renumeroteId(Node * n, int & id);
void annoteType(Node * n, TreeTemplate<Node> * S);
Node * ajouteNoeudPerte(Node * x_S, Node * f_G,  int &id, TreeTemplate<Node> * S);
void ajoutePerte(Node * n, int &id, TreeTemplate<Node> S);
void reconcil(TreeTemplate<Node > * S, TreeTemplate<Node > * G);
void reconciliation(TreeTemplate<Node > * S, vector<Tree *> &Arbres, map<string,string> &gene_species_EXT);
//   METHODS Step1 -> Step4   //
void StoreSpeciesTree(TreeTemplate<Node> * S);
void StoreReconciledTrees(vector<Tree *> &Arbres);
void AssociateGeneWithReconciledTreeNb(vector<Tree *> &Arbres, map<string,vector<int> > &gene_GeneTreeID, map<string,string> &gene_species_EXT);
int GeneSpeciesClean(map<string,vector<int> > &gene_GeneTreeID, map<string,string> &gene_species_EXT);

//   Step3(_proba) -> Step4   //
int countNodesWithBranchPropertyBySpecies(Node * n, string prop, string val, map<int,int>& especes_nb_genes);
void AdjClassVector(map<string,vector<int> > &classes_adjacences);
void AdjVector_Esp_AdjExtant(TreeTemplate<Node> * S, map<string,string> &gene_species_EXT, vector<adjacence> &adjacencies, map<int,vector<adjacence> > &Adj_actuelles);

#endif //_GENERAL_H_
