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

File: Step4_statistics.h                       Last modified on: 16/06/2015
Created by: Yoann Anselmetti/Sèverine Bérard   Created on: 10/03/2014
=========================================================================*/
#ifndef _STATS_H_
#define _STATS_H_

#include "General.h"

string E(Node * n);
bool isLinear(vector<string> v, map<string,int> deg);
vector<adjacence> composantesConnexes(vector<adjacence> DUPLI, map<int,int> &distrib_taille_compo_connexe, int i);
void File_DOT_for_Adj_Graph(map<int,vector<adjacence> > Adj_actuelles, map<int,vector<adjacence> > Adj_nouvelles, TreeTemplate<Node> * S);
void File_TULIP_for_Adj_Graph(map<int,vector<adjacence> > Adj_actuelles, map<int,vector<adjacence> > Adj_nouvelles, TreeTemplate<Node> * S);
string nomGene(string nom);
void retourneFeuillesNomsGenes(Node * n, vector<string> & feuilles);
void ecritGenes(Node * n, int num_arbre, ofstream & Offic_SORTIE_genes);
int countNodesWithBranchProperty(Node * n, string prop, string val);
void lectureFicRelBin(string fic,map<int,vector<adjacence> > & miva);

int main(int argc, char* argv[]);

#endif //_STATS_H_
