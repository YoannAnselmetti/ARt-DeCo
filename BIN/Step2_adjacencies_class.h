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

File: Step2_adjacencies_class.h                Last modified on: 01/04/2015
Created by: Yoann Anselmetti/Sèverine Bérard   Created on: 28/02/2014
--------------------------------------------------------------------------
Specification:
Class file for Step2_adjacencies_class.cpp
=========================================================================*/
#include "General.h"

bool isAncestor(int id1, int id2, TreeTemplate<Node> *A);
int trouveNoeudSpeSAuDessusDeN(int s, int id_n,TreeTemplate<Node> * G);
void retourneIdSpeSAuDessousDeN(vector<int> & Ids, int s, Node * n,TreeTemplate<Node> * G);

void AdjClassPreConstruction(TreeTemplate<Node> *S, vector<Tree*> &Arbres, map<string,vector<int> > &gene_GeneTreeID, map<string,string> &gene_species_EXT, map<string,vector<int> > &classes_adjacences, vector<adjacence> &adjacencies);
void AdjClassConstruction(TreeTemplate<Node> *S, vector<Tree*> &Arbres, map<string,string> &gene_species_EXT, map<string,vector<int> > &classes_adjacences, vector<adjacence> &adjacencies);
void WriteAdjClassFile(map<string,vector<int> > &classes_adjacences);
void WriteAdjClassFileHumanRead(map<string,vector<int> > &classes_adjacences, vector<adjacence> &adjacencies, map<string,string> &gene_species_EXT);

int main(int argc, char* argv[]);
