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

File: Step3_DECO.h                     Last modified on: 21/11/2014
Created by: Sèverine Bérard               Created on: 10/01/2011
=========================================================================*/
#ifndef _DECO_H_
#define _DECO_H_

#include <Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Phyl/Tree.h>
#include <Phyl/TreeTemplate.h>
#include <Phyl/TreeTemplateTools.h>
#include <Phyl/Node.h>

using namespace bpp;

// The STL library:
#include <iostream>
#include <string>
#include <map>
#include <vector>

using namespace std;

//#include "General.h"
//#include "CostFunctions.h"
#include "Backtracking.h"

void afficheInfoNoeud(Node * n);
void afficheInfoNoeuds(Node * n);
void ordonnePostfixe(Node * n, vector<Node *> & v);
int compteFeuilles(Node * n);
void ecritAdjAncestrale(Node * n, int num1, int num2, TreeTemplate<Node> * A1, TreeTemplate<Node> * A2, ofstream & OfficAdj2);
void ecritAdjAncestrales(Node * n, int num1, int num2, TreeTemplate<Node> * A1, TreeTemplate<Node> * A2, ofstream & OfficAdj2);
void ecritSORTIEAdj(Node * n, int num1, int num2, ofstream & OfficAdj2);
void ecritSORTIEDup(vector<adjacence> DUP, int num1, TreeTemplate<Node>* Arbre1, int num2, ofstream & OfficSORTIE_dup);
void DECO(TreeTemplate<Node>* S, TreeTemplate<Node>* Arbre1, TreeTemplate<Node>* Arbre2, vector<adjacence> * Adj_classe, ofstream & OfficAdj2, ofstream & OfficArbresAdj, ofstream & OfficSORTIE_dup, string nom_classe);

int main(int argc, char* argv[]);

#endif //_DECO_H_
