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

File: Backtracking.h               Last modified on: 20/11/2014
Created by: Sèverine Bérard            Created on: 10/01/2011
=========================================================================*/
#ifndef _PARCOURS_ARRIERE_H_
#define _PARCOURS_ARRIERE_H_

#include <Phyl/TreeTemplate.h>
using namespace bpp;

// The STL library:
#include <vector>

using namespace std;

#include "CostFunctions.h"

//VARIABLES GLOBALES 
Node * parcoursArriereCassure(Node * v1, Node * v2, bool v1dansA1);
Node * parcoursArriereC1ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe, vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
Node * parcoursArriereC1GLosAny(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
Node * parcoursArriereC1GDupExtantOrSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
Node * parcoursArriereC1GLosGLos(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
Node * parcoursArriereC1SpecSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
Node * parcoursArriereC1GDupGDup(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
void parcoursArriereC0ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe, vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
void parcoursArriereC0GLosAny(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences);
void parcoursArriereC0GDupExtantOrSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
void parcoursArriereC0GLosGLos(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
void parcoursArriereC0SpecSpec(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
void parcoursArriereC0GDupGDup(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
Node * parcoursArriereC1(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);
void parcoursArriereC0(Node * v1, Node * v2, vector<adjacence> * Adj_classe,vector<Tree *> * ArbresDAdjacences, bool v1dansA1);

#endif //_PARCOURS_ARRIERE_H_
