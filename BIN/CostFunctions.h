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

File: CostFunctions.h                  Last modified on: 20/11/2014
Created by: Sèverine Bérard               Created on: 10/01/2011
=========================================================================*/
#ifndef _FCOUT_H_
#define _FCOUT_H_

#include <Phyl/TreeTemplate.h>
using namespace bpp;

// The STL library:
#include <iostream>
#include <vector>

using namespace std;

#include "General.h"

//VARIABLES GLOBALES

//renvoie l'événement associé au noeud n
string E(Node * n);
int Espece(Node * n);
int EvtToInt(string e);

//min à plusieurs éléments
float min(float x, float y);
float min(float x, float y, float z);
float min(float x, float y, float z, float a, float b);
int min_pos(float x, float y);
int min_pos(float x, float y, float z);
int min_pos(float x, float y, float z, float a, float b);

//Calcul pour C0 et C1
bool appAdj(adjacence a, vector<adjacence> * Adjacences);
float recupCoutC1(Node * v1, Node * v2);
float recupCoutC0(Node * v1, Node * v2);
float C1ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe);
float C1GLosAny(Node * v1, Node * v2);
float C1GDupExtantOrSpec(Node * v1, Node * v2);
float C1GLosGLos(Node * v1, Node * v2);
float C1SpecSpec(Node * v1, Node * v2);
float C1GDupGDup(Node * v1, Node * v2);
float C0ExtantExtant(Node * v1, Node * v2, vector<adjacence> * Adj_classe);
float C0GLosAny(Node * v1, Node * v2);
float C0GDupExtantOrSpec(Node * v1, Node * v2);
float C0GLosGLos(Node * v1, Node * v2);
float C0SpecSpec(Node * v1, Node * v2);
float C0GDupGDup(Node * v1, Node * v2);
float calculeC0(Node * v1, Node * v2, vector<adjacence> * Adj_classe);
float calculeC1(Node * v1, Node * v2, vector<adjacence> * Adj_classe);

#endif //_FCOUT_H_
