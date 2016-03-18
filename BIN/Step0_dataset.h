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

File: Step0_dataset.h                      Last modified on: 29/09/2015
Created by: Yoann Anselmetti               Created on: 08/04/2014
--------------------------------------------------------------------------
Specification:
Create Gene_Specie association and adjacences from Genomicus file on
syntenies from Ensembl Database.
(Future => File use to create produce INPUT data from Ensembl Database )
=========================================================================*/   
#include "General.h"

//Pour création des fichiers en entrée
#include <algorithm>
#include <sys/types.h>
//#include <regex>
//#include <boost/regex.hpp>

vector<string> read_directory(string& path);
void minuscule(string &chaine);
set<string> extantGenes(vector<Tree*> Arbres);
string getcwd_string();

int main(int argc, char* argv[]);
