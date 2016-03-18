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

File: Step1_reconciliation.h                    Last modified on: 19/03/2015
Created by: Yoann Anselmetti/Sèverine Bérard    Created on: 28/02/2014
--------------------------------------------------------------------------
=========================================================================*/
#include "General.h"

//Methods also used in Step0_dataset.cpp & Step0_dataset_emf_only.cpp

void RewriteSpeciesGeneFile(map<string,string> &);
void RewriteAdjFile(map<string,string> &);

int main(int argc, char* argv[]);
