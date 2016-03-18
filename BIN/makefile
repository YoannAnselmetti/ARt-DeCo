# ========================================================================
# This file is part of DeCo.

# Copyright (c) 2012, UM2, UL1, INRIA

# DeCo is free software. This software is governed by the CeCILL license
# under French law and abiding by the rules of distribution of free
# software

# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided
# only with a limited warranty and the software's author, the holder of
# the economic rights, and the successive licensors have only limited
# liability. See the file license_CeCILL_V2-en.txt for more details.

# File: makefile                                Last modified on: 2016/03/17
# Created by: Sèverine Bérard/Yoann Anselmetti  Created on: 2011/01/10
# --------------------------------------------------------------------------
# WARNING: You may have to modify in CFLAGS1 the "-I/usr/include/Bpp/" part
# depending on your installation of Bio++
# =========================================================================

CC = g++

CFLAGS1 = -pg -g -Wall -O2 -I/usr/include/Bpp/ -L/usr/lib
CFLAGS1_2 = -pg -g -Wall -O2 -static -I/usr/include/Bpp/ -L/usr/lib
CFLAGS2 = -lbpp-phyl -lbpp-seq -lbpp-core



all : Step0 Step1 Step2 Step3 Step3_proba Step4
#$(CC) $(CFLAGS1) $(OBJS) $(CFLAGS2) -o DeCo

static : Step0_static Step1_static Step2_static Step3_static Step3_proba_static Step4_static

Step0 : General.o Step0_dataset_emf_only.o
	$(CC) $(CFLAGS1) General.o Step0_dataset_emf_only.o $(CFLAGS2) -o Step0

Step1 : General.o Step1_reconciliation.o
	$(CC) $(CFLAGS1) General.o Step1_reconciliation.o  $(CFLAGS2) -o Step1

Step2 : General.o Step2_adjacencies_class.o
	$(CC) $(CFLAGS1) General.o Step2_adjacencies_class.o $(CFLAGS2) -o Step2

Step3 : General.o Step3_DECO.o CostFunctions.o Backtracking.o
	$(CC) $(CFLAGS1) General.o Step3_DECO.o CostFunctions.o Backtracking.o $(CFLAGS2) -o Step3

Step3_proba : General.o Step3_DECOProba.o CostFunctionsProba.o BacktrackingProba.o
	$(CC) $(CFLAGS1) General.o Step3_DECOProba.o CostFunctionsProba.o BacktrackingProba.o $(CFLAGS2) -o Step3_proba

Step4 : General.o Step4_statistics.o
	$(CC) $(CFLAGS1) General.o Step4_statistics.o $(CFLAGS2) -o Step4



General.o : General.cpp General.h
	$(CC) $(CFLAGS1) -c General.cpp $(CFLAGS2) -o General.o

Step0_dataset_emf_only.o : Step0_dataset_emf_only.cpp Step0_dataset.h
	$(CC) $(CFLAGS1) -c Step0_dataset_emf_only.cpp $(CFLAGS2) -o Step0_dataset_emf_only.o

Step1_reconciliation.o : Step1_reconciliation.cpp Step1_reconciliation.h
	$(CC) $(CFLAGS1) -c Step1_reconciliation.cpp $(CFLAGS2) -o Step1_reconciliation.o

Step2_adjacencies_class.o : Step2_adjacencies_class.cpp Step2_adjacencies_class.h
	$(CC) $(CFLAGS1) -c Step2_adjacencies_class.cpp $(CFLAGS2) -o Step2_adjacencies_class.o

Step3_DECO.o : Step3_DECO.cpp Step3_DECO.h
	$(CC) $(CFLAGS1) -c Step3_DECO.cpp $(CFLAGS2) -o Step3_DECO.o

Step3_DECOProba.o : Step3_DECOProba.cpp Step3_DECOProba.h
	$(CC) $(CFLAGS1) -c Step3_DECOProba.cpp $(CFLAGS2) -o Step3_DECOProba.o

CostFunctions.o : CostFunctions.cpp CostFunctions.h
	$(CC) $(CFLAGS1) -c CostFunctions.cpp $(CFLAGS2) -o CostFunctions.o

CostFunctionsProba.o : CostFunctionsProba.cpp CostFunctionsProba.h
	$(CC) $(CFLAGS1) -c CostFunctionsProba.cpp $(CFLAGS2) -o CostFunctionsProba.o

Backtracking.o : Backtracking.cpp Backtracking.h
	$(CC) $(CFLAGS1) -c Backtracking.cpp $(CFLAGS2) -o Backtracking.o

BacktrackingProba.o : BacktrackingProba.cpp BacktrackingProba.h
	$(CC) $(CFLAGS1) -c BacktrackingProba.cpp $(CFLAGS2) -o BacktrackingProba.o

Step4_statistics.o : Step4_statistics.cpp Step4_statistics.h
	$(CC) $(CFLAGS1) -c Step4_statistics.cpp $(CFLAGS2) -o Step4_statistics.o



Step0_static : General.cpp General.h Step0_dataset_emf_only.cpp Step0_dataset.h
	$(CC) $(CFLAGS1_2) General.cpp Step0_dataset_emf_only.cpp $(CFLAGS2) -o Step0

Step1_static : General.cpp General.h Step1_reconciliation.cpp Step1_reconciliation.h
	$(CC) $(CFLAGS1_2) General.cpp Step1_reconciliation.cpp  $(CFLAGS2) -o Step1

Step2_static : General.cpp General.h Step2_adjacencies_class.cpp Step2_adjacencies_class.h
	$(CC) $(CFLAGS1_2) General.cpp Step2_adjacencies_class.cpp $(CFLAGS2) -o Step2

Step3_static : General.cpp General.h Step3_DECO.cpp Step3_DECO.h
	$(CC) $(CFLAGS1_2) General.cpp Step3_DECO.cpp CostFunctions.cpp Backtracking.cpp $(CFLAGS2) -o Step3

Step3_proba_static : General.cpp General.h Step3_DECOProba.cpp Step3_DECOProba.h
	$(CC) $(CFLAGS1_2) General.cpp Step3_DECOProba.cpp CostFunctionsProba.cpp BacktrackingProba.cpp $(CFLAGS2) -o Step3_proba

Step4_static : General.cpp General.h Step4_statistics.cpp Step4_statistics.h
	$(CC) $(CFLAGS1_2) General.cpp Step4_statistics.cpp $(CFLAGS2) -o Step4



depend :
	$(CC) $(CFLAGS1) -MM General.cpp General.h General.o Step0_dataset.cpp Step0_dataset.h Step0_dataset.o Step1_reconciliation.cpp Step1_reconciliation.h Step1_reconciliation.o Step2_adjacencies_class.cpp Step2_adjacencies_class.h Step2_adjacencies_class.o Step3_DECO.cpp Step3_DECO.h Step3_DECO.o Step3_DECOProba.cpp Step3_DECOProba.h Step3_DECOProba.o CostFunctions.cpp CostFunctions.h CostFunctions.o CostFunctionsProba.cpp CostFunctionsProba.h CostFunctionsProba.o Backtracking.cpp Backtracking.h Backtracking.o BacktrackingProba.cpp BacktrackingProba.h BacktrackingProba.o Step4_statistics.cpp Step4_statistics.h Step4_statistics.o > depend_file

clean :
	rm -f Step0 Step1 Step2 Step3 Step3_proba Step4 core *.o 
	touch depend_file

extremiste :
	rm -f *~

include depend_file
