#!/bin/bash
###########################################################################
###########################################################################
###                                                                     ###
### Goal:                                                               ###
###    Script permettant de lancer 1 run de DeCo                        ###
### INPUT:                                                              ###
###    1- Config file of DeCo/ARt-DeCo                                  ###
###    2- 0 -> Generate data from Ensembl | 1 -> Begin to Step1         ###
###    3- 0 -> DeCo algorithm | 1 -> ARt-DeCo algorithm                 ###
### OUTPUT:                                                             ###
###    OUTPUT files of DeCo/ARt-DeCo                                    ###
###                                                                     ###
### Name: deco_run.sh              Author: Yoann Anselmetti             ###
### Creation: 2014/11/21           Last modification: 2015/03/26        ###
###                                                                     ###
###########################################################################
###########################################################################


if [[ $2 == 1 ]]
then
   echo -en "\nFrom deco_run.sh: Run Step0 + Step2 + "
   if [[ $3 == 1 ]]
   then
      echo -en "Step3_proba + Step4:\n\n"
      echo -e "\tFrom deco_run.sh: Run ./Step0 "$1
      ./Step0 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step0 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step2 "$1
      ./Step2 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step2 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step3_proba "$1
      ./Step3_proba $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step3_proba '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step4 "$1
      ./Step4 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step4 '$1' execution!'
			exit
      fi
   elif [[ $3 == 0 ]]
   then
      echo -en "Step3 + Step4:\n\n"
      echo "From deco_run.sh: Run ./Step0 "$1
      ./Step0 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step0 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step2 "$1
      ./Step2 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step2 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step3 "$1
      ./Step3 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step3 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step4 "$1
      ./Step4 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step4 '$1' execution!'
			exit
      fi
   else
      echo -en "\nFrom deco_run.sh: Wrong value for 3rd parameter: '$3' Should be 0 (DeCo) or 1 (DeCo++)"
   fi      
elif [[ $2 == 0 ]]
then
   echo -en "\nFrom deco_run.sh: Run Step0 + Step2 + "
   if [[ $3 == 1 ]]
   then
      echo -en "Step3_proba + Step4:\n\n"
      echo -e "\tFrom deco_run.sh: Run ./Step1 "$1
      ./Step1 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step1 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step2 "$1
      ./Step2 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step2 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step3_proba "$1
      ./Step3_proba $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step3_proba '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step4 "$1
      ./Step4 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step4 '$1' execution!'
			exit
      fi
   elif [[ $3 == 0 ]]
   then
      echo "Step3 + Step4:"
      echo -e "\tFrom deco_run.sh: Run ./Step1 "$1
      ./Step1 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step1 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step2 "$1
      ./Step2 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step2 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step3 "$1
      ./Step3 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step3 '$1' execution!'
			exit
      fi
      echo -e "\tFrom deco_run.sh: Run ./Step4 "$1
      ./Step4 $1
		if [[ $? != 0 ]]
      then
			echo -e '\tFrom deco_run.sh: ERROR during ./Step4 '$1' execution!'
			exit
      fi
   else
      echo -en "\nFrom deco_run.sh: Wrong value for 3rd parameter: '$3' Should be 0 (DeCo) or 1 (DeCo++)"
   fi   
else
   echo -en "\nFrom deco_run.sh: Wrong value for 2nd parameter: '$2' Should be 0 (DeCo) or 1 (DeCo++)"
fi
