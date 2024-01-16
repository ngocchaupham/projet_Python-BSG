#!/bin/bash

# Demande à l'utilisateur de saisir des arguments
echo -e " Ce script permer d'excécuter l'ensemble des scripts Python dans ce projet\n" 

echo -e "Voulez-vous comparer entre deux organismes ? Pour cela, il faut disposer deux fichiers rapport.txt de vos deux organismes !
0 : Non, je ne dispose pas ces deux fichiers. Je dois faire d'abord une analyse des séquences 
1 : Oui, et je veux obtenir un fichier 'comparaison_report.txt'
2 : Oui, et je veux visualiser les graphes de comparaison \n"
read arg2

if [ "$arg2" -eq 1 ] || [ "$arg2" -eq 2 ]; then
	arg1=1 
	if [ "$arg2" -eq 2 ]; then
		if [ -d ./figures ]; then
	    	
	    	echo "Votre analyse est en cours..."
		else
			mkdir ./figures
		fi

	fi
else
	echo -e  "Veuillez saisir l'organisme de votre analyse :
	Propos : 
	0 : Homo_sapiens (chromosome 1)
	1 : Escherichia coli K12\n"
	read arg1
fi



# Exécute le script Python en lui passant l'argument
python ./main.py $arg1 $arg2 