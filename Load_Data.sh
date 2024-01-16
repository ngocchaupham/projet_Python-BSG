#!/bin/bash
#Créer un répertoire pour ce projet : 
current_dir=$(pwd)
echo "Le répertoire actuel est : $current_dir"
read -p "Entrez le chemin absolu ou relatif du répertoire où vous voulez travailler: " dir
cd $dir
mkdir Projet_de_programmation_Python_PHAM_2023
cd ./Projet_de_programmation_Python_PHAM_2023


#Créer  un fichier texte contenant le nom des organismes traitées :
FILE=Nom_des_organismes.txt
if [ -e $FILE ]; then # Vérifier l'existance de ce fichier.
    echo "Le fichier existe déjà"
else
    echo 'Escherichia coli str. K-12 substr. MG1655 str. K12 (GCA_000005845) and Homo_sapiens (chromosome 1)' > $FILE
    exec 1>&- # Fermer le fichier
fi

# Créer un fichier qui contient des liens d'Ensembl pour télécharger des données :
LIENS=urls.txt
if [ -e $LIENS ]; then # Vérifier si le fichier existe
    echo "Le fichier existe déjà"
else
    echo -e 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-56/bacteria//gtf/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.56.gtf.gz \n' > $LIENS

    echo -e 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/release-56/bacteria//fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/dna/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa.gz \n' >> $LIENS

    echo -e 'https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz \n' >> $LIENS

    echo 'https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz' >> $LIENS
    exec 1>&-
fi



#Télécharger des données depuis le site Ensembl: 
mkdir data 
while read url; do
    wget "$url" -P ./data
done < urls.txt

#Extraire des fichiers .gz :
cd data 
for file in *.gz; do
    gunzip "$file"
done

#Télécharger le code génétique: 
cd data
curl -o "codon_table_human.txt" "https://zenodo.org/record/7728908/files/codon_table_human.txt?download=1" 
curl -o "codon_table_ecoli.txt" "https://zenodo.org/record/7728908/files/codon_table_ecoli.txt?download=1"
cd ..


