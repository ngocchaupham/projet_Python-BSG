#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd
import extraction as ex
import analyse as anl
import matplotlib.pyplot as plt

# Recevoir des arguments indiquant le fichier FASTA contenant les séquences et le fichier GTF contenant les informations des régions d'intérêt :
# Des données initiales pour l'organisme (Ex : Homo_sapien et E.Coli)
Organisme = int(sys.argv[1])
tache = int(sys.argv[2])
if tache != (1, 2):
    if Organisme == 0:
        fichier_gtf = './data/Homo_sapiens.GRCh38.109.chr.gtf'
        fichier_fasta = './data/Homo_sapiens.GRCh38.dna.chromosome.1.fa'
        num_chromosome = 1
        table_codon = './data/codon_table_human.txt'
        file_out_name = './data/sequences_genes_H.fa'
        nom_fichier_traduction = './data/sequences_traduites_genes_H.fa'
        rapport_nom = './data/rapport_H.txt'
    if Organisme == 1:
        fichier_gtf = './data/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.56.gtf'
        fichier_fasta = './data/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.chromosome.Chromosome.fa'
        num_chromosome = 'Chromosome'
        table_codon = './data/codon_table_ecoli.txt'
        file_out_name = './data/sequences_genes_E.fa'
        nom_fichier_traduction = './data/sequences_traduites_genes_E.fa'
        rapport_nom = './data/rapport_E.txt'


    ########## Extraction des données ##########

        # Séquence génomique du chromosome
    genome_seq = ex.Sequence_genome(fichier_fasta)
        # Convertir gtf en dataframe
    gtf_df = ex.LireGTF(fichier_gtf)
        # Extraire le dataframe gtf entier pour un chromosome spécifique (ex: chr1)
    gtf_chr = ex.ExtraireGTF(gtf_df, num_chromosome)
        # Créer un dataframe gtf qui ne contient des annotations pour gène
    genes_df = ex.Genes_gtf(gtf_chr)
        # Sélectionner 50 gènes aléatoires
    list_50genes = ex.Selection_genes(genes_df)
        # Extraire des coordonnées des gènes sélectionnés
    gene_positions = ex.Position_genes(list_50genes, genes_df)
        # Extraire la séquence génomique du gène
    genes_sequences = ex.Extraire_sequences_genes(gene_positions, genome_seq, file_out_name)
        # Générer un dictionnaire contenant des transcrits du gène (Ex : genes_transcrits = {'gene_id':[liste des transcrits_id]})
    genes_transcrits = ex.Gene_transcrits(gtf_chr, list_50genes)
        # Générer un dictionnaire contenant la position des CDS du transcrit
    transcrit_cds = ex.Extraire_transcrit_cds(genes_transcrits, gtf_chr)
        # Ajouter dans le dictionnaire transcrit_cds, le séquence codante du transcrit
    transcrit_cds_exp = ex.Sequence_codante(transcrit_cds, genes_transcrits, genome_seq)
        # Lire le fichier du code génétique et générer la table des codons sous forme d'un dictionnaire
    codon_table = ex.Lire_table_codon(table_codon)
        # Générer un dictionnaire contenant la séquence traduite des transcrits
    genes_traductions = ex.Traduire_sequence(transcrit_cds_exp, genes_transcrits, codon_table, nom_fichier_traduction)

    ########## Analyse des séquences ##########

        # Calculer le pourcentage GC, résultant sous forme d'un dictionnaire {'gene_id':'pourcentage GC
    GC_percent = anl.Pourcentage_gc(gene_positions, genome_seq)
        # Calculer la taille de la séquence du gène, résultant sous forme d'un dictionnaire
    taille_gene = anl.Taille_nu(genes_sequences)
        # Calculer la taille moyenne des gènes
    taille_moyenne_genes = anl.Taille_nu_moyenne(taille_gene)
        # Calculer la taille de la séquence protéique
    taille_prot = anl.Taille_prot(genes_traductions)
        # Calculer la taille moyenne des séquences protéiques
    taille_moyenne_prot = anl.Taille_prot_moyenne(taille_prot)
        # Calculer le nombre moyen des nucléotides non codants
    Nombre_nu_non_codante = anl.Moyenne_non_codante(genes_sequences, genes_transcrits, transcrit_cds)


    # Créer un dataframe contenant des résultats d'analyses:
    Gene_ID = list_50genes
    Positions = list(gene_positions.values())
    Transcrits_id = list(genes_transcrits.values())
    Nb_transcrits = []
    # Parcourir les éléments de la série
    for sublist in Transcrits_id:
        # Calculer la longueur de chaque sous-liste et l'ajouter à la liste
        Nb_transcrits.append(len(sublist))
    Gene_taille = list(taille_gene.values())
    Prot_taille = list(taille_prot.values())
    Pourcentage_GC = list(GC_percent.values())

    data = {
        "Gene_ID": Gene_ID,
        "Positions": Positions,
        "Nombre_transcrits": Nb_transcrits,
        "Transcrits_id": Transcrits_id,
        "Gene_taille": Gene_taille,
        "Prot_taille": Prot_taille,
        "Pourcentage_GC": Pourcentage_GC
    }
    rapport_df = pd.DataFrame(data)
    # Exporter le dataframe résultant dans un fichier .txt
    rapport_df.to_csv(rapport_nom, sep='\t', index=False)


########## Comparer entre deux organismes  ##########
"""Résultats se retrouvent dans le fichier './data/comparaison_report.txt' """
if tache == 1:
    if not os.path.exists('./data/comparaison_report.txt'):
        # Lire le fichier rapport des deux organismes
        rapport_nom_euka = './data/rapport_H.txt'
        rapport_nom_proka = './data/rapport_E.txt'
        rapport_df_euka = pd.read_csv(rapport_nom_euka, sep="\t")
        rapport_df_proka = pd.read_csv(rapport_nom_proka, sep="\t")
            # Définir la colonne 'Gene_ID' comme index
        rapport_df_euka = rapport_df_euka.set_index('Gene_ID')
        rapport_df_proka = rapport_df_proka.set_index('Gene_ID')
            # Comparer le pourcentage GC moyen
                # Convertir les deux colonnes en un dictionnaire
        GC_percent_euka = rapport_df_euka.to_dict()['Pourcentage_GC']
        GC_percent_proka = rapport_df_proka.to_dict()['Pourcentage_GC']
        anl.Comparer_gc_moyen(GC_percent_euka, GC_percent_proka)

            # Comparer la taille moyenne des gènes
        gene_taille_euka = rapport_df_euka.to_dict()['Gene_taille']
        gene_taille_proka = rapport_df_proka.to_dict()['Gene_taille']
        gene_taille_moyenne_euka = anl.Taille_nu_moyenne(gene_taille_euka)
        gene_taille_moyenne_proka = anl.Taille_nu_moyenne(gene_taille_proka)
        anl.Comparer_taille_nu_moyenne(gene_taille_moyenne_euka, gene_taille_moyenne_proka)

            # Comparer la taille moyenne des séquences protéiques
        prot_taille_euka = rapport_df_euka.to_dict()['Prot_taille']
        prot_taille_proka = rapport_df_proka.to_dict()['Prot_taille']  # Problème : les valeurs du dictionnaire sont sous forme str
            # Convertir en dictionnaire avec les valeurs étant sous forme int
        prot_taille_euka_int = {}
        for cle, valeur in prot_taille_euka.items():
            liste_str = valeur.strip('[]').split(',')  # enlever les crochets et spliter par la virgule
            liste_int = [int(chiffre) for chiffre in liste_str]  # convertir chaque chiffre en int
            prot_taille_euka_int[cle] = liste_int
        prot_taille_proka_int = {}
        for cle, valeur in prot_taille_proka.items():
            liste_str = valeur.strip('[]').split(',')  # enlever les crochets et spliter par la virgule
            liste_int = [int(chiffre) for chiffre in liste_str]  # convertir chaque chiffre en int
            prot_taille_proka_int[cle] = liste_int

        prot_taille_moyenne_euka = anl.Taille_prot_moyenne(prot_taille_euka_int)
        prot_taille_moyenne_proka = anl.Taille_prot_moyenne(prot_taille_proka_int)
        anl.Comparer_taille_prot_moyen(prot_taille_moyenne_euka, prot_taille_moyenne_proka)

            # Comparer les moyennes d'acides nucléiques par acide aminé
        anl.Comparer_rapport_GeneProt(gene_taille_moyenne_euka, gene_taille_moyenne_proka, prot_taille_moyenne_euka, prot_taille_moyenne_proka)
    else:
        print("L'analyse est déjà faite. Veillez retrouver le fichier comparaison_report.txt dans dossier data")

########## Graphiques  ##########
if tache == 2:
    # Lire le fichier rapport des deux organismes et convertir en dataframe
    rapport_nom_euka = './data/rapport_H.txt'
    rapport_nom_proka = './data/rapport_E.txt'
    rapport_df_euka = pd.read_csv(rapport_nom_euka, sep="\t")
    rapport_df_proka = pd.read_csv(rapport_nom_proka, sep="\t")

    ##### Plotter pourcentage GC selon organisme
    if not os.path.exists('./figures/Pourcentage_GC.png'):   # Vérifier l'existence de cette figure
        # Extraire les colonnes Pourcentage_GC de chaque dataframe
        pourcentage_gc_H = rapport_df_euka['Pourcentage_GC']
        pourcentage_gc_E = rapport_df_proka['Pourcentage_GC']

        # Créer un tableau avec les pourcentages GC des deux organismes
        data = [pourcentage_gc_H, pourcentage_gc_E]

        # Créer un plot avec deux boxplots pour les deux organismes
        fig, ax = plt.subplots()
        ax.boxplot(data)

        # Ajouter des étiquettes d'axe et un titre
        ax.set_xticklabels(['Homo sapiens', 'E.Coli'])
        ax.set_ylabel('Pourcentage GC')
        ax.set_title('Comparaison des pourcentages GC des deux organismes')

        # Enregistrement de la figure dans un dossier
        fig.savefig("./figures/Pourcentage_GC.png")
    else:
        print("L'analyse du pourcentage GC a été faite. Veillez retrouver votre figure dans dossier 'figures'")

    ###### Plotter la longueur du gène en fonction de l'organisme
    if not os.path.exists('./figures/Gene_taille.png'):   # Vérifier l'existence de cette figure
        # Extraire les colonnes Pourcentage_GC de chaque dataframe
        gene_taille_H = rapport_df_euka['Gene_taille']
        gene_taille_E = rapport_df_proka['Gene_taille']

        # Créer un tableau avec les pourcentages GC des deux organismes
        data = [gene_taille_H, gene_taille_E]

        # Créer un plot avec deux boxplots pour les deux organismes
        fig, ax = plt.subplots()
        ax.boxplot(data)

        # Ajouter des étiquettes d'axe et un titre
        ax.set_xticklabels(['Homo sapiens', 'E.Coli'])
        ax.set_ylabel('La longueur du gène ')
        ax.set_title('Comparaison la longueur du gènes des deux organismes')

        # Enregistrement de la figure dans un dossier
        fig.savefig("./figures/Gene_taille.png")
    else:
        print("L'analyse de la longueur du gène a été faite. Veillez retrouver votre figure dans dossier 'figures'")

    ##### Plotter la longueur protéique en fonction de l'organisme
    if not os.path.exists('./figures/Prot_taille.png'):  # Vérifier l'existence de cette figure
        # Extraire les données de chaque dataframe
            # Comparer la taille moyenne des séquences protéiques
        prot_taille_H = rapport_df_euka['Prot_taille'].tolist()
        prot_taille_E = rapport_df_proka['Prot_taille'].tolist()  # Problème : les valeurs sont sous forme str
            # Convertir en dictionnaire avec les valeurs étant sous forme int
        prot_taille_H_int = []
        for sublist in prot_taille_H:
            liste_str = sublist.strip('[]').split(',')  # enlever les crochets et spliter par la virgule
            liste_int = [int(chiffre) for chiffre in liste_str]  # convertir chaque chiffre en int
            prot_taille_H_int.extend(liste_int)
        prot_taille_E_int = []
        for sublist in prot_taille_E:
            liste_str = sublist.strip('[]').split(',')  # enlever les crochets et spliter par la virgule
            liste_int = [int(chiffre) for chiffre in liste_str]  # convertir chaque chiffre en int
            prot_taille_E_int.extend(liste_int)
        # Trier chaque série de longueurs
        prot_taille_H_sorted = sorted(prot_taille_H_int)
        prot_taille_E_sorted = sorted(prot_taille_E_int)

        # Créer un axe y pour représenter les longueurs triées
        y1 = prot_taille_H_sorted
        y2 = prot_taille_E_sorted
        # Tracer les histogrammes
        plt.hist(y1, alpha=0.5, label='Homo sapiens')
        plt.hist(y2, alpha=0.5, label='E.Coli')


        # Ajouter une légende et un titre
        plt.legend()
        plt.ylabel('Le nombre de séquences protéiques')
        plt.xlabel('La longueur protéique ')
        plt.title('Distribution des longueurs protéiques pour 50 gènes')

        # Enregistrement de la figure dans un dossier
        plt.savefig("./figures/Prot_taille.png")
    else:
        print("L'analyse de la longueur protéique a été faite. Veillez retrouver votre figure dans dossier 'figures'")