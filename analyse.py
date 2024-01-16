#!/usr/bin/python3
# -*- coding: utf-8 -*-
import extraction as ex
# Importer la librairie statistics
import statistics as stat
import warnings
# Ignorer les warnings
warnings.filterwarnings('ignore')


def Pourcentage_gc(gene_positions, genome_seq):
    """Calculer le pourcentage GC d'une séquence d'ADN."""
    # Dictionnaire pour stocker le gene_id et son pourcentage GC
    GC_percent = {}
    for gene_id, (strand, start, end) in gene_positions.items():
        sequence = genome_seq[start - 1:end]
        # La longueur totale de la séquence
        length_total = len(sequence)
        # Le nombre de paires 'GC' contenu dans la séquence
        nb_GC = sequence.count('G') + sequence.count('C')
        # Pourcentage GC :
        GC_percent[gene_id] = (nb_GC / length_total) * 100
    return GC_percent


def Comparer_gc_moyen(GC_percent_euka, GC_percent_proka):
    """ Comparer le pourcentage GC moyen entre l'organisme procaryote et eucaryote."""
    GC_moyen_euka = stat.mean(GC_percent_euka.values())
    GC_moyen_proka = stat.mean(GC_percent_proka.values())
    with open("./data/comparaison_report.txt", "a") as f:
        print("Voici des comparaisons entre l'organisme procaryote et eucaryote:", '\n\n', file=f)
        print("Le pourcentage GC moyen de Homo sapiens est de {:.2f}%.".format(GC_moyen_euka), '\n', file=f)
        print("Le pourcentage GC moyen de E.coli est de {:.2f}%.".format(GC_moyen_proka), '\n', file=f)
        if GC_moyen_euka > GC_moyen_proka:
            print("Le pourcentage GC moyen de l'organisme eucaryote est supérieur à celui de l'organisme procaryote.", '\n\n', file=f)
        elif GC_moyen_euka < GC_moyen_proka:
            print("Le pourcentage GC moyen de l'organisme eucaryote est inférieur à celui de l'organisme procaryote.", '\n\n', file=f)
        else:
            print("Le pourcentage GC moyen des deux organismes sont égaux.", '\n\n', file=f)


def Taille_nu(genes_sequences):
    """Déterminer la taille de la séquence nucléique pour chaque gène"""
    # Dictionnaire pour stocker le gene_id et la longueur de sa séquence nucléique
    gene_taille = {}
    for gene_id, sequence in genes_sequences.items():
        seq_length = len(sequence)
        gene_taille[gene_id] = seq_length
    return gene_taille


def Taille_nu_moyenne(gene_taille):
    """Calculer la taille moyenne des séquences nucléiques entre deux organismes"""
    taille_moyenne_nu = stat.mean(gene_taille.values())
    return taille_moyenne_nu

def Comparer_taille_nu_moyenne(taille_moyenne_nu_euka, taille_moyenne_nu_proka):
    """Comparer la taille moyenne des séquences nucléiques entre deux organismes"""
    with open("./data/comparaison_report.txt", "a") as f:
        print("La taille moyenne des séquences nucléiques de l'organisme eucaryote est de", taille_moyenne_nu_euka, 'pb', '\n', file=f)
        print("La taille moyenne des séquences nucléiques de l'organisme procaryote est de", taille_moyenne_nu_proka, 'pb', '\n', file=f)
        if taille_moyenne_nu_euka > taille_moyenne_nu_proka:
            print("La taille moyenne des séquences nucléiques est plus grande chez l'organisme eucaryotes que chez l'organisme procaryote.", '\n\n', file=f)
        elif taille_moyenne_nu_euka < taille_moyenne_nu_proka:
            print("La taille moyenne des séquences nucléiques est plus petite chez l'organisme eucaryotes que chez l'organisme procaryote.", '\n\n', file=f)
        else:
            print("La taille moyenne des séquences nucléiques est égale chez les deux organismes.", '\n\n', file=f)

def Taille_prot(genes_traduction):
    """Déterminer la taille de la séquence protéique pour chaque gène"""
    # Dictionnaire pour stocker le gene_id et la longueur de sa séquence protéique
    prot_taille = {}
    for gene_id, transcrits in genes_traduction.items():
        for transcrit_id, seq_prot in transcrits.items():
            seq_prot_length = len(seq_prot)
            if gene_id not in prot_taille:
                prot_taille[gene_id] = [seq_prot_length]
            else:
                prot_taille[gene_id].append(seq_prot_length)
    return prot_taille

def Taille_prot_moyenne(prot_taille):
    """Calculer la taille moyenne des séquences nucléiques entre deux organismes"""
    liste_prot_taille = []
    for gene_id, prot_length in prot_taille.items():
        liste_prot_taille.append(prot_length)  # Cette liste contient des sous-listes
        # Concaténer toutes les sous-listes
        liste_con = []
        for sublist in liste_prot_taille:
            liste_con.extend(sublist)
    taille_moyenne_prot = stat.mean(liste_con)
    return taille_moyenne_prot

def Comparer_taille_prot_moyen(taille_moyenne_prot_euka, taille_moyenne_prot_proka):
    """Comparer la taille moyenne des séquences protéiques entre deux organismes"""
    with open("./data/comparaison_report.txt", "a") as f:
        print("La taille moyenne des séquences protéiques de l'organisme eucaryote est de", taille_moyenne_prot_euka, '\n', file=f)
        print("La taille moyenne des séquences protéiques de l'organisme procaryote est de", taille_moyenne_prot_proka, '\n', file=f)
        if taille_moyenne_prot_euka > taille_moyenne_prot_proka:
            print("La taille moyenne des séquences protéiques est plus grande chez l'organisme eucaryotes que chez l'organisme procaryote.", '\n\n', file=f)
        elif taille_moyenne_prot_euka < taille_moyenne_prot_proka:
            print("La taille moyenne des séquences protéiques est plus petite chez l'organisme eucaryotes que chez l'organisme procaryote.", '\n\n', file=f)
        else:
            print("La taille moyenne des séquences protéiques est égale chez les deux organismes.", '\n\n', file=f)


def Comparer_rapport_GeneProt(taille_moyenne_nu_euka, taille_moyenne_nu_proka, taille_moyenne_prot_euka, taille_moyenne_prot_proka):
    """Comparer les moyennes d'acides nucléiques par acide aminé des deux organismes"""
    nucl_par_prot_euka = taille_moyenne_nu_euka/taille_moyenne_prot_euka
    nucl_par_prot_proka = taille_moyenne_nu_proka/taille_moyenne_prot_proka
    with open("./data/comparaison_report.txt", "a") as f:
        print("Le rapport moyen d'acides nucléiques par acide aminé de l'organisme eucaryote est de", nucl_par_prot_euka, '\n', file=f)
        print("Le rapport moyen d'acides nucléiques par acide aminé de l'organisme procaryote est de", nucl_par_prot_proka, '\n', file=f)
        # Comparer les moyennes d'acides nucléiques par acide aminé des deux organismes
        if nucl_par_prot_euka > nucl_par_prot_proka:
            print("Organisme eucaryote a un rapport moyenne d'acides nucléiques par acide aminé plus élevé.", '\n\n', file=f)
        elif nucl_par_prot_euka > nucl_par_prot_proka:
            print("Organisme procaryote a un rapport moyenne d'acides nucléiques par acide aminé plus élevé.", '\n\n', file=f)
        else:
            print("Les deux organismes ont le même rapport moyenne d'acides nucléiques par acide aminé.", '\n\n', file=f)


def Moyenne_non_codante(genes_sequence, genes_transcrits, transcrit_cds):
    """Calculer la moyenne du nombre d'acides nucléiques non codants par gène"""
    #Liste pour stocker le nombre de nucléotides non codants pour chaque gène
    non_codant_counts = []
    for gene_id, sequence in genes_sequence.items():
        #Calculer la taille de la séquence génomique du gène
        gene_seq = genes_sequence[gene_id]
        #Récupérer la liste des transcrits du gène
        transcrits = genes_transcrits[gene_id]
        for transcrit_id in transcrits:
            if transcrit_id in transcrit_cds.keys():
                #Calculer la taille de la séquence codante du gène
                cds_seq = transcrit_cds[transcrit_id]['seq_codante']
                cds_seq_length = len(cds_seq)
                non_codant_seq_length = len(gene_seq) - cds_seq_length
                non_codant_counts.append(non_codant_seq_length)
    # Calcule la moyenne
        if non_codant_counts:
            moyen = stat.mean(non_codant_counts)
            print("La moyenne du nombre d'acides nucléiques non codants par gène chez l'organisme eucaryote est de", moyen)
        else:
            moyen = 0.0
        return moyen