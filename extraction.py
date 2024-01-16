#!/usr/bin/python3
# -*- coding: utf-8 -*-
import pandas as pd
import random
import warnings
# Ignorer les warnings
warnings.filterwarnings('ignore')


def LireGTF(fichier_gtf):
    """ Lire le fichier gtf et convertir en dataframe de pandas """
    columns_name = ["chrname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gtf_df = pd.read_csv(fichier_gtf, sep="\t", comment="#", header=None, names=columns_name)
    return gtf_df


def ExtraireGTF(gtf_df, num_chromosome):
    """Créer un nouveau gtf qui ne contient que des annotations sur un chromosome"""
    gtf_chr = gtf_df[gtf_df['chrname'] == num_chromosome]
    gtf_chr.to_csv('./data/Homo_sapiens.GRCh38.109.chr.1.gtf', sep='\t', index=False)
    return gtf_chr


def Genes_gtf(gtf_chr):
    """Créer un nouveau dataframe qui ne contient que les gènes"""
    genes_df = gtf_chr[(gtf_chr['feature'] == 'gene') & (gtf_chr['attribute'].str.contains('gene_biotype "protein_coding"') == True)]
    genes_df.reset_index(drop=True)
    # Ajouter d'une colonne 'gene_id' à partir des infos dans la colonne 'attribue'
    genes_df['gene_id'] = genes_df['attribute'].str.extract('gene_id "(.*?)";', expand=False)
    return genes_df


def Sequence_genome(fasta_file):
    """Récupérer la séquence génomique à partir du fichier Fasta"""
    with open(fasta_file, 'r') as file:
        # Variables pour stocker l'identifiant et la séquence du génome
        genome_id = ''
        genome_seq = ''
        # Boucle pour parcourir les lignes du fichier fasta
        for line in file:
            # Si la ligne commence par le caractère ">", c'est l'identifiant du génome
            if line.startswith('>'):
                genome_id = line.strip()[1:]
            # Sinon, c'est une ligne de séquence
            else:
                genome_seq += line.strip()
    return genome_seq


def Selection_genes(genes_df):
    """Sélectionner 50 gènes aléatoires"""
    genes_random = genes_df.sample(50)
    # Convertir une série de pandas en une liste
    list_50genes_random = genes_random['gene_id'].tolist()
    return list_50genes_random


def Position_genes(list_id_genes, genes_df):
    """Extraire les coordonnées des gènes dans la liste"""
    # Dictionnaire pour stocker le gene_id et ses positions
    gene_positions = {}
    # Pour chaque gene_id dans la liste, extraire ses positions dans dataframe
    for i in list_id_genes:
        # Extraire toutes les colonnes de la ligne correspondant au gene_id
        gene_info = genes_df.loc[genes_df['gene_id'] == i]
        # Créer des variables pour stocker des infos sur la coordonnée
        strand = str(gene_info['strand'].iloc[0])
        start = int(gene_info['start'].iloc[0])
        end = int(gene_info['end'].iloc[0])
        # Stocker ces infos dans un dictionnaire
        gene_positions[i] = [strand, start, end]
    return gene_positions


def Reverse_complement(sequence):
    """Convertir la séquence lue start-end sur le brin '+' en une séquence reversée
    et complémentaire (lue sur le brin '-') """
    # Dictionnaire de complémentarité
    base_complement = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G', 'N': 'N'}
    # inverser la séquence sur le brin '+'
    sequence_rev = reversed(sequence)
    sequence_rev_com = ''
    # Liste de bases complémentaires inversées
    for base in sequence_rev:
        sequence_rev_com += base_complement[base]
    return sequence_rev_com


def Extraire_sequences_genes(gene_positions, genome_seq, file_out_name):
    """Parcourir le dictionnaire contenant les coordonnées du gène,
    puis extraire sa séquence et écrire dans un fichier fasta"""
    #Dictionnaire pour stocker le gene_id et sa séquence
    genes_sequence = {}
    with open(file_out_name, 'w') as f:
        for gene_id, (strand, start, end) in gene_positions.items():
            # Soustraire 1 pour aligner avec les index de Python
            sequence = genome_seq[start - 1:end]
            if strand == '+':
                f.write(f">{gene_id}:{start}-{end} strand = '{strand}'\n{sequence}\n")
                genes_sequence[gene_id] = sequence

            else:
                sequence_rev_com = Reverse_complement(sequence)
                genes_sequence[gene_id] = sequence_rev_com
                f.write(f">{gene_id}:{start}-{end} strand ='{strand}'\n{sequence_rev_com}\n")
        return genes_sequence

def Gene_transcrits(gtf_chr, list_id_genes):
    """Extraire tous les transcrits que possède d'un gène"""
    # Dictionnaire pour stocker le gene_id et le transcrit_id
    genes_transcrits = {}
    # Extraire les transcrit_id du gène sélectionné
    for gene_id in list_id_genes:
        transcrits_genes_df = gtf_chr[(gtf_chr['feature'] == 'transcript') & (gtf_chr['attribute'].str.contains(gene_id) == True)]
        transcrit_id = transcrits_genes_df['attribute'].str.extract('transcript_id "(.*?)";', expand=False)
        genes_transcrits[gene_id] = transcrit_id.tolist()
    return genes_transcrits


def Extraire_transcrit_cds(genes_transcrits, gtf_chr):
    """Extraire les positions des CDS d'un transcrit à partir des informations d'annotation du gène"""
    # Dictionnaire pour stocker transcrit_id et ses CDS
    transcrit_cds = {}
    for gene_id in genes_transcrits.keys():
        for transcrit_id in genes_transcrits[gene_id]:
            # Sélectionner les lignes du dataframe correspondant au transcrit spécifié
            transcrits_df = gtf_chr[gtf_chr["attribute"].str.contains(f"transcript_id \"{transcrit_id}\"")]
            # Parcourir ligne par ligne le nouveau dataframe généré
            for i, row in transcrits_df.iterrows():
                if row["feature"] == 'CDS':
                    # Récupérer la position de start, stop et l'orientation
                    cds_start = row['start']
                    cds_end = row['end']
                    cds_strand = row['strand']
                    if transcrit_id not in transcrit_cds:
                        # Ajouter un nouveau transcrit au dictionnaire si nécessaire
                        transcrit_cds[transcrit_id] = {'cds': [(cds_start, cds_end, cds_strand)]}
                    else:
                        transcrit_cds[transcrit_id]['cds'].append((cds_start, cds_end, cds_strand))
    return transcrit_cds

def Sequence_codante(transcrit_cds, genes_transcrit, genome_seq):
    """Extraire la séquence codante d'un transcrit"""
    # Pour chaque gène, récupérer ses transcrits
    for gene_id, transcrit_id_s in genes_transcrit.items():
        for transcrit_id in transcrit_id_s:
            if transcrit_id in transcrit_cds.keys():
                # Pour chaque transcrit, s'il existe dans le dictionnaire, récupérer la position de ses CDS
                cds_positions = transcrit_cds[transcrit_id]['cds']
                seq_codante = ''
                # Pour chaque CDS, récupérer sa séquence génomique et concaténer la séquence des CDS d'un transcrit
                for (cds_start, cds_end, cds_strand) in cds_positions:
                    # Si le CDS est sur le brin "+", sa séquence est la même que la séquence génomique
                    if cds_strand == '+':
                        seq_codante += genome_seq[cds_start-1:cds_end]
                    # Si le CDS est sur le brin "-", sa séquence est la séquence reversée et complémentaire de la séquence génomique
                    else:
                        seq_codante += Reverse_complement(genome_seq[cds_start-1:cds_end])
                    # Ajouter dans le dictionnaire, une valeur contenant la séquence codante du transcrit
                transcrit_cds[transcrit_id]['seq_codante'] = seq_codante
                # Le dictionnaire transcrit_cds adopte la structure suivante : {'transcrit_id': {'cds' : [(start,end], 'seq_codante':'ATGC...'}}
    return transcrit_cds


def Lire_table_codon(table_codon):
    """Définir la table des codons. Le fichier doit avoir un format spécifique
    où chaque ligne contient un triplet de nucléotides suivi d'un espace et d'un acide aminé correspondant"""
    # Créer un dictionnaire pour stocker le code génétique ('codon' : 'acide aminé')
    codon_table = {}
    with open(table_codon, "r") as file:
        for line in file:
            line = line.strip().split()
            codon = line[0]
            amino_acid = line[1]
            #Le code génétique pour ch
            codon_table[codon] = amino_acid
    return codon_table


def Traduire_sequence(transcrit_cds, genes_transcrits, codon_table, nom_fichier_traduction):
    """Traduire la séquence nucléique en une séquence d'acides aminées.
    Sauvegarder dans un fichier multifasta"""
    # Créer un dictionnaire pour stocker la séquence protéique
    genes_traduction = {}
    with open(nom_fichier_traduction, 'w') as f:
        #Pour chaque gene_id dans le dictionnaire
        for gene_id, transcrit_id_s in genes_transcrits.items():
            genes_traduction[gene_id] = {}
            for transcrit_id in transcrit_id_s:
                if transcrit_id in transcrit_cds.keys():
                    proteine_seq = ""
                    cds_seq = transcrit_cds[transcrit_id]['seq_codante']
                    #Lire successivement la séquence chaque 3 nucléotides et attribue un acide aminé correspondant
                    for i in range(0, len(cds_seq), 3):
                        codon = cds_seq[i:i + 3]
                        if codon in codon_table:
                            proteine_seq += codon_table[codon]
                        else:
                            proteine_seq += "X"
                    # Écrire la séquence traduite sous format fasta
                    f.write(f">{transcrit_id}:{gene_id}\n{proteine_seq}\n")
                    # Ajouter dans le dictionnaire la séquence protéique pour chaque transcrit du gène
                    genes_traduction[gene_id][transcrit_id] = proteine_seq
                    #Lee dictionnaire genes_traduction a la structure suivante : {'gene_id':{'transcrit_id':'MLQT...'}}
    return genes_traduction