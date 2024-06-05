import os
import sys
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map
import scanpy as sc
import pickle
import re
import csv
import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp
from gseapy import barplot, dotplot, enrichment_map, Msigdb
from gseapy.plot import gseaplot2
from sklearn.cluster import AgglomerativeClustering
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram, linkage
import networkx as nx
import seaborn as sns
from upsetplot import UpSet, generate_counts
import time
import tkinter as tk
from tkinter import ttk

# Récupération de la base de données MSigDB
msig = Msigdb()
gmt = msig.get_gmt(category = "c6.all")
class DataPreprocessor:
    def __init__(self, fichier_compt, fichier_design):
        self.fichier_compt = fichier_compt
        self.fichier_design = fichier_design

    def load_data_compt(self):
        # Load compt data
        fichier_compt = self.fichier_compt

        if os.path.exists(fichier_compt):
            print("Fichier compt chargé")
            return pd.read_csv(fichier_compt, sep="\t", index_col=0)
        else:
            raise FileNotFoundError(fichier_compt)

    def load_data_design(self):
        # Load design data
        fichier_design = self.fichier_design

        if os.path.exists(fichier_design):
            print("Fichier design chargé")
            return pd.read_csv(fichier_design, sep="\t", index_col="SampleName")
        else:
            raise FileNotFoundError(fichier_design)

def create_dds(fichier_compt,fichier_design):
    # Créer et prétraiter les données
    preprocessor = DataPreprocessor(fichier_compt, fichier_design)
    compt_data = preprocessor.load_data_compt()
    compt_data = compt_data[compt_data.sum(axis=1) > 0].T # On sélectionne les gènes qui ont plus d'une apparition, et on transverse le tout
    design_data = preprocessor.load_data_design()
    design_data.index = compt_data.index # On matche les deux index pour l'objet dds
    # Créer l'objet DeseqDataSet
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=compt_data,
        metadata=design_data,
        design_factors="Condition",
        refit_cooks=True,
        inference=inference,
    ) # Se rééférer à la documentation de Gseapy
    dds.fit_size_factors()
    dds.obsm["size_factors"]
    dds.fit_genewise_dispersions()
    dds.varm["genewise_dispersions"]
    dds.fit_dispersion_trend()
    dds.uns["trend_coeffs"]
    dds.varm["fitted_dispersions"]
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.varm["MAP_dispersions"]
    dds.varm["dispersions"]
    dds.fit_LFC()
    dds.varm["LFC"]
    dds.calculate_cooks()
    if dds.refit_cooks:
        # Replace outlier counts
        dds.refit()
    return dds

    # Tous les paramètres ont été faits à la main, car cela fonctionnait mieux, documentation Gseapy
########################################################################################

def create_volcano_plot(results_df, ligne, t,de_dir, xlim=None):
    volcano_plot_dir = os.path.join(de_dir, "volcano_plot")
    if not os.path.exists(volcano_plot_dir):
        os.makedirs(volcano_plot_dir)
    # Définir les seuils pour l'expression différentielle
    fold_change_threshold = 1.5 # Valeur arbitraire, à modifier si besoin
    p_value_threshold = 0.05

    # Filtrer les gènes significativement régulés
    sig_genes_up = results_df[(results_df['log2FoldChange'] > fold_change_threshold) & (results_df['padj'] < p_value_threshold)]
    sig_genes_down = results_df[(results_df['log2FoldChange'] < -fold_change_threshold) & (results_df['padj'] < p_value_threshold)]

    # Trier les gènes par -log10(padj) pour récupérer les gènes les plus significatifs
    top_sig_genes_up = sig_genes_up.sort_values('-log10(padj)',ascending = False)
    top_sig_genes_down = sig_genes_down.sort_values('-log10(padj)',ascending = False)

    # Créer le volcano plot
    plt.figure(figsize=(10, 6))
    plt.scatter(results_df['log2FoldChange'], (results_df['-log10(padj)']), color='grey', alpha=0.5, label='Non-significant')
    plt.scatter(sig_genes_up['log2FoldChange'], (sig_genes_up['-log10(padj)']), color='green', label='Over-expressed')
    plt.scatter(sig_genes_down['log2FoldChange'], (sig_genes_down['-log10(padj)']), color='red', label='Under-expressed')


    plt.title(f"Volcano_plot_{ligne}_at_{t}_VS_{ligne}_ctrl")
    plt.xlabel('log2 Fold Change')
    plt.ylabel('-log10 p-value')
    plt.axvline(x=fold_change_threshold, color='black', linestyle='--', linewidth=1)
    plt.axvline(x=-fold_change_threshold, color='black', linestyle='--', linewidth=1)
    plt.axhline(y=-np.log10(p_value_threshold), color='black', linestyle='--', linewidth=1)

    if xlim is not None:
        plt.xlim(xlim) # Permet de créer un csv avec une échelle des abscisses arbitraire

    # Ajouter les annotations pour les gènes les plus significatifs
    for index, row in top_sig_genes_up.iterrows():
        plt.annotate(index, (row['log2FoldChange'],(row['-log10(padj)'])), textcoords="offset points", xytext=(0,5), ha='center', fontsize=8, color='green')
    for index, row in top_sig_genes_down.iterrows():
        plt.annotate(index, (row['log2FoldChange'],(row['-log10(padj)'])), textcoords="offset points", xytext=(0,5), ha='center', fontsize=8, color='red')

    # Enregistrer le plot
    filename = f"Volcano_plot_{ligne}_at_{t}_VS_{ligne}_ctrl"
    if xlim is not None:
        filename += "_adjusted_axis"
    filename += ".png"
    filepath = os.path.join(volcano_plot_dir, filename)
    plt.savefig(filepath)
    print(f"Volcano Plot enregistré : {filepath}")
    plt.close()

######################################################################################

def create_volcano_plot_gsea(results_df, ligne, t,gsea_dir):
# Léger problème ici, la fonction GSEA retourne des p-value = 0, ainsi les termes en questions ne sont pas affichés sur le plot, donc je leur attribue une valeur aléatoire pour les faire apparaitre si p-value = 0
# Donc une fois ce problème de p-valeu fixé, il n'y aura aucun problème pour les plots.
# Même philosophie que la fonction précédente

    # Définir les seuils pour l'expression différentielle
    fold_change_threshold = 0
    p_value_threshold = 0.05

    # Filtrer les gènes significativement régulés
    sig_pathway_up = results_df[(results_df['nes'] > fold_change_threshold) & (results_df['p_value'] < p_value_threshold)]
    sig_pathway_down = results_df[(results_df['nes'] < fold_change_threshold) & (results_df['p_value'] < p_value_threshold)]

    # Trier les termes par NES pour récupérer les 5 termes les plus significatifs
    top_sig_pathway_up = sig_pathway_up.nlargest(5, '-log10(pvalue)')
    top_sig_pathway_down = sig_pathway_down.nlargest(5, '-log10(pvalue)')

    # Créer deux dictionnaires séparés pour les termes sur-exprimés et sous-exprimés
    term_to_number_up = {row["Term"]: i + 1 for i, (_, row) in enumerate(top_sig_pathway_up.iterrows())}
    term_to_number_down = {row["Term"]: i + 1 for i, (_, row) in enumerate(top_sig_pathway_down.iterrows())}

    # Créer le volcano plot
    plt.figure(figsize=(10, 6))
    plt.scatter(results_df['nes'], (results_df['-log10(pvalue)']), color='grey', alpha=0.5, label='Non-significant')
    plt.scatter(sig_pathway_up['nes'], (sig_pathway_up['-log10(pvalue)']), color='green', label='Over-expressed')
    plt.scatter(sig_pathway_down['nes'], (sig_pathway_down['-log10(pvalue)']), color='red', label='Under-expressed')
    plt.title(f"Volcano_plot_{ligne}_at_{t}_VS_{ligne}_ctrl")
    plt.xlabel('NES')
    plt.ylabel('-log10 p-value')
    plt.axvline(x=fold_change_threshold, color='black', linestyle='--', linewidth=1)
    plt.axvline(x=-fold_change_threshold, color='black', linestyle='--', linewidth=1)
    plt.axhline(y=-np.log10(p_value_threshold), color='black', linestyle='--', linewidth=1)

    # Ajouter les annotations pour les termes les plus significatifs
    for term, number in term_to_number_up.items():
        row = top_sig_pathway_up.loc[top_sig_pathway_up['Term'] == term]
        plt.annotate(number, (row['nes'].iloc[0], (row['-log10(pvalue)'].iloc[0])), textcoords="offset points", xytext=(0,5), ha='center', fontsize=8, color='green')

    for term, number in term_to_number_down.items():
        row = top_sig_pathway_down.loc[top_sig_pathway_down['Term'] == term]
        plt.annotate(number, (row['nes'].iloc[0],(row['-log10(pvalue)'].iloc[0])), textcoords="offset points", xytext=(0,5), ha='center', fontsize=8, color='red')

    # Créer une légende avec les chiffres associés aux termes
    handles_up = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label=f"{number}: {term}") for term, number in term_to_number_up.items()]
    handles_down = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label=f"{number}: {term}") for term, number in term_to_number_down.items()]
    handles_legend = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='Over-expressed'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Under-expressed')]
    plt.legend(handles=handles_up + handles_down + handles_legend, loc="lower right", bbox_to_anchor=(2.1, 0))
    # Enregistrer le plot
    filename = f"Volcano_plot_{ligne}_at_{t}_VS_{ligne}_ctrl.png"
    filepath = os.path.join(gsea_dir, filename)
    plt.savefig(filepath, bbox_inches='tight')
    print(f"Volcano Plot enregistré : {filepath}")
    plt.close()
######################################################################################

def enr_function(df_expressed,background,ligne,t,enrch_dir,mode):
    df_expressed = df_expressed.reset_index()
    df_expressed["Gene"] = df_expressed["Gene"].str.upper() # Les gènes sont en majuscules dans les BDD
    list_genes = df_expressed["Gene"].tolist() # On a seulement besoin de la liste de gènes
    liste_background = background.index.str.upper()
    liste_background.tolist() # Le background est l'ensemble des gènes reconnus par le séquencage, dans la condition d'intérêt, TRES IMPORTANT
    df_expressed.drop_duplicates(subset ="Gene",inplace = True)
    if mode == "time":
        enr_res = gp.enrichr(gene_list =list_genes,background=liste_background, gene_sets = ["GO_Biological_Process_2023", "KEGG_2019_Mouse"], organism="Mouse")
    else:
        enr_res = gp.enrich(gene_list =list_genes,background=liste_background, gene_sets = gmt)

    try:
        ax = dotplot(enr_res.res2d, title=f"{ligne}_{t}_VS_{ligne}_ctrl", ofname=os.path.join(enrch_dir, f"{ligne}_{t}_VS_{ligne}_ctrl_ORA.png"),top_term = 20,size = 3, cmap = "viridis")
        enr_res.res2d = enr_res.res2d.sort_values(by='Adjusted P-value')
        bx = barplot(enr_res.res2d, title=f"{ligne}_{t}_VS_{ligne}_ctrl", ofname=os.path.join(enrch_dir, f"{ligne}_{t}_VS_{ligne}_ctrl_ORA_barplot.png"),top_term = 20)
    except Exception as e:
        print("No pathway found")
        print(e)
        pass
    enr_res.res2d = enr_res.res2d.set_index("Term")
    save_results_to_csv_df(ligne, t, enrch_dir, enr_res.res2d)
######################################################################################


def gsea_function(df,ligne,t,gsea_dir,mode):
    no_pathway = False # Permet de vérifier si des pathways ont été trouvés
    df = df.reset_index()
    df.sort_values(by ="Rank",ascending = False,inplace = True) # On utilise une valeur "Rank" qui prend en compte la p-value et le L2FC

    df_rank = df.loc[:,["Gene","Rank"]].reset_index(drop = True) # On ne garde que les Genes et le Rank
    df_rank["Gene"] = df_rank["Gene"].str.upper()
    df_rank = df_rank.set_index("Gene")
    try :
        if mode == "time":
            gsea_res = gp.prerank(rnk = df_rank, gene_sets = ["GO_Biological_Process_2023", "KEGG_2019_Mouse"],permutation_num = 1000,organism="Mouse")
        else:
            gsea_res = gp.prerank(rnk = df_rank, gene_sets =gmt,permutation_num = 1000, min_size=3, max_size=1000)
        out_list = []
        for term in list(gsea_res.results):
            out_list.append([term,
                            gsea_res.results[term]["pval"],
                            gsea_res.results[term]["fdr"],
                            gsea_res.results[term]["nes"],
                            gsea_res.results[term]["es"],
                            gsea_res.results[term]["tag %"],
                            gsea_res.results[term]["gene %"],
                            gsea_res.results[term]["lead_genes"]])
        df_out = pd.DataFrame(out_list,columns = ["Term","p_value","fdr","nes","es","tag %", "gene %", "Genes"]).sort_values('p_value').reset_index(drop = True)
        non_zero_p_values = df_out[df_out["p_value"] != 0]["p_value"]
        best_p_value = np.min(non_zero_p_values)
        range_best_p_value = best_p_value - best_p_value * 0.05
        df_out.loc[df_out['p_value'] == 0, 'p_value'] = np.random.uniform(range_best_p_value, best_p_value, df_out['p_value'].eq(0).sum()) # ici, la fonction retourne des p-value à 0 pour certains termes, et cela pose problème pour les volcano-plots. Ainsi on leur attribue une valeur si p-value = 0
        # Lorsque cela sera fixé, alors ce bout de code ne sera plus utile
        df_out["-log10(pvalue)"] = -np.log10(df_out["p_value"])
        df_out.sort_values(by="p_value")
        create_volcano_plot_gsea(df_out, ligne, t,gsea_dir)
        terms = df_out["Term"]
        # Calculer les hits pour les 5 premiers termes
        hits = [gsea_res.results[t]['hits'] for t in terms[:5]]

        # Calculer les runes pour les 5 premiers termes
        runes = [gsea_res.results[t]['RES'] for t in terms[:5]]

        # Créer le graphique avec gseaplot2
        fig = gseaplot2(terms=terms[:5], RESs=runes, hits=hits,
                        rank_metric=gsea_res.ranking,
                        legend_kws={'loc': (1.2, 0)},
                        figsize=(8, 6), ofname=os.path.join(gsea_dir, f"{ligne}_{t}_VS_{ligne}_ctrl_GSEA.png"))

        df_out_filtrate = df_out[df_out["p_value"] < 0.05]
        save_results_to_csv_df(ligne, t, gsea_dir, df_out)
        save_results_to_csv_df(ligne, t, gsea_dir, df_out_filtrate, text="filtered")
    except Exception as e:
        print(f"Error: {e}")
        print("No pathway found")
        no_pathway = True
        pass


    # Créer un graphe vide
    if no_pathway:
        pass
    else:
        G = nx.Graph()

        # Ajouter les nœuds (termes)
        terms_subset = df_out_filtrate['Term'].head(22) # On affiche les 22 termes les plus significatifs, choix arbitraire

        # Créer un sous-ensemble de votre DataFrame contenant uniquement ces 22 termes
        df_subset = df_out_filtrate[df_out_filtrate['Term'].isin(terms_subset)]
        G.add_nodes_from(terms_subset)

        # Ajouter les arêtes basées sur les relations entre les termes
        for i in range(len(terms_subset)):
            for j in range(i+1, len(terms_subset)):
                # Vérifier si la colonne 'Genes' existe
                if 'Genes' in df_subset.columns:
                    # Vérifier s'il n'y a qu'un seul gène dans la ligne
                    if ';' not in df_subset.loc[i, 'Genes'] and ';' not in df_subset.loc[j, 'Genes']:
                        gene_i = df_subset.loc[i, 'Genes']
                        gene_j = df_subset.loc[j, 'Genes']
                        if gene_i == gene_j:
                            G.add_edge(terms_subset[i], terms_subset[j], genes={gene_i})
                    else:
                        # Séparer les gènes par le point-virgule
                        genes_i = set(df_subset.loc[i, 'Genes'].split(';'))
                        genes_j = set(df_subset.loc[j, 'Genes'].split(';'))
                        common_genes = genes_i.intersection(genes_j)
                        # Vérifier s'il y a des gènes communs
                        if common_genes:
                            maximum = 6 # Car sinon les traits sont trop épais ensuite
                            # Calculer l'épaisseur du trait en fonction du nombre de gènes communs
                            edge_width = 0.5 + (len(common_genes) - 1) * 2
                            if edge_width > maximum:
                                edge_width=maximum
                            G.add_edge(terms_subset[i], terms_subset[j], genes=common_genes, width=edge_width)
                else:
                    print("La colonne 'Genes' n'existe pas dans le DataFrame.")
        gene_count_colors = {}
        for node in G.nodes():
            edge_genes = [data['genes'] for _, _, data in G.edges(node, data=True)]
            total_gene_count = sum(len(genes) for genes in edge_genes)
            if total_gene_count > 9:
                gene_count_colors[node] = 'red'  # Couleur rouge pour plus de 5 gènes
            elif total_gene_count > 6:
                gene_count_colors[node] = 'orange'  # Couleur orange pour 4 à 5 gènes
            elif total_gene_count > 3:
                gene_count_colors[node] = 'yellow'
            else:
                gene_count_colors[node] = 'skyblue'  # Couleur bleu ciel pour moins de 3 gènes


        # Tracer le graphe
        plt.figure(figsize=(16, 10))
        pos = nx.spring_layout(G, k=2.5)  # Positionne les nœuds en utilisant un algorithme de mise en page spring
        if all('width' in data for _, _, data in G.edges(data=True)):
            # Si toutes les arêtes ont l'attribut 'width' défini
            nx.draw(G, pos, with_labels=True, node_size=500, node_color=[gene_count_colors[node] for node in G.nodes()], font_size=10, width=[data['width'] for _, _, data in G.edges(data=True)])
        else:
            # Si certaines arêtes n'ont pas 'width' défini, utiliser une largeur par défaut
            nx.draw(G, pos, with_labels=True, node_size=1000, node_color=[gene_count_colors[node] for node in G.nodes()], font_size=10)
        legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='More than 9 genes'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=10, label='From 7 to 9 genes'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='yellow', markersize=10, label='From 4 to 6 genes'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='skyblue', markersize=10, label='Less than 4 genes')
    ]
        plt.legend(handles=legend_elements, loc='upper right')
        plt.title("Graph of relationships between terms")
        plt.savefig(os.path.join(gsea_dir, f"{ligne}_{t}_VS_{ligne}_ctrl_GSEA_graph.png"))
        plt.close()

######################################################################################

def dictionary_to_csv(de_dir, filename, dict_genes):
    with open(os.path.join(de_dir, filename), mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(["Condition", "Gene", "p_value", "Log2FoldChange"])  # En-têtes des colonnes
        # Parcours de chaque condition dans le dictionnaire
        for condition, gene_info_list in dict_genes.items():
            # Triage des gènes par padj
            sorted_gene_info = sorted(gene_info_list, key=lambda x: x[1])
            # Utilisation d'un ensemble pour stocker les gènes déjà écrits
            written_genes = set()
            # Parcours de chaque triplet (gene, p_value, log2FoldChange) trié par padj
            for gene_info in sorted_gene_info:
                gene, p_value, log2FoldChange = gene_info
                # Vérification si le gène n'a pas déjà été écrit
                if gene not in written_genes:
                    writer.writerow([condition, gene, p_value, log2FoldChange])
                    # Ajout du gène à l'ensemble des gènes déjà écrits
                    written_genes.add(gene)
    print(f"Fichier CSV enregistré : {filename}")



########################################################################################


def save_results_to_csv(ligne, t, de_dir, stat_res):
    stat_res.results_df.sort_values(by='padj', inplace=True)
    filename = f"{ligne}_{t}_VS_{ligne}_ctrl.csv"
    filepath = os.path.join(de_dir, filename)
    with open(filepath, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(['Gene'] + list(stat_res.results_df.columns))
        # Écrire les données
        for index, row in stat_res.results_df.iterrows():
            writer.writerow([index] + list(row))
    print(f"Fichier CSV enregistré : {filepath}")

########################################################################################

def create_upset_plot(dict,nb,dossier,name,de_dir):
    all_elements = set()
    for elements in dict.values():
        all_elements.update(elements)

    # Créer une liste de dictionnaires pour stocker les comptages de chaque élément dans chaque ensemble
    counts = []
    for ensemble, elements in dict.items():
        counts.append({element: elements.count(element) for element in all_elements})

    # Créer une DataFrame à partir de la liste de dictionnaires
    df_upset = pd.DataFrame(counts, index=dict.keys())

    # Trier les colonnes par ordre alphabétique
    df_upset = df_upset.reindex(sorted(df_upset.columns), axis=1)
    df_bool = df_upset.astype(bool)
    df_bool = df_bool.T
    # Créer un MultiIndex à partir de la liste de tuples
    multi_index = pd.MultiIndex.from_tuples(df_bool.apply(tuple, axis=1),names=df_bool.columns)
    # Créer une série avec un MultiIndex
    series = pd.Series(index=multi_index)
    upset = UpSet(series, show_counts=True, subset_size="count",min_subset_size=nb)
    upset.plot()
    plt.savefig(os.path.join(dossier, name))
    plt.close()
########################################################################################

def save_results_to_csv_df(ligne, t, dir, df, text=""):
    filename = f"{ligne}_{t}_VS_{ligne}_ctrl_{text}.csv"
    filepath = os.path.join(dir, filename)
    with open(filepath, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        writer.writerow(['Pathway'] + list(df.columns))
        for index, row in df.iterrows():
            writer.writerow([index] + list(row))
    print(f"Fichier CSV enregistré : {filepath}")

########################################################################################


def create_box_plot(df, nb, ligne, t,de_dir):
    best_fold = df.loc[df['Rank'].sort_values(ascending=False).index] # On trie les gènes par Rank
    best_fold = best_fold.head(nb)
    mean_base_mean = df['baseMean'].mean()
    median_base_mean = df['baseMean'].median()
    plt.figure(figsize=(12, 8))
    sns.boxplot(y=best_fold.index, x="baseMean", data=best_fold, palette='Set2', whis=1.5, showfliers=False)
    for index, row in best_fold.iterrows():
        plt.text(row['baseMean'], index, f"LFC: {row['log2FoldChange']:.2f}", color='black', ha='center', va='center')
    plt.axvline(x=mean_base_mean, color='r', linestyle='--', label=f'Moyenne baseMean: {mean_base_mean:.2f}')
    plt.axvline(x=median_base_mean, color='b', linestyle='--', label=f'Médiane baseMean: {median_base_mean:.2f}')
    plt.legend()
    plt.title(f'Boxplots per gene for line {ligne} at {t} VS {ligne}_ctrl')
    plt.xlabel('baseMean')
    plt.ylabel('Gène')
    plt.tight_layout()  # Optimiser le placement des sous-graphiques
    plt.savefig(os.path.join(de_dir, f"Boxplot_{ligne}_at_{t}_VS_{ligne}_ctrl.png"))
    plt.close()


########################################################################################

def box_plot_from_dds(dds,ligne, temps,sig_genes_all, genes_of_interest,de_dir):
    boxplot_dir = os.path.join(de_dir, "boxplot")
    if not os.path.exists(boxplot_dir):
        os.makedirs(boxplot_dir)
    if isinstance(genes_of_interest,int): # Si l'utilisateur a rentré un entier, alors on prend les n premiers gènes selon le Rank
        sig_genes_all = sig_genes_all.loc[sig_genes_all['Rank'].sort_values(ascending=False).index]
        genes_index = sig_genes_all.index.tolist()
        genes_index = genes_index[:80]
        unique_genes = [gene for i, gene in enumerate(genes_index) if gene not in genes_index[:i]]
        list_genes = unique_genes[:genes_of_interest]
    else:
        list_genes = genes_of_interest # Sinon on prend les gènes d'intérêt rentrés par l'utilisateur
        print(f"Liste des gènes d'intérêt : {list_genes}")
    for gene_interest in list_genes:
        all_gene_data = pd.DataFrame()
        condition_data = []
        subset_data_ctrl = dds.obs['Condition'] == f"{ligne}-ctrl"
        gene_expression_ctrl = dds[subset_data_ctrl, :].X[:, dds.var_names == gene_interest] # On récupère l'expression du gène d'intérêt pour la condition contrôle
        if gene_expression_ctrl.size == 0:
            print(f"Le gène {gene_interest} n'est pas présent dans la lignée {ligne}.") # Si on ne trouve pas alors le gène n'est pas présent dans la lignée
            continue
        sample_names_ctrl = dds.obs.index[subset_data_ctrl] # On récupère les noms des échantillons
        temp_data_ctrl = pd.Series(gene_expression_ctrl.flatten(), index=sample_names_ctrl, name=f"{ligne}-ctrl") # On crée une série avec les données
        all_gene_data = pd.concat([all_gene_data, temp_data_ctrl], axis=1) # On concatène les données
        condition_data.append(gene_expression_ctrl.flatten()) # On ajoute les données à la liste
        for t in temps: # On fait la même chose pour les autres conditions, selon le temps, dans l'ordre croissant
            condition_name = f"{ligne}-{t}"
            subset_data = dds.obs['Condition'] == condition_name
            gene_expression = dds[subset_data, :].X[:, dds.var_names == gene_interest]
            sample_names = dds.obs.index[subset_data]
            temp_data = pd.Series(gene_expression.flatten(), index=sample_names, name=condition_name)
            all_gene_data = pd.concat([all_gene_data, temp_data], axis=1)
            condition_data.append(gene_expression.flatten())
        all_gene_data.fillna(np.nan, inplace=True)
        csv_filename = os.path.join(boxplot_dir, f"{ligne}_{gene_interest}_expression.csv")
        all_gene_data.to_csv(csv_filename)
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=condition_data, palette='Set2', whis=1.5, showfliers=False)
        sns.stripplot(data=condition_data, color='black', size=4)
        plt.title(f'Boxplot for the gene {gene_interest}  under different conditions')
        plt.xlabel("Conditions")
        plt.ylabel("Expression values")
        plt.xticks(range(len(temps) + 1), [f"{ligne}-ctrl"] + [f"{ligne}-{t}" for t in temps], rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(boxplot_dir, f"Boxplot_{ligne}_{gene_interest}.png"))
        plt.close()


########################################################################################

def pca_from_dds(dds,ligne,de_dir):
    pca = PCA(.90)
    pcs = pca.fit_transform(dds.X)
    condition_name = f"{ligne}"
    subset_data = dds[dds.obs['Condition'].str.contains(condition_name), :]
    pcs_filtered = pcs[dds.obs['Condition'].str.contains(condition_name), :]
    labels_filtered = subset_data.obs['Condition']
    plt.figure(figsize=(10, 6))
    sns.scatterplot(x=pcs_filtered[:, 0], y=pcs_filtered[:, 1], hue=labels_filtered, palette="viridis", s=400)
    plt.title(f"PCA for elements of condition {condition_name}")
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.2%})')
    plt.ylabel('PC2 ({:.2%})'.format(pca.explained_variance_ratio_[1]))
    plt.savefig(os.path.join(de_dir, f"PCA_{ligne}.png"))
    plt.close()

########################################################################################

def process_differential_analysis(dds, lignees, temps, output_dir,mode):
    de_dir = os.path.join(output_dir, "DE")
    gsea_dir = os.path.join(output_dir, "GSEA")
    enrch_dir = os.path.join(output_dir, "ORA")
    dict_genes = {}
    dict_genes_all = {}
    dict_intersections = {}
    condition_name = []
    expressed_genes = []
    not_expressed_genes = []
    for ligne in lignees:
        sig_genes_df = [] # Créer une liste vide pour stocker les DataFrames des gènes significatifs
        genes_ligne = []  # Créer une liste vide pour stocker les gènes significatifs pour cette ligne
        for t in temps:
            print(f"Comparaison pour la lignée {ligne} entre la condition {t} et contrôle ")
            # Définir les contrastes pour la comparaison
            stat_res = DeseqStats(dds, alpha=0.05, cooks_filter=True, contrast=("Condition",f"{ligne}-{t}", f"{ligne}-ctrl"), independent_filter=True) # Comparaison entre la condition et le contrôle
            # Exécuter le test statistique
            stat_res.summary()
            stat_res.p_values
            if stat_res.independent_filter:
                stat_res._independent_filtering()
            else:
                stat_res._p_value_adjustment()

            stat_res.padj
            # Réinitialiser l'index pour faire de 'Gene' une colonne
            stat_res.results_df.reset_index(inplace=True)

            # Renommer la colonne 'index' en 'Gene'

            stat_res.results_df.rename(columns={'index': 'Gene'}, inplace=True)
            expr_genes = stat_res.results_df[stat_res.results_df['padj'] < 0.05]["Gene"].count() # Compter le nombre de gènes exprimés significativement différemment
            not_expr_genes = stat_res.results_df[stat_res.results_df['padj'] >= 0.05]["Gene"].count() # Compter le nombre de gènes non exprimés significativement différemment
            expressed_genes.append(expr_genes) # Ajouter le nombre de gènes exprimés à la liste
            not_expressed_genes.append(not_expr_genes) # Ajouter le nombre de gènes non exprimés à la liste
            condition_name.append(f"{ligne}_{t}") # Ajouter le nom de la condition à la liste
            stat_res.results_df.set_index('Gene', inplace=True)
            # Ajout de nouvelles colonnes

            stat_res.results_df["Rank"] = -np.log10(stat_res.results_df["padj"]) *stat_res.results_df["log2FoldChange"]
            stat_res.results_df["-log10(padj)"] = -np.log10(stat_res.results_df["padj"])
            # Filtrer les gènes significatifs

            sig_genes = stat_res.results_df[(stat_res.results_df['padj'] < 0.05)]
            sig_genes_df.append(sig_genes)
            dict_genes[f"{ligne}_{t}"] = sig_genes.index.tolist() # Ajouter les gènes significatifs de la condition (clé) en tant que liste (valeurs)
            gene_info = [(gene, p_value, log2FoldChange) for gene, p_value, log2FoldChange
                         in zip(sig_genes.index.tolist(), sig_genes["padj"], sig_genes["log2FoldChange"])]

            # Ajoutez cette liste au dictionnaire avec la clé correspondant à la condition
            dict_intersections[f"{ligne}_{t}"] = gene_info
            genes_ligne.extend(gene_info)
            create_box_plot(sig_genes,30, ligne, t,de_dir)
            stat_res.results_df = stat_res.results_df.dropna()
            save_results_to_csv(ligne, t, de_dir, stat_res)
            create_volcano_plot(stat_res.results_df, ligne, t,de_dir = de_dir)
            create_volcano_plot(stat_res.results_df, ligne, t,de_dir = de_dir, xlim=(-9,9))
            enr_function(sig_genes,stat_res.results_df,ligne,t,enrch_dir,mode)

            stat_res.results_df = stat_res.results_df[(stat_res.results_df['padj'] < 0.90)]
            gsea_function(stat_res.results_df,ligne,t,gsea_dir,mode)

        dict_genes_all[f"{ligne}"] = genes_ligne
        sig_genes = pd.DataFrame()
        sig_genes_all = pd.concat(sig_genes_df)
        pca_from_dds(dds,ligne,de_dir)
    create_upset_plot(dict_genes,8,de_dir,"upset_plot_all.png",de_dir)
    create_upset_plot(dict_genes_all,0,de_dir,"upset_plot_lignees.png",de_dir)

    dictionary_to_csv(de_dir, "all_genes_significatifs.csv", dict_intersections)
    dictionary_to_csv(de_dir, "all_genes_significatifs_all_lignes.csv", dict_genes_all)
    create_bar_chart(de_dir, expressed_genes, not_expressed_genes, condition_name)

########################################################################################

def create_bar_chart(de_dir,expr_genes,non_expr_genes,condition_name):
    # Créer le diagramme à barres
    plt.figure(figsize=(10, 6))
    x = range(len(condition_name))
    plt.bar(x, expr_genes, color='blue', label='Expressed significantly differently')
    plt.bar(x, non_expr_genes, color='grey', bottom=expr_genes, label='Not expressed significantly differently')
    plt.title("Number of genes expressed differently for each condition")
    plt.xlabel('Condition')
    plt.ylabel('Number of genes')
    plt.xticks(x, condition_name, rotation=60, ha = 'right')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(de_dir, "expr_genes_all_conditions.png"))
    plt.close()

############################################################

class ConditionSelector:
    # Interface graphique permettant de sélectionner les conditions pour créer un fichier CSV (lors du second appel du programme)
    def __init__(self, master, dict_genes, de_dir):
        self.master = master
        self.dict_genes = dict_genes
        self.de_dir = de_dir
        self.selected_conditions = []  # Liste pour stocker les conditions sélectionnées

        self.master.title("Sélecteur de conditions")

        self.frame = ttk.Frame(self.master)
        self.frame.pack(fill=tk.BOTH, expand=True)

        # Définir les styles pour les boutons
        self.style = ttk.Style()
        self.style.configure("TButton", background="#f0f0f0")
        self.style.map("TButton", background=[('active', '#d9d9d9')])

        self.style.configure("Selected.TButton", background="#add8e6")
        self.style.map("Selected.TButton", background=[('active', '#87ceeb')])

        self.condition_buttons = {}
        for i, condition in enumerate(self.dict_genes.keys()):
            button = ttk.Button(self.frame, text=condition, command=lambda cond=condition: self.toggle_condition(cond))
            button.grid(row=i, column=0, sticky="w")
            self.condition_buttons[condition] = button

        ttk.Label(self.frame, text="Nom du fichier CSV :").grid(row=len(self.dict_genes), column=0, sticky="w")
        self.filename_entry = ttk.Entry(self.frame)
        self.filename_entry.grid(row=len(self.dict_genes), column=1, pady=5)

        self.launch_button = ttk.Button(self.frame, text="Lancer", command=self.launch)
        self.launch_button.grid(row=len(self.dict_genes)+1, column=0, columnspan=2, pady=5)
        self.finish_button = ttk.Button(self.frame, text="Fin", command=self.finish)
        self.finish_button.grid(row=len(self.dict_genes)+2, column=0, columnspan=2, pady=5)
        self.after_ids = []
    def toggle_condition(self, condition):
        button = self.condition_buttons[condition]
        if condition in self.selected_conditions:
            self.selected_conditions.remove(condition)  # Déselectionne la condition si elle est déjà sélectionnée
            button.configure(style="TButton")  # Réinitialiser le style par défaut
        else:
            self.selected_conditions.append(condition)  # Sélectionne la condition si elle n'est pas déjà sélectionnée
            button.configure(style="Selected.TButton")  # Appliquer le style sélectionné
    def cancel_after_callbacks(self):
        for after_id in self.after_ids:
            self.master.after_cancel(after_id)
        self.after_ids.clear()

    def finish(self):
        self.master.destroy()

    def launch(self):
        filename = self.filename_entry.get() + ".csv"

        if not filename:
            tk.messagebox.showerror("Erreur", "Veuillez entrer un nom de fichier CSV.")
            return

        if not self.selected_conditions:
            tk.messagebox.showerror("Erreur", "Veuillez sélectionner au moins une condition.")
            return

        print("Conditions sélectionnées :", self.selected_conditions)
        self.write_to_csv(filename, self.selected_conditions)
        self.selected_conditions = []
        self.filename_entry.delete(0, tk.END)
        self.cancel_after_callbacks()
    # Réinitialiser la liste des conditions sélectionnées
    def write_to_csv(self, filename, selected_conditions):
        csv_path = os.path.join(self.de_dir,filename)
        with open(csv_path, mode='w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')

            # Écrire l'en-tête avec les noms des conditions, p-values et Log2FC
            header = ["Gene"]
            for condition in selected_conditions:
                header.extend([f"{condition}_p-value", f"{condition}_Log2FC"])
            writer.writerow(header)

            # Récupérer tous les gènes présents dans la première condition sélectionnée
            all_genes = set(gene_info[0] for gene_info in self.dict_genes[selected_conditions[0]])

            # Parcourir tous les gènes pour vérifier leur présence dans les autres conditions
            for gene in all_genes:
                if all(any(gene == info[0] for info in self.dict_genes[condition]) for condition in selected_conditions[1:]):
                    row = [gene]
                    for condition in selected_conditions:
                        # Trouver les informations sur le gène pour chaque condition
                        gene_info = next(info for info in self.dict_genes[condition] if info[0] == gene)
                        row.extend([gene_info[1], gene_info[2]])
                    writer.writerow(row)
        print(f"Fichier CSV '{filename}' créé avec succès.")