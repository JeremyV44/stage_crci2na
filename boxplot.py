import os
import sys
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import customtkinter as ctk
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import pickle
import re
import csv
import numpy as np
import matplotlib.pyplot as plt
import utils
from gseapy import gseaplot2
import anndata as ad
import time

### Ce code sert à lancer l'analyse une seconde fois, après etude des premiers résultats. Il permet de sélectionner les gènes à étudier plus en détail.
### Ainsi, on doit de nouveau refaire toutes les analyses différentielles entre chaque condition et le contrôle.

def load_files():
    global fichier_pickle_entry, dossier_output_entry
    fichier_pickle = ctk.filedialog.askopenfilename(title="Select the pickle file")
    dossier_output = ctk.filedialog.askdirectory(title="Select the output directory")
    dossier_output_entry.configure(state=tk.NORMAL)
    fichier_pickle_entry.configure(state=tk.NORMAL)
    fichier_pickle_entry.delete(0, tk.END)
    dossier_output_entry.delete(0, tk.END)
    fichier_pickle_entry.insert(0, fichier_pickle)
    dossier_output_entry.insert(0, dossier_output)
    fichier_pickle_entry.configure(state=tk.DISABLED)
    dossier_output_entry.configure(state=tk.DISABLED)

def run_analysis():
    fichier_pickle = fichier_pickle_entry.get()
    try:
        with open(fichier_pickle, 'rb') as f:
            dds = pickle.load(f)
        liste_genes = liste_genes_entry.get("1.0", tk.END).strip() # Etant donné que c'est un Textbox, on doit récupérer le texte avec get("1.0", tk.END)
        mode = mode_var.get()
        dossier_output = dossier_output_entry.get()
        main(dds, liste_genes,mode,dossier_output)
    except Exception as e:
        ctk.messagebox.showerror("Erreur", f"Erreur lors du chargement du fichier pickle ou de l'exécution de l'analyse:\n{str(e)}")
    finally:
        root.destroy()  # Fermer la fenêtre après l'exécution de l'analyse
def main(dds, liste_genes, mode, dossier_output):
    de_dir = os.path.join(dossier_output, 'DE')
    lignee = dds.obs["Condition"].unique()
    temps = []
    lignees = []

    for element in lignee:
        if "-ctrl" in element:
            continue
        if mode == "time":
            temps_match = re.search(r'(\d+h)', element,re.IGNORECASE)
            if temps_match:
                temps_str = temps_match.group(1)
            temps.append(temps_str)
            lignee_sans_temps = element.replace(temps_match.group(1), '')
            lignee_sans_dash = lignee_sans_temps.rstrip('-')
            lignees.append(lignee_sans_dash)
        elif mode == "treated":
            temps_match = re.search(r'ttt', element, re.IGNORECASE)
            temps.append(temps_match.group(0))
            lignee_sans_temps = element.replace(temps_match.group(0), '')
            lignee_sans_dash = lignee_sans_temps.rstrip('-')
            lignees.append(lignee_sans_dash)
        else:
            print("Mode non reconnu, veuillez choisir entre 'time' et 'treated'")
            sys.exit(1)

    if mode == "time":
        temps = list(set(filter(None, temps)))
        temps = sorted(temps, key=lambda x: int(x[:-1]))
    else:
        temps = list(set(temps))
    lignees = list(set(lignees))
    dict_genes = {}
    dict_genes_all = {}
    if liste_genes.isdigit():
        liste_genes = int(liste_genes)
    else:
        liste_genes = liste_genes.split(",")
    dict_intersections = {}
    for ligne in lignees:
        sig_genes_df = []
        genes_ligne = []  # Créer une liste vide pour stocker les gènes significatifs pour cette ligne
        for t in temps:
            print(f"Comparaison pour la lignée {ligne} entre la condition {t} et contrôle ")
            # Définir les contrastes pour la comparaison
            stat_res = DeseqStats(dds, alpha=0.05, cooks_filter=True, contrast=("Condition",f"{ligne}-{t}", f"{ligne}-ctrl"), independent_filter=True)
            # Exécuter le test statistique
            stat_res.summary()
            stat_res.p_values
            if stat_res.independent_filter:
                stat_res._independent_filtering()
            else:
                stat_res._p_value_adjustment()

            stat_res.padj
            # Renommer la première colonne du DataFrame pour avoir le nom "Gene"
            # Réinitialiser l'index pour faire de 'Gene' une colonne
            stat_res.results_df.reset_index(inplace=True)

            # Renommer la colonne 'index' en 'Gene'

            stat_res.results_df.rename(columns={'index': 'Gene'}, inplace=True)
            stat_res.results_df.set_index('Gene', inplace=True)
            # Ajout de nouvelles colonnes

            stat_res.results_df["Rank"] = -np.log10(stat_res.results_df["padj"]) *stat_res.results_df["log2FoldChange"]
            stat_res.results_df["-log10(padj)"] = -np.log10(stat_res.results_df["padj"])
            # Filtrer les gènes significatifs

            sig_genes = stat_res.results_df[(stat_res.results_df['padj'] < 0.05)]
            sig_genes_df.append(sig_genes)
            stat_res.results_df = stat_res.results_df.dropna()


            stat_res.results_df = stat_res.results_df[(stat_res.results_df['padj'] < 0.90)]
            gene_info = [(gene, p_value, log2FoldChange) for gene, p_value, log2FoldChange
                         in zip(sig_genes.index.tolist(), sig_genes["padj"], sig_genes["log2FoldChange"])]
            dict_intersections[f"{ligne}_{t}"] = gene_info

        sig_genes_all = pd.concat(sig_genes_df)
        utils.box_plot_from_dds(dds,ligne, temps,sig_genes_all, liste_genes,de_dir)
    while True:
        root = tk.Tk()
        app = utils.ConditionSelector(root, dict_intersections,de_dir)
        root.mainloop()
        break

root = ctk.CTk()
root.geometry("800x600")
root.title("Paramètres de l'analyse")
# Interface graphique pour lancer l'analyse

# Frame pour les paramètres de l'analyse
parameter_frame = ctk.CTkFrame(root)
parameter_frame.pack(padx=20, pady=20, expand=True, fill="both")

fichier_pickle_label = ctk.CTkLabel(parameter_frame, text="Pickle file:")
fichier_pickle_label.grid(row=1, column=0, sticky="w")

fichier_pickle_entry = ctk.CTkEntry(parameter_frame, state=tk.DISABLED)
fichier_pickle_entry.grid(row=1, column=1, padx=5, pady=5, sticky="ew")

fichier_pickle_button = ctk.CTkButton(parameter_frame, text="Parse...", command=load_files)
fichier_pickle_button.grid(row=1, column=2, padx=5, pady=5)

dossier_output_label = ctk.CTkLabel(parameter_frame, text="Name of the output folder:")
dossier_output_label.grid(row=2, column=0, sticky="w")

dossier_output_entry = ctk.CTkEntry(parameter_frame, state=tk.DISABLED)
dossier_output_entry.grid(row=2, column=1, padx=5, pady=5, sticky="ew")

dossier_output_button = ctk.CTkButton(parameter_frame, text="Parse...", command=load_files)
dossier_output_button.grid(row=2, column=2, padx=5, pady=5)

liste_genes_label = ctk.CTkLabel(parameter_frame, text="List of genes (or number) :")
liste_genes_label.grid(row=3, column=0, sticky="w")

liste_genes_entry = ctk.CTkTextbox(parameter_frame, width=300, height=100)
liste_genes_entry.grid(row=3, column=1, padx=5, pady=5, sticky="ew")

mode_label = ctk.CTkLabel(parameter_frame, text="Mode:")
mode_label.grid(row=4, column=0, sticky="w")

mode_var = ctk.StringVar()
mode_var.set("time")  # Valeur par défaut

mode_time_radio = ctk.CTkRadioButton(parameter_frame, text="Time", variable=mode_var, value="time")
mode_time_radio.grid(row=4, column=1, padx=5, pady=5, sticky="w")

mode_vaccine_radio = ctk.CTkRadioButton(parameter_frame, text="Treated", variable=mode_var, value="treated")
mode_vaccine_radio.grid(row=4, column=2, padx=5, pady=5, sticky="w")

# Bouton pour exécuter l'analyse
run_button = ctk.CTkButton(root, text="Run", command=run_analysis)
run_button.pack(padx=10, pady=10)

root.mainloop()
