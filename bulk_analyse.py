import os
import sys
import tkinter as tk
from tkinter import filedialog
import customtkinter as ctk
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
import pickle
import re
import csv
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
import utils
from gseapy import gseaplot2, Msigdb
import anndata as ad
from sklearn.decomposition import PCA
import seaborn as sns
import time
import scanpy as sc
# Fonction pour charger les fichiers de comptage et de design
def load_files():
    global fichier_compt_entry, fichier_design_entry
    fichier_compt = ctk.filedialog.askopenfilename(title="Select the count file")
    fichier_design = ctk.filedialog.askopenfilename(title="Select the design file")
    fichier_compt_entry.configure(state=tk.NORMAL)
    fichier_design_entry.configure(state=tk.NORMAL)
    fichier_compt_entry.delete(0, tk.END)
    fichier_design_entry.delete(0, tk.END)
    fichier_compt_entry.insert(0, fichier_compt)
    fichier_design_entry.insert(0, fichier_design)
    fichier_compt_entry.configure(state=tk.DISABLED)
    fichier_design_entry.configure(state=tk.DISABLED)

# Fonction pour exécuter l'analyse avec les fichiers et paramètres spécifiés
def run_analysis():
    fichier_compt = fichier_compt_entry.get()
    fichier_design = fichier_design_entry.get()
    pickle_name = pickle_entry.get() + ".pkl"
    mode = mode_var.get()
    output_dir = output_dir_entry.get() + time.strftime("_%d_%m_%Y")  # Récupérer le nom du dossier de sortie
    organism = organism_var.get()
    uniprot_number = int(uniprot_entry.get())
    dds = utils.create_dds(fichier_compt, fichier_design)
    os.makedirs(output_dir, exist_ok=True)

    # Appel de votre fonction main avec les paramètres spécifiés
    save_pickle(dds, pickle_name)
    main(dds, mode, output_dir, organism, uniprot_number)
    root.destroy() # Permet de fermer la fenêtre principale après l'exécution de l'analyse

def save_pickle(dds, pickle_name):
    with open(pickle_name, "wb") as f:
        pickle.dump(dds, f)

# Fonction principale
def main(dds, mode, output_dir, organism, uniprot_number):
    de_dir = os.path.join(output_dir, 'DE')
    os.makedirs(de_dir, exist_ok=True)
    lignee = dds.obs["Condition"].unique() # Ici on récupère toutes les conditions, puis on va chercher un extraire le nom des lignées et des temps
    temps = []
    lignees = []

    for element in lignee:
        if "-ctrl" in element:
            continue
        if mode == "time":
            temps_match = re.search(r'(\d+h)', element,re.IGNORECASE) # On cherche les temps sous la forme de chiffres suivis de h
            if temps_match:
                temps_str = temps_match.group(1) # On récupère le temps
            temps.append(temps_str) # On ajoute le temps à la liste des temps
            lignee_sans_temps = element.replace(temps_match.group(1), '') # On retire le temps de la condition
            lignee_sans_dash = lignee_sans_temps.rstrip('-') # On retire le tiret final
            lignees.append(lignee_sans_dash) # On ajoute la condition à la liste des lignées
        elif mode == "treated": # Si le mode est treated, on cherche les conditions contenant "ttt"
            temps_match = re.search(r'ttt', element, re.IGNORECASE)
            temps.append(temps_match.group(0))
            lignee_sans_temps = element.replace(temps_match.group(0), '')
            lignee_sans_dash = lignee_sans_temps.rstrip('-')
            lignees.append(lignee_sans_dash)
        else:
            print("Mode non reconnu, veuillez choisir entre 'time' et 'treated'")
            sys.exit(1)

    if mode == "time":
        temps = list(set(filter(None, temps))) # On retire les valeurs nulles et on récupère les temps uniques
        temps = sorted(temps, key=lambda x: int(x[:-1])) # On trie les temps par ordre croissant
    else:
        temps = list(set(temps))
    lignees = list(set(lignees)) # On retire les valeurs nulles et on récupère les lignées uniques
    de_dir = os.path.join(output_dir, 'DE')
    gsea_dir = os.path.join(output_dir, 'GSEA')
    enrch_dir = os.path.join(output_dir, 'ORA')
    os.makedirs(de_dir, exist_ok=True)
    os.makedirs(gsea_dir, exist_ok=True)
    os.makedirs(enrch_dir, exist_ok=True)
    # Processus d'analyse différentielle
    utils.process_differential_analysis(dds, lignees, temps, output_dir,mode, organism, uniprot_number)

    # PCA global en utilisant dds
    sc.tl.pca(dds)
    sc.pl.pca(dds, color=("Condition"), size=400,annotate_var_explained=True, show=False)
    pca_save = os.path.join(de_dir, 'pca_global.png')
    plt.savefig(pca_save, bbox_inches='tight')
    plt.close()

# Créer la fenêtre principale de l'interface utilisateur
root = ctk.CTk()
root.geometry("800x600")
root.title("Paramètres de l'analyse")

# Frame pour les paramètres de l'analyse
parameter_frame = ctk.CTkFrame(root)
parameter_frame.pack(padx=10, pady=10)

# Entrée pour spécifier le fichier de comptage
fichier_compt_label = ctk.CTkLabel(parameter_frame, text="Count file :") 
fichier_compt_label.grid(row=0, column=0, sticky="w")

fichier_compt_entry = ctk.CTkEntry(parameter_frame, state=tk.DISABLED)
fichier_compt_entry.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

fichier_compt_button = ctk.CTkButton(parameter_frame, text="Parse... ", command=load_files)
fichier_compt_button.grid(row=0, column=2, padx=5, pady=5)

# Entrée pour spécifier le fichier de design
fichier_design_label = ctk.CTkLabel(parameter_frame, text="Design file:")
fichier_design_label.grid(row=1, column=0, sticky="w")

fichier_design_entry = ctk.CTkEntry(parameter_frame, state=tk.DISABLED)
fichier_design_entry.grid(row=1, column=1, padx=5, pady=5, sticky="ew")

fichier_design_button = ctk.CTkButton(parameter_frame, text="Parse...", command=load_files)
fichier_design_button.grid(row=1, column=2, padx=5, pady=5)

# Entrée pour spécifier le nom du dossier de sortie
output_dir_label = ctk.CTkLabel(parameter_frame, text="Name of the output folder:")
output_dir_label.grid(row=2, column=0, sticky="w")

output_dir_entry = ctk.CTkEntry(parameter_frame)
output_dir_entry.grid(row=2, column=1, padx=5, pady=5, sticky="ew")

# Entrée pour spécifier la liste des gènes, ou le nombre de gènes à afficher
pickle_label = ctk.CTkLabel(parameter_frame, text="Name of the pickle file:")
pickle_label.grid(row=3, column=0, sticky="w")

pickle_entry = ctk.CTkEntry(parameter_frame)
pickle_entry.grid(row=3, column=1, padx=5, pady=5, sticky="ew")

uniprot_label = ctk.CTkLabel(parameter_frame, text="Number of Uniprot IDs:")
uniprot_label.grid(row=4, column=0, sticky="w")

uniprot_entry = ctk.CTkEntry(parameter_frame)
uniprot_entry.grid(row=4, column=1, padx=5, pady=5, sticky="ew")

# Choix du mode
mode_label = ctk.CTkLabel(parameter_frame, text="Mode:")
mode_label.grid(row=5, column=0, sticky="w")

mode_var = ctk.StringVar()
mode_var.set("time")  # Valeur par défaut

mode_time_radio = ctk.CTkRadioButton(parameter_frame, text="Time", variable=mode_var, value="time")
mode_time_radio.grid(row=5, column=1, padx=5, pady=5, sticky="w")

mode_vaccine_radio = ctk.CTkRadioButton(parameter_frame, text="Treated", variable=mode_var, value="treated")
mode_vaccine_radio.grid(row=5, column=2, padx=5, pady=5, sticky="w")


organism_label = ctk.CTkLabel(parameter_frame, text="Organism:")
organism_label.grid(row=6, column=0, sticky="w")

organism_var = ctk.StringVar()
organism_var.set("Human")  # Valeur par défaut

human_radio = ctk.CTkRadioButton(parameter_frame, text="Human", variable=organism_var, value="Human")
human_radio.grid(row=6, column=1, padx=5, pady=5, sticky="w")

mouse_radio = ctk.CTkRadioButton(parameter_frame, text="Mouse", variable=organism_var, value="Mouse")
mouse_radio.grid(row=6, column=2, padx=5, pady=5, sticky="w")
# Bouton pour exécuter l'analyse
run_button = ctk.CTkButton(root, text="Run", command=run_analysis)
run_button.pack(padx=10, pady=10)

root.mainloop()
