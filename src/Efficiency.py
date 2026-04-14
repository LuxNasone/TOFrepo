import argparse 
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import argparse
import pandas as pd

#Programma per calcolare l'efficienza di ogni piano di PMT di ogni telescopio per ogni valore della tensione è stato testato


Data_path="" #inizializzo data path, da modificare a seconda di dove si trovano i file da analizzare


# parser per input del file txt con cui fare l'efficienza 

#da lanciare come python3 Efficiency.py -f "file.txt"

parser = argparse.ArgumentParser()
parser.add_argument("file", help="Percorso del file")
args = parser.parse_args()

percorso = args.file
print("Percorso:", percorso)

df = pd.read_csv(percorso, sep=r"\s+")
V, S_1, S_2, S_3, D, T = np.loadtxt(percorso,unpack=True)

#leggi il nome del file e a seconda se trova 01,02,03 cambia il tipo di piano che si sta analizzando
#di conseguenza cambia il modo in cui calcolare le doppie fake 

w = (50. + 50. - 2 * 3.5) *1e-9 # tempo di coincidenza per doppie fake, 50 ns è il tempo di segnale, 3.5 ns è il tempo di sovrapposizione, moltiplico per 1e-9 per convertire in secondi
dt = 100 # tempo totale di acquisizione 100s 
R_1 = S_1/dt # rate di singole per il piano 1
R_2 = S_2/dt # rate di singole per il piano 2
R_3 = S_3/dt # rate di singole per il piano 3

# se definiamo tre piani di ogni telescopio come i,j,k 
# allora l'efficienza di ogni piano sarà eps_i=#triple(i&j&k)/#doppie(j&k)-#doppie_fake 
# le doppie fake si calcolano come 
# l'errore associato all'efficienza così trovata sarà... radice di(eps*(1-eps)/#doppie(j&k)-#doppie_fake)
# 
#  

if "01" in percorso:
    D_fake = dt * w * (R_2 * R_3) #calcolo delle doppie fake per il piano 1, che sono date dal prodotto dei rate di singole degli altri due piani moltiplicato per il tempo di coincidenza e il tempo totale di acquisizione
    Eff = T/(D-D_fake) #calcolo dell'efficienza del piano 1, che è data dal numero di triple  diviso  doppie 
    sig_eps = np.sqrt((Eff * (1-Eff))/(D-D_fake)) #calcolo dell'errore associato all'efficienza del piano formula dell'errore binomiale
    contamination = D_fake/D *100
    piano="i"
    k=1

elif "02" in percorso:
    D_fake = dt * w * (R_1 * R_3) #calcolo delle doppie fake per il piano 2, che sono date dal prodotto dei rate di singole degli altri due piani moltiplicato per il tempo di coincidenza e il tempo totale di acquisizione
    Eff = T/(D-D_fake) #calcolo dell'efficienza del piano 2, che è data dal numero di triple coincidences diviso il numero di doppie coincidences, 
    sig_eps = np.sqrt((Eff * (1-Eff))/(D-D_fake)) #calcolo dell'errore associato all'efficienza del piano formula dell'errore binomiale
    contamination = D_fake/D *100
    piano="j"
    k=2
elif "03" in percorso:
    D_fake = dt * w * (R_1 * R_2) #calcolo delle doppie fake per il piano 3, che sono date dal prodotto dei rate di singole degli altri due piani moltiplicato per il tempo di coincidenza e il tempo totale di acquisizione
    Eff = T/(D-D_fake) #calcolo dell'efficienza del piano 3, che è data dal numero di triple coincidences diviso il numero di doppie coincidences,
    contamination = D_fake/D *100 
    sig_eps = np.sqrt((Eff * (1-Eff))/(D-D_fake)) #calcolo dell'errore associato all'efficienza del piano formula dell'errore binomiale
    piano="k"
    k=3

# Grafico
plt.figure(figsize=(8, 5))
plt.errorbar(V, Eff,xerr=1, yerr=sig_eps, fmt='o')
plt.xlabel("Tensione (V)")
plt.ylabel("Efficienza")
plt.title(f"Rapporto triple/doppie in funzione della tensione PMT0{k}".format(k=piano))
plt.grid(True, linestyle='--', alpha=0.6)
#plt.legend()
plt.tight_layout()
plt.savefig(f"../Dati/Eff/grafico_PMT{k}.png", dpi=300, bbox_inches="tight")
plt.show()


print("Efficienza del piano",Eff)


# aggiunge le colonne

df["    Efficienza"] = Eff
df["    Sig_eps"] = sig_eps
df["    Contamination (%)"] = contamination

df["    Efficienza"] = df["    Efficienza"].map(lambda x: f"{x:.2g}")
df["    Sig_eps"] = df["    Sig_eps"].map(lambda x: f"{x:.2g}")
df["    Contamination (%)"] = df["    Contamination (%)"].map(lambda x: f"{x:.2g}")

# salva
if piano=="i":
    df.to_csv("../Dati/Dati_tabelle/FinalPMT01.txt", sep="\t", index=False)
elif piano=="j":
    df.to_csv("../Dati/Dati_tabelle/FinalPMT02.txt", sep="\t", index=False)
elif piano=="k":
    df.to_csv("../Dati/Dati_tabelle/FinalPMT03.txt", sep="\t", index=False)
