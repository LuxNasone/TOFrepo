
# Istruzioni per l'utilizzo dei programmi di TOF

## I programmi di TOF sono:
### - TOF_AcceptanceMC_v2.cpp
### - TOF_Calibration_v9.cpp
### - TOF_Analysis_v6.cpp

Sono necessari i risultati di tutti e tre i programmi per avere l'analisi completa.

# ==============================================

#                            TOF_AcceptanceMC_v2.cpp
# ==============================================


### Programa per il montecarlo di simulazione dell'accettanza del PMT03 aggiornato con la funzione per simulare un numero preimpostato $N$  di raggi cosmici generati a posizione definita sulla superficie superiore della barra scintillatrice, vengono poi generati l'angolo $\theta$ secondo la distribuzione $\cos^n{\theta}$  con $n$ impostato  da linea di comando e poi viene generato in maniera uniforme $\phi$ tra $[0,2\pi]$, a questo punto si ricostruisce la traiettoria del cosmico e si verifica se intercetta il piano dello scintillatore 3 (la cui posizione è modificabile da terminale)


### - Caricamento:  .L TOF_AcceptanceMC_v2.cpp
### ParallaxScan(1e7, 2.0, +4.5, -135.0, +135.0, 2.0,true, "parallax_GuidaB.root");

- #### 1e7 è il numero $N$ di cosmici simulati per ogni poszione 
- #### 2 è l'esponente della distribuzione dei cosmici 
- #### +4.5 è la distanza tra lo scintillatore del PMT03 e la superficie superiore della barra scintillatrice  
- #### -135 è l'estremo inferiore della barra da simulare
- #### +135 è l'estremo superiore della barra da simulare
- #### 2 è il passo in centimetri tra le posizioni impostate così vengono esplorate -135,-133,-131,...,+135



##### Una volta conclusa l'analisi  verrà fornito un file ROOT con tutti i grafici e i valori da caricare nel programma di calibrazione
Il nome di default del file root in uscita è : parallax_GuidaB.root



# ==============================================

#                            TOF_Calibration_v9.cpp
# ==============================================


#### Codice utilizzato per le analisi della calibrazione della barra scintillatrice
### Scopi principali:

- #### Ricostruire la retta di calibrazione $\Delta t_{1,2}$ in funzione di $x_k$ (poszione del PMT03), in questo modo otteniamo la relazione
	$$
		\Delta t_{1,2} = \frac{2\cdot x_k}{v_{bar}} + (\Delta_1 + \Delta_2)
					$$
		in cui $v_{bar}$ è la velocità effettiva della luce all'interno della barra scintillatrice
	Per fare questo è necessario avere diversi punti di calibrazione, sono utilizzati quindi diversi dataset che serviranno alla realizzazione della retta e del suo fit.
	In particolare noi abbiamo utlizzato 10 punti $[-130,-112,-84,-56,-28,0,+28,+56,+84,+112,+130]$
	Questi compongono i dataset che sono indicizzati con N per i negativi e X per i positivi 

	*esempio* : posizione $x_k=-130 \rightarrow$ N130.xml 
	lo 0 è indicizzato come X0.xml !!

	Questi dataset devono avere questo specifico nome per essere utilizzati dal programma.
	
- #### Calcolare la Costante di calibrazione C(x_k) necessaria per la ricostruzione dei TOF:
	 Il tempo misurato dalla DRS è legato al TOF e un ritardo intrinseco di tutta l'elettronica e la propagazione dei segnali pari a 
	 $$
		T_{mis} = TOF + C(x_k)
		$$
		con $T_{mis}$ dato da 
		$$ 
			T_{mis} = t_3 - \frac{t_1+t_2}{2}
			$$
	Considerando il TOF nella configurazione di misura utlizzata (scintillatore del PMT3 a 5cm dalla superficie superiore della barra) TOF=5/30 ns = 0.166ns circa zero ( DA RIFLETTERCI)
	 Otteniamo quindi 
	 $$
		T_{mis}=  t_3 - \frac{t_1+t_2}{2} = C(x_k)
		$$
	si realizzano quindi gli istogrammi di $T_{mis}$ per ogni dataset, si esegue un fit con una gaussiana nel picco nell'intervallo del 90% degli eventi e si utilizzano i valori così ottenuti per realizzare il grafico C in funzione di $x_k$


- #### Stima della risoluzione temporale della barra:
	Osservare la larghezza della distribuzione $\Delta t_{1,2}$ per ogni posizione $x_k$ del PMT03
	Si realizza l'istogramma che poi viene fittato nel picco on intervallo 90% degli eventi  in maniera da trovare la media $\mu$ e la larghezza $\sigma$	
	Questi costituiscono anche i punti del grafico per la retta di calibrazione del punto 1
		
 

Nel programma sono presenti istruzioni per scartare i segnali oscillanti ed effettuare altri controlli ma ora non ho voglia di spiegare tutto.


# Come utilizzare il programma

### - Caricare con: .L TOF_Calibration_v9.cpp
### - Eseguire l'analisi con: 
###        TOF_Calibration(const char* folder, const char* outname = "TOF_Calibration_output.root",const char* parallax_file = nullptr)

##### Gli argomenti in ordine sono:
- ###### Percorso della CARTELLA contenente i dataset di calibrazione (*N130,N112,X0, ecc*)  
- ###### Nome del file root di output desiderato, esempio : Calibration.root
- ##### File (con eventuale percorso relativo) del risultato dell'analisi montecarlo effettuata in precedenza (vedi primo programma)

### Il risultato sarà un file root (ad esempio "Calibration.root") che servirà per utilizzare correttamente il programma di analisi dei TOF 



# ==============================================

#                            TOF_Analysis_v6.cpp
# ==============================================


#### INSERISCI SPIEGAZIONE SOMMARIA






## Istruzioni del programma:

- ### LoadCalibrationFromFile():
- ### Set Geometry(): 
- ### TOF_Analysis(): 

## Il risultato finale sarà un file .root con tutti gli istogrammi di analisi realizzati come: 
- #### istogramma dei $\beta$ 
- ### istogramma dei TOF 
- ### Istogramma delle posizioni di impatto ricostruite $x_{imp}$
- ### Istogramma dei $\theta$ ricostruiti 
ecc ecc