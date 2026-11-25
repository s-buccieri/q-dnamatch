
# Analytics: Benchmark di Scalabilità e Visualizzazione NISQ

---

## 1. Spiegazione del Problema

Per dimostrare in modo rigoroso il valore dell'algoritmo di Grover nella ricerca su sequenze di DNA, non è sufficiente eseguire una singola ricerca: servono strumenti oggettivi per rispondere a tre domande chiave:

1. **Vantaggio asintotico $(O(N)$ vs $O(\sqrt{N}))$:** sui simulatori locali possiamo eseguire circuiti completi solo per sequenze piccole $(N = 4, 8, 16)$, dove lo speedup non è ancora evidente. Serve quindi combinare i punti simulati con un'estrapolazione analitica fino a $N = 10^6$ basi per rendere visibile il vantaggio quadratico.
2. **Leggibilità probabilistica:** i dizionari grezzi restituiti dal runtime quantistico devono essere tradotti in visualizzazioni immediate che distinguano i picchi amplificati da Grover dal rumore di fondo.
3. **Impatto reale del rumore NISQ:** è necessario confrontare sullo stesso piano il simulatore ideale (`AerSimulator`) e il modello hardware reale (`FakeBrisbane`), affiancando alla distribuzione delle probabilità le metriche fisiche di transpilazione (`depth` e gate a 2 qubit).

---

## 2. Formulazione della Soluzione

Il modulo `analytics.py` separa la raccolta dei dati dalla visualizzazione, seguendo un'architettura **Interface-First**:

- **Benchmark riproducibile (`build_scaling_dataframe`):** genera un `pd.DataFrame` standardizzato con `seed` fissato (default `42`), combinando punti simulati (collegabili alla pipeline reale dei Rami 1–4 tramite callback `pipeline_runner`) ed estrapolazione analitica fino a $N = 10^6$.
- **Formula esatta con padding quantistico (`compute_grover_iterations`):** calcola il numero ottimale di iterazioni di Grover tenendo conto dell'arrotondamento alla potenza di 2 superiore:

```math
k = \left\lfloor \frac{\pi}{4} \sqrt{\frac{2^{n_{\text{idx}}}}{S}} \right\rfloor, \quad n_{\text{idx}} = \max\left(1,\; \lceil \log_2(\text{valid windows}) \rceil\right)
```

- **Tre grafici standardizzati (`matplotlib.figure.Figure`):** restituiscono oggetti `Figure` pronti sia per il rendering diretto in Streamlit (`st.pyplot(fig)`) sia per l'esportazione ad alta risoluzione tramite `export_figure`.

---

## 3. Soluzione: Componenti Principali e Contratto API

### 1. `compute_grover_iterations(num_windows, num_solutions=1) -> int`
Calcola le iterazioni ottimali di Grover gestendo i casi limite (input nulli/negativi, singola finestra, o numero di soluzioni $S \ge 2^{n_{\text{idx}}}$).

### 2. `build_scaling_dataframe(...) -> pd.DataFrame`
Costruisce la tabella di benchmark ordinata per $N$ crescente con colonne:
- `N`: lunghezza della sequenza DNA in basi.
- `valid_windows`: numero di finestre valide ($N - M + 1$).
- `classical_queries`: numero di confronti richiesti dalla ricerca classica $O(N)$.
- `grover_iterations`: numero di iterazioni quantistiche $O(\sqrt{N})$.
- `speedup`: rapporto `classical_queries / grover_iterations`.
- `mode`: `"simulated"` (zona eseguita su AerSimulator) oppure `"analytical"` (estrapolazione fino a $10^6$).

### 3. I 3 Grafici Standardizzati

- **`plot_scaling(df) -> Figure` (`fig_scaling`)**
  Grafico in scala **log-log** che confronta la curva classica $O(N)$ con la curva di Grover $O(\sqrt{N})$, evidenziando la fascia simulata su Qiskit Aer e annotando lo speedup raggiunto a $N = 10^6$ (~1243x per $M=2, S=1$).
- **`plot_histogram(exec_result, target_indices=None) -> Figure` (`fig_histogram`)**
  Istogramma delle probabilità per ciascun indice DNA a partire da `ExecutionResult`. Evidenzia in verde gli indici sopra soglia (o i target) e in grigio il rumore di fondo.
- **`plot_noise_vs_ideal(ideal_result, noisy_result, target_indices=None) -> Figure` (`fig_noise_vs_ideal`)**
  Grafico a barre affiancate sull'unione degli indici osservati tra simulazione ideale e `FakeBrisbane`. Include un riquadro riassuntivo con `depth`, conteggio dei gate a 2 qubit (`cx` / `ecr`) e probabilità aggregata sui target.

### 4. `export_figure(fig, output_path, dpi=300) -> Path`
Crea automaticamente le directory di destinazione se assenti e salva la figura in formato PNG o SVG con `bbox_inches="tight"`.
