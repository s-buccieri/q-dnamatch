# Q_Runtime: Motore di Esecuzione, Transpilazione e Simulazione NISQ

---

## 1. Spiegazione del Problema

Una volta costruito il circuito astratto dell'algoritmo di Grover (`q_grover.py`), non è possibile eseguirlo direttamente su un processore quantistico reale né valutarne le prestazioni pratiche senza affrontare tre problemi tipici dell'era **NISQ (*Noisy Intermediate-Scale Quantum*)**:

1. **Astrazione circuitale vs Hardware reale:** il circuito generato da `q_grover.py` utilizza gate logici ad alto livello (come `MCX` o il blocco composito `QROM_DNA`) che non esistono fisicamente sul silicio. Un chip quantistico reale supporta solo un insieme ristretto di *gate nativi* e una *topologia di connessione* limitata tra i qubit.
2. **Degradazione dello stato e rumore:** sui dispositivi attuali, lo stato quantistico si degrada principalmente per due fattori:
   - **Profondità del circuito (`depth`):** più il circuito è lungo temporalmente, più i qubit perdono coerenza (decoerenza).
   - **Gate a 2 qubit (`cx` / `ecr`):** le operazioni che generano entanglement tra coppie di qubit sono le più lente e le più soggette a errori fisici.
3. **Decodifica e filtraggio statistico:** l'esecuzione quantistica è probabilistica e viene ripetuta per un numero finito di campionamenti (`shots`). Qiskit restituisce i conteggi grezzi come dizionario di stringhe binarie in formato **little-endian** (dove il qubit $q_0$ è il bit meno significativo, cioè quello più a destra). In presenza di rumore, oltre ai picchi amplificati da Grover compaiono stati spuri a bassa probabilità che devono essere filtrati.

---

## 2. Formulazione della Soluzione

Il modulo `q_runtime.py` funge da livello di esecuzione e astrazione hardware per la pipeline Q-DNAMatch. La strategia si articola in tre meccanismi:

- **Doppia modalità di simulazione (`AerSimulator`):**
  - **Ideale (`use_noise=False`):** simulazione matematica priva di errori per verificare la correttezza teorica della ricerca di Grover.
  - **Rumorosa (`use_noise=True`):** tramite `AerSimulator.from_backend(FakeBrisbane())` (da `qiskit_ibm_runtime.fake_provider`), il simulatore locale importa il modello di rumore reale, la topologia e i gate nativi del processore **IBM Brisbane** (127 qubit, architettura Eagle), emulando errori realistici di gate e di lettura.
- **Transpilazione ottimizzata (`optimization_level=3`):** tramite `generate_preset_pass_manager(optimization_level=3)` il circuito viene compilato e ottimizzato in modo aggressivo per il backend selezionato (analogamente al flag `-O3` di un compilatore classico), eliminando gate ridondanti e minimizzando `depth` e gate a 2 qubit.
- **Decodifica Little-Endian e soglia statistica:** converte le bitstring in indici interi del DNA e filtra il rumore di fondo tramite una soglia adattiva, preservando sia i singoli match sia le occorrenze multiple ($S > 1$).

---

## 3. Soluzione: Architettura del Modulo e Contratto API

### Componenti Principali

#### 1. `get_backend(use_noise: bool = False) -> AerSimulator`
Funzione *factory* che istanzia il simulatore: restituisce `AerSimulator()` in modalità ideale oppure `AerSimulator.from_backend(FakeBrisbane())` quando `use_noise=True`.

#### 2. `count_two_qubit_gates(qc: QuantumCircuit) -> int`
Scorre le istruzioni del circuito transpilato (`qc.data`) e conta tutti i gate che agiscono su esattamente 2 qubit, escludendo le istruzioni non fisiche (`barrier` e `measure`).
- **Scelta progettuale:** i processori IBM Eagle come `FakeBrisbane` utilizzano il gate nativo **`ecr`** (*Echoed Cross-Resonance*) al posto del `cx` (`CNOT`). Contare genericamente `num_qubits == 2` garantisce che la metrica `cx_count` sia sempre accurata su entrambi i backend, evitando di restituire erroneamente `0` in modalità rumorosa.

#### 3. `decode_results(counts: dict[str, int], threshold: float | None = None)`
Riceve il dizionario grezzo `counts` restituito da Qiskit (dove la somma dei valori è pari a `shots`) ed esegue:
1. **Conversione little-endian:** ripulisce la bitstring e la converte in base 2 tramite `int(clean_bits, 2)` (es. `'001'` $\rightarrow$ indice `1`).
2. **Normalizzazione:** calcola la probabilità di ciascun indice dividendo le occorrenze per `total_shots`.
3. **Filtro del rumore:** se non viene fornita una soglia esplicita `threshold`, applica una soglia statistica adattiva:

$$\text{threshold} = \max\left(0.10,\; \frac{2}{2^{n_{\text{idx}}}}\right)$$

Poiché una distribuzione di puro rumore uniforme assegnerebbe probabilità $1 / 2^{n_{\text{idx}}}$ a ogni stato, richiedere almeno il doppio della media uniforme (con minimo $10\%$) isola i picchi reali di Grover dal rumore di fondo.

> **Gestione di indici multipli e verifica dei falsi positivi:** in presenza di occorrenze multiple $(S > 1)$ o di rumore elevato, `found_indices` può contenere più indici con probabilità differenti. Dal punto di vista computazionale, mentre la **ricerca** quantistica costa $O(\sqrt{N})$, la **verifica** classica di un singolo indice candidato sulla sequenza DNA costa $O(1)$, permettendo di validare istantaneamente i candidati senza compromettere il vantaggio quantistico.

#### 4. `run_circuit(grover_qc, shots=1024, use_noise=False, threshold=None, seed=42)`
Orchestratore pubblico del modulo:
1. Inizializza il backend tramite `get_backend(use_noise)`.
2. Compila `grover_qc` con `generate_preset_pass_manager(optimization_level=3, backend=backend, seed_transpiler=seed)`.
3. Estrae le metriche hardware `depth` e `cx_count` dal circuito transpilato.
4. Esegue la simulazione misurando il tempo in millisecondi con `time.perf_counter()` e `seed_simulator=seed` per garantire la riproducibilità.
5. Decodifica `result.get_counts()` tramite `decode_results` e restituisce il dizionario `ExecutionResult`.

---

### Contratto di Interfaccia (`ExecutionResult`)

Seguendo l'approccio **Interface-First**, `run_circuit` espone un dizionario standardizzato `ExecutionResult` destinato a `analytics.py` e `app.py`:

```python
{
    "counts": dict[str, int],          # Conteggi grezzi {bitstring: occorrenze}
    "found_indices": list[int],        # Lista ordinata degli indici DNA sopra soglia
    "probabilities": dict[int, float], # Distribuzione {indice_intero: probabilità}
    "depth": int,                      # Profondità del circuito transpilato
    "cx_count": int,                   # Numero di gate a 2 qubit (CX / ECR)
    "exec_time_ms": float              # Tempo di esecuzione della simulazione (ms)
}