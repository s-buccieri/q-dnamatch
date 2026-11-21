# Q-DNAMatch — Quantum DNA Pattern Matching

Prototipo sviluppato in **Qiskit** per il **Qiskit Fall Fest 2026 (Catania)** — traccia *Quantum per le Scienze della Vita*.  
L'obiettivo è individuare rapidamente una sottostringa genetica (es. una mutazione o un *motif*) all'interno di una sequenza di DNA utilizzando una ricerca di **Grover basata su QROM**, dimostrando il vantaggio quantistico quadratico $O(\sqrt{N})$ rispetto alla ricerca lineare classica $O(N)$.

---

## Struttura della Repository

```text
q-dnamatch/
├── src/
│   ├── __init__.py
│   ├── dna_core.py        # [Coppia 3] Codifica DNA 2-bit, sliding window, baseline classica
│   ├── q_loader.py        # [Coppia 1] Circuito QROM per caricamento finestre DNA + inverse()
│   ├── q_grover.py        # [Coppia 1] Oracolo comparatore MCZ, uncomputation e diffusore
│   ├── q_runtime.py       # [Coppia 2] Esecuzione AerSimulator (ideale & FakeBackend), transpilazione
│   └── analytics.py       # [Coppia 2] Benchmark O(N) vs O(sqrt(N)) e generazione grafici
├── tests/
│   ├── test_dna_core.py
│   ├── test_quantum_core.py
│   └── test_runtime.py
├── docs/
│   └── presentation.pdf   # [Coppia 3] Slide per il pitch da 5 minuti
├── app.py                 # [Coppia 3] Web app interattiva in Streamlit
├── requirements.txt
└── README.md
```

---

## Strategia dei Branch e Specifiche per Coppia

Il progetto utilizza un modello a **5 branch**: `main` (produzione/demo finale), `dev` (integrazione graduale e test parziali) e **3 branch di lavoro** assegnati alle coppie.

### 1. `feat/bio-product` — Coppia 3: Bio-Product & Pitch
* **File di competenza:** `src/dna_core.py`, `app.py`, `tests/test_dna_core.py`, `docs/`
* **Obiettivi:**
  * Validare le sequenze DNA (`A, C, G, T`), codificarle a 2 bit (`A=00, C=01, G=10, T=11`), generare le finestre (*sliding window*) con padding a potenza di 2 e implementare la ricerca classica con conteggio esatto delle query.
  * Sviluppare l'interfaccia Streamlit (3 tab: *Demo Live*, *Classico vs Quantistico*, *Hardware & Noise*) e coordinare le slide del pitch da 5 minuti.
* **Contratto in uscita (verso `feat/quantum-core` e `feat/runtime-analytics`):**
  * Dizionario `DNAPayload`:
    ```python
    {
        "windows_bin": list[str],
        "pattern_bin": str,
        "n_idx": int,
        "n_data": int,
        "valid_len": int,
        "classical_res": dict  # {"indices": list[int], "queries": int}
    }
    ```

### 2. `feat/quantum-core` — Coppia 1: Quantum Core
* **File di competenza:** `src/q_loader.py`, `src/q_grover.py`, `tests/test_quantum_core.py`
* **Obiettivi:**
  * Implementare `build_dna_loader(windows_bin, n_idx, n_data)` come `qiskit.circuit.Gate` (`"QROM_DNA"`) invertibile tramite `.inverse()`.
  * Implementare l'oracolo di fase (`MCZ` + *uncomputation* del registro dati), il diffusore sui soli qubit di indice e il calcolo del numero ottimale di iterazioni $k = \lfloor \frac{\pi}{4}\sqrt{2^{n_{\text{idx}}} / S} \rfloor$.
* **Contratto in ingresso:** `DNAPayload` (da `feat/bio-product`).
* **Contratto in uscita (verso `feat/runtime-analytics`):**
  * Circuito `grover_qc: qiskit.QuantumCircuit` con misura applicata solo sui `n_idx` qubit di indice + metadato `iterations: int`.

### 3. `feat/runtime-analytics` — Coppia 2: Runtime & Advantage
* **File di competenza:** `src/q_runtime.py`, `src/analytics.py`, `tests/test_runtime.py`
* **Obiettivi:**
  * Eseguire `grover_qc` su `AerSimulator` in modalità ideale e rumorosa (`FakeBackend`, es. `FakeBrisbane`), applicando transpilazione ottimizzata (`optimization_level=3`) e decodifica *little-endian* con filtro del rumore.
  * Produttore dei 3 grafici chiave: scalabilità $O(N)$ vs $O(\sqrt{N})$, istogramma delle probabilità e confronto Ideale vs Rumoroso.
* **Contratto in ingresso:** `grover_qc` (da `feat/quantum-core`) e `DNAPayload` (da `feat/bio-product`).
* **Contratto in uscita (verso `feat/bio-product` / `app.py`):**
  * Dizionario `ExecutionResult`:
    ```python
    {
        "counts": dict[str, int],
        "found_indices": list[int],
        "probabilities": dict[int, float],
        "depth": int,
        "cx_count": int,
        "exec_time_ms": float
    }
    ```
  * Oggetti grafici `fig_scaling`, `fig_histogram`, `fig_noise_vs_ideal` (`matplotlib`/`plotly`).

---

## Workflow di Merge Graduale e Test Parziali

Per evitare conflitti dell'ultimo minuto, **nessuna coppia fa merge diretto su `main`**. L'integrazione avviene progressivamente sul branch `dev` tramite Pull Request (PR):

1. **Fase 1 — Sviluppo isolato (Giorni 1–2):**
   * Ogni coppia lavora sul proprio branch `feat/*` utilizzando dati *mock* conformi ai contratti API sopra definiti.
   * Prima di aprire una PR, tutti gli unit test locali (`pytest`) del proprio modulo devono passare.
2. **Fase 2 — Merge graduali su `dev` per test parziali (Giorno 3):**
   * **Step A:** `feat/bio-product` fa merge di `dna_core.py` su `dev`. Le altre coppie fanno `git pull origin dev` per usare il vero encoding DNA al posto dei mock.
   * **Step B:** `feat/quantum-core` fa merge su `dev`. **Test parziale 1:** verifica `dna_core` $\to$ `q_loader` $\to$ `q_grover` su `Statevector`.
   * **Step C:** `feat/runtime-analytics` fa merge su `dev`. **Test parziale 2:** verifica `grover_qc` $\to$ `q_runtime` (ideale e rumoroso) $\to$ `analytics`.
   * **Step D:** `feat/bio-product` collega `app.py` ai moduli reali presenti su `dev` e apre PR verso `dev`.
3. **Fase 3 — Riunione Finale e Merge su `main` (Giorno 4, Code Freeze):**
   * Tutto il team si riunisce per eseguire il **test end-to-end completo** su `dev` (UI Streamlit, casi limite, grafici e benchmark).
   * Solo quando l'intera pipeline funziona senza errori, si apre la PR finale `dev` $\to$ `main`, si effettua il merge e si blocca il codice per la presentazione del 25 Novembre.

---

## Installazione e Avvio Rapido

```bash
# 1. Crea e attiva un ambiente virtuale
python -m venv .venv
source .venv/bin/activate  # Su Windows: .venv\Scripts\activate

# 2. Installa le dipendenze
pip install -r requirements.txt

# 3. Esegui i test unitari
pytest

# 4. Avvia la web app
streamlit run app.py
```