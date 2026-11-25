# Q-DNAMatch — Quantum DNA Pattern Matching

Prototipo sviluppato in **Qiskit** per il **Qiskit Fall Fest 2026 (Catania)** — traccia *Quantum per le Scienze della Vita*.  
L'obiettivo è individuare rapidamente una sottostringa genetica (es. una mutazione o un *motif*) all'interno di una sequenza di DNA utilizzando una ricerca di **Grover basata su QROM**, dimostrando il vantaggio quantistico quadratico $O(\sqrt{N})$ rispetto alla ricerca lineare classica $O(N)$.

---

## Struttura della Repository

```text
q-dnamatch/
├── src/
│   ├── __init__.py
│   ├── dna_core.py        # Codifica DNA 2-bit, sliding window, baseline classica
│   ├── q_loader.py        # Circuito QROM per caricamento finestre DNA + inverse()
│   ├── q_grover.py        # Oracolo comparatore MCZ, uncomputation e diffusore
│   ├── q_runtime.py       # Esecuzione AerSimulator (ideale & FakeBackend), transpilazione
│   └── analytics.py       # Benchmark O(N) vs O(sqrt(N)) e generazione grafici
├── tests/
│   ├── dna_core/
|       └── test_dna_core.py
│   ├── quantum_core/
|       ├── test_q_loader.py
|       └── test_q_grover.py
│   └── runtime/
|       └── test_q_runtime.py
|   └── partials/
|       ├── test_end_to_end.py
|       ├── test_quantum_runtime.py
|       └── test_dna_quantum.py
├── docs/
|   ├── explanations/      # Documentazione dettagliata delle implementazioni
│   └── presentation.pdf   # Slide per il pitch
├── app.py                 # Web app interattiva in Streamlit
├── requirements.txt
└── README.md
```

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

---

## Pagina dedicata
```bash
# Accedi alla pagina dedicata
https://q-dnamatch.streamlit.app/
```