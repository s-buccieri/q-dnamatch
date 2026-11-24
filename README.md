# Q-DNAMatch — Quantum DNA Pattern Matching

**Q-DNAMatch** è un prototipo sviluppato in **Qiskit** per il **Qiskit Fall Fest 2026 (Catania)** nell'ambito della traccia *Quantum per le Scienze della Vita*. Il progetto affronta uno dei problemi fondamentali della bioinformatica — l'individuazione di specifiche sottostringhe genetiche (come mutazioni o *motif*) all'interno di sequenze di DNA — sfruttando l'algoritmo di **Grover** per dimostrare il vantaggio quantistico quadratico $O(\sqrt{N})$ rispetto alla ricerca lineare classica $O(N)$.

L'architettura traduce le basi azotate in una codifica binaria a 2 bit (`A=00`, `C=01`, `G=10`, `T=11`) e scompone la sequenza tramite *sliding window*. Un modulo dedicato costruisce una **QROM (Quantum Read-Only Memory)** reversibile che carica tutte le finestre di DNA in sovrapposizione quantistica, proteggendo gli stati di padding tramite un qubit di validità dedicato per evitare falsi positivi (come il pattern `"AA"`). L'oracolo di fase confronta in parallelo le finestre con il pattern cercato tramite un gate `MCZ`, esegue l'*uncomputation* dei registri ausiliari per rimuovere l'entanglement residuo e lascia al diffusore di Grover il compito di amplificare la probabilità degli indici corretti nel numero ottimale di iterazioni.

Il circuito viene infine transpilato al massimo livello di ottimizzazione (`optimization_level=3`) ed eseguito tramite `AerSimulator` sia in regime ideale sia con il modello di rumore realistico **IBM `FakeBrisbane`**, permettendo di valutare l'impatto dell'hardware NISQ su fedeltà, profondità e numero di gate a due qubit. L'intera pipeline è accessibile attraverso una web app interattiva in **Streamlit** che integra la demo live con evidenziazione dei match sul DNA, l'ispezione del circuito, i grafici di scalabilità fino a $10^6$ basi e il confronto diretto tra simulazione ideale e rumorosa.

---

## Struttura della Repository

```text
q-dnamatch/
├── src/
│   ├── __init__.py
│   ├── dna_core.py        # Codifica DNA 2-bit, sliding window e baseline classica O(N)
│   ├── q_loader.py        # Circuito QROM invertibile per il caricamento delle finestre DNA
│   ├── q_grover.py        # Oracolo di fase MCZ, uncomputation, diffusore e circuito Grover
│   ├── q_runtime.py       # Transpilazione, esecuzione AerSimulator (ideale e FakeBrisbane) e decodifica
│   └── analytics.py       # Benchmark di scalabilità O(N) vs O(sqrt(N)) e generazione grafici
├── tests/                 # Test unitari e di integrazione end-to-end (pytest)
├── docs/                  # Documentazione tecnica e slide del pitch
├── app.py                 # Web app interattiva in Streamlit
├── requirements.txt
└── README.md