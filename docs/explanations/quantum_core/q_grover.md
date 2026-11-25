# Q_Grover: Implementazione dell'algoritmo di Grover per la ricerca quantistica del pattern

## 1. Spiegazione del Problema

Dopo aver caricato tutte le finestre di DNA nella sovrapposizione quantistica tramite la QROM (`q_loader.py`), ci troviamo di fronte al problema centrale della ricerca quantistica: **trovare il pattern** riducendo la complessità da $O(N)$ a $O(\sqrt{N})$, dove $N = 2^{n_{\text{idx}}}$ è la dimensione dello spazio degli indici.

Il solo caricamento dei dati in memoria quantistica non è sufficiente per risolvere il problema, a causa di due limiti fondamentali della meccanica quantistica:

1. **Distribuzione di probabilità uniforme:** se misurassimo i qubit subito dopo l'esecuzione della QROM, ogni possibile posizione $i \in \{0, \dots, N-1\}$ avrebbe la stessa identica probabilità di essere osservata:

$$P(i) = \left|\frac{1}{\sqrt{N}}\right|^2 = \frac{1}{N}$$

   Il risultato della misura sarebbe quindi un indice completamente casuale.
2. **Impossibilità di ispezione diretta:** in un computer quantistico non è possibile "sbirciare" all'interno della sovrapposizione per controllare quale ramo contiene il pattern cercato. Qualsiasi operazione di lettura diretta provoca il collasso della funzione d'onda su un singolo stato scelto casualmente.

È quindi necessario un meccanismo che manipoli le ampiezze di probabilità **prima della misura**, concentrando l'intera massa di probabilità sugli indici in cui compare il pattern ed eliminando la probabilità sugli indici errati.

---

## 2. Formulazione della Soluzione

Il modulo `q_grover.py` risolve il problema implementando l'**algoritmo di Grover**, utilizzando la QROM di `q_loader.py` come oracolo di accesso ai dati. La strategia si fonda su due principi:

* **Distinguibilità tramite fase relativa:** poiché la probabilità dipende dal modulo quadro dell'ampiezza $(P(i) = |\alpha_i|^2)$, cambiare il segno di un'ampiezza da $+\alpha$ a $-\alpha$ non altera ancora le probabilità osservabili, ma rende gli stati soluzione $(i^{\ast})$ matematicamente distinguibili da tutti gli altri.
* **Amplificazione di ampiezza (Inversione rispetto alla media):** applicando ripetutamente una riflessione attorno alla media delle ampiezze, gli stati con ampiezza negativa vengono spinti verso valori positivi molto elevati (interferenza costruttiva), mentre gli stati con ampiezza positiva scendono verso zero (interferenza distruttiva).

### Geometria e numero ottimale di iterazioni

Geometricamente, lo stato del sistema ruota in un piano bidimensionale definito dalla sovrapposizione degli stati errati e dalla sovrapposizione degli $S$ stati soluzione. Ogni iterazione completa (Oracolo + Diffusore) ruota il vettore di stato di un angolo costante $2\theta$ verso il sottospazio delle soluzioni, dove:

$$\theta = \arcsin\left(\sqrt{\frac{S}{N}}\right)$$

Per massimizzare la probabilità di successo ed evitare il fenomeno dell'**over-rotation** (superare i $90^\circ$ facendo calare nuovamente la probabilità della soluzione), il numero esatto di iterazioni $k$ è calcolato dalla funzione `compute_optimal_iterations`:

$$k = \left\lfloor \frac{\pi}{4} \sqrt{\frac{2^{n_{\text{idx}}}}{S}} \right\rfloor$$

| Spazio degli indici $(N)$ | Occorrenze del pattern $(S)$ | Iterazioni ottimali $(k)$ | Probabilità di successo |
| :---: | :---: | :---: | :---: |
| $8$ | $1$ | $2$ | $94.53\%$ |
| $8$ | $2$ | $1$ | $100.00\%$ (divisa sui $2$ indici) |

---

## 3. Soluzione: Architettura del Circuito in `q_grover.py`

La funzione `build_grover_circuit` costruisce il circuito quantistico completo coordinando tre registri quantistici (`idx`, `data`, `valid`) e un registro classico (`c_idx`).

### Fase A: Inizializzazione in sovrapposizione uniforme
All'avvio del circuito, viene applicato un gate di Hadamard $(H)$ su ciascuno degli $n_{\text{idx}}$ qubit del registro `idx`, mentre i registri `data` e `valid` restano inizializzati a $|0\rangle$. Tutti gli $N$ indici partono così con la medesima ampiezza positiva $\alpha = \frac{1}{\sqrt{N}}$:

$$|\psi_0\rangle = H^{\otimes n_{\text{idx}}} |0\rangle_{\text{idx}} \otimes |0\dots 0\rangle_{\text{data}} \otimes |0\rangle_{\text{valid}} = \frac{1}{\sqrt{N}} \sum_{i=0}^{N-1} |i\rangle_{\text{idx}} \otimes |0\dots 0\rangle_{\text{data}} \otimes |0\rangle_{\text{valid}}$$

---

### Fase B: Il Ciclo di Grover $(k$ iterazioni)
Il circuito esegue in sequenza per $k$ volte due blocchi unitari: l'**Oracolo di Fase** (`build_pattern_oracle`) e il **Diffusore** (`build_diffuser`).

#### 1. L'Oracolo di Fase (`build_pattern_oracle`)
L'oracolo lascia invariati gli stati errati e moltiplica per $-1$ esclusivamente gli stati $|i^{\ast}\rangle$ in cui la finestra di DNA coincide con `pattern_bin` **e** il flag di validità è attivo (`valid == 1`, escludendo il padding). Opera in **5 passaggi simmetrici**:

1. **Caricamento QROM:** invoca `loader_gate` per legare in entanglement ogni indice $|i\rangle$ alla propria finestra di DNA $|w_i\rangle$ e al bit di validità $|v_i\rangle$.
2. **Matching del pattern:** applica gate `X` sui qubit del registro `data` in corrispondenza dei bit in cui `pattern_bin` vale `'0'`. In questo modo, esclusivamente i rami in cui `data == pattern_bin` assumono la configurazione $|11\dots 1\rangle_{\text{data}}$.
3. **Marcatura di fase (`MCZ` senza ancilla):** applica in sequenza `H`, `MCX` (controllato da tutti i qubit di `data` e `valid`) e `H` sul qubit `valid_reg[0]`. Sfruttando l'identità quantistica $HXH = Z$, questa sequenza realizza un gate **Multi-Controlled Z (`MCZ`)** senza richiedere qubit ancilla aggiuntivi, imprimendo una fase $-1$ solo se `data == pattern_bin` e `valid == 1`.
4. **Ripristino del registro `data`:** riapplica gli stessi gate `X` del passaggio 2 per riportare i qubit di `data` al loro valore binario originale $|w_i\rangle$.
5. **Uncomputation:** applica `loader_gate.inverse()`, riavvolgendo la QROM per riportare `data` e `valid` allo stato pulito $|0\dots 0\rangle|0\rangle$ e sciogliendo l'entanglement, mentre la fase $-1$ rimane impressa sull'indice soluzione:

$$|\psi_{\text{oracle}}\rangle = \frac{1}{\sqrt{N}} \left(\sum_{i \neq i^{\ast}} |i\rangle_{\text{idx}} - |i^{\ast}\rangle_{\text{idx}} \right) \otimes |0\dots 0\rangle_{\text{data}} \otimes |0\rangle_{\text{valid}}$$

#### 2. Il Diffusore (`build_diffuser`)
Liberato `idx` dall'entanglement, il Diffusore agisce esclusivamente sugli $n_{\text{idx}}$ qubit di indice implementando l'operatore di inversione rispetto alla media:

$$D = 2|s\rangle\langle s| - I$$

Matematicamente, trasforma ogni ampiezza $\alpha_i$ in:

$$\alpha_i' = 2\mu - \alpha_i$$

dove $\mu = \frac{1}{N} \sum_{j=0}^{N-1} \alpha_j$ è la media aritmetica delle ampiezze correnti. Poiché gli stati marcati hanno ampiezza negativa molto inferiore a $\mu$, la riflessione rispetto a $\mu$ ne amplifica drasticamente il valore positivo, riducendo contemporaneamente le ampiezze degli stati errati.

A livello circuitale, il blocco esegue la sequenza:
1. `H` su tutti i qubit di `idx`
2. `X` su tutti i qubit di `idx`
3. `MCZ` su `idx` (realizzato tramite `H`, `MCX`, `H` sull'ultimo qubit di indice, oppure un singolo gate `Z` nel caso limite $n_{\text{idx}} = 1$)
4. `X` su tutti i qubit di `idx`
5. `H` su tutti i qubit di `idx`

Questa sequenza implementa $-D$, equivalente a $D$ a meno di una fase globale $-1$ fisicamente inosservabile.

---

### Fase C: Misura e Contratto di Interfaccia
Al termine delle $k$ iterazioni, `build_grover_circuit` applica l'operazione di misura (`qc.measure`) **esclusivamente sul registro `idx`** verso il registro classico `c_idx`:

* I registri `data` e `valid` non vengono misurati: grazie all'uncomputation, sono già garantiti nello stato $|0\dots 0\rangle|0\rangle$.
* Escluderli dalla misura evita di aggiungere bit superflui alle bitstring restituite, rispettando rigorosamente il contratto di interfaccia richiesto dal modulo `q_runtime.py` per l'estrazione finale dell'indice vincente.
