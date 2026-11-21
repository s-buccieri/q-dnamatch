# Q_Loader: Costruzione della QROM per la corretta applicazione dell'algoritmo di Grover

---

## 1. Spiegazione del Problema

Nello *string/pattern matching* applicato alla bioinformatica, l'obiettivo è individuare l'indice esatto in cui una sottosequenza (pattern) di lunghezza $M$ basi compare all'interno di una sequenza di DNA più lunga composta da $N$ basi azotate $(\text{A, C, G, T})$.

Applicando la tecnica dello *sliding window* (finestra scorrevole), la sequenza principale viene suddivisa in $K$ finestre sovrapposte, dove:

$$K = N - M + 1$$

Dal punto di vista della complessità computazionale:
* **Approccio Classico:** richiede un controllo sequenziale, con una complessità temporale di $O(K)$ operazioni nel caso peggiore.
* **Approccio Quantistico (Algoritmo di Grover):** permette di ridurre il numero di interrogazioni a $O(\sqrt{K})$, garantendo un **vantaggio quadratico**.

### Il collo di bottiglia hardware/quantistico
Affinché l'algoritmo di Grover possa operare, l'oracolo quantistico deve valutare **contemporaneamente** tutte le $K$ finestre in sovrapposizione. Tuttavia, un circuito quantistico presenta due limiti fondamentali rispetto al calcolo classico:
1. Non può accedere dinamicamente a una variabile Python o a un array salvato nella RAM durante l'esecuzione del circuito.
2. Non può eseguire costrutti condizionali classici come `if windows[i] == pattern`.

Inoltre, il registro di ricerca di Grover contiene inizialmente solo la sovrapposizione degli **indici** $(|0\rangle, |1\rangle, |2\rangle, \dots, |K-1\rangle)$ e non le informazioni relative alle basi azotate corrispondenti. È quindi necessario creare un ponte che associ in modo coerente ogni indice ai rispettivi dati.

---

## 2. Formulazione della Soluzione

Per risolvere il problema del caricamento dei dati, il modulo `q_loader.py` implementa una **QROM (*Quantum Read-Only Memory*)**, un circuito quantistico progettato per caricare l'intero database delle finestre di DNA direttamente nello stato quantistico.

### Architettura dei Registri Quantistici
La QROM opera su tre registri quantistici distinti:

1. **`idx` (Registro degli Indici):** composto da $n_{\text{idx}} = \lceil \log_2(K) \rceil$ qubit, utilizzati per rappresentare l'indirizzo binario $|i\rangle$ di ciascuna finestra.
2. **`data` (Registro dei Dati):** composto da $n_{\text{data}} = 2M$ qubit. Poiché le basi azotate sono quattro, ogni base viene codificata utilizzando 2 bit:
   * $\text{A} = 00$
   * $\text{C} = 01$
   * $\text{G} = 10$
   * $\text{T} = 11$
3. **`valid` (Registro di Validità):** composto da $1$ singolo qubit che funge da "timbro di autenticità" per distinguere i dati reali dal *padding*.

### Gestione del Padding (Falsi Positivi)
Un registro di $n_{\text{idx}}$ qubit spazia sempre su una dimensione totale pari a $2^{n_{\text{idx}}}$ stati. Poiché quasi sempre $K < 2^{n_{\text{idx}}}$, gli slot in eccesso (da $K$ a $2^{n_{\text{idx}}}-1$) rimangono vuoti, mantenendo il registro `data` inizializzato allo stato $|00\dots 0\rangle$.

Nella nostra codifica binaria, la sequenza di zeri `00...0` corrisponde a una serie di adenine (`"AA...A"`). Senza un meccanismo di controllo, la ricerca del pattern `"AA"` amplificherebbe erroneamente anche gli stati di *padding*, generando **falsi positivi**. Impostando il qubit `valid` a $|1\rangle$ esclusivamente per le $K$ finestre reali (e lasciandolo a $|0\rangle$ sul *padding*), l'oracolo può filtrare i risultati validi imponendo la doppia condizione:

$$\text{data} == \text{pattern} \quad \land \quad \text{valid} == 1$$

---

## 3. Soluzione: Funzionamento del Circuito e Uncomputation

### Inizializzazione
Prima di invocare il caricatore QROM, il sistema viene preparato nello stato iniziale:
* Il registro `idx` viene posto in una **sovrapposizione uniforme** applicando porte di Hadamard $(H^{\otimes n_{\text{idx}}})$ a tutti i suoi qubit.
* I registri `data` e `valid` partono inizializzati nello stato fondamentale $|0\dots 0\rangle_{\text{data}} \otimes |0\rangle_{\text{valid}}$.

$$|\psi_0\rangle = \frac{1}{\sqrt{2^{n_{\text{idx}}}}} \sum_{i=0}^{2^{n_{\text{idx}}}-1} |i\rangle_{\text{idx}} \otimes |0\dots 0\rangle_{\text{data}} \otimes |0\rangle_{\text{valid}}$$

### Correlazione indice-finestra tramite entanglement
Per associare ogni indice $|i\rangle$ alla rispettiva finestra $|w_i\rangle$ senza alterare gli altri rami della sovrapposizione, il loader esegue un ciclo iterativo per ogni finestra reale $i \in \{0, \dots, K-1\}$ articolato in tre fasi:

1. **Fase di Apertura:** si applicano porte NOT (`qc.x`) sui qubit del registro `idx` che corrispondono ai bit uguali a `0` nella rappresentazione binaria di $i$. In questo modo, **solo** il ramo corrispondente allo stato $|i\rangle$ viene temporaneamente trasformato in $|11\dots 1\rangle$, mentre tutti gli altri stati della sovrapposizione mantengono almeno uno `0`.
2. **Fase di Scrittura:** si utilizzano porte *Multi-Controlled NOT* (`qc.mcx`) controllate dall'intero registro `idx`. Poiché il controllo si attiva solo quando tutti i qubit di controllo sono a $|1\rangle$, i gate copiano i bit a `1` della finestra $w_i$ nel registro `data` e portano il qubit `valid` allo stato $|1\rangle$ **esclusivamente per il ramo attivo**.
3. **Fase di Chiusura:** si riapplicano le stesse porte `qc.x` della fase di apertura sul registro `idx`, riportando i qubit allo stato originale $|i\rangle$ e lasciando intatti tutti gli altri rami.

Al termine del ciclo, per la linearità della meccanica quantistica, si ottiene uno **stato entangled** in cui ogni indice è correlato in parallelo alla propria finestra di DNA:

$$|\psi_{\text{QROM}}\rangle = \frac{1}{\sqrt{2^{n_{\text{idx}}}}} \left( \sum_{i=0}^{K-1} |i\rangle_{\text{idx}} \otimes |w_i\rangle_{\text{data}} \otimes |1\rangle_{\text{valid}} + \sum_{i=K}^{2^{n_{\text{idx}}}-1} |i\rangle_{\text{idx}} \otimes |0\dots 0\rangle_{\text{data}} \otimes |0\rangle_{\text{valid}} \right)$$

### Reversibilità e Uncomputation nell'Oracolo di Grover
Nel calcolo quantistico non è possibile effettuare un *reset* irreversibile su un registro entangled (`data` e `valid`) senza far collassare la funzione d'onda, distruggendo la coerenza quantistica e impedendo l'interferenza costruttiva di Grover.

Per ovviare a questo problema, la funzione `build_dna_loader` incapsula l'intero circuito in un **Gate unitario** chiamato `QROM_DNA`, dotato di un inverso esatto richiamabile tramite il metodo `.inverse()`.

L'esecuzione completa dell'oracolo di Grover segue quindi tre passaggi fondamentali:

1. **Caricamento (`QROM_DNA`):** crea l'entanglement tra indici e finestre di DNA.
2. **Marcatura di Fase (Comparatore):** controlla i registri `data` e `valid` e applica un *phase kickback* (fase negativa $-1$) esclusivamente allo stato $|i^*\rangle$ che soddisfa il match con il pattern cercato.
3. **Pulizia (`QROM_DNA.inverse()` - *Uncomputation*):** riavvolge le operazioni della QROM riportando i registri `data` e `valid` allo stato pulito $|0\dots 0\rangle |0\rangle$, ma **preservando la fase $-1$** impressa sull'indice vincente $|i^*\rangle$:

$$|\psi_{\text{oracle}}\rangle = \frac{1}{\sqrt{2^{n_{\text{idx}}}}} \left( \sum_{i \neq i^*} |i\rangle_{\text{idx}} - |i^*\rangle_{\text{idx}} \right) \otimes |0\dots 0\rangle_{\text{data}} \otimes |0\rangle_{\text{valid}}$$

Una volta che il registro `idx` è stato completamente liberato dall'entanglement con `data` e `valid`, il **Diffusore di Grover** può operare correttamente su `idx`, amplificando l'ampiezza di probabilità dell'indice corretto $|i^*\rangle$ fino a portarla prossima al $100\%$.