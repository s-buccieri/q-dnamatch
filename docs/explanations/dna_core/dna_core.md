# Spiegazione Modulo: DNA Core (Baseline Classica ed Encoding)

## 1. Cosa fa questo modulo
Il modulo `dna_core.py` rappresenta il punto di ingresso della pipeline di **Q-DNAMatch**. 

Il suo scopo principale è:
1. Tradurre il problema biologico del pattern matching in un formato binario ottimizzato e compatibile con i registri quantistici (**QROM**).
2. Fornire una **baseline classica lineare O(N)** con conteggio esatto delle query per dimostrare empiricamente il vantaggio quadratico dell'algoritmo di **Grover O(√N)**.

---

## 2. Come funziona (Le tre fasi dell'algoritmo)

### Fase 1: Validazione e Codifica a 2-bit (Encoding)
Le sequenze biologiche ammettono esclusivamente le 4 basi azotate (A, C, G, T). Ciascuna base viene mappata in modo univoco su 2 bit:
* **A** → `00`
* **C** → `01`
* **G** → `10`
* **T** → `11`

> **Efficienza sui Qubit:** La codifica a 2 bit è la più compatta possibile per 4 simboli. Un pattern di M basi richiederà esattamente `n_data = 2 * M` qubit per il registro dati del circuito quantistico.

### Fase 2: Sliding Window e Padding Quantistico
Data una sequenza genomica di lunghezza N e un pattern target di lunghezza M (con M ≤ N):
1. **Estrazione finestre:** Estraiamo tutte le possibili sottostringhe consecutive di lunghezza M con passo unitario (1 base = 2 bit). Il numero di finestre valide è:
   `valid_len = N - M + 1`
2. **Padding a potenza di 2:** Il registro quantistico di indirizzamento della QROM deve essere posto in sovrapposizione uniforme tramite porte di Hadamard, richiedendo una dimensione pari a una potenza di due (`2^n_idx`):
   `n_idx = max(1, ceil(log2(valid_len)))`
3. **Finestre dummy:** Se `valid_len < 2^n_idx`, completiamo la lista con finestre fittizie (*dummy padding* di soli zeri) per non alterare la ricerca e garantire la corretta dimensione del registro quantistico.

### Fase 3: Baseline Classica e Metriche
Eseguiamo una ricerca lineare esaustiva (*classical linear search*) confrontando il pattern con le sole finestre valide (ignorando il padding):
* Ogni confronto effettuato incrementa di 1 il contatore delle `queries`.
* Il numero totale di query classiche sarà sempre esattamente pari a `valid_len`, confermando la complessità teorica:
  `T_classico = O(N)`
* Questo dato viene passato al `runtime-analytics` per tracciare il grafico comparativo di scalabilità rispetto alle iterazioni di Grover:
  `T_quantistico = O(√N)`

---

## 3. Perché è strutturato così (Approccio Interface-First)
Seguendo la metodologia *Interface-First*, il modulo espone un'unica interfaccia di alto livello `process_dna(seq, pattern)`. 

La funzione restituisce il dizionario standardizzato **`DNAPayload`**, che rispetta rigorosamente il contratto API concordato tra le coppie:

```python
{
    "windows_bin": list[str],  # Lista di 2^n_idx stringhe binarie (ciascuna di 2*M bit)
    "pattern_bin": str,        # Pattern target codificato a 2*M bit
    "n_idx": int,              # Numero di qubit di indice necessari (log2 delle finestre totali)
    "n_data": int,             # Numero di qubit dati necessari (2 * M)
    "valid_len": int,          # Finestre biologiche reali escluso il padding
    "classical_res": {
        "indices": list[int],  # Indici in cui è stato trovato il match
        "queries": int         # Numero esatto di controlli eseguiti (pari a valid_len)
    }
}
```

Questo payload alimenta direttamente:
* **`q_loader.py`:** per caricare le finestre nella porta composita `QROM_DNA`.
* **`q_grover.py`:** per impostare la condizione di fase dell'oracolo MCZ su `pattern_bin`.
* **`analytics.py`:** per confrontare le `queries` classiche con il benchmark quantistico.
* **`app.py`:** per mostrare i risultati e l'allineamento nella web app Streamlit.
