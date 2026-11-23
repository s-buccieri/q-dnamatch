import math

# Dizionario per la codifica a 2 bit
DNA_MAPPING = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}

def validate_and_encode(sequence: str) -> str:
    """Valida la sequenza e la converte in una stringa binaria a 2 bit."""
    if not sequence:
        raise ValueError("La sequenza non può essere vuota.")
    sequence = sequence.upper()
    encoded = ""
    for base in sequence:
        if base not in DNA_MAPPING:
            raise ValueError(f"Base non valida trovata: {base}. Sono ammesse solo A, C, G, T.")
        encoded += DNA_MAPPING[base]
    return encoded

def classical_linear_search(windows_bin: list[str], pattern_bin: str, valid_len: int) -> dict:
    """Esegue la ricerca classica lineare per stabilire la baseline O(N)."""
    indices = []
    queries = 0
    # Iteriamo solo sulle finestre valide, ignorando il padding
    for i in range(valid_len):
        queries += 1     # Ogni confronto conta come una query
        if windows_bin[i] == pattern_bin:
            indices.append(i)
    
    return {"indices": indices, "queries": queries}

def process_dna(seq: str, pattern: str) -> dict:
    """
    Funzione principale che genera il DNAPayload richiesto dai contratti API.
    """
    # 1. Validazione e codifica
    seq_bin_full = validate_and_encode(seq)
    pattern_bin = validate_and_encode(pattern)
    
    N = len(seq)
    M = len(pattern)
    
    if M > N:
        raise ValueError("Il pattern non può essere più lungo della sequenza di DNA.")
    
    # 2. Sliding window (Estrazione finestre)
    valid_len = N - M + 1
    windows_bin = []
    
    # Estraiamo finestre di lunghezza M (che in binario sono 2*M bit)
    bit_window_size = 2 * M
    for i in range(valid_len):
        # L'indice di partenza in bit è i * 2
        start_bit = i * 2
        window = seq_bin_full[start_bit : start_bit + bit_window_size]
        windows_bin.append(window)
        
    # 3. Calcolo qubit indice e padding
    n_idx = max(1, math.ceil(math.log2(valid_len)))
    total_windows_needed = 2 ** n_idx
    
    # Applichiamo il padding fino a 2^n_idx
    padding_needed = total_windows_needed - valid_len
    # Usiamo una stringa di zeri come padding (non influirà grazie a valid_len)
    dummy_padding = '0' * bit_window_size
    windows_bin.extend([dummy_padding] * padding_needed)
    
    # 4. Ricerca classica per baseline
    n_data = bit_window_size
    classical_res = classical_linear_search(windows_bin, pattern_bin, valid_len)
    
    # 5. Restituzione del contratto API (DNAPayload)
    return {
        "windows_bin": windows_bin,
        "pattern_bin": pattern_bin,
        "n_idx": n_idx,
        "n_data": n_data,
        "valid_len": valid_len,
        "classical_res": classical_res
    }