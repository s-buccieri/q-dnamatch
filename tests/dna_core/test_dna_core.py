import pytest
from src.dna_core import process_dna, validate_and_encode

def test_validation():
    # Test validazione corretta e case insensitive
    assert validate_and_encode("aCgT") == "00011011"
    
    # Test errore su stringa vuota
    with pytest.raises(ValueError):
        validate_and_encode("")
    # Test errore su stringa vuota
    with pytest.raises(ValueError):
        validate_and_encode("ACGTX")
        

def test_normal_match():
    # Sequenza: ACGTAC (6 basi) - Pattern: GT (2 basi)
    # N-M+1 = 6 - 2 + 1 = 5 finestre valide. 
    # Qubit indice: ceil(log2(5)) = 3 qubit -> 8 finestre totali (3 di padding)
    payload = process_dna("ACGTAC", "GT")
    
    assert payload["n_idx"] == 3
    assert payload["valid_len"] == 5
    assert len(payload["windows_bin"]) == 8 # 2^3
    assert payload["pattern_bin"] == "1011" # G=10, T=11
    assert payload["n_data"] == 4 # 2 basi * 2 bit
    
    # GT si trova all'indice 2 (A=0, C=1, G=2)
    assert payload["classical_res"]["indices"] == [2]
    assert payload["classical_res"]["queries"] == 5

def test_multiple_occurrences():
    # Sequenza con "A" ripetute, cerchiamo "AA"
    payload = process_dna("AAAA", "AA")
    assert payload["valid_len"] == 3
    assert payload["classical_res"]["indices"] == [0, 1, 2] # Trovato in tutte le finestre

def test_pattern_not_found():
    payload = process_dna("ACGT", "TT")
    assert payload["classical_res"]["indices"] == [] # Lista vuota
    assert payload["classical_res"]["queries"] == 3 # Ha comunque fatto 3 confronti

def test_pattern_longer_than_seq():
    # Pattern più lungo della sequenza principale deve sollevare errore
    with pytest.raises(ValueError):
        process_dna("AC", "ACGT")
def test_single_window():
    # Caso limite N == M: una sola finestra valida
    payload = process_dna("ACGT", "ACGT")
    assert payload["valid_len"] == 1
    assert payload["n_idx"] == 1 # Almeno 1 qubit per il circuito quantistico
    assert len(payload["windows_bin"]) == 2 # 2^1
    assert payload["classical_res"]["indices"] == [0]