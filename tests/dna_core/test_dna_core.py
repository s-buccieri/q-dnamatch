import pytest
from src.dna_core import process_dna, validate_and_encode

def test_validation():
    
    assert validate_and_encode("aCgT") == "00011011"
    
    with pytest.raises(ValueError):
        validate_and_encode("")
    
    with pytest.raises(ValueError):
        validate_and_encode("ACGTX")
        

def test_normal_match():
    payload = process_dna("ACGTAC", "GT")
    
    assert payload["n_idx"] == 3
    assert payload["valid_len"] == 5
    assert len(payload["windows_bin"]) == 8
    assert payload["pattern_bin"] == "1011"
    assert payload["n_data"] == 4
    
    assert payload["classical_res"]["indices"] == [2]
    assert payload["classical_res"]["queries"] == 5

def test_multiple_occurrences():

    payload = process_dna("AAAA", "AA")
    assert payload["valid_len"] == 3
    assert payload["classical_res"]["indices"] == [0, 1, 2]

def test_pattern_not_found():
    payload = process_dna("ACGT", "TT")
    assert payload["classical_res"]["indices"] == []
    assert payload["classical_res"]["queries"] == 3

def test_pattern_longer_than_seq():
    with pytest.raises(ValueError):
        process_dna("AC", "ACGT")
def test_single_window():
    
    payload = process_dna("ACGT", "ACGT")
    assert payload["valid_len"] == 1
    assert payload["n_idx"] == 1
    assert len(payload["windows_bin"]) == 2
    assert payload["classical_res"]["indices"] == [0]