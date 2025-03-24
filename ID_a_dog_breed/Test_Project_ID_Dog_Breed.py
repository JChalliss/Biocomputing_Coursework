""" IT DOESNT WORK YET. NEED TO READ UP ON PYTEST"""


import pytest
from unittest.mock import patch, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from dog_breeds import read_mystery_seq, read_dog_breeds, find_most_similar_breed  # Import from your script

# Mock for reading the mystery sequence
@pytest.fixture
def mock_mystery_seq():
    return SeqRecord(Seq("ATCGTACGATCG"), id="mystery_sequence", description="Mystery Sequence")

# Mock for reading the dog breeds
@pytest.fixture
def mock_dog_breeds():
    return {
        "Breed1": "ATCGTACGATCG",
        "Breed2": "ATCGGACGATCG",
        "Breed3": "GTCAGTACGATC"
    }

# Test for the read_mystery_seq function
def test_read_mystery_seq(mock_mystery_seq):
    with patch("Bio.SeqIO.read", return_value=mock_mystery_seq):
        mystery_sequence = read_mystery_seq("fake_path.fasta")
        assert mystery_sequence.id == "mystery_sequence"
        assert str(mystery_sequence.seq) == "ATCGTACGATCG"
        assert mystery_sequence.description == "Mystery Sequence"

# Test for the read_dog_breeds function
def test_read_dog_breeds(mock_dog_breeds):
    with patch("Bio.SeqIO.parse", return_value=[SeqRecord(Seq(seq), id=breed, description=breed) for breed, seq in mock_dog_breeds.items()]):
        dog_breeds = read_dog_breeds("fake_dog_breeds.fasta")
        assert len(dog_breeds) == 3
        assert dog_breeds["Breed1"] == "ATCGTACGATCG"
        assert dog_breeds["Breed2"] == "ATCGGACGATCG"
        assert dog_breeds["Breed3"] == "GTCAGTACGATC"

# Test for the find_most_similar_breed function
def test_find_most_similar_breed(mock_mystery_seq, mock_dog_breeds):
    # Mock PairwiseAligner and its methods
    aligner_mock = MagicMock()
    aligner_mock.align.return_value = [MagicMock(score=1.0, seqs=[str(mock_mystery_seq.seq), mock_dog_breeds["Breed1"]])]

    with patch("Bio.Align.PairwiseAligner", return_value=aligner_mock):
        most_similar_description, most_similar_sequence, p_value = find_most_similar_breed(mock_mystery_seq, mock_dog_breeds)
        assert most_similar_description == "Breed1"
        assert most_similar_sequence == "ATCGTACGATCG"
        assert isinstance(p_value, float)  # Ensure p-value is a float
        assert p_value >= 0  # p-value should be >= 0

    # Check for the case where there is no alignment
    aligner_mock.align.return_value = []
    with patch("Bio.Align.PairwiseAligner", return_value=aligner_mock):
        most_similar_description, most_similar_sequence, p_value = find_most_similar_breed(mock_mystery_seq, mock_dog_breeds)
        assert most_similar_description is None
        assert most_similar_sequence is None
        assert p_value is None
