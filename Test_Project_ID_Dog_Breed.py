import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator
import sys
import os
from Project_ID_Dog_Breed import read_mystery_seq, read_dog_breeds, find_most_similar_breed, compute_distance_matrix

class TestDogBreedIdentification(unittest.TestCase):

    def test_read_mystery_seq(self):
        """Test if mystery sequence is read correctly"""
        mystery_seq = read_mystery_seq("test_mystery.fa")
        self.assertIsInstance(mystery_seq, SeqRecord)  # Must return a SeqRecord
        self.assertGreater(len(mystery_seq.seq), 0)  # Ensure sequence is not empty

    def test_read_dog_breeds(self):
        """Test if dog breeds sequences are read correctly"""
        dog_breeds = read_dog_breeds("test_dog_breeds.fa")
        self.assertIsInstance(dog_breeds, dict)  # Must return a dictionary
        self.assertGreater(len(dog_breeds), 0)  # Ensure it's not empty

    def test_find_most_similar_breed(self):
        """Test similarity function with mock sequences"""
        mystery_seq = SeqRecord(Seq("ATCGATCG"))
        dog_breeds = {
            "Breed 1": "ATCGATCG",
            "Breed 2": "ATCGTTTT",
        }
        most_similar, _, p_value = find_most_similar_breed(mystery_seq, dog_breeds)
        self.assertEqual(most_similar, "Breed 1")  # Must match the most similar sequence
        self.assertLess(p_value, 0.05)  # p-value should be significant

    def test_compute_distance_matrix(self):
        """Test distance matrix calculation"""
        alignment = MultipleSeqAlignment([
            SeqRecord(Seq("ATCG"), id="Breed 1"),
            SeqRecord(Seq("ATGG"), id="Breed 2")
        ])
        distance_matrix = compute_distance_matrix(alignment)
        self.assertAlmostEqual(distance_matrix[0, 1], 0.25)  # Example expected distance

if __name__ == "__main__":
python -m unittest Test_Project_ID_Dog_Breed.py
