# ID_a_dog_bread

## Project Overview
This project identifies the closest matching dog breed for a given **mystery DNA sequence** by performing **global sequence alignment**.  
It uses the **Biopython** library to read and compare DNA sequences from **FASTA files**. The p-value is calculated to show the similarity of the closest matching dog breed is statistically significant. It goes on to construct a phylogenetic tree using the **Neighbor Joining method**


## FIles in this Project
- `dog_breeds.fa` - Contains DNA sequences of different dog breeds.
- `mystery.fa` - Contains an unknown DNA sequence to compare.
- `Project_ID_Breed.py` - The main Python script.
- `README.md` - This file, explaining the project.

##  How it Works
1. **Reads** the mystery DNA sequence from a FASTA file.
2. **Loads** multiple dog breed DNA sequences from another FASTA file.
3. **Compares** the mystery sequence to each breed using **global sequence alignment**.
4. **Finds** the best-matching breed based on the alignment score.
5. **Calculates** the p-value to assess whether the similarity of the best matcing breed is statistically significant relative to other scores.
6. **Aligns** the dog breed sequences
7. **Computes** the distance matrix, using **the identity model** to calculate the distances.
8. **Constructs** a phylogenetic tree using the **Neighbor Joining method**

9. **Outputs**  highlights the sequence differences, the closest match sequence and description, the p-value and simlarity score and a phylogenetic tree

## Requirements
Ensure you have the following libraries: 
pip install biopython 
pip install scipy
pip intsall matplotlib

## Function descriptions 
**read_mystery_seq(mystery_fa)** Reads the mystery sequence from a FASTA file and returns it as a SeqRecord object.
**read_dog_breeds(dog_breeds_fa)**Reads the dog breed sequences from a FASTA file and returns them as a dictionary.
**find_most_similar_breed(mystery_sequence, dog_breeds)**Finds the most similar breed using pairwise sequence alignment and statistical analysis.
**align_sequences(dog_breeds)**Performs global alignment on the dog breed sequences and returns a MultipleSeqAlignment object.
**compute_distance_matrix(alignment)**Computes the pairwise distance matrix using an identity-based model.
**construct_phylogenetic_tree(distance_matrix)**Builds a phylogenetic tree using the Neighbor-Joining method and saves it in Newick format.
**visualize_tree(tree)**Generates a visualization of the phylogenetic tree using Matplotlib.
**main()**Executes the entire workflow from reading data to visualizing results.

## Expected Output

1. Identifies the most similar dog breed to the mystery sequence and prints:
2. The closest matching breed.
3. The alignment score.
4. A statistical p-value to assess significance.
5. Differences between the sequences.
6. Constructs and saves a phylogenetic tree.
7. Displays a visual representation of the phylogenetic tree.