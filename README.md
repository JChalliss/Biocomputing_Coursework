# ID_a_dog_bread

## Project Overview
This project identifies the closest matching dog breed for a given **mystery DNA sequence** by performing **global sequence alignment**.  
It uses the **Biopython** library to read and compare DNA sequences from **FASTA files**.

## FIles in this Project
- `dog_breeds.fa` - Contains DNA sequences of different dog breeds.
- `mystery.fa` - Contains an unknown DNA sequence to compare.
- `alignment_script.py` - The main Python script for sequence alignment.
- `README.md` - This file, explaining the project.

##  How it Works
1. **Reads** the mystery DNA sequence from a FASTA file.
2. **Loads** multiple dog breed DNA sequences from another FASTA file.
3. **Compares** the mystery sequence to each breed using **global sequence alignment**.
4. **Finds** the best-matching breed based on the alignment score.
5. **Calculates** the  p-value to assess whether the similarity of the best matcing breed is statistically significant relative to other scores.

5. **Outputs** the closest match, highlights the sequence differences, the p-value and simlarity score 

## Requirements
Ensure you have the following libraries: 
'pip install biopython' 
'pip install scipy'
