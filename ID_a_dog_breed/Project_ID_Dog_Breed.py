from Bio import SeqIO
from Bio.Align import PairwiseAligner
from scipy import stats
from Bio.Phylo import TreeConstruction
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
import matplotlib.pyplot as plt
"""
    Reads the mystery sequence from the provided FASTA file.
    
    Parameters:
    mystery_fa (str): File path to the mystery sequence in FASTA format.
    
    Returns:
    SeqRecord: The sequence record of the mystery sequence.
    """
dog_breeds_fa= r"C:\Users\jadec\OneDrive\biocomputing\Coursework Project\project_dog_dna\dog_breeds.fa"
mystery_fa= r"C:\Users\jadec\OneDrive\biocomputing\Coursework Project\project_dog_dna\mystery.fa"

def read_mystery_seq(mystery_fa):
    mystery_sequence= SeqIO.read(mystery_fa, "fasta") 
    return mystery_sequence


"""
    Reads the dog breeds from the provided FASTA file and builds a dictionary
    where the keys are the descriptions (headers) and the values are the sequences.
    
    Parameters:
    dog_breeds_fa (str): File path to the dog breeds sequences in FASTA format.
    
    Returns:
    dict: A dictionary with descriptions as keys and sequences as values.
    """
def read_dog_breeds(dog_breeds_fa):
    dog_breeds = {}
    for record in SeqIO.parse(dog_breeds_fa, "fasta"):
        description = record.description
        dog_breeds[description] = str(record.seq)
    return dog_breeds

"""
    Finds the most similar dog breed to the mystery sequence by performing a
    global sequence alignment and calculating the alignment score. Additionally,
    it compares the sequences, prints the differences, and calculates a p-value
    to assess the statistical significance of the similarity.
    
    Parameters:
    mystery_sequence (SeqRecord): The mystery DNA sequence.
    dog_breeds (dict): A dictionary containing dog breed names (as keys) and DNA sequences (as values).
    
    Returns:
    tuple: A tuple containing:
        - The description of the most similar dog breed (str)
        - The sequence of the most similar dog breed (str)
        - The p-value from the statistical test (float)
    """

# itialises 'PairwaiseAligner' with specific scoring parameters
# compares the mysetery sequence agaiast each breed's sequence using global alignment amd calculates alignment score
# tracks highest alignment score and corresponding sequence
# if alignment is found the most similar breed, it compares the sequences and prints the differences of the two.
# Performs a one-sample t-test (ttest_1samp) on the alignment scores with the top score as the sample.
# Computes a p-value to assess whether the similarity of the top match is statistically significant relative to other scores.

def find_most_similar_breed(mystery_sequence, dog_breeds):
    aligner = PairwiseAligner()
    aligner.mode = "global"  
    aligner.match_score = 1  
    aligner.mismatch_score = 0   
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    
    top_score = float('-inf')
    most_similar_description = None
    most_similar_sequence = None
    highest_alignment = None
    similarity_scores = []
    
    for description, breed_sequence in dog_breeds.items():
        alignments = aligner.align(mystery_sequence, breed_sequence)
        score = alignments[0].score  
        similarity_scores.append(score)
        
        if score > top_score:
            top_score = score
            most_similar_description = description
            most_similar_sequence = breed_sequence
            highest_alignment = alignments[0]
    if highest_alignment:
         aligned_seq1 = highest_alignment[0]  
         aligned_seq2 = highest_alignment[1]
    
    print("\nDifferences between 'mystery sequence' and most similar sequence match found:")
    for i, (char1, char2) in enumerate(zip(aligned_seq1, aligned_seq2)):
        if char1 != char2:
             if char1 == "-":
                print(f"Position {i+1}: gap in mystery sequence -> {char2}")
             elif char2 == "-":
                print(f"Position {i+1}: {char1} -> gap in breed sequence")
             else:
                print(f"Position {i+1}: {char1} -> {char2}")
                print(f"Position {i+1}: {char1} -> {char2}")  

    p_value = stats.ttest_1samp(similarity_scores, top_score).pvalue
    
    print(f"\nThe closest sequence is: {most_similar_description,most_similar_sequence}")
    print(f"Similarity Score: {top_score}")
    print(f"P-value: {p_value}")

    return most_similar_description, most_similar_sequence, p_value

"""
Aligns the dog breed sequences using global pairwise alignment and returns the
multiple sequence alignment.

Parameters:
dog_breeds (dict): A dictionary containing dog breed names (as keys) and DNA sequences (as values).

Returns:
MultipleSeqAlignment: An alignment of the dog breed sequences.
"""
# Align each breed sequence with itself (this could be replaced with a multiple sequence alignment)
# Wrap the aligned sequence in a SeqRecord object
def align_sequences(dog_breeds):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    
    sequences = []
    
    for description, breed_sequence in dog_breeds.items():
        
        aligned_sequence = aligner.align(breed_sequence, breed_sequence)
       
        seq_record = SeqRecord(Seq(str(aligned_sequence[0])), id=description)
        sequences.append(seq_record)
    
    alignment = MultipleSeqAlignment(sequences)
    return alignment
"""
Computes the distance matrix for the given multiple sequence alignment using
the identity model to calculate pairwise distances between sequences.

Parameters:
alignment (MultipleSeqAlignment): The multiple sequence alignment to compute distances for.

Returns:
DistanceMatrix: A matrix representing the pairwise distances between sequences.
"""
# Compute the distance matrix using a pairwise distance calculation
# Using the identity distance model
def compute_distance_matrix(alignment):
    calculator = DistanceCalculator('identity')  
    distance_matrix = calculator.get_distance(alignment)
    return distance_matrix
"""
Constructs a phylogenetic tree using the Neighbor Joining method from a distance matrix.

Parameters:
distance_matrix (DistanceMatrix): The matrix of pairwise distances between sequences.

Returns:
Phylo.BaseTree.Tree: The constructed phylogenetic tree.
"""
# Construct the phylogenetic tree from the distance matrix
# Using Neighbor Joining algorithm
def construct_phylogenetic_tree(distance_matrix):
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)  
    Phylo.write(tree, "phylogenetic_tree.nwk", "newick")
    return tree
"""
Visualizes the phylogenetic tree using Matplotlib and Bio.Phylo, adjusting the size
and collapsing branches with small lengths for better readability.

Parameters:
tree (Phylo.BaseTree.Tree): The phylogenetic tree to visualize.

Returns:
None
"""
# Visualize the phylogenetic tree
# Increase figure size
# Only collapse clades if they have valid branch lengths
# Ignore collapse errors if they occur
def visualize_tree(tree):
    fig = plt.figure(figsize=(20, 10))  
    ax = fig.add_subplot(1, 1, 1)
    for clade in tree.find_clades():  
        if clade.branch_length is not None and clade.branch_length < 0.0001:
            try:
                clade.collapse()
            except ValueError:
                pass 

    Phylo.draw(tree, axes=ax)
    return plt.show()


"""
Main function that orchestrates the process:
- Reads in the dog breed and mystery sequences.
- Finds the most similar breed to the mystery sequence.
- Aligns the dog breed sequences and generates the distance matrix.
- Constructs and visualizes the phylogenetic tree.

Returns:
None
"""
def main():
    
    dog_breeds = read_dog_breeds(dog_breeds_fa)
    mystery_sequence = read_mystery_seq(mystery_fa)

    find_most_similar_breed(mystery_sequence, dog_breeds)

    alignment = align_sequences(dog_breeds)
    distance_matrix = compute_distance_matrix(alignment)

    phylogenetic_tree = construct_phylogenetic_tree(distance_matrix)


    visualize_tree(phylogenetic_tree)


if __name__ == "__main__":
    main()
 