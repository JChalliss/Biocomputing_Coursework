
# PairwiseAlinger is imported to align and compare the sequences and scipy imported to provide statistical functions
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from scipy import stats
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
    
    print(f"\nThe closest sequence is: {most_similar_description}")
    print(f"Similarity Score: {top_score}")
    print(f"P-value: {p_value}")

    return most_similar_description, most_similar_sequence, p_value
    
find_most_similar_breed(read_mystery_seq(mystery_fa), read_dog_breeds(dog_breeds_fa)) 