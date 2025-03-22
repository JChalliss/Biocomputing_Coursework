#variables to hold my file pathways
#using libraru biopython, importing SeqIO module to read the files 
from Bio import SeqIO
from Bio.Align import PairwiseAligner

dog_breeds_fa= r"C:\Users\jadec\OneDrive\biocomputing\Coursework Project\project_dog_dna\dog_breeds.fa"
mystery_fa= r"C:\Users\jadec\OneDrive\biocomputing\Coursework Project\project_dog_dna\mystery.fa"

def read_mystery_seq(mystery_fa):
    SeqIO.read(mystery_fa, "fasta") 
    return()

# processes a FASTA file containing information about dog breeds in its headers and builds a dictionary 
# where the keys are dog breed names and the values are the corresponding DNA/protein sequences.
def read_dog_breeds(dog_breeds_fa):
    dog_breeds = {}
    for record in SeqIO.parse(dog_breeds_fa, "fasta"):
        description_parts = record.description.split()
        breed = description_parts[9] if len(description_parts) > 7 else "Unknown"
        dog_breeds[breed] = str(record.seq)
    return dog_breeds


def find_most_similar_breed(mystery_sequence, dog_breeds):
    aligner = PairwiseAligner()
    aligner.mode = "global"  
    aligner.match_score = 1  
    aligner.mismatch_score = 0   
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -1
    
    top_score = float('-inf')
    most_similar_breed = None
    most_similar_sequence = None
    
    for breed, breed_sequence in dog_breeds.items():
        alignments = aligner.align(mystery_sequence, breed_sequence)
        score = alignments[0].score  
        
        if score > max_score:
            max_score = score
            most_similar_breed = breed
            most_similar_sequence = breed_sequence
    return most_similar_breed, most_similar_sequence

print(find_most_similar_breed(read_mystery_seq(mystery_fa), read_dog_breeds(dog_breeds_fa)))