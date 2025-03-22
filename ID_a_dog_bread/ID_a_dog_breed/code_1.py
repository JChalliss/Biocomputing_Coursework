<<<<<<< HEAD
dog_breeds_fa= r"C:\Users\jadec\OneDrive\biocomputing\Coursework Project\project_dog_dna\dog_breeds.fa"
mystery_fa= r"C:\Users\jadec\OneDrive\biocomputing\Coursework Project\project_dog_dna\mystery.fa"


from Bio import SeqIO

mystery_seq = SeqIO.read(mystery_fa, "fasta")

def dbf_descritpitons(dog_breeds_fa):
    dog_breeds= SeqIO.parse(dog_breeds_fa, "fasta")
    for record in dog_breeds:
        breed = record
        print(f"Breed: {breed}")

      #  desc= record.description.split("[")
       # for part in desc:
       # if part.startswith("breed="):
        #    return(part)

print(dbf_descritpitons)
         #if str(mystery_seq.seq) == str(record.seq):
             #breed_name = record.description.split()["breed"]
    #return(breed_name)


    def dbf_descriptions(dog_breed_data):
    """Parses the FASTA file and extracts breed names."""
    dog_breeds = SeqIO.parse(dog_breed_data, "fasta")
    for record in dog_breeds:
        breed = extract_breed(record)
        print(f"Breed: {breed}")




        def dog_breed_compare(dog_breeds_fa):
    SeqIO.parse(dog_breeds_fa, "fasta")
    
def find_most_similar_breed(mystery_sequence, dog_breeds):
   max_score = float('-inf')
   most_similar_breed = None
   most_similar_sequence = None
   alignment = None
=======
dog_breeds_fa= r"C:\Users\jadec\OneDrive\biocomputing\Coursework Project\project_dog_dna\dog_breeds.fa"
mystery_fa= r"C:\Users\jadec\OneDrive\biocomputing\Coursework Project\project_dog_dna\mystery.fa"


from Bio import SeqIO

mystery_seq = SeqIO.read(mystery_fa, "fasta")

def dbf_descritpitons(dog_breeds_fa):
    dog_breeds= SeqIO.parse(dog_breeds_fa, "fasta")
    for record in dog_breeds:
        breed = record
        print(f"Breed: {breed}")

      #  desc= record.description.split("[")
       # for part in desc:
       # if part.startswith("breed="):
        #    return(part)

print(dbf_descritpitons)
         #if str(mystery_seq.seq) == str(record.seq):
             #breed_name = record.description.split()["breed"]
    #return(breed_name)


    def dbf_descriptions(dog_breed_data):
    """Parses the FASTA file and extracts breed names."""
    dog_breeds = SeqIO.parse(dog_breed_data, "fasta")
    for record in dog_breeds:
        breed = extract_breed(record)
        print(f"Breed: {breed}")




        def dog_breed_compare(dog_breeds_fa):
    SeqIO.parse(dog_breeds_fa, "fasta")
    
def find_most_similar_breed(mystery_sequence, dog_breeds):
   max_score = float('-inf')
   most_similar_breed = None
   most_similar_sequence = None
   alignment = None
>>>>>>> a7d1fc05acf6fcff6358aad5dee319f2f97d5ec2
