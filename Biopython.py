#importing package & basics of DNA Sequence
from Bio.Seq import Seq
dna= Seq("ATGCATGCATGCATGC")
print(dna)
print(dna.complement())

print(dna.reverse_complement())

len(dna)

dna.startswith("ATG")

dna.endswith("TAG")

# Count nucleotide

from Bio.Seq import Seq

dna = Seq("ATGCGCGTTAA")

print("A =", dna.count("A"))
print("T =", dna.count("T"))
print("G =", dna.count("G"))
print("C =", dna.count("C"))

# Generate complementary DNA
from Bio.Seq import Seq

dna = Seq("ATGCCGTA")

complement = dna.complement()

print("Complement :", complement)

# Molecular Weight of DNA
from Bio.SeqUtils import molecular_weight

sequence = "ATGCGTA"

weight = molecular_weight(sequence, "DNA")

print("Molecular Weight :", weight)



# Finding GC Content

dna = '''ATGTCTGATAATGGACCCCAAAATCAGCGAAATGCACCCCGCATTACGTTTGGTGGACCCTCAGATTCAA
CTGGCAGTAACCAGAATGGAGAACGCAGTGGGGCGCGATCAAAACAACGTCGGCCCCAAGGTTTACCCAA
TAATACTGCGTCTTGGTTCACCGCTCTCACTCAACATGGCAAGGAAGACCTTAAATTCCCTCGAGGACAA
GGCGTTCCAATTAACACCAATAGCAGTCCAGATGACCAAATTGGCTACTACCGAAGAGCTACCAGACGAA
TTCGTGGTGGTGACGGTAAAATGAAAGATCTCAGTCCAAGATGGTATTTCTACTACCTAGGAACTGGGCC
AGAAGCTGGACTTCCCTATGGTGCTAACAAAGACGGCATCATATGGGTTGCAACTGAGGGAGCCTTGAAT
ACACCAAAAGATCACATTGGCACCCGCAATCCTGCTAACAATGCTGCAATCGTGCTACAACTTCCTCAAG
GAACAACATTGCCAAAAGGCTTCTACGCAGAAGGGAGCAGAGGCGGCAGTCAAGCCTCTTCTCGTTCCTC
ATCACGTAGTCGCAACAGTTCAAGAAATTCAACTCCAGGCAGCAGTAGGGGAACTTCTCCTGCTAGAATG
GCTGGCAATGGCGGTGATGCTGCTCTTGCTTTGCTGCTGCTTGACAGATTGAACCAGCTTGAGAGCAAAA
TGTCTGGTAAAGGCCAACAACAACAAGGCCAAACTGTCACTAAGAAATCTGCTGCTGAGGCTTCTAAGAA
GCCTCGGCAAAAACGTACTGCCACTAAAGCATACAATGTAACACAAGCTTTCGGCAGACGTGGTCCAGAA
CAAACCCAAGGAAATTTTGGGGACCAGGAACTAATCAGACAAGGAACTGATTACAAACATTGGCCGCAAA
TTGCACAATTTGCCCCCAGCGCTTCAGCGTTCTTCGGAATGTCGCGCATTGGCATGGAAGTCACACCTTC
GGGAACGTGGTTGACCTACACAGGTGCCATCAAATTGGATGACAAAGATCCAAATTTCAAAGATCAAGTC
ATTTTGCTGAATAAGCATATTGACGCATACAAAACATTCCCACCAACAGAGCCTAAAAAGGACAAAAAGA
AGAAGGCTGATGAAACTCAAGCCTTACCGCAGAGACAGAAGAAACAGCAAACTGTGACTCTTCTTCCTGC
TGCAGATTTGGATGATTTCTCCAAACAATTGCAACAATCCATGAGCAGTGCTGACTCAACTCAGGCCTAA'''

# (g+c)/total nucleotide *100

g_count=dna.count("G")
c_count=dna.count("C")
total_gc=g_count+ c_count
total_nucleotides= len(dna)
gc_percent=(total_gc/total_nucleotides)*100
print(f"The GC percent of our dna sequence is (round(gc_percent,2))%")

# Transcription
dna 
transcript= dna.maketrans("ATGC","UACG")

print(dna.translate(transcript))

# Codon Search
sequence= ("ATGCACGCTACGTCGACGGCTAGCGGAGCTACGATAGCCGTCGATCGACGGCTAGCGAGACGCGTCGATCGATAGAGCGGGAGTAGCGATCAGCAGCGACTAGCGCCGTTCGAGTCGAGTCGCCAGACTAGCAGGCTATGGCATGTGACGG")
codon_list=["ATG","TAC"]
cdn_position=[]
seq_length=len(sequence)

# loop through sequence
for i in range(0, seq_length - 2):

    # extract 3-base codon
    possible_cdn = sequence[i:i+3]

    # check if codon present in codon_list
    if possible_cdn in codon_list:

        cdn_position.append((possible_cdn, i))

# print output
print("Codon Positions:\n")

for codon, position in cdn_position:
    print(codon, "found at position", position)
    

#Sequence alignment

from Bio import SeqIO
for sequence in SeqIO.parse("Chicken.fasta", "fasta"):
    print(sequence.id)
    print(repr(sequence.seq))
    print(len(sequence))

from Bio import Align
aligner = Align.PairwiseAligner(match_score=1.0)
Seq1="ATGC" 
Seq2="GATC" 
score = aligner.score(Seq1,Seq2)
print ("THe Score of the alignment is {score}")   
alignments=aligner.align(Seq1,Seq2)
for alignment in alignments:
    print(alignment)
    aligner.mode='local'
    Seq3="ACGTAT"
    Seq4="GTACT"
    scorelocal = aligner.align(Seq3,Seq4)
    print("The local alignment score is {scorelocal} ")

    alignLocal=aligner.align(Seq3,Seq4)
    for x in alignLocal:
        print(x)


# Pairwise Sequence alignment
from Bio import pairwise2

seq1 = "ATGCTAG"
seq2 = "ATGCGAG"

alignments = pairwise2.align.globalxx(seq1, seq2)

for alignment in alignments:
    print(alignment)


# Local Sequence Alignment

    from Bio import pairwise2

seq1 = "ATGCGTA"
seq2 = "GCG"

alignment = pairwise2.align.localxx(seq1, seq2)

print(alignment)


# Protein Sequence Analysis

from Bio.SeqUtils.ProtParam import ProteinAnalysis

protein = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLR"

analysis = ProteinAnalysis(protein)

print("Molecular Weight :", analysis.molecular_weight())
print("Aromaticity :", analysis.aromaticity())
print("Isoelectric Point :", analysis.isoelectric_point())

#NCBI Sequence downloading
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "shindevaishnavi965@gmail.com"

handle = Entrez.efetch(db="nucleotide", id="AC004528.1", rettype="fasta", retmode= "text")                    

record = SeqIO.read(handle, "fasta")

print(record.id)
print(record.seq)


