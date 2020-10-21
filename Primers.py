# Processes input. First makes it lowercase, then removes the spaces to get a pure sequence.
rawSequence = "gcggcgcgcc tgggcgctaa gatggcggcg gcgtgagttg catgttgtgt gaggatcccggggccgccgc gtcgctcggg ccccgcc atg gccgtcacca tcacgctcaa aacgctgcagcagcagacct tcaagatccg catggagcct gacgagacgg tgaaggtgct aaaggagaagatagaagctg agaagggtcg tgatgccttc cccgtggctg gacagaaact catctatgccggcaagatct tgagtgacga tgtccctatc agggactatc gcatcgatga gaagaactttgtggtcgtca tggtgaccaa gaccaaagcc ggccagggta cctcagcacc cccagaggcctcacccacag ctgccccaga gtcctctaca tccttcccgc ctgcccccac ctcaggcatgtcccatcccc cacctgccgc cagagaggac aagagcccat cagaggaatc cgcccccacgacgtccccag agtctgtgtc aggctctgtt ccctcttcag gtagcagcgg gcgagaggaagacgcggcct ccacgctagt gacgggctct gagtatgaga cgatgctgac ggagatcatgtccatgggct atgagcgaga gcgggtcgtg gccgccctga gagccagcta caacaacccccaccgagccg tggagtatct gctcacggga attcctggga gccccgagcc ggaacacggttctgtccagg agagccaggt atcggagcag ccggccacgg aagcagcagg agagaaccccctggagttcc tgcgggacca gccccagttc cagaacatgc ggcaggtgat tcagcagaaccctgcgctgc tgcccgccct gctccagcag ctgggccagg agaaccctca gcttttacagcaaatcagcc ggcaccagga gcagttcatc cagatgctga acgagccccc tggggagctggcggacatct cagatgtgga gggggaggtg ggcgccatag gagaggaggc cccgcagatgaactacatcc aggtgacgcc gcaggagaaa gaagctatag agaggttgaa ggccctgggcttcccagaga gcctggtcat ccaggcctat ttcgcgtgtg aaaaaaatga gaacttggctgccaacttcc tcctgagtca gaactttgat gacgagtgat gccaggaagc caggccaccgaagcccccac cctaccctta ttccatgaaa gttttataaa agaaaaaata tatatatattcatgtttatt taagaaatgg aaaaaaaaat caaaaatctt aaaaaaacaa gcaaacagtccagcttcctg tcctcctaaa gtggcccctg ttcccatctc ccgggccaga cagctgtccccccgtcctcc tccccagccc agcctgctca gagaagctgg caggactggg aggcgacagatgggcccctc ttggcctctg tcccagctct ctgcagccag acggaaaggc ggctgcttgcctctccatcc tccgaaaaac ccctgaggac ccccccccat cctcttctag gatgaggggaagctggagcc ccaactttga tcctccattg gagtggccca aatctttcca tctagggcaagtcctgaaag gcccaaggcc ccctccccag tctggccttg gcctccagcc tggagaagggctaacatcag ctcattgtca aggccacccc caccccagaa cagaaccgtg tctctgataaaggttttgaa gtgaataaag ttttaaaaac ta"
lowerSequence = rawSequence.lower()
noSpace = rawSequence.replace(" ", "")

nucleotideCounter = 0
primer_min_len = 18
primer_max_len = 30
seq_len = len(noSpace)

# Defines primer_dict as a list to hold libraries with potential primer information.
primer_dict = []

for i in range(seq_len):
    # Counts at which nucleotide the printed primer starts
    nucleotideCounter += 1
    for j in range(primer_min_len, primer_max_len + 1):

        if i + j > seq_len:
            break

        primerSeq = noSpace[i:i + j]
        AT = sum([1 for x in primerSeq if x in ['a', 't']]) # Counts and calculates AT and GC%
        GC = j - AT

        # Calculates melting temperature according Marmur and Doty, 1962:
        tm = 64.9 + 41 * (GC - 16.4) / (AT + GC)            # Tm = 64 + 41 * (yG + zC - 16.4) / (wA + xT + yG + zC)

        percentage = 100 / j * GC                           # Calculates the GC percentage

        # Filters all sequence bites and filters them, then appends to, and stores them in the primer_dict dictionary.
        if 40 <= percentage <= 60 and 50 <= tm <= 60:
            d = {"Begin": nucleotideCounter,
                 "End": nucleotideCounter + len(primerSeq),
                 "Length:": len(primerSeq),
                 "Sequence:": primerSeq,
                 "GC%": percentage,
                 "Tm": tm}
            primer_dict.append(d)

#def product_length:                 # function for fetching primers for desired length
#    for primer in primer_dict:      # Searching through every primer stored in dictionary


print(primer_dict)

# We have every single suitable primer saved into primer_dict.
# Now the following tasks remain:
#
#   1. Request product length from user and fetch two primers, roughly the desired length apart
#      This should be relatively easy, considering the primer start/end nt are saved in primer_dict.
#        https://medium.com/@rob3hr/searching-through-a-list-of-dictionaries-in-python-618fe77b2799
#
#   2. Make the reverse primer reverse complimentary. Also easy enough, turn A into T, G into C.
#
#   3. Check for self-annealing. This one's a bit more tricky but I'm sure we'll figure it out
#
