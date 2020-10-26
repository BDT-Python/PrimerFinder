# Processes input. First makes it lowercase, then removes the spaces to get a pure sequence.
rawSequence = "gcggcgcgcc tgggcgctaa gatggcggcg gcgtgagttg catgttgtgt gaggatcccggggccgccgc gtcgctcggg ccccgcc atg gccgtcacca tcacgctcaa aacgctgcagcagcagacct tcaagatccg catggagcct gacgagacgg tgaaggtgct aaaggagaagatagaagctg agaagggtcg tgatgccttc cccgtggctg gacagaaact catctatgccggcaagatct tgagtgacga tgtccctatc agggactatc gcatcgatga gaagaactttgtggtcgtca tggtgaccaa gaccaaagcc ggccagggta cctcagcacc cccagaggcctcacccacag ctgccccaga gtcctctaca tccttcccgc ctgcccccac ctcaggcatgtcccatcccc cacctgccgc cagagaggac aagagcccat cagaggaatc cgcccccacgacgtccccag agtctgtgtc aggctctgtt ccctcttcag gtagcagcgg gcgagaggaagacgcggcct ccacgctagt gacgggctct gagtatgaga cgatgctgac ggagatcatgtccatgggct atgagcgaga gcgggtcgtg gccgccctga gagccagcta caacaacccccaccgagccg tggagtatct gctcacggga attcctggga gccccgagcc ggaacacggttctgtccagg agagccaggt atcggagcag ccggccacgg aagcagcagg agagaaccccctggagttcc tgcgggacca gccccagttc cagaacatgc ggcaggtgat tcagcagaaccctgcgctgc tgcccgccct gctccagcag ctgggccagg agaaccctca gcttttacagcaaatcagcc ggcaccagga gcagttcatc cagatgctga acgagccccc tggggagctggcggacatct cagatgtgga gggggaggtg ggcgccatag gagaggaggc cccgcagatgaactacatcc aggtgacgcc gcaggagaaa gaagctatag agaggttgaa ggccctgggcttcccagaga gcctggtcat ccaggcctat ttcgcgtgtg aaaaaaatga gaacttggctgccaacttcc tcctgagtca gaactttgat gacgagtgat gccaggaagc caggccaccgaagcccccac cctaccctta ttccatgaaa gttttataaa agaaaaaata tatatatattcatgtttatt taagaaatgg aaaaaaaaat caaaaatctt aaaaaaacaa gcaaacagtccagcttcctg tcctcctaaa gtggcccctg ttcccatctc ccgggccaga cagctgtccccccgtcctcc tccccagccc agcctgctca gagaagctgg caggactggg aggcgacagatgggcccctc ttggcctctg tcccagctct ctgcagccag acggaaaggc ggctgcttgcctctccatcc tccgaaaaac ccctgaggac ccccccccat cctcttctag gatgaggggaagctggagcc ccaactttga tcctccattg gagtggccca aatctttcca tctagggcaagtcctgaaag gcccaaggcc ccctccccag tctggccttg gcctccagcc tggagaagggctaacatcag ctcattgtca aggccacccc caccccagaa cagaaccgtg tctctgataaaggttttgaa gtgaataaag ttttaaaaac ta"
lowerSequence = rawSequence.lower()
noSpace = rawSequence.replace(" ", "")

nucleotideCounter = 0
primer_min_len = 15
primer_max_len = 19
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
        AT = sum([1 for x in primerSeq if x in ['a', 't']])  # Counts and calculates AT and GC%
        GC = j - AT

        # Calculates melting temperature according Marmur and Doty, 1962:
        tm = 64.9 + 41 * (GC - 16.4) / (AT + GC)            # Tm = 64 + 41 * (yG + zC - 16.4) / (wA + xT + yG + zC)

        percentage = 100 / j * GC                           # Calculates the GC percentage

        # Filters all sequence bites and filters them, then appends to, and stores them in the primer_dict dictionary.
        if 40 <= percentage <= 60 and 50 <= tm <= 60:
            d = {"Begin": nucleotideCounter,                # Creates dictionary 'd' with primer information
                 "End": nucleotideCounter + len(primerSeq),
                 "Length:": len(primerSeq),
                 "Sequence:": primerSeq,
                 "GC%": (round(percentage, 2)),
                 "Tm": (round(tm, 1)), }
            primer_dict.append(d)                           # Appends and stores each primer dictionary in the list

# Finding the primer pairs in the primer dictionary (a reverse that matches the forward) which creates a product of minimum 30nt and has a Tm difference of 5 or less 
min_distance = 30
pairs = []


def make_reverse(toReverse):               # Function to fetch the reverse compliment of the reverse primer.
    normal = "acgt"
    compliment = "tgca"
    reversePrimer = str.maketrans(normal, compliment)
    primer_dict[x]["Sequence:"] = toReverse.translate(reversePrimer)[::-1]
    return


for i in range(len(primer_dict)):
    options = []
    for x in range(i, (len(primer_dict))):
        if abs(primer_dict[i]["Tm"] - primer_dict[x]["Tm"]) <= 5:
            if primer_dict[x]["Begin"] - primer_dict[i]["End"] >= min_distance:
                make_reverse(primer_dict[x]["Sequence:"])
                options.append(primer_dict[x])
    if options:
        pairs.append([primer_dict[i], options])

print(pairs)

