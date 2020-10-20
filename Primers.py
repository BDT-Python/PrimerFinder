#processes input. First makes it lowercase, then removes the spaces to get a pure sequence.
rawSequence = "gcggcgcgcc tgggcgctaa gatggcggcg gcgtgagttg catgttgtgt gaggatcccggggccgccgc gtcgctcggg ccccgcc atg gccgtcacca tcacgctcaa aacgctgcagcagcagacct tcaagatccg catggagcct gacgagacgg tgaaggtgct aaaggagaagatagaagctg agaagggtcg tgatgccttc cccgtggctg gacagaaact catctatgccggcaagatct tgagtgacga tgtccctatc agggactatc gcatcgatga gaagaactttgtggtcgtca tggtgaccaa gaccaaagcc ggccagggta cctcagcacc cccagaggcctcacccacag ctgccccaga gtcctctaca tccttcccgc ctgcccccac ctcaggcatgtcccatcccc cacctgccgc cagagaggac aagagcccat cagaggaatc cgcccccacgacgtccccag agtctgtgtc aggctctgtt ccctcttcag gtagcagcgg gcgagaggaagacgcggcct ccacgctagt gacgggctct gagtatgaga cgatgctgac ggagatcatgtccatgggct atgagcgaga gcgggtcgtg gccgccctga gagccagcta caacaacccccaccgagccg tggagtatct gctcacggga attcctggga gccccgagcc ggaacacggttctgtccagg agagccaggt atcggagcag ccggccacgg aagcagcagg agagaaccccctggagttcc tgcgggacca gccccagttc cagaacatgc ggcaggtgat tcagcagaaccctgcgctgc tgcccgccct gctccagcag ctgggccagg agaaccctca gcttttacagcaaatcagcc ggcaccagga gcagttcatc cagatgctga acgagccccc tggggagctggcggacatct cagatgtgga gggggaggtg ggcgccatag gagaggaggc cccgcagatgaactacatcc aggtgacgcc gcaggagaaa gaagctatag agaggttgaa ggccctgggcttcccagaga gcctggtcat ccaggcctat ttcgcgtgtg aaaaaaatga gaacttggctgccaacttcc tcctgagtca gaactttgat gacgagtgat gccaggaagc caggccaccgaagcccccac cctaccctta ttccatgaaa gttttataaa agaaaaaata tatatatattcatgtttatt taagaaatgg aaaaaaaaat caaaaatctt aaaaaaacaa gcaaacagtccagcttcctg tcctcctaaa gtggcccctg ttcccatctc ccgggccaga cagctgtccccccgtcctcc tccccagccc agcctgctca gagaagctgg caggactggg aggcgacagatgggcccctc ttggcctctg tcccagctct ctgcagccag acggaaaggc ggctgcttgcctctccatcc tccgaaaaac ccctgaggac ccccccccat cctcttctag gatgaggggaagctggagcc ccaactttga tcctccattg gagtggccca aatctttcca tctagggcaagtcctgaaag gcccaaggcc ccctccccag tctggccttg gcctccagcc tggagaagggctaacatcag ctcattgtca aggccacccc caccccagaa cagaaccgtg tctctgataaaggttttgaa gtgaataaag ttttaaaaac ta"
#lowerSequence = rawSequence.lower()
noSpace = rawSequence.replace(" ", "")

primer_min_len = 18
primer_max_len = 30
seq_len = len(rawSequence)

for i in range(seq_len):
    
    for j in range(primer_min_len,primer_max_len+1):
        
        if i + j > seq_len:
            break
        
        primer = rawSequence[i:i+j]
        AT = sum([1 for x in primer if x in ['a','t']])
        GC = j - AT
        tm = AT * 2 + GC * 4
        percentage = 100 / j * GC
        if percentage > 39 or percentage < 61 and tm > 64 or tm < 76:
                        print(tm)
                        print(primer)
                        print(round(percentage,2))

