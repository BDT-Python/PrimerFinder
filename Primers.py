# processes input. First makes it lowercase, then removes the spaces to get a pure sequence.
rawSequence = "gcggcgcgcc tgggcgctaa gatggcggcg gcgtgagttg catgttgtgt gaggatcccggggccgccgc gtcgctcggg ccccgcc atg gccgtcacca tcacgctcaa aacgctgcagcagcagacct tcaagatccg catggagcct gacgagacgg tgaaggtgct aaaggagaagatagaagctg agaagggtcg tgatgccttc cccgtggctg gacagaaact catctatgccggcaagatct tgagtgacga tgtccctatc agggactatc gcatcgatga gaagaactttgtggtcgtca tggtgaccaa gaccaaagcc ggccagggta cctcagcacc cccagaggcctcacccacag ctgccccaga gtcctctaca tccttcccgc ctgcccccac ctcaggcatgtcccatcccc cacctgccgc cagagaggac aagagcccat cagaggaatc cgcccccacgacgtccccag agtctgtgtc aggctctgtt ccctcttcag gtagcagcgg gcgagaggaagacgcggcct ccacgctagt gacgggctct gagtatgaga cgatgctgac ggagatcatgtccatgggct atgagcgaga gcgggtcgtg gccgccctga gagccagcta caacaacccccaccgagccg tggagtatct gctcacggga attcctggga gccccgagcc ggaacacggttctgtccagg agagccaggt atcggagcag ccggccacgg aagcagcagg agagaaccccctggagttcc tgcgggacca gccccagttc cagaacatgc ggcaggtgat tcagcagaaccctgcgctgc tgcccgccct gctccagcag ctgggccagg agaaccctca gcttttacagcaaatcagcc ggcaccagga gcagttcatc cagatgctga acgagccccc tggggagctggcggacatct cagatgtgga gggggaggtg ggcgccatag gagaggaggc cccgcagatgaactacatcc aggtgacgcc gcaggagaaa gaagctatag agaggttgaa ggccctgggcttcccagaga gcctggtcat ccaggcctat ttcgcgtgtg aaaaaaatga gaacttggctgccaacttcc tcctgagtca gaactttgat gacgagtgat gccaggaagc caggccaccgaagcccccac cctaccctta ttccatgaaa gttttataaa agaaaaaata tatatatattcatgtttatt taagaaatgg aaaaaaaaat caaaaatctt aaaaaaacaa gcaaacagtccagcttcctg tcctcctaaa gtggcccctg ttcccatctc ccgggccaga cagctgtccccccgtcctcc tccccagccc agcctgctca gagaagctgg caggactggg aggcgacagatgggcccctc ttggcctctg tcccagctct ctgcagccag acggaaaggc ggctgcttgcctctccatcc tccgaaaaac ccctgaggac ccccccccat cctcttctag gatgaggggaagctggagcc ccaactttga tcctccattg gagtggccca aatctttcca tctagggcaagtcctgaaag gcccaaggcc ccctccccag tctggccttg gcctccagcc tggagaagggctaacatcag ctcattgtca aggccacccc caccccagaa cagaaccgtg tctctgataaaggttttgaa gtgaataaag ttttaaaaac ta"
lowerSequence = rawSequence.lower()
noSpace = lowerSequence.replace(" ", "")

# defines counters used to loop chunks
x = 0
y = 10
z = 0
AT = 0
GC = 0

# starts at x0 - y10 in noSpace, and adds 1 to x and y every loop, breaking it into chunks of 10.
for i in range(len(noSpace)):
    # defines variable "primer" length (by default: 10)
    primer = noSpace[x:y]

    for j in range(len(primer)):
        if len(primer) < 10:
            break
        else:
            if primer[z] == "a" or "t":
                AT = AT + 1
            elif primer[z] == "g" or "c":
                GC = GC + 1

            print(primer[z])

            if z == 9:
                z = 0
            else:
                z = z + 1

    print("Primer sequence\n" + primer + "\n\n AT count\n" + str(AT) + "\n\n GC Count\n" + str(GC) +"\n --------------------------")
    x = x + 1
    y = y + 1


#for i in range(len(primer)):
    #count += 1
#print(count/len(primer)*100)

#percentage = 100*(sum([1 for [x,y] in zip(data1,data2) if x == y])/len(data1))
#print('the two sequences are',round(percentage,2),' % the same')
