# Processes input. First makes it lowercase, then removes the spaces to get a pure sequence.
rawSequence = input("Please enter your sequence: ")
lowerSequence = rawSequence.lower()
noSpace = rawSequence.replace(" ", "")

nucleotideCounter = 0
primer_min_len = 28     # Defines minimum primer length
primer_max_len = 30     # Defines maximum primer length
seq_len = len(noSpace)  # Stores length of the processed sequence as an integer.


def make_reverse(toReverse):               # Function to fetch the reverse compliment of the reverse primer.
    normal = "acgt"
    compliment = "tgca"
    reversePrimer = str.maketrans(normal, compliment)
    # print("toReverse: " + toReverse)
    # print("Original primer_dict[x]: " + primer_dict[x]["Sequence:"])
    tmp = toReverse.translate(reversePrimer)[::-1]
    return tmp
    # print("Revcomp primer_dict[x]: " + primer_dict[x]["Sequence:"])
    # print("-----------------------------------------")
    # return



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
min_distance = 100
pairs = []
cnt = 0

f = open("output.txt", "w")                         # Opens output file
i = 0
for primer1 in primer_dict:
    for primer2 in primer_dict[i:]:                 # Ensures it doesn't check for reverse primers behind the fwd.
        if abs( primer1["Tm"] - primer2["Tm"]) <= 5 and (primer1["Sequence:"].endswith("g") or primer1["Sequence:"].endswith("c")):
            if primer2["Begin"] - primer1["End"] >= min_distance:
                if primer2["Sequence:"].startswith("g") or primer2["Sequence:"].startswith("c"):
                    temp = primer2
                    temp["Sequence:"] = make_reverse(temp["Sequence:"])
                    f.write("Pair #: " + str(cnt) + "\n")
                    f.write("Forward Primer: " + str(primer1) + "\n")
                    f.write("Reverse Primer: " + str(temp) + "\n")
                    f.write("Product Length: " + str(primer2["End"] - primer1["Begin"]) + "\n")
                    f.write("\n")
                    pairs.append([primer1, primer2])
                    cnt = cnt + 1                       # Counts valid pairs
    i = i+1

f.close()                                           # Closes output file

num = 0

print("Potential Primers: " + str(len(primer_dict)))
print("Valid Pairs:" + str(cnt))
