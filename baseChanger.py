def reverseComplemented(sequence):
    reverseComplement = ""
    sequence_upper = sequence[::-1].upper()
    for nucleotide in sequence_upper:
        if nucleotide == "A":
            reverseComplement += "T"
        elif nucleotide == "T":
            reverseComplement += "A"
        elif nucleotide == "C":
            reverseComplement += "G"
        elif nucleotide == "G":
            reverseComplement += "C"
    return reverseComplement

def SaPAMFinder(sequence):
    #SaCas9 PAM = NNGRRT
    R = ["A","G"]
    fragment_indices = []
    reverse_fragment_indices = []
    sequence_upper = sequence.upper()
    reverse_sequence_upper = reverseComplemented(sequence).upper()
    for increment in list(range((len(sequence_upper)-5))):
        PAM_indices = [increment, increment + 6]
        fragment = sequence_upper[increment:increment + 6]
        if fragment[2] == "G" and fragment[3] in R and fragment[4] in R and fragment[5] == "T":
            fragment_indices.append(PAM_indices)

    for increment in list(range((len(reverse_sequence_upper)-5))):
        PAM_indices = [increment, increment + 6]
        fragment = reverse_sequence_upper[increment:increment + 6]
        if fragment[2] == "G" and fragment[3] in R and fragment[4] in R and fragment[5] == "T":
            reverse_fragment_indices.append(PAM_indices)
    return fragment_indices, reverse_fragment_indices

def SpPAMFinder(sequence):
    #SpCas9 PAM = NGG
    fragment_indices = []
    sequence_upper = sequence.upper()
    reverse_sequence_upper = reverseComplemented(sequence).upper()
    for increment in list(range((len(sequence_upper)-2))):
        PAM_indices = [increment, increment + 3]
        fragment = sequence_upper[increment:increment + 3]
        if fragment[2] == "G" and fragment[3] == "G":
            fragment_indices.append(PAM_indices)
    for increment in list(range((len(reverse_sequence_upper)-2))):
        PAM_indices = [increment, increment + 3]
        fragment = reverse_sequence_upper[increment:increment + 3]
        if fragment[2] == "G" and fragment[3] == "G":
            reverse_fragment_indices.append(PAM_indices)
    return fragment_indices, reverse_fragment_indices

def generalPAMFinder(PAM, sequence):
    R = ["A","G"]
    Y = ["C","T"]
    N = ["A","T","C","G"]
    fragment_indices = []
    reverse_fragment_indices = []
    PAM_upper = PAM.upper()
    sequence_upper = sequence.upper()
    reverse_sequence_upper = reverseComplemented(sequence).upper()

    for increment in list(range((len(sequence_upper)-(len(PAM_upper)-1)))):
        PAM_indices = [increment, increment + len(PAM_upper)]
        fragment = sequence_upper[increment:increment + len(PAM_upper)]

        n = 0
        onOff = 0
        for nucleotide in PAM_upper:
            if nucleotide == "N":
                n += 1
            elif nucleotide == "G" and fragment[n] == "G":
                n += 1
            elif nucleotide == "R" and fragment[n] in R:
                n += 1
            elif nucleotide == "T" and fragment[n] == "T":
                n += 1
            elif nucleotide == "Y" and fragment[n] in Y:
                n += 1
            elif nucleotide == "C" and fragment[n] == "C":
                n += 1
            else:
                n += 1
                onOff = 1
        if onOff != 1:
            fragment_indices.append(PAM_indices)
        n += 1

    for increment in list(range((len(reverse_sequence_upper)-(len(PAM_upper)-1)))):
        PAM_indices = [increment, increment + len(PAM_upper)]
        fragment = reverse_sequence_upper[increment:increment + len(PAM_upper)]

        n = 0
        onOff = 0
        for nucleotide in PAM_upper:
            if nucleotide == "N":
                n += 1
            elif nucleotide == "G" and fragment[n] == "G":
                n += 1
            elif nucleotide == "R" and fragment[n] in R:
                n += 1
            elif nucleotide == "T" and fragment[n] == "T":
                n += 1
            elif nucleotide == "Y" and fragment[n] in Y:
                n += 1
            elif nucleotide == "C" and fragment[n] == "C":
                n += 1
            else:
                n += 1
                onOff = 1
        if onOff != 1:
            reverse_fragment_indices.append(PAM_indices)
        n += 1

    return fragment_indices, reverse_fragment_indices

def guideIdentifier(fragment_indices, reverse_fragment_indices, sequence, length_of_guide):
    sequence_upper = sequence.upper()
    reverse_sequence_upper = reverseComplemented(sequence).upper()
    guides = []
    for fragment in fragment_indices:
        if (fragment[0] - length_of_guide) >= 0:
            start_of_guide = fragment[0] - length_of_guide
            end_of_guide = fragment[0]
            guide = sequence_upper[start_of_guide:end_of_guide]
            guides.append(guide)

    guides_reverse = []
    for fragment in reverse_fragment_indices:
        if (fragment[0] - length_of_guide) >= 0:
            start_of_guide = fragment[0] - length_of_guide
            end_of_guide = fragment[0]
            guide = reverse_sequence_upper[start_of_guide:end_of_guide]
            guides_reverse.append(guide)

    if guides != [] and guides_reverse != []:
        return guides, guides_reverse
    elif guides != [] and guides_reverse == []:
        return guides
    elif guides == [] and guides_reverse != []:
        return guides_reverse
    else:
        return "No guides found"

print(reverseComplemented("GTGGGAAAGCAGAGTCGGGGGGTACCGGCGAGTCCGACGAATCGGC"))
print(SaPAMFinder("GTGGGAAAGCAGAGTCGGGGGGTACCGGCGAGTCCGACGAATCGGC"))
print(*generalPAMFinder("NNGRRT","GTGGGAAAGCAGAGTCGGGGGGTACCGGCGAGTCCGACGAATCGGC"))
#print(guideIdentifier(*generalPAMFinder("NNGRRT","GTGGGAAAGCAGAGTCGGGGGGTACCCAAACCGGCGAGTCCGACGAATCGGC"),"GTGGGAAAGCAGAGTCGGGGGGTACCCAAACCGGCGAGTCCGACGAATCGGC",20))

print(guideIdentifier(*generalPAMFinder("NNGRRT","GTGGGAAAGCAGAGTCGGGGGGTACCGGCGAGTCCGACGAATCGGC"),"GTGGGAAAGCAGAGTCGGGGGGTACCGGCGAGTCCGACGAATCGGC",20))
