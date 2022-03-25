"""
Library for finding potential genes in a strand of DNA.
"""
import helpers


def get_complement(nucleotide):
    """
    A function that returns a complementary nucleotide
    correlated to the arguement.

    Args:
        nucleotide: a string representing a nucleotide

    Returns:
        A string representing the inputted nucleotide's
        complement
    """
    if nucleotide == "A":
        return "T"
    if nucleotide == "T":
        return "A"
    if nucleotide == "G":
        return "C"
    if nucleotide == "C":
        return "G"
    return None


def get_reverse_complement(strand):
    """
    A function that returns the complement of a strand
    of DNA.

    Args:
        strand: a string representing a strand of DNA.

    Returns:
        A string representing the complement of the
        inputted strand of DNA.
    """
    comp_strand = []
    for character in strand:
        comp_strand.append(get_complement(character))
    comp_strand.reverse()
    return "".join(comp_strand)


def rest_of_orf(strand):
    """
    A function that scans a strand of DNA beginning with
    a start codon for a stop codon and returns the ORF.

    Args:
        strand: a string representing the strand of DNA
        to be analyzed.

    Returns:
        A string representing the found ORF (from start
        codon up until end codon) or the entire strand
        if not end codon is detected.
    """
    for i, _ in enumerate(strand):
        # Enumerate through the strand.
        if (i + 1) % 3 == 0:
            if strand[(i - 2):(i + 1)] == "TAA" or \
                    strand[(i - 2):(i + 1)] == "TAG" or \
                    strand[(i - 2):(i + 1)] == "TGA":
                return f"{strand[0:(i - 2)]}"
    # If stop codon is detected, we stop enumerating
    # and return what we found.
    return strand
    # If stop codon is never detected, we return whole
    # string.


def find_all_orfs_one_frame(strand):
    """
    A function that scans a strand of DNA for ORFs and
    returns a list of all ORFs found.

    Args:
       strand: a string representing the strand of DNA
       to be analyzed.

    Returns:
        A list representing all ORFs found in the analyzed
        strand of DNA.
    """
    orf_list = []
    i = 0
    while i <= len(strand):
        # While loop with "counter" variable to analyze the
        # string.
        i += 1
        if (i + 1) % 3 == 0 and strand[(i - 2):(i + 1)] == "ATG":
            # Search for valid start codon.
            end = i + 1
            while strand[end:(end + 3)] != "TAA" and \
                    strand[end:(end + 3)] != "TGA" and \
                    strand[end:(end + 3)] != "TAG":
                # Search for valid stop codon.
                end += 3
                if end + 3 > len(strand):
                    # Make sure we don't go out-of-bounds to prevent
                    # infinite while loop.
                    if strand[((len(strand)) - 3):(len(strand))] == "TAA" or \
                            strand[((len(strand)) - 3):(len(strand))] == "TGA" or \
                            strand[((len(strand)) - 3):(len(strand))] == "TAG":
                        orf_list.append(
                            strand[(i - 2):((len(strand)) - 3)])
                        return orf_list
    # If strand ends in end codon, add last ORF and
    # return ORF list.
                    orf_list.append(strand[(i - 2):(len(strand))])
                    return orf_list
    # Otherwise, add rest of strand and return ORF
    # list.
            end += 3
            orf_list.append(strand[(i - 2):(end - 3)])
    # Append the valid ORF to our list.
            i = end - 1
    return orf_list
    # Return completed list at the end of the while loop.


def find_all_orfs(strand):
    """
    A function that finds ORFs in a strand of DNA, as well
    as ORFs in the same strand shifted by one and two
    nucleotides.

    Args:
       strand: a string representing the strand of DNA
       to be analyzed.

    Returns:
        A list containing all ORFs found in the strand of
        DNA as well as all ORFs found in relevant shifts.
    """
    shift_1 = strand[len(strand) - 1:] + \
        strand[:len(strand) - 1]
    # Create a once shifted strand.
    shift_2 = strand[len(strand) - 2:] + \
        strand[:len(strand) - 2]
    # Create a twice shifted strand.
    list_noshift = find_all_orfs_one_frame(strand)
    list_oneshift = find_all_orfs_one_frame(shift_1)
    list_twoshift = find_all_orfs_one_frame(shift_2)
    all_orf_shifts = list_noshift + list_oneshift + list_twoshift
    # Compile and append all ORFs to one list
    return all_orf_shifts


def find_all_orfs_both_strands(strand):
    """
    A function that finds ORFs in a strand of DNA as well as
    its complementary strand.

    Args:
       strand: a string representing the strand of DNA
       to be analyzed.

    Returns:
        A list containing all ORFs found in the analyzed
        strand of DNA and its complement.
    """
    strand_comp = get_reverse_complement(strand)
    left_orfs = find_all_orfs(strand_comp)
    right_orfs = find_all_orfs(strand)
    all_orfs = left_orfs + right_orfs
    return all_orfs


def find_longest_orf(strand):
    """
    A function that finds the longest ORF in a strand of
    DNA.

    Args:
       strand: a string representing the strand of DNA
       to be analyzed.

    Returns:
        A string representing the longest ORF found.
    """
    all_orfs_to_sort = find_all_orfs_both_strands(strand)
    sorted_orfs = sorted(all_orfs_to_sort, key=len)
    if len(sorted_orfs) > 0:
        return sorted_orfs[-1]
    return None


def noncoding_orf_threshold(strand, num_trials):
    """
    A function that randomizes a strand of DNA and
    finds the shortest ORF among the longest ORFs
    found in the randomization process.

    Args:
       strand: a string representing the strand of DNA
       to be analyzed.
       num_trials: a positive integer representing the
       number of randomization trials to hold.

    Returns:
        min_orf_len: an integer representing the length
        of the shortest ORF from all the trials maximum
        ORFs.
    """
    threshold_orfs = []
    for _ in range(num_trials):
        current_strand = helpers.shuffle(strand)
        current_long_orf = find_longest_orf(current_strand)
        if isinstance(current_long_orf, str) is True:
            # If statement makes sure list doesn't get corrupted
            # by adding None type.
            threshold_orfs.append(current_long_orf)
        else:
            continue
    # Hold the randomization trials.
    longest_orf = min(threshold_orfs, key=len)
    # Find the shortest threshold ORF.
    min_orf_len = int(len(longest_orf))
    # Find the length of the minimum ORF as int.
    return min_orf_len


def encode_amino_acids(orf):
    """
    A function that converts ORFs into the animo acids
    they encode.

    Args:
       orf: a string representing the orf to be analyzed.

    Returns:
        A string representing the amino acids the
        analyzed ORF encodes.
    """
    amino_acids = []
    i = 0
    if len(orf) % 3 != 0:
        orf = orf[:-1]
    if len(orf) % 3 != 0:
        orf = orf[:-1]
    # Make length of orf divisable by three. We repeat
    # this process twice to ensure this occurs.
    while i <= len(orf):
        i += 1
        if (i + 1) % 3 == 0:
            current_acid = helpers.amino_acid(orf[(i - 2):(i + 1)])
            amino_acids.append(current_acid)
            continue
    # Turn ORF into relevant animo acids.
    return "".join(amino_acids)


def find_genes(path):
    """
    A function to load a FASTA file of nucleotide and
    find all amino acids of significance it can encode.

    Args:
        path: a string representing the file path where
        the FASTA file to be analyzed is located.

    Returns:
        all_acids: a list with all of the relevant amino
        acids that can be encoded by the FASTA file.
    """
    strand = helpers.load_fasta_file(path)
    min_length = int(noncoding_orf_threshold(strand, 1500))
    print(min_length)
    all_orfs = find_all_orfs_both_strands(strand)
    all_acids_temp = []
    all_acids = []
    for item in all_orfs:
        acid_curr = encode_amino_acids(item)
        all_acids_temp.append(acid_curr)
# Encode amino acids and add them to our list
    for item in all_acids_temp:
        if len(item) < (min_length // 3):
            del item
        else:
            all_acids.append(item)
# Remove amino-acids from ORFs smaller than minimum
# size threshold.
    return all_acids
