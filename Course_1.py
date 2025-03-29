import random
import copy

def pattern_matching(sequence, pattern):
    # returns all indices where a pattern is found in a sequence
    main_sequence = sequence
    search_pattern = pattern
    sequence_length = len(main_sequence)
    pattern_length = len(search_pattern)
    match_indices = []
    for i in range(0, sequence_length - pattern_length + 1):
        is_match = 1
        for j in range(i, i + pattern_length):
            if main_sequence[j] != search_pattern[j - i]:
                is_match = 0
                break
        if is_match: match_indices.append(i)
    return match_indices

def reverse_complement(dna_sequence):
    # returns reverse complement of a sequence
    complement_map = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }
    dna_sequence = dna_sequence.upper()
    reversed_complement = ''
    for i in range(len(dna_sequence) - 1, -1, -1):
        reversed_complement += complement_map[dna_sequence[i]]
    return reversed_complement

def frequency_table(sequence, k):
    # returns frequencies of all k-mers in a sequence given sequence and k
    dna_sequence = sequence
    sequence_length = len(dna_sequence)
    kmer_counts = {}
    for i in range(0, sequence_length - k + 1):
        current_kmer = dna_sequence[i:i + k]
        if current_kmer not in kmer_counts: kmer_counts[current_kmer] = 0
        kmer_counts[current_kmer] += 1
    return kmer_counts

def find_clump_pattern(sequence, k, window_size, min_occurrences):
    # returns all k-mers occurring at least min_occurrences times in a window
    dna_sequence = sequence
    sequence_length = len(dna_sequence)
    clump_kmers = set()
    for i in range(0, sequence_length - window_size + 1):
        window = dna_sequence[i:i + window_size]
        window_counts = frequency_table(window, k)
        for kmer, count in window_counts.items():
            if count >= min_occurrences: clump_kmers.add(kmer)
    return list(clump_kmers)

def skew(dna_sequence):
    # returns skew list for a sequence
    sequence = dna_sequence
    current_skew = 0
    skew_values = []
    for i in range(len(sequence)):
        if sequence[i] == 'G': current_skew += 1
        elif sequence[i] == 'C': current_skew -= 1
        skew_values.append(current_skew)
    return skew_values

def hamming_distance(sequence1, sequence2):
    distance = 0
    for i in range(len(sequence1)):
        if sequence1[i] != sequence2[i]: distance += 1
    return distance

def approximate_pattern_matching(sequence, pattern, max_distance):
    # returns indices of k-mers with hamming distance to pattern at most max_distance
    main_sequence = sequence
    search_pattern = pattern
    sequence_length = len(main_sequence)
    pattern_length = len(search_pattern)
    match_indices = []
    for i in range(sequence_length - pattern_length + 1):
        current_kmer = main_sequence[i:i + pattern_length]
        if hamming_distance(current_kmer, search_pattern) <= max_distance: 
            match_indices.append(i)
    return match_indices

def neighbors(pattern, max_distance):
    # returns all strings with hamming distance to pattern at most max_distance
    original_pattern = pattern
    pattern_length = len(original_pattern)
    neighbor_list = []
    if pattern_length < 1: return neighbor_list
    elif max_distance == 0: return [pattern,]
    elif pattern_length == 1: return ['A','C','G','T']
    prefix = pattern[:pattern_length - 1]
    prefix_neighbors = neighbors(prefix, max_distance)
    bases = ('A','C','G','T')
    for neighbor in prefix_neighbors:
        if hamming_distance(neighbor, original_pattern) == max_distance:
            neighbor_list.append(neighbor + original_pattern[pattern_length - 1])
        else:
            for base in bases:
                neighbor_list.append(neighbor + base)
    return neighbor_list

def better_frequent_kmers(sequence, k, max_distance):
    # returns max frequency kmers with hamming distance to sequence at most max_distance
    dna_sequence = sequence
    sequence_length = len(dna_sequence)
    kmer_frequencies = {}
    for i in range(sequence_length - k + 1):
        current_kmer = dna_sequence[i:i + k]
        kmer_neighbors = neighbors(current_kmer, max_distance)
        for neighbor in kmer_neighbors:
            if neighbor not in kmer_frequencies: kmer_frequencies[neighbor] = 0
            kmer_frequencies[neighbor] += 1
    max_frequency = max(kmer_frequencies.values())
    top_kmers = []
    for kmer, frequency in kmer_frequencies.items():
        if frequency == max_frequency: top_kmers.append(kmer)
    return top_kmers

def better_frequent_kmers_with_rc(sequence, k, max_distance):
    # returns max frequency kmers including reverse complements
    dna_sequence = sequence
    sequence_length = len(dna_sequence)
    kmer_frequencies = {}
    for i in range(sequence_length - k + 1):
        current_kmer = dna_sequence[i:i + k]
        reverse_kmer = reverse_complement(current_kmer)
        kmer_neighbors = neighbors(current_kmer, max_distance)
        reverse_neighbors = neighbors(reverse_kmer, max_distance)
        for neighbor in kmer_neighbors:
            if neighbor not in kmer_frequencies: kmer_frequencies[neighbor] = 0
            kmer_frequencies[neighbor] += 1
        for neighbor in reverse_neighbors:
            if neighbor not in kmer_frequencies: kmer_frequencies[neighbor] = 0
            kmer_frequencies[neighbor] += 1
    max_frequency = max(kmer_frequencies.values())
    top_kmers = []
    for kmer, frequency in kmer_frequencies.items():
        if frequency == max_frequency: top_kmers.append(kmer)
    return top_kmers

def motif_enumeration(sequence_list, k, max_distance):
    # returns all kmers with hamming distance at most max_distance to each sequence
    initial_patterns = set()
    for sequence in sequence_list:
        for i in range(len(sequence) - k + 1):
            initial_patterns.add(sequence[i:i + k])
    expanded_patterns = set()
    for pattern in initial_patterns:
        for variant in neighbors(pattern, max_distance): 
            expanded_patterns.add(variant)
    valid_motifs = set()
    for candidate in expanded_patterns:
        is_valid = True
        for sequence in sequence_list:
            if len(approximate_pattern_matching(sequence, candidate, max_distance)) == 0:
                is_valid = False
                break
        if is_valid: valid_motifs.add(candidate)
    return list(valid_motifs)

def pattern_and_list_distance(sequence_list, pattern):
    # returns sum of hamming distances from pattern to all sequences
    search_pattern = pattern
    total_distance = 0
    pattern_length = len(search_pattern)
    for sequence in sequence_list:
        sequence_length = len(sequence)
        min_distance = sequence_length
        for i in range(sequence_length - pattern_length + 1):
            current_distance = hamming_distance(search_pattern, sequence[i:i + pattern_length])
            if current_distance < min_distance: 
                min_distance = current_distance
        total_distance += min_distance
    return total_distance

def nucleotide_bases_permutations(length):
    # returns all possible dna sequences of given length
    bases = ['A','C','G','T']
    current_strings = ['',]
    next_strings = []
    for i in range(length):
        for s in current_strings:
            for base in bases:
                next_strings.append(s + base)
        current_strings = next_strings.copy()
        next_strings = []
    return current_strings

def median_string(sequence_list, k):
    # returns kmers with least sum of hamming distances to all sequences
    all_possible_kmers = nucleotide_bases_permutations(k)
    best_kmers = []
    min_total_distance = len(sequence_list) * k
    distance_map = {}
    for kmer in all_possible_kmers:
        current_distance = pattern_and_list_distance(sequence_list, kmer)
        distance_map[kmer] = current_distance
    min_distance = min(distance_map.values())
    for kmer, distance in distance_map.items():
        if distance == min_distance: best_kmers.append(kmer)
    return best_kmers

def profile_without_pseudocounts(motifs):
    motif_collection = motifs
    num_motifs = len(motif_collection)
    motif_length = len(motif_collection[0])
    profile = {'A':[0]*motif_length, 'C':[0]*motif_length, 
               'G':[0]*motif_length, 'T':[0]*motif_length}
    for pos in range(motif_length):
        for motif in motif_collection:
            profile[motif[pos]][pos] += 1
    for base in 'ACGT':
        for pos in range(motif_length):
            profile[base][pos] /= num_motifs
    profile_matrix = [profile[base] for base in 'ACGT']
    return profile_matrix

def profile_with_pseudocounts(motifs):
    motif_collection = motifs
    num_motifs = len(motif_collection)
    motif_length = len(motif_collection[0])
    profile = {'A':[1]*motif_length, 'C':[1]*motif_length,
               'G':[1]*motif_length, 'T':[1]*motif_length}
    for pos in range(motif_length):
        for motif in motif_collection:
            profile[motif[pos]][pos] += 1
    for base in 'ACGT':
        for pos in range(motif_length):
            profile[base][pos] /= (num_motifs + 4)
    profile_matrix = [profile[base] for base in 'ACGT']
    return profile_matrix

def profile_most_probable_kmer(profile, sequence, k):
    dna_sequence = sequence
    sequence_length = len(dna_sequence)
    max_probability = 0
    best_kmer = dna_sequence[:k]
    base_index = {'A':0, 'C':1, 'G':2, 'T':3}
    for i in range(sequence_length - k + 1):
        current_kmer = dna_sequence[i:i + k]
        probability = 1
        for pos in range(k):
            probability *= profile[base_index[current_kmer[pos]]][pos]
        if probability > max_probability:
            max_probability = probability
            best_kmer = current_kmer
    return best_kmer

def profile_kmer_probabilities(profile, sequence, k):
    dna_sequence = sequence
    sequence_length = len(dna_sequence)
    probabilities = []
    base_index = {'A':0, 'C':1, 'G':2, 'T':3}
    for i in range(sequence_length - k + 1):
        current_kmer = dna_sequence[i:i + k]
        probability = 1
        for pos in range(k):
            probability *= profile[base_index[current_kmer[pos]]][pos]
        probabilities.append(probability)
    return probabilities

def consensus(motifs):
    # returns consensus for a given list of motifs
    num_motifs = len(motifs)
    motif_length = len(motifs[0])
    consensus_sequence = ''
    for pos in range(motif_length):
        a_count = t_count = g_count = c_count = 0
        for motif in motifs:
            if motif[pos] == 'A': a_count += 1
            elif motif[pos] == 'T': t_count += 1
            elif motif[pos] == 'G': g_count += 1
            else: c_count += 1
        max_count = max(a_count, t_count, g_count, c_count)
        if max_count == a_count: consensus_sequence += 'A'
        elif max_count == c_count: consensus_sequence += 'C'
        elif max_count == g_count: consensus_sequence += 'G'
        else: consensus_sequence += 'T'
    return consensus_sequence

def score_of_motifs(motifs, consensus_seq):
    # returns score for a given list of motifs and consensus
    score = 0
    num_motifs = len(motifs)
    motif_length = len(motifs[0])
    for pos in range(motif_length):
        for motif in motifs:
            if motif[pos] != consensus_seq[pos]: score += 1
    return score

def greedy_motif_search_with_pseudocounts(sequence_list, k):
    num_sequences = len(sequence_list)
    sequences = sequence_list
    best_motifs = []
    for seq in sequences: best_motifs.append(seq[:k])
    best_score = score_of_motifs(best_motifs, consensus(best_motifs))
    for i in range(len(sequence_list[0]) - k + 1):
        current_motifs = []
        counts = {'A':[0]*k, 'C':[0]*k, 'G':[0]*k, 'T':[0]*k}
        initial_kmer = sequences[0][i:i + k]
        for pos in range(k):
            counts[initial_kmer[pos]][pos] += 1
        current_motifs.append(initial_kmer)
        for t in range(1, num_sequences):
            profile = [counts[base] for base in 'ACGT']
            for row in range(len(profile)):
                for col in range(len(profile[0])):
                    profile[row][col] += 1
            next_motif = profile_most_probable_kmer(profile, sequences[t], k)
            current_motifs.append(next_motif)
            for pos in range(k):
                counts[next_motif[pos]][pos] += 1
        current_score = score_of_motifs(current_motifs, consensus(current_motifs))
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs.copy()
    return best_motifs

def randomized_motif_search(sequence_list, k):
    # returns best_motifs and its score
    sequences = sequence_list
    num_sequences = len(sequence_list)
    sequence_length = len(sequences[0])
    best_motifs = []
    for i in range(num_sequences):
        start = random.randint(0, sequence_length - k)
        best_motifs.append(sequences[i][start:start + k])
    best_profile = profile_with_pseudocounts(best_motifs)
    while True:
        current_motifs = []
        for i in range(num_sequences):
            current_motifs.append(profile_most_probable_kmer(best_profile, sequences[i], k))
        current_profile = profile_with_pseudocounts(current_motifs)
        best_motif_score = score_of_motifs(best_motifs, consensus(best_motifs))
        current_score = score_of_motifs(current_motifs, consensus(current_motifs))
        if current_score < best_motif_score:
            best_motifs = current_motifs.copy()
            best_profile = copy.deepcopy(current_profile)
        else:
            return best_motifs, best_motif_score

def gibbs_sampler(sequence_list, k, iterations):
    sequences = sequence_list
    num_sequences = len(sequence_list)
    sequence_length = len(sequences[0])
    best_motifs = []
    for i in range(num_sequences):
        start = random.randint(0, sequence_length - k)
        best_motifs.append(sequences[i][start:start + k])
    best_score = score_of_motifs(best_motifs, consensus(best_motifs))
    current_motifs = best_motifs.copy()
    for _ in range(iterations):
        excluded_index = random.randint(0, num_sequences - 1)
        del current_motifs[excluded_index]
        profile = profile_with_pseudocounts(current_motifs)
        positions = []
        probabilities = profile_kmer_probabilities(profile, sequences[excluded_index], k)
        chosen_position = random.choices(range(len(probabilities)), weights=probabilities, k=1)[0]
        current_motifs.insert(excluded_index, sequences[excluded_index][chosen_position:chosen_position + k])
        current_score = score_of_motifs(current_motifs, consensus(current_motifs))
        if current_score < best_score:
            best_motifs = current_motifs.copy()
            best_score = current_score
    return best_motifs, best_score
