import networkx as nx

codon_to_aa = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S',
    'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I',
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GGA': 'G', 'GGC': 'G', 'GGA': 'G', 'GGU': 'G',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'UAA': 'Stop', 'UAC': 'Y', 'UAG': 'Stop', 'UAU': 'Y',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UGA': 'Stop', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C',
    'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F'
}

aa_to_codon = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAU', 'AAC'],
    'D': ['GAU', 'GAC'],
    'C': ['UGU', 'UGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'K': ['AAA', 'AAG'],
    'M': ['AUG'],
    'F': ['UUU', 'UUC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'Stop': ['UAA', 'UAG', 'UGA']
}

amino_acid_mass = {
    'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 
    'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 
    'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 
    'F': 147, 'R': 156, 'Y': 163, 'W': 186
}

amino_weights = {57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186}

mass_to_aa = {
    57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 
    103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'K', 
    129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 
    163: 'Y', 186: 'W'
}

def output_file(output_list):
    with open("output.txt", "a") as file:
        for item in output_list:
            file.write(f"{item} ")

def input_file(file_path):
    with open(file_path, "r") as file:
        content = file.read()
    return content.split()

def sequence_composition(sequence, k):
    kmer_list = []
    for i in range(len(sequence) - k + 1):
        kmer_list.append(sequence[i:i+k])
    return kmer_list

def ordered_kmers_to_sequence(kmer_list):
    assembled_sequence = kmer_list[0]
    k = len(kmer_list[0])
    for i in range(1, len(kmer_list)):
        assembled_sequence += kmer_list[i][k-1]
    return assembled_sequence

def overlap_graph(kmer_list):
    prefix_map = {}
    k = len(kmer_list[0])
    for kmer in kmer_list:
        prefix = kmer[:k-1]
        if prefix not in prefix_map:
            prefix_map[prefix] = []
        prefix_map[prefix].append(kmer)
    
    graph_edges = {}
    for kmer in kmer_list:
        suffix = kmer[1:]
        if suffix not in prefix_map:
            continue
        if suffix == kmer[:k-1]:
            prefix_map[suffix].remove(kmer)
        graph_edges[kmer] = prefix_map[suffix]
    return graph_edges

def debruijn_graph(kmer_list):
    adjacency_map = {}
    k = len(kmer_list[0])
    for kmer in kmer_list:
        prefix = kmer[:k-1]
        if prefix not in adjacency_map:
            adjacency_map[prefix] = []
        adjacency_map[prefix].append(kmer[1:])
    
    edge_list = []
    for prefix, suffixes in adjacency_map.items():
        for suffix in suffixes:
            edge_list.append((prefix, suffix))
    return edge_list

def paired_debruijn_graph(paired_kmer_list):
    adjacency_map = {}
    k = len(paired_kmer_list[0][0])
    for pair in paired_kmer_list:
        prefix_a = pair[0][:k-1]
        prefix_b = pair[1][:k-1]
        if (prefix_a, prefix_b) not in adjacency_map:
            adjacency_map[(prefix_a, prefix_b)] = []
        adjacency_map[(prefix_a, prefix_b)].append((pair[0][1:], pair[1][1:]))
    
    edge_list = []
    for prefix_pair, suffix_pairs in adjacency_map.items():
        for suffix_pair in suffix_pairs:
            edge_list.append((prefix_pair, suffix_pair))
    return edge_list

def eulerian_path(edge_list):
    graph = nx.DiGraph()
    graph.add_edges_from(edge_list)
    return list(nx.eulerian_path(graph))

def path_to_sequence(path):
    sequence = path[0][0]
    k = len(path[0][0])
    for i in range(1, len(path)):
        sequence += path[i][0][k-1]
    sequence += path[-1][1][k-1]
    return sequence

def paired_path_to_sequence(path, d):
    n = len(path)
    k = len(path[0][0][0])
    size = 2*(k+1) + n + d - 1
    forward = [-1] * size
    backward = [-1] * size
    
    forward_ptr = k
    backward_ptr = size - k - 1
    
    for i in range(k):
        forward[i] = path[0][0][0][i]
    for i in range(k-1, -1, -1):
        backward[size + i - k] = path[n-1][1][1][i]
    
    for i in range(1, n):
        forward[forward_ptr] = path[i][0][0][k-1]
        forward_ptr += 1
    forward[forward_ptr] = path[n-1][1][0][k-1]
    
    for i in range(n-2, -1, -1):
        backward[backward_ptr] = path[i][1][1][0]
        backward_ptr -= 1
    backward[backward_ptr] = path[0][0][1][0]
    
    result = ''
    for i in range(size):
        if forward[i] == -1 and backward[i] == -1:
            print("both -1")
        elif forward[i] == -1:
            result += backward[i]
        elif backward[i] == -1:
            result += forward[i]
        elif forward[i] != backward[i]:
            print("not -1 but not equal")
        else:
            result += forward[i]
    return result

def contigs_from_kmers(kmer_list):
    k = len(kmer_list[0])
    contigs = []
    adjacency_map = {}
    
    for kmer in kmer_list:
        prefix = kmer[:k-1]
        if prefix not in adjacency_map:
            adjacency_map[prefix] = []
        adjacency_map[prefix].append(kmer[1:])
    
    in_degree = {}
    out_degree = {}
    for prefix, suffixes in adjacency_map.items():
        out_degree[prefix] = len(suffixes)
        for suffix in suffixes:
            if suffix not in in_degree:
                in_degree[suffix] = 0
            in_degree[suffix] += 1
    
    for node in adjacency_map:
        if in_degree.get(node, 0) == out_degree.get(node, 0) == 1:
            continue
        
        current_node = node
        while adjacency_map.get(current_node, []):
            contig = current_node
            next_node = adjacency_map[current_node].pop()
            
            while in_degree.get(next_node, 0) == out_degree.get(next_node, 0) == 1:
                contig += next_node[-1]
                next_node = adjacency_map[next_node][0]
            
            contig += next_node[-1]
            contigs.append(contig)
    
    return contigs

def rna_to_peptide(rna_sequence):
    peptide = ''
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        amino_acid = codon_to_aa[codon]
        if amino_acid == 'Stop':
            continue
        peptide += amino_acid
    return peptide

def peptide_to_rna(peptide):
    possible_rna = []
    for amino_acid in peptide:
        codons = aa_to_codon[amino_acid]
        if not possible_rna:
            possible_rna = codons.copy()
        else:
            new_rna = []
            for rna in possible_rna:
                for codon in codons:
                    new_rna.append(rna + codon)
            possible_rna = new_rna
    return possible_rna

def cycle_spectrum_of_peptide(peptide):
    n = len(peptide)
    spectrum = [0]
    prefix_mass = [0] * (n + 1)
    
    for i in range(1, n + 1):
        prefix_mass[i] = prefix_mass[i-1] + amino_acid_mass[peptide[i-1]]
        spectrum.append(prefix_mass[i])
    
    for length in range(1, n):
        for start in range(n):
            if start == 0:
                current_mass = prefix_mass[length]
            else:
                current_mass += amino_acid_mass[peptide[(start + length - 1) % n]]
                current_mass -= amino_acid_mass[peptide[start - 1]]
            spectrum.append(current_mass)
    
    spectrum.sort()
    return spectrum

def linear_spectrum_of_peptide(peptide):
    n = len(peptide)
    spectrum = [0]
    prefix_mass = [0] * (n + 1)
    
    for i in range(1, n + 1):
        prefix_mass[i] = prefix_mass[i-1] + amino_acid_mass[peptide[i-1]]
        spectrum.append(prefix_mass[i])
    
    for length in range(1, n):
        for start in range(n):
            if start + length - 1 >= n:
                break
            if start == 0:
                current_mass = prefix_mass[length]
            else:
                current_mass += amino_acid_mass[peptide[start + length - 1]]
                current_mass -= amino_acid_mass[peptide[start - 1]]
            spectrum.append(current_mass)
    
    spectrum.sort()
    return spectrum

def cycle_spectrum_of_peptide_weights(weights):
    n = len(weights)
    spectrum = [0]
    prefix_mass = [0] * (n + 1)
    
    for i in range(1, n + 1):
        prefix_mass[i] = prefix_mass[i-1] + weights[i-1]
        spectrum.append(prefix_mass[i])
    
    for length in range(1, n):
        for start in range(n):
            if start == 0:
                current_mass = prefix_mass[length]
            else:
                current_mass += weights[(start + length - 1) % n]
                current_mass -= weights[start - 1]
            spectrum.append(current_mass)
    
    spectrum.sort()
    return spectrum

def linear_spectrum_of_peptide_weights(weights):
    n = len(weights)
    spectrum = [0]
    prefix_mass = [0] * (n + 1)
    
    for i in range(1, n + 1):
        prefix_mass[i] = prefix_mass[i-1] + weights[i-1]
        spectrum.append(prefix_mass[i])
    
    for length in range(1, n):
        for start in range(n):
            if start + length - 1 >= n:
                break
            if start == 0:
                current_mass = prefix_mass[length]
            else:
                current_mass += weights[start + length - 1]
                current_mass -= weights[start - 1]
            spectrum.append(current_mass)
    
    spectrum.sort()
    return spectrum

def expand_peptides(peptide_list):
    expanded_peptides = []
    for peptide in peptide_list:
        for mass in amino_weights:
            expanded_peptides.append(peptide + [mass])
    expanded_peptides.sort()
    return expanded_peptides

def cyclopeptide_sequencing(spectrum):
    candidate_peptides = [[]]
    final_peptides = []
    target_spectrum = set(spectrum)
    
    while candidate_peptides:
        candidate_peptides = expand_peptides(candidate_peptides)
        new_candidates = []
        
        for peptide in candidate_peptides:
            peptide_spectrum = set(cycle_spectrum_of_peptide_weights(peptide))
            if peptide_spectrum == target_spectrum:
                if peptide not in final_peptides:
                    final_peptides.append(peptide)
            else:
                linear_spectrum = linear_spectrum_of_peptide_weights(peptide)
                if all(mass in target_spectrum for mass in linear_spectrum):
                    new_candidates.append(peptide)
        
        candidate_peptides = new_candidates
    
    return final_peptides

def peptide_spectrum_score(peptide, spectrum):
    theoretical_spectrum = cycle_spectrum_of_peptide(peptide)
    spectrum_counts = {}
    theoretical_counts = {}
    
    for mass in theoretical_spectrum:
        if mass > max(spectrum):
            return 0
        theoretical_counts[mass] = theoretical_counts.get(mass, 0) + 1
    
    for mass in spectrum:
        spectrum_counts[mass] = spectrum_counts.get(mass, 0) + 1
    
    score = 0
    for mass, count in theoretical_counts.items():
        if mass in spectrum_counts:
            score += min(count, spectrum_counts[mass])
    
    return score

def peptide_spectrum_weights_score(peptide_weights, spectrum):
    theoretical_spectrum = cycle_spectrum_of_peptide_weights(peptide_weights)
    spectrum_counts = {}
    theoretical_counts = {}
    
    for mass in theoretical_spectrum:
        if mass > max(spectrum):
            return 0
        theoretical_counts[mass] = theoretical_counts.get(mass, 0) + 1
    
    for mass in spectrum:
        spectrum_counts[mass] = spectrum_counts.get(mass, 0) + 1
    
    score = 0
    for mass, count in theoretical_counts.items():
        if mass in spectrum_counts:
            score += min(count, spectrum_counts[mass])
    
    return score

def trim_leaderboard(leaderboard, spectrum, N):
    if not leaderboard:
        return leaderboard
    
    scored_peptides = []
    for peptide in leaderboard:
        score = peptide_spectrum_weights_score(peptide, spectrum)
        scored_peptides.append((score, peptide))
    
    scored_peptides.sort(reverse=True)
    
    if len(scored_peptides) <= N:
        return [peptide for (score, peptide) in scored_peptides]
    
    cutoff_score = scored_peptides[N-1][0]
    trimmed_leaderboard = []
    
    for score, peptide in scored_peptides:
        if score >= cutoff_score and score > 0:
            trimmed_leaderboard.append(peptide)
        else:
            break
    
    return trimmed_leaderboard

def leaderboard_cyclopeptide_sequencing(spectrum, N):
    leaderboard = [[]]
    best_peptide = []
    best_score = 0
    target_mass = max(spectrum)
    
    while leaderboard:
        leaderboard = expand_peptides(leaderboard)
        new_leaderboard = []
        
        for peptide in leaderboard:
            peptide_mass = sum(peptide)
            
            if peptide_mass == target_mass:
                current_score = peptide_spectrum_weights_score(peptide, spectrum)
                if current_score > best_score:
                    best_peptide = peptide.copy()
                    best_score = current_score
            elif peptide_mass > target_mass:
                continue
            
            new_leaderboard.append(peptide)
        
        leaderboard = trim_leaderboard(new_leaderboard, spectrum, N)
    
    return best_peptide

def spectral_convolution(spectrum):
    spectrum = sorted(spectrum)
    convolution = []
    for i in range(len(spectrum)):
        for j in range(i):
            mass_diff = spectrum[i] - spectrum[j]
            if mass_diff >= 57 and mass_diff <= 200:
                convolution.append(mass_diff)
    return convolution

def convolution_cyclopeptide_sequencing(spectrum, M, N):
    convolution = spectral_convolution(spectrum)
    mass_counts = {}
    
    for mass in convolution:
        mass_counts[mass] = mass_counts.get(mass, 0) + 1
    
    sorted_masses = sorted(mass_counts.items(), key=lambda x: x[1])
    cutoff_count = sorted_masses[-M][1]
    
    original_weights = amino_weights.copy()
    amino_weights.clear()
    
    for mass, count in mass_counts.items():
        if count >= cutoff_count:
            amino_weights.add(mass)
    
    result = leaderboard_cyclopeptide_sequencing(spectrum, N)
    amino_weights = original_weights
    
    return result
