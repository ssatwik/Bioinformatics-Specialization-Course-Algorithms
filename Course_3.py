import sys
sys.setrecursionlimit(10000)

def output_file(output_list):
    n = len(output_list)
    with open("output.txt", "a") as file:
        for i in range(n):
            file.write(f"{output_list[i]} ")  

def input_file(path):
    with open(path, "r") as file:
        content = file.read()
    input_list = content.split()
    return input_list

def dp_change_problem(target_sum, coins):
    n = len(coins)
    min_coins = [-1] * (target_sum + 1)
    min_coins[0] = 0
    for current_sum in range(target_sum + 1):
        for coin in coins:
            if coin + current_sum <= target_sum:
                if min_coins[current_sum] != -1: 
                    if min_coins[coin + current_sum] > min_coins[current_sum] + 1 or min_coins[coin + current_sum] == -1:
                        min_coins[coin + current_sum] = min_coins[current_sum] + 1
    return min_coins[target_sum]

def manhattan_tourist_problem(grid, dp_table, x, y):
    if dp_table[x][y] != -1: 
        return dp_table[x][y]
    dp_table[x][y] = max(manhattan_tourist_problem(grid, dp_table, x+1, y) + grid[x][y][1],
                         manhattan_tourist_problem(grid, dp_table, x, y+1) + grid[x][y][0])
    return dp_table[x][y]

def longest_common_subsequence(seq1, seq2, idx1, idx2, dp_table):
    if dp_table[idx1][idx2][0] != -1: 
        return dp_table[idx1][idx2]
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    if idx1 >= len_seq1 or idx2 >= len_seq2:
        dp_table[idx1][idx2][0] = 0
        dp_table[idx1][idx2][1] = ''
        return dp_table[idx1][idx2]
    
    match_len, match_str = 0, ''
    if seq1[idx1] == seq2[idx2]:
        match_result = longest_common_subsequence(seq1, seq2, idx1+1, idx2+1, dp_table)
        match_len = match_result[0] + 1
        match_str = seq1[idx1] + match_result[1]
    
    skip_seq1 = longest_common_subsequence(seq1, seq2, idx1+1, idx2, dp_table)
    skip_seq2 = longest_common_subsequence(seq1, seq2, idx1, idx2+1, dp_table)
    
    max_len = max(match_len, skip_seq1[0], skip_seq2[0])
    
    if max_len == match_len:
        dp_table[idx1][idx2] = [match_len, match_str]
    elif max_len == skip_seq1[0]:
        dp_table[idx1][idx2] = [skip_seq1[0], skip_seq1[1]]
    else:
        dp_table[idx1][idx2] = [skip_seq2[0], skip_seq2[1]]
    return dp_table[idx1][idx2]

def longest_path_in_DAG(graph, current_node, end_node, dp_table, visited):
    if dp_table[current_node][0] != -1: 
        return dp_table[current_node]
    visited[current_node] = 1
    num_nodes = len(graph)
    max_length = 0
    path = []
    
    for neighbor in graph[current_node]:
        if visited[neighbor[0]] == 1:
            print("cycle detected at node", neighbor[0])
        
        neighbor_result = longest_path_in_DAG(graph, neighbor[0], end_node, dp_table, visited)
        if neighbor[1] + neighbor_result[0] > max_length:
            max_length = neighbor[1] + neighbor_result[0]
            path = neighbor_result[1].copy()
            path.reverse()
            path.append(current_node)
            path.reverse()
    
    dp_table[current_node] = [max_length, path]
    visited[current_node] = 0
    return dp_table[current_node]

def global_alignment_problem(seq1, seq2, idx1, idx2, match_score, mismatch_penalty, indel_penalty, dp_table):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    
    if dp_table[idx1][idx2][0] != 0.5: 
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1 and idx2 >= len_seq2:
        dp_table[idx1][idx2] = [0, '', '']
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1:
        dp_table[idx1][idx2] = [(len_seq2 - idx2) * indel_penalty * (-1), 
                               '-' * (len_seq2 - idx2), 
                               seq2[idx2:]]
        return dp_table[idx1][idx2]
    
    if idx2 >= len_seq2:
        dp_table[idx1][idx2] = [(len_seq1 - idx1) * indel_penalty * (-1), 
                               seq1[idx1:], 
                               '-' * (len_seq1 - idx1)]
        return dp_table[idx1][idx2]
    
    # Match or mismatch
    if seq1[idx1] == seq2[idx2]:
        align_result = global_alignment_problem(seq1, seq2, idx1+1, idx2+1, match_score, mismatch_penalty, indel_penalty, dp_table)
        score = align_result[0] + match_score
    else:
        align_result = global_alignment_problem(seq1, seq2, idx1+1, idx2+1, match_score, mismatch_penalty, indel_penalty, dp_table)
        score = align_result[0] - mismatch_penalty
    
    aligned_seq1 = seq1[idx1] + align_result[1]
    aligned_seq2 = seq2[idx2] + align_result[2]
    
    # Insertion in seq1 (gap in seq2)
    insert_seq1_result = global_alignment_problem(seq1, seq2, idx1, idx2+1, match_score, mismatch_penalty, indel_penalty, dp_table)
    insert_seq1_score = insert_seq1_result[0] - indel_penalty
    insert_seq1_aligned1 = '-' + insert_seq1_result[1]
    insert_seq1_aligned2 = seq2[idx2] + insert_seq1_result[2]
    
    # Insertion in seq2 (gap in seq1)
    insert_seq2_result = global_alignment_problem(seq1, seq2, idx1+1, idx2, match_score, mismatch_penalty, indel_penalty, dp_table)
    insert_seq2_score = insert_seq2_result[0] - indel_penalty
    insert_seq2_aligned1 = seq1[idx1] + insert_seq2_result[1]
    insert_seq2_aligned2 = '-' + insert_seq2_result[2]
    
    max_score = max(score, insert_seq1_score, insert_seq2_score)
    
    if max_score == score:
        dp_table[idx1][idx2] = [score, aligned_seq1, aligned_seq2]
    elif max_score == insert_seq1_score:
        dp_table[idx1][idx2] = [insert_seq1_score, insert_seq1_aligned1, insert_seq1_aligned2]
    else:
        dp_table[idx1][idx2] = [insert_seq2_score, insert_seq2_aligned1, insert_seq2_aligned2]
    
    return dp_table[idx1][idx2]

def global_alignment_problem_with_scoring_matrix(seq1, seq2, idx1, idx2, score_matrix, dp_table):
    nucleotide_map = {'A':0, 'G':1, 'C':2, 'T':3, '-':4}
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    
    if dp_table[idx1][idx2][0] != 'inf': 
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1 and idx2 >= len_seq2:
        dp_table[idx1][idx2] = [0, '', '']
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1:
        gap_penalty = 0
        for i in range(idx2, len_seq2):
            gap_penalty += score_matrix[nucleotide_map['-']][nucleotide_map[seq2[i]]]
        dp_table[idx1][idx2] = [gap_penalty, '-' * (len_seq2 - idx2), seq2[idx2:]]
        return dp_table[idx1][idx2]
    
    if idx2 >= len_seq2:
        gap_penalty = 0
        for i in range(idx1, len_seq1):
            gap_penalty += score_matrix[nucleotide_map[seq1[i]]][nucleotide_map['-']]
        dp_table[idx1][idx2] = [gap_penalty, seq1[idx1:], '-' * (len_seq1 - idx1)]
        return dp_table[idx1][idx2]
    
    # Match or mismatch
    align_result = global_alignment_problem_with_scoring_matrix(seq1, seq2, idx1+1, idx2+1, score_matrix, dp_table)
    score = align_result[0] + score_matrix[nucleotide_map[seq1[idx1]]][nucleotide_map[seq2[idx2]]]
    aligned_seq1 = seq1[idx1] + align_result[1]
    aligned_seq2 = seq2[idx2] + align_result[2]
    
    # Insertion in seq1 (gap in seq2)
    insert_seq1_result = global_alignment_problem_with_scoring_matrix(seq1, seq2, idx1+1, idx2, score_matrix, dp_table)
    insert_seq1_score = insert_seq1_result[0] + score_matrix[nucleotide_map[seq1[idx1]]][nucleotide_map['-']]
    insert_seq1_aligned1 = seq1[idx1] + insert_seq1_result[1]
    insert_seq1_aligned2 = '-' + insert_seq1_result[2]
    
    # Insertion in seq2 (gap in seq1)
    insert_seq2_result = global_alignment_problem_with_scoring_matrix(seq1, seq2, idx1, idx2+1, score_matrix, dp_table)
    insert_seq2_score = insert_seq2_result[0] + score_matrix[nucleotide_map['-']][nucleotide_map[seq2[idx2]]]
    insert_seq2_aligned1 = '-' + insert_seq2_result[1]
    insert_seq2_aligned2 = seq2[idx2] + insert_seq2_result[2]
    
    max_score = max(score, insert_seq1_score, insert_seq2_score)
    
    if max_score == score:
        dp_table[idx1][idx2] = [score, aligned_seq1, aligned_seq2]
    elif max_score == insert_seq1_score:
        dp_table[idx1][idx2] = [insert_seq1_score, insert_seq1_aligned1, insert_seq1_aligned2]
    else:
        dp_table[idx1][idx2] = [insert_seq2_score, insert_seq2_aligned1, insert_seq2_aligned2]
    
    return dp_table[idx1][idx2]

def local_alignment_problem_with_scoring_matrix(seq1, seq2, idx1, idx2, score_matrix, indel_penalty, dp_table):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    
    if dp_table[idx1][idx2][0] != 'inf': 
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1 and idx2 >= len_seq2:
        dp_table[idx1][idx2] = [0, '', '']
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1:
        dp_table[idx1][idx2] = [(len_seq2 - idx2) * indel_penalty * (-1), 
                               '-' * (len_seq2 - idx2), 
                               seq2[idx2:]]
        return dp_table[idx1][idx2]
    
    if idx2 >= len_seq2:
        dp_table[idx1][idx2] = [(len_seq1 - idx1) * indel_penalty * (-1), 
                               seq1[idx1:], 
                               '-' * (len_seq1 - idx1)]
        return dp_table[idx1][idx2]
    
    # Match or mismatch
    align_result = local_alignment_problem_with_scoring_matrix(seq1, seq2, idx1+1, idx2+1, score_matrix, indel_penalty, dp_table)
    score = align_result[0] + score_matrix[seq1[idx1]][seq2[idx2]]
    aligned_seq1 = seq1[idx1] + align_result[1]
    aligned_seq2 = seq2[idx2] + align_result[2]
    
    # Insertion in seq1 (gap in seq2)
    insert_seq1_result = local_alignment_problem_with_scoring_matrix(seq1, seq2, idx1+1, idx2, score_matrix, indel_penalty, dp_table)
    insert_seq1_score = insert_seq1_result[0] - indel_penalty
    insert_seq1_aligned1 = seq1[idx1] + insert_seq1_result[1]
    insert_seq1_aligned2 = '-' + insert_seq1_result[2]
    
    # Insertion in seq2 (gap in seq1)
    insert_seq2_result = local_alignment_problem_with_scoring_matrix(seq1, seq2, idx1, idx2+1, score_matrix, indel_penalty, dp_table)
    insert_seq2_score = insert_seq2_result[0] - indel_penalty
    insert_seq2_aligned1 = '-' + insert_seq2_result[1]
    insert_seq2_aligned2 = seq2[idx2] + insert_seq2_result[2]
    
    # Local alignment option (start new alignment)
    local_score = 0
    local_aligned1 = ''
    local_aligned2 = ''
    
    if idx1 == 0 and idx2 == 0:
        for i in range(len_seq1):
            for j in range(len_seq2):
                if i == 0 and j == 0: 
                    continue
                local_result = local_alignment_problem_with_scoring_matrix(seq1, seq2, i, j, score_matrix, indel_penalty, dp_table)
                if local_result[0] > local_score:
                    local_score = local_result[0]
                    local_aligned1 = local_result[1]
                    local_aligned2 = local_result[2]
    else:
        local_score = score_matrix[seq1[idx1]][seq2[idx2]]
        local_aligned1 = seq1[idx1]
        local_aligned2 = seq2[idx2]
    
    max_score = max(score, insert_seq1_score, insert_seq2_score, local_score)
    
    if max_score == score:
        dp_table[idx1][idx2] = [score, aligned_seq1, aligned_seq2]
    elif max_score == insert_seq1_score:
        dp_table[idx1][idx2] = [insert_seq1_score, insert_seq1_aligned1, insert_seq1_aligned2]
    elif max_score == insert_seq2_score:
        dp_table[idx1][idx2] = [insert_seq2_score, insert_seq2_aligned1, insert_seq2_aligned2]
    else:
        dp_table[idx1][idx2] = [local_score, local_aligned1, local_aligned2]
    
    return dp_table[idx1][idx2]

def fitting_alignment_problem(seq1, seq2, idx1, idx2, score_matrix, indel_penalty, dp_table):
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    
    if dp_table[idx1][idx2][0] != 'inf': 
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1 and idx2 >= len_seq2:
        dp_table[idx1][idx2] = [0, '', '']
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1:
        dp_table[idx1][idx2] = [(len_seq2 - idx2) * indel_penalty * (-1), 
                               '-' * (len_seq2 - idx2), 
                               seq2[idx2:]]
        return dp_table[idx1][idx2]
    
    if idx2 >= len_seq2:
        dp_table[idx1][idx2] = [(len_seq1 - idx1) * indel_penalty * (-1), 
                               seq1[idx1:], 
                               '-' * (len_seq1 - idx1)]
        return dp_table[idx1][idx2]
    
    # Match or mismatch
    align_result = fitting_alignment_problem(seq1, seq2, idx1+1, idx2+1, score_matrix, indel_penalty, dp_table)
    score = align_result[0] + score_matrix[seq1[idx1]][seq2[idx2]]
    aligned_seq1 = seq1[idx1] + align_result[1]
    aligned_seq2 = seq2[idx2] + align_result[2]
    
    # Insertion in seq1 (gap in seq2)
    insert_seq1_result = fitting_alignment_problem(seq1, seq2, idx1+1, idx2, score_matrix, indel_penalty, dp_table)
    insert_seq1_score = insert_seq1_result[0] - indel_penalty
    insert_seq1_aligned1 = seq1[idx1] + insert_seq1_result[1]
    insert_seq1_aligned2 = '-' + insert_seq1_result[2]
    
    # Insertion in seq2 (gap in seq1)
    insert_seq2_result = fitting_alignment_problem(seq1, seq2, idx1, idx2+1, score_matrix, indel_penalty, dp_table)
    insert_seq2_score = insert_seq2_result[0] - indel_penalty
    insert_seq2_aligned1 = '-' + insert_seq2_result[1]
    insert_seq2_aligned2 = seq2[idx2] + insert_seq2_result[2]
    
    # Special cases for fitting alignment
    fitting_score = 0
    fitting_aligned1 = ''
    fitting_aligned2 = ''
    
    if idx1 == 0 and idx2 == 0:
        for i in range(1, len_seq1):
            local_result = fitting_alignment_problem(seq1, seq2, i, 0, score_matrix, indel_penalty, dp_table)
            if local_result[0] > fitting_score:
                fitting_score = local_result[0]
                fitting_aligned1 = local_result[1]
                fitting_aligned2 = local_result[2]
    elif idx2 == len_seq2 - 1 and idx1 != len_seq1 - 1:
        fitting_score = score_matrix[seq1[idx1]][seq2[idx2]]
        fitting_aligned1 = seq1[idx1]
        fitting_aligned2 = seq2[idx2]
    else:
        max_score = max(score, insert_seq1_score, insert_seq2_score)
        if max_score == score:
            dp_table[idx1][idx2] = [score, aligned_seq1, aligned_seq2]
        elif max_score == insert_seq1_score:
            dp_table[idx1][idx2] = [insert_seq1_score, insert_seq1_aligned1, insert_seq1_aligned2]
        else:
            dp_table[idx1][idx2] = [insert_seq2_score, insert_seq2_aligned1, insert_seq2_aligned2]
        return dp_table[idx1][idx2]
    
    max_score = max(score, insert_seq1_score, insert_seq2_score, fitting_score)
    
    if max_score == score:
        dp_table[idx1][idx2] = [score, aligned_seq1, aligned_seq2]
    elif max_score == insert_seq1_score:
        dp_table[idx1][idx2] = [insert_seq1_score, insert_seq1_aligned1, insert_seq1_aligned2]
    elif max_score == insert_seq2_score:
        dp_table[idx1][idx2] = [insert_seq2_score, insert_seq2_aligned1, insert_seq2_aligned2]
    else:
        dp_table[idx1][idx2] = [fitting_score, fitting_aligned1, fitting_aligned2]
    
    return dp_table[idx1][idx2]

def gap_penalties_alignment_problem(seq1, seq2, idx1, idx2, match_score, mismatch_penalty, gap_open, gap_extend, in_gap, dp_table):
    current_gap_penalty = gap_open
    if in_gap == 1: 
        current_gap_penalty = gap_extend
    
    len_seq1 = len(seq1)
    len_seq2 = len(seq2)
    
    if dp_table[idx1][idx2][0] != 'inf': 
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1 and idx2 >= len_seq2:
        dp_table[idx1][idx2] = [0, '', '']
        return dp_table[idx1][idx2]
    
    if idx1 >= len_seq1:
        dp_table[idx1][idx2] = [(len_seq2 - idx2 - 1) * in_gap * (-1) - current_gap_penalty, 
                               '-' * (len_seq2 - idx2), 
                               seq2[idx2:]]
        return dp_table[idx1][idx2]
    
    if idx2 >= len_seq2:
        dp_table[idx1][idx2] = [(len_seq1 - idx1 - 1) * current_gap_penalty * (-1) - current_gap_penalty, 
                               seq1[idx1:], 
                               '-' * (len_seq1 - idx1)]
        return dp_table[idx1][idx2]
    
    # Match or mismatch
    if seq1[idx1] == seq2[idx2]:
        align_result = gap_penalties_alignment_problem(seq1, seq2, idx1+1, idx2+1, match_score, mismatch_penalty, gap_open, gap_extend, 0, dp_table)
        score = align_result[0] + match_score
    else:
        align_result = gap_penalties_alignment_problem(seq1, seq2, idx1+1, idx2+1, match_score, mismatch_penalty, gap_open, gap_extend, 0, dp_table)
        score = align_result[0] - mismatch_penalty
    
    aligned_seq1 = seq1[idx1] + align_result[1]
    aligned_seq2 = seq2[idx2] + align_result[2]
    
    # Insertion in seq1 (gap in seq2)
    insert_seq1_result = gap_penalties_alignment_problem(seq1, seq2, idx1+1, idx2, match_score, mismatch_penalty, gap_open, gap_extend, 1, dp_table)
    insert_seq1_score = insert_seq1_result[0] - current_gap_penalty
    insert_seq1_aligned1 = seq1[idx1] + insert_seq1_result[1]
    insert_seq1_aligned2 = '-' + insert_seq1_result[2]
    
    # Insertion in seq2 (gap in seq1)
    insert_seq2_result = gap_penalties_alignment_problem(seq1, seq2, idx1, idx2+1, match_score, mismatch_penalty, gap_open, gap_extend, 1, dp_table)
    insert_seq2_score = insert_seq2_result[0] - current_gap_penalty
    insert_seq2_aligned1 = '-' + insert_seq2_result[1]
    insert_seq2_aligned2 = seq2[idx2] + insert_seq2_result[2]
    
    max_score = max(score, insert_seq1_score, insert_seq2_score)
    
    if max_score == score:
        dp_table[idx1][idx2] = [score, aligned_seq1, aligned_seq2]
    elif max_score == insert_seq1_score:
        dp_table[idx1][idx2] = [insert_seq1_score, insert_seq1_aligned1, insert_seq1_aligned2]
    else:
        dp_table[idx1][idx2] = [insert_seq2_score, insert_seq2_aligned1, insert_seq2_aligned2]
    
    return dp_table[idx1][idx2]

def middle_edge_in_alignment_graph(seq1, seq2, x_start, x_end, y_start, y_end, match_score, mismatch_penalty, indel_penalty):
    if x_start == x_end and y_start == y_end: 
        return None
    
    if y_start == y_end:
        mid_x = (x_start + x_end) // 2
        return (mid_x, y_start), (mid_x + 1, y_start)
    
    if x_start == x_end:
        mid_y = (y_start + y_end) // 2
        return (x_start, mid_y), (x_start, mid_y + 1)
    
    # Convert penalties to positive values for easier calculation
    indel_penalty *= -1
    mismatch_penalty *= -1
    
    mid_x = (x_start + x_end) // 2
    num_cols = y_end - y_start + 1
    
    # Forward pass
    forward_scores = []
    prev_forward = [i * indel_penalty for i in range(num_cols)]
    
    if x_start == mid_x: 
        forward_scores = prev_forward
    
    for i in range(x_start, mid_x):
        current_forward = ['-'] * num_cols
        for j in range(num_cols):
            if j != num_cols - 1:
                if seq1[i] == seq2[y_start + j]:
                    if current_forward[j + 1] == '-' or prev_forward[j] + match_score > current_forward[j + 1]:
                        current_forward[j + 1] = prev_forward[j] + match_score
                else:
                    if current_forward[j + 1] == '-' or prev_forward[j] + mismatch_penalty > current_forward[j + 1]:
                        current_forward[j + 1] = prev_forward[j] + mismatch_penalty
                
                prev_forward[j + 1] = max(prev_forward[j] + indel_penalty, prev_forward[j + 1])
            
            if current_forward[j] == '-' or prev_forward[j] + indel_penalty > current_forward[j]:
                current_forward[j] = prev_forward[j] + indel_penalty
        
        if i != mid_x - 1: 
            prev_forward = current_forward.copy()
    
    # Backward pass
    backward_scores = []
    prev_backward = [(num_cols - i - 1) * indel_penalty for i in range(num_cols)]
    
    if x_end == mid_x: 
        backward_scores = prev_backward
    
    for i in range(x_end, mid_x, -1):
        current_backward = ['-'] * num_cols
        for j in range(num_cols - 1, -1, -1):
            if j != 0:
                if seq1[i - 1] == seq2[y_start + j - 1]:
                    if current_backward[j - 1] == '-' or prev_backward[j] + match_score > current_backward[j - 1]:
                        current_backward[j - 1] = prev_backward[j] + match_score
                else:
                    if current_backward[j - 1] == '-' or prev_backward[j] + mismatch_penalty > current_backward[j - 1]:
                        current_backward[j - 1] = prev_backward[j] + mismatch_penalty
                
                prev_backward[j - 1] = max(prev_backward[j] + indel_penalty, prev_backward[j - 1])
            
            if current_backward[j] == '-' or prev_backward[j] + indel_penalty > current_backward[j]:
                current_backward[j] = prev_backward[j] + indel_penalty
        
        if i != mid_x + 1: 
            prev_backward = current_backward.copy()
    
    # Combine scores
    combined_scores = []
    for i in range(num_cols):
        try:
            combined_scores.append(backward_scores[i] + forward_scores[i])
        except:
            print("Error in score calculation at index", i)
            quit()
    
    max_score = 0
    if len(combined_scores):
        max_score = max(combined_scores)
    
    max_index = 0
    for i in range(num_cols):
        if combined_scores[i] == max_score:
            max_index = i
            break
    
    if max_index == num_cols - 1:
        return (mid_x, max_index + y_start), (mid_x + 1, max_index + y_start)
    
    if seq1[mid_x] == seq2[max_index + y_start]:
        return (mid_x, max_index + y_start), (mid_x + 1, max_index + 1 + y_start)
    
    # Determine the best edge
    diagonal_score = prev_backward[max_index + 1] + mismatch_penalty
    down_score = prev_backward[max_index] + indel_penalty
    right_score = backward_scores[max_index + 1] + indel_penalty
    
    if diagonal_score >= down_score and diagonal_score >= right_score:
        return (mid_x, max_index + y_start), (mid_x + 1, max_index + y_start + 1)
    elif down_score >= right_score and down_score >= diagonal_score:
        return (mid_x, max_index + y_start), (mid_x + 1, max_index + y_start)
    else:
        return (mid_x, max_index + y_start), (mid_x, max_index + 1 + y_start)

def linear_space_alignment(seq1, seq2, x_start, x_end, y_start, y_end, match_score, mismatch_penalty, indel_penalty, edges):
    if x_start == x_end and y_start == y_end: 
        return
    
    start_edge, end_edge = middle_edge_in_alignment_graph(seq1, seq2, x_start, x_end, y_start, y_end, match_score, mismatch_penalty, indel_penalty)
    edges[start_edge] = end_edge
    
    # Recurse on the left part
    linear_space_alignment(seq1, seq2, x_start, start_edge[0], y_start, start_edge[1], match_score, mismatch_penalty, indel_penalty, edges)
    
    # Recurse on the right part
    linear_space_alignment(seq1, seq2, end_edge[0], x_end, end_edge[1], y_end, match_score, mismatch_penalty, indel_penalty, edges)

def two_break_sorting(permutation, target):
    target_edges = set()
    permutation_map = {}
    
    # Create edges from target permutation
    for i in range(len(target)):
        if i == len(target) - 1:
            target_edges.add((target[i], -target[0]))
            target_edges.add((-target[0], target[i]))
        else:
            target_edges.add((target[i], -target[i + 1]))
            target_edges.add((-target[i + 1], target[i]))
    
    num_blocks = len(target_edges) // 2
    
    # Create edges from initial permutation
    for i in range(len(permutation)):
        if i == len(permutation) - 1:
            permutation_map[permutation[i]] = -permutation[0]
            permutation_map[-permutation[0]] = permutation[i]
        else:
            permutation_map[permutation[i]] = -permutation[i + 1]
            permutation_map[-permutation[i + 1]] = permutation[i]
    
    breakpoint_graph = []
    nodes = set()
    visited = set()
    
    for i in range(1, num_blocks + 1):
        nodes.add(i)
        nodes.add(-i)
    
    # Find cycles in initial permutation
    for node in nodes:
        if node in visited: 
            continue
        
        cycle = [node]
        visited.add(node)
        visited.add(-node)
        
        next_node = permutation_map[node]
        next_node *= -1
        
        while next_node != node:
            cycle.append(next_node)
            visited.add(next_node)
            visited.add(-next_node)
            next_node = permutation_map[next_node]
            next_node *= -1
        
        breakpoint_graph.append([cycle])
    
    # Perform 2-breaks to transform permutation into target
    for edge in target_edges:
        x = edge[0]
        y = edge[1]
        x_next = permutation_map[x]
        y_next = permutation_map[y]
        
        if x_next == y and y_next == x: 
            continue
        
        # Perform the 2-break
        permutation_map[x] = y
        permutation_map[y] = x
        permutation_map[x_next] = y_next
        permutation_map[y_next] = x_next
        
        # Find the new cycles
        nodes = set()
        visited = set()
        current_breakpoint_graph = []
        
        for i in range(1, num_blocks + 1):
            nodes.add(i)
            nodes.add(-i)
        
        for node in nodes:
            if node in visited: 
                continue
            
            cycle = [node]
            visited.add(node)
            visited.add(-node)
            
            next_node = permutation_map[node]
            next_node *= -1
            
            while next_node != node:
                cycle.append(next_node)
                visited.add(next_node)
                visited.add(-next_node)
                next_node = permutation_map[next_node]
                next_node *= -1
            
            current_breakpoint_graph.append(cycle)
        
        breakpoint_graph.append(current_breakpoint_graph)
    
    return breakpoint_graph
