# Using Hamming's distance to determine if we can accept the amount of missmatches

from bm_preproc import BoyerMoore


def boyer_moore(p, p_bm, t):  # p: pattern , p_bm: Boyer Moore class object, t: index
    i = 0
    occurences = []
    while i < len(t) - len(p) + 1:
        shift = 1  # how many chars should be skipped
        missmatched = False
        for j in range(len(p) - 1, -1, -1):  # right-to-left search (last to first)
            if p[j] != t[i + j]:  # missmatch -> bad_char_rule -> good_sufix_rule
                skip_bc = p_bm.bad_character_rule(j, t[i + j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(skip_bc, skip_gs)
                missmatched = True
                break
        if missmatched is False:
            occurences.append(i)
            skip_gs = p_bm.match_skip()  # skipp pattern in text
            shift = max(shift, skip_gs)
        i += shift
    return occurences


def naive_hamming(p, t, max_distance):
    occurences = []
    for i in range(len(t) - len(p) + 1):
        missmatches = 0
        for j in range(len(p)):
            if p[j] != t[i + j]:
                missmatches += 1
                if missmatches > max_distance:
                    break
        if missmatches <= max_distance:
            occurences.append(i)
    return occurences


G = 'ACTGCTGACATG'
P = 'TGGT'
naive_hamming(P, G, 2)


# Trying to find where our pattern matches with the text with least edits
# If we know there are k edits in the align, split the pattern into k+1 partitions ->
# WORST CASE scenario - there will be AT LEAST 1 partition which we can match with EXACT MATCHING ALGORITHM
#   and every partition needs one edit (and one exact)...this is not usually the case
# after we split the pattern, we run any of the EXACT MATCHING algorithms on all partitions
# once the partition is aligned to the text,its an indicator that there might be an approximate alignment
# therefore we should do a verification

def approximate_match(p, t, n):  # n: max distance
    segment_length = int(round(len(p) / (n + 1)))
    # we use a set because if the whole pattern exact matches,every partition will match,
    # verification for every partition will succeed ,and every partition would return
    # a begining to the pattern,so to avoid repetition we use a set
    all_matches = set()
    for i in range(n + 1):
        start = i * segment_length
        end = min((i + 1) * segment_length, len(p))
        # Make BoyerMoore object which preprocesses the pattern to use good suffix and bad character rules
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = boyer_moore(p[start:end], p_bm, t)
        for m in matches:
            # check if pattern goes out of bounds
            if m - start + len(p) > len(t) or m < start:
                continue
            missmatches = 0
            for j in range(0, start):  # verification for pattern before the begining of partition
                if p[j] != t[m - start + j]:
                    missmatches += 1
                if missmatches > n:
                    break

            for j in range(end, len(p)):  # verification for pattern after end of partition
                if p[j] != t[m - start + j]:
                    missmatches += 1
                if missmatches > n:
                    break
            if missmatches <= n:
                all_matches.add(m - start)
    return list(all_matches)


p1 = 'AACTTG'
t1 = 'CACTTAATTTG'
print(approximate_match(p1, t1, 2))


def hamming_distance(x, y):
    number_missmatches = 0
    for i in range(len(x)):
        if x[i] != y[i]:
            number_missmatches += 1
    return number_missmatches


# Edit / Levensteins distance
def edit_distance_recursive(a, b):
    # If one string is empty,distance is length(other)
    if len(a) == 0:
        return len(b)
    if len(b) == 0:
        return len(a)
    # if last letters in string are equal ,delta=0,else delta =1
    if a[-1] == b[-1]:
        delt = 0
    else:
        delt = 1
    # minimum of 3 distances recursively
    return min(edit_distance_recursive(a[:-1], b[:-1]) + delt,
               edit_distance_recursive(a[:-1], b) + 1,
               edit_distance_recursive(a, b[:-1]) + 1)


def edit_distance(x, y):
    # Matrix of dynamic programming (x+1) * (y+1) , +1 due to row and collumn of empty words
    D = []
    for i in range(len(x) + 1):
        D.append([0] * (len(y) + 1))
    # Fill first row and column
    for i in range(len(x) + 1):
        D[i][0] = i
    for i in range(len(y) + 1):
        D[0][i] = i
    # Fill the rest of matrix, fields in matrix are calculated based on the field left from it
    # or field above or diagonal up-left
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            horizontal_distance = D[i][j - 1] + 1  # field left
            vertical_distance = D[i - 1][j] + 1  # field above
            if x[i - 1] == y[j - 1]:  # diagonal field
                diagonal_distance = D[i - 1][j - 1]
            else:
                diagonal_distance = D[i - 1][j - 1] + 1
            D[i][j] = min(horizontal_distance, vertical_distance, diagonal_distance)
    return D[-1][-1]  # Lower right corner is returned
