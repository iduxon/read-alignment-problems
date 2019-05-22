# Similar to naive-matching , but a bit smarter,skipping some already checked positions to increase performance

# Logic: T= base DNA string P = pattern which needs to be found

# Rule of 'bad character'
# If a missmatch happens on an 'i' position in T,start checking backwards in P untill
#   that character is reached -> then realign P and T,and check again
#   if the end of P is reached without reaching the T(i) character ,move the whole strand

# Rule of good sufix - same logic but instead of using the 1 missmatch character which we look for in P string
#   ,we will use the MATCHED pattern 's' (s is a string match in P and T)
#   and look for the whole string in P, if found,realign
#   if not found,move to the right untill matching

# B-M algorithm : mix of these two rules


#                       *****************************
#                       ********* IMPORTANT *********
#                       *****************************
# for this code to work we need a Boyer Moore class,
# methods : good_suffix_rule(self, i) , match_skip(self), bad_character_rule(self,i ,c )
# -match_skip is used when the pattern is found
# the code for this has been downloaded from the web and can be found in the 'bm_preproc.py' file
# therefore,the BM algorithm uses preprocessing ,unlike naive-matching algorithm which does not


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


t = 'GCTAGCTCTACGAGTCTA'
p = 'TCTA'
p_bm = BoyerMoore(p, alphabet='ACGT')
print(boyer_moore(p, p_bm, t))
