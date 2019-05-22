# 2 types of indexing (preprocessing a genome) : binary search (multimaps) and hash maps (dictionary)

# **************************************************************************
# sorted table of k-mers :
# making a table(key : value) where key = k-mer (k characters) and value = position(s) of this k-mer in the read
# if we were looking for (k+n)-mer in a read,we can align our patern to the k-mer (from value position) and check the remaining n-mer
# **************************************************************************

# enables binary search
import bisect


# constructor does the preprocessing of the text
class Index(object):
    def __init__(self, t, k):  # t:text , k:size of kmer
        self.k = k
        self.index = []  # index table of all kmers in text,elements are touples(kmer,position)
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i + k], i))
        self.index.sort()  # sorting the index table makes an efficient binary search

    def query(self, p):  # pattern p
        kmer = p[:self.k]  # p will be the pattern from begining
        i = bisect.bisect_left(self.index, (kmer, -1))  # 1st position of kmer
        # kmer gets the list which we are searching, -1 returns first iteration
        hits = []
        # there can be multiple kmers which are the same ,but on different positions
        # but the list is sorted , so if the first kmer has been found, we advance in the list
        # untill end is reached, or the searched kmer has been found
        # if we get a hit,we put it in hits
        while i < len(self.index):
            if self.index[i][0] != kmer:  # if the kmer is different from curent,break
                break
            hits.append(self.index[i][1])  # else add its index
            i += 1
        return hits  # return list of offsets in text where kmer appears

    # p: pattern ,t: text, index: object


def query_index(p, t, index):
    k = index.k
    offsets = []  # list of offsets of found patterns
    for i in index.query(p):  # for every HIT we do a verification
        if p[k:] == t[i + k:i + len(p)]:
            offsets.append(i)
    return offsets


t = 'GCTACGATCTAGAATCTA'
p = 'TCTA'

index = Index(t, 2)
print(query_index(p, t, index))
