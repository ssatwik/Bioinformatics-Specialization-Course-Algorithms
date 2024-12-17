import random
import copy

def pattern_matching(seq,pattern):
    # returns all indices where a pattern is found in a seq
    s1=seq
    s2=pattern
    n=len(s1)
    k=len(s2)
    lis=[]
    for i in range(0,n-k+1):
        f=1
        for j in range(i,i+k):
            if s1[j]!=s2[j-i]:
                f=0
                break
        if f: lis.append(i)
    return lis

def reverse_complement(s):
    # returns reverse complement of a sequence
    m={
        'A':'T',
        'T':'A',
        'G':'C',
        'C':'G'
    }
    s=s.upper()
    ss=''
    for i in range(len(s)-1,-1,-1): ss+=m[s[i]]
    return ss

def frequency_table(seq,k):
    # returns frequencies of all k-mers in a sequence given seq and k
    s=seq
    n=len(s)
    m={}
    for i in range(0,n-k+1):
        ss=s[i:i+k]
        if ss not in m: m[ss]=0
        m[ss]+=1
    return m

def find_clump_pattern(seq,k,L,t):
    # returns all k-mers occuring atleast t times in a window of length L in the sequence
    s=seq
    n=len(s)
    lis=set()
    for i in range(0,n-L+1):
        ss=s[i:i+L]
        m=frequency_table(ss,k)
        for p,v in m.items():
            if v>=t: lis.add(p)
    return list(lis)

def skew(seq):
    # returns skew list for a sequence
    s=seq
    skew=0
    lis=[]
    for i in range(len(seq)):
        if s[i]=='G': skew+=1
        elif s[i]=='C': skew-=1
        lis.append(skew)
    return lis

def hamming_distance(p,q):
    # returns hamming distance of 2 sequences of equal length
    h=0
    for i in range(len(p)):
        if p[i]!=q[i]: h+=1
    return h

def approximate_pattern_matching(seq,pattern,d):
    # returns indices of k-mers with hamming distance to pattern atmost d
    s1=seq
    s2=pattern
    n=len(s1)
    k=len(s2)
    lis=[]
    for i in range(n-k+1):
        ss=s1[i:i+k]
        if hamming_distance(ss,s2)<=d: lis.append(i)
    return lis

def neighbors(pattern,d):
    # returns all strings with hamming distance to pattern atmost d
    s=pattern
    n=len(s)
    lis=[]
    if n<1: return lis
    elif d==0: return [pattern,]
    elif n==1: return ['A','C','G','T']
    ss=pattern[:n-1]
    ll=neighbors(ss,d)
    b=('A','C','G','T')
    for st in ll:
        if hamming_distance(st,s)==d:
            lis.append(st+s[n-1])
        else:
            for bb in b:
                lis.append(st+bb)
    return lis

def better_frequent_kmers(seq,k,d):
    # returns max frequency kmers with hamming distance to sequence atmost d
    s=seq
    n=len(s)
    m={}
    for i in range(n-k+1):
        ss=s[i:i+k]
        ll=neighbors(ss,d)
        for x in ll:
            if x not in m: m[x]=0
            m[x]+=1
    mm=max(m.values())
    lis=[]
    for p,v in m.items():
        if v==mm: lis.append(p)
    return lis

def better_frequent_kmers_with_rc(seq,k,d):
    # returns max frequency kmers with hamming distance to sequence atmost d including reverse compliments
    s=seq
    n=len(s)
    m={}
    for i in range(n-k+1):
        ss=s[i:i+k]
        rss=reverse_complement(ss)
        ll=neighbors(ss,d)
        rll=neighbors(rss,d)
        for x in ll:
            if x not in m: m[x]=0
            m[x]+=1
        for x in rll:
            if x not in m: m[x]=0
            m[x]+=1
    mm=max(m.values())
    lis=[]
    for p,v in m.items():
        if v==mm: lis.append(p)
    return lis

def motif_enumeration(seq_list,k,d):
    # returns all kmeres with hamming distance atmost d to each sequence in the list
    pat=set()
    for s in seq_list:
        for i in range(len(s)-k+1):
            pat.add(s[i:i+k])
    lis=set()
    for s in pat:
        for ss in neighbors(s,d): lis.add(ss)
    ans=set()
    for s in lis:
        for seq in seq_list:
            f=1
            if len(approximate_pattern_matching(seq,s,d))==0:
                f=0
                break
        if f==1: ans.add(s)
    return list(ans)
        
def pattern_and_list_distance(seq_list,pat):
    # returns sum of hamming distances from pattern to all sequences in the list
    ss=pat
    ans=0
    k=len(ss)
    for s in seq_list:
        n=len(s)
        dist=n
        for i in range(n-k+1):
            a=hamming_distance(ss,s[i:i+k])
            if dist>a: dist=a
        ans+=dist
    return ans

def nucleotide_bases_permutations(n):
    # returns all possible dna sequences of length n
    b=['A','C','G','T']
    l1=['',]
    l2=[]
    for i in range(n):
        for s in l1:
            for bb in b:
                l2.append(s+bb)
        l1=l2.copy()
        l2=[]
    return l1

def median_string(seq_list,k):
    # returns list of kmers with least sum of hamming distances to all sequences in the list
    ll=nucleotide_bases_permutations(k)
    ans=[]
    d=len(seq_list)*k
    m={}
    for ss in ll:
        dd=pattern_and_list_distance(seq_list,ss)
        m[ss]=dd
    mi=min(m.values())
    for p,v in m.items():
        if v==mi: ans.append(p)
    return ans

def profile_without_pseudocounts(motifs):
    s=motifs
    n=len(s)
    k=len(s[0])
    m={'A':[0]*k,'C':[0]*k,'G':[0]*k,'T':[0]*k}
    for i in range(k):
        for j in range(n):
            m[s[j][i]][i]+=1
    for i in 'ACGT':
        for j in range(k):
            m[i][j]/=n
    lis=[m[key] for key in 'ACGT']
    return lis
    

def profile_with_pseudocounts(motifs):
    s=motifs
    n=len(s)
    k=len(s[0])
    m={'A':[1]*k,'C':[1]*k,'G':[1]*k,'T':[1]*k}
    for i in range(k):
        for j in range(n):
            m[s[j][i]][i]+=1
    for i in 'ACGT':
        for j in range(k):
            m[i][j]/=(n+4)
    lis=[m[key] for key in 'ACGT']
    return lis

def profile_most_probable_kmer(profile,seq,k):
    s=seq
    n=len(s)
    mx=0
    ans=s[:k]
    m={'A':0,'C':1,'G':2,'T':3}
    for i in range(n-k+1):
        ss=s[i:i+k]
        d=1
        for j in range(0,k):
            d*=profile[m[ss[j]]][j]
        if d>mx:
            mx=d
            ans=ss
    return ans

def profile_kmer_probabilities(profile,seq,k):
    s=seq
    n=len(s)
    lis=[]
    m={'A':0,'C':1,'G':2,'T':3}
    for i in range(n-k+1):
        ss=s[i:i+k]
        d=1
        for j in range(0,k):
            d*=profile[m[ss[j]]][j]
        lis.append(d)
    return lis

def consensus(l):
    # returns consensus for a given list of motifs
    n=len(l)
    m=len(l[0])
    con=''
    for i in range(m):
        a,t,g,c=0,0,0,0
        for j in range(n):
            if l[j][i]=='A': a+=1
            elif l[j][i]=='T': t+=1
            elif l[j][i]=='G': g+=1
            else: c+=1
        e=max(a,t,g,c)
        if e==a: con+='A'
        elif e==c: con+='C'
        elif e==g: con+='G'
        else: con+='T'
    return con

def score_of_motifs(l,con):
    # returns score for a given list of motifs and consensus
    ans=0
    n=len(l)
    m=len(l[0])
    for i in range(m):
        for j in range(n):
            if l[j][i]!=con[i]: ans+=1
    return ans


def greedy_motif_search_with_pseudocounts(seq_list,k):
    n=len(seq_list)
    s=seq_list
    best_motifs=[]
    for ss in s: best_motifs.append(ss[:k])
    best_score=score_of_motifs(best_motifs,consensus(best_motifs))
    for i in range(len(seq_list[0])-k+1):
        motifs=[]
        count={'A':[0]*k,'C':[0]*k,'G':[0]*k,'T':[0]*k}
        ss=s[0][i:i+k]
        for j in range(k):
            count[ss[j]][j]+=1
        motifs.append(ss)
        for t in range(1,n):
            lis=[count[key] for key in 'ACGT']
            for xx in range(len(lis)):
                for yy in range(len(lis[0])):
                    lis[xx][yy]+=1
            p=profile_most_probable_kmer(lis,s[t],k)
            motifs.append(p)
            for j in range(k):
                count[p[j]][j]+=1
        sc=score_of_motifs(motifs,consensus(motifs))
        if sc<best_score:
            best_score=sc
            best_motifs=motifs.copy()
    return best_motifs

def randomized_motif_search(seq_list,k):
    # returns best_motifs and its score
    s=seq_list
    t=len(seq_list)
    n=len(s[0])
    best_motifs=[]
    for i in range(t):
        a=random.randint(0,n-k)
        best_motifs.append(s[i][a:a+k])
    best_profile=profile_with_pseudocounts(best_motifs)
    while True:
        motifs=[]
        for i in range(t):
            motifs.append(profile_most_probable_kmer(best_profile,s[i],k))
        profile=profile_with_pseudocounts(motifs)
        best_score=score_of_motifs(best_motifs,consensus(best_motifs))
        score=score_of_motifs(motifs,consensus(motifs))
        if score<best_score:
            best_motifs=motifs.copy()
            best_profile=copy.deepcopy(profile)
        else:
            return best_motifs,best_score

def gibbs_sampler(seq_list,k,N):
    s=seq_list
    t=len(seq_list)
    n=len(s[0])
    best_motifs=[]
    for i in range(t):
        a=random.randint(0,n-k)
        best_motifs.append(s[i][a:a+k])
    best_score=score_of_motifs(best_motifs,consensus(best_motifs))
    motifs=[]
    motifs=best_motifs.copy()
    for q in range(N):
        a=random.randint(0,t-1)
        del motifs[a]
        profile=profile_with_pseudocounts(motifs)
        ll=[]
        for i in range(n-k+1): ll.append(i)
        x=random.choices(ll,weights=profile_kmer_probabilities(profile,s[a],k),k=1)
        x=x[0]
        motifs.insert(a,s[a][x:x+k])
        score=score_of_motifs(motifs,consensus(motifs))
        if score<best_score:
            best_motifs=motifs.copy()
            best_score=score
    return best_motifs,best_score

