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
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
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

amino_acid_mass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

amino_weights = {57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186}


rm = {
    57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'K', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'
}


def output_file(lis):
    n=len(lis)
    with open("output.txt", "a") as file:
        for i in range(n):
            file.write(f"{lis[i]} ")  

def input_file(path):
    with open(path, "r") as file:
        content = file.read()
    lis = content.split()
    return lis

def sequence_composition(seq,k):
    lis=[]
    for i in range(len(seq)-k+1): lis.append(seq[i:i+k])
    return lis

def ordered_kmers_to_sequence(seq_list):
    n=len(seq_list)
    k=len(seq_list[0])
    s=seq_list[0]
    for i in range(1,n): s+=seq_list[i][k-1]
    return s

def overlap_graph(seq_list):
    n=len(seq_list)
    k=len(seq_list[0])
    m={}
    for s in seq_list:
        a=s[:k-1]
        if a not in m: m[a]=[]
        m[a].append(s)
    mm={}
    for s in seq_list:
        a=s[1:]
        if a not in m: continue
        if s[1:]==s[:k-1]: m[a].remove(s)
        mm[s]=m[a]
    return mm

def debruijn_graph(seq_list):
    m={}
    n=len(seq_list)
    k=len(seq_list[0])
    for s in seq_list:
        a=s[:k-1]
        if a not in m: m[a]=[]
        m[a].append(s[1:])
    # return m
    # to return edge list instead of adj map
    e=[]
    for k,v in m.items():
        for x in v:
            e.append((k,x))
    return e

def paired_debruijn_graph(seq_list):
    m={}
    n=len(seq_list)
    k=len(seq_list[0][0])
    for s in seq_list:
        a=s[0][:k-1]
        b=s[1][:k-1]
        if (a,b) not in m: m[(a,b)]=[]
        m[(a,b)].append((s[0][1:],s[1][1:]))
    # return m
    # to return edge list instead of adj map
    e=[]
    for k,v in m.items():
        for x in v:
            e.append((k,x))
    return e

def eulerian_path(edge_list):
    # edge list is a list of tuples. eg: [(0, 1), (1, 2), (2, 0), (0, 3), (3, 1)]
    G = nx.DiGraph()
    G.add_edges_from(edge_list)
    path = list(nx.eulerian_path(G))
    return path

def path_to_seq(path):
    s=path[0][0]
    n=len(path)
    k=len(path[0][0])
    for i in range(1,n): s+=path[i][0][k-1]
    s+=path[n-1][1][k-1]
    return s

def paired_path_to_seq(path,d):
    n=len(path)
    k=len(path[0][0][0])
    size=2*(k+1)+n+d-1
    f=[-1 for i in range(size)]
    b=[-1 for i in range(size)]
    ff=k
    bb=size-k-1
    for i in range(k): f[i]=path[0][0][0][i]
    for i in range(k-1,-1,-1): b[size+i-k]=path[n-1][1][1][i]
    for i in range(1,n):
        f[ff]=path[i][0][0][k-1]
        ff+=1
    f[ff]=path[n-1][1][0][k-1]
    for i in range(n-2,-1,-1):
        b[bb]=path[i][1][1][0]
        bb-=1
    b[bb]=path[0][0][1][0]
    ans=''
    for i in range(size):
        if f[i]==-1 and b[i]==-1:
            print("both -1")
        elif f[i]==-1:
            ans+=b[i]
        elif b[i]==-1:
            ans+=f[i]
        elif f[i]!=b[i]:
            print("not -1 but not equal")
        else:
            ans+=f[i]
    return ans

def contigs_from_kmers(ll):
    t=len(ll)
    p=len(ll[0])
    ans=[]
    m={}
    for s in ll:
        if s[:p-1] not in m: m[s[:p-1]]=[]
        m[s[:p-1]].append(s[1:])
    n = len(m.items())
    indeg={}
    outdeg={}
    nice={}
    for k,v in m.items():
        outdeg[k]=len(v)
        for x in v:
            if x not in indeg: indeg[x]=0
            indeg[x]+=1
    for k,v in indeg.items():
        if k in outdeg:
            if v==outdeg[k] and v==1:
                nice[k]=1
    for k,v in m.items():
        if k in nice: continue
        kk=k
        while len(v):
            s=kk
            k=v.pop()
            while k in nice:
                s+=k[p-2]
                k=m[k][0]
            s+=k[p-2]
            ans.append(s)
    return ans

def rna_to_peptide(s):
    n=len(s)
    ans=''
    for i in range(0,n,3):
        x=codon_to_aa[s[i:i+3]]
        if x=='Stop': continue
        ans+=x
    return ans

def peptide_to_rna(s):
    ans=[]
    f=0
    for a in s:
        temp=[]
        if f==0:
            f=1
            l=aa_to_codon[a]
            for ss in l:
                temp.append(ss)
        else:
            for i in range(len(ans)):
                l=aa_to_codon[a]
                for ss in l:
                    temp.append(ans[i]+ss)
        ans=temp.copy()
    return ans

def cycle_spectrum_of_peptide(s):
    n=len(s)
    ans=[]
    ans.append(0)
    lis=[0 for i in range(n+1)]
    for i in range(1,n+1):
        lis[i]=lis[i-1]+amino_acid_mass[s[i-1]]
        ans.append(lis[i])
    for k in range(1,n):
        for b in range(0,n):
            if b==0:
                sum=lis[k]
                continue
            sum+=amino_acid_mass[s[(b+k-1)%n]]
            sum-=amino_acid_mass[s[b-1]]
            ans.append(sum)
    ans.sort()
    return ans

def linear_spectrum_of_peptide(s):
    n=len(s)
    ans=[]
    ans.append(0)
    lis=[0 for i in range(n+1)]
    for i in range(1,n+1):
        lis[i]=lis[i-1]+amino_acid_mass[s[i-1]]
        ans.append(lis[i])
    for k in range(1,n):
        for b in range(0,n):
            if b==0:
                sum=lis[k]
                continue
            if b+k-1>=n: break
            sum+=amino_acid_mass[s[b+k-1]]
            sum-=amino_acid_mass[s[b-1]]
            ans.append(sum)
    ans.sort()
    return ans

def cycle_spectrum_of_peptide_weights(ll):
    n=len(ll)
    ans=[]
    ans.append(0)
    lis=[0 for i in range(n+1)]
    for i in range(1,n+1):
        lis[i]=lis[i-1]+ll[i-1]
        ans.append(lis[i])
    for k in range(1,n):
        for b in range(0,n):
            if b==0:
                sum=lis[k]
                continue
            sum+=ll[(b+k-1)%n]
            sum-=ll[b-1]
            ans.append(sum)
    ans.sort()
    return ans

def linear_spectrum_of_peptide_weights(ll):
    n=len(ll)
    ans=[]
    ans.append(0)
    lis=[0 for i in range(n+1)]
    for i in range(1,n+1):
        lis[i]=lis[i-1]+ll[i-1]
        ans.append(lis[i])
    for k in range(1,n):
        for b in range(0,n):
            if b==0:
                sum=lis[k]
                continue
            if b+k-1>=n: break
            sum+=ll[b+k-1]
            sum-=ll[b-1]
            ans.append(sum)
    ans.sort()
    return ans

def expand_peptides(ll):
    weights = amino_weights
    ans=[]
    for new in weights:
        for lis in ll:
            ans.append(lis+[new])
    ans.sort()
    return ans

def cyclopeptide_sequencing(lis):
    # handling duplicates might be wrong
    pep=[]
    pep.append([])
    ans=[]
    lis.sort()
    s=set(lis)
    while len(pep):
        pepp = expand_peptides(pep)
        pep=[]
        # print(pepp)
        for p in pepp:
            pep.append(p)
            if set(cycle_spectrum_of_peptide_weights(p))==s and p not in ans:
                ans.append(p)
                pep.remove(p)
            else:
                f=0
                for x in linear_spectrum_of_peptide_weights(p):
                    if x not in s:
                        f=1
                if f==1: pep.remove(p)
    return ans

def peptide_spectrum_score(pep,spectrum):
    # ll = linear_spectrum_of_peptide(pep)
    ll=cycle_spectrum_of_peptide(pep)
    m1={}
    m2={}
    for x in ll:
        if x not in m1: m1[x]=0
        m1[x]+=1
    for x in spectrum:
        if x not in m2: m2[x]=0
        m2[x]+=1
    score=0
    maxx=max(spectrum)
    for k,v in m1.items():
        if k>maxx: return 0
        if k in m2:
            score+=min(v,m2[k])
    return score

def peptide_spectrum_weights_score(pep,spectrum):
    ll = cycle_spectrum_of_peptide_weights(pep)  
    # ll = linear_spectrum_of_peptide_weights(pep)   
    m1={}
    m2={}
    for x in ll:
        if x not in m1: m1[x]=0
        m1[x]+=1
    for x in spectrum:
        if x not in m2: m2[x]=0
        m2[x]+=1
    score=0
    maxx=max(spectrum)
    for k,v in m1.items():
        if k>maxx: return 0
        if k in m2:
            score+=min(v,m2[k])
    return score

def trim_leaderboard(leaderboard,spectrum,N):
    lb=leaderboard
    if len(lb)==0: return lb
    l=[]
    for p in lb:
        l.append((peptide_spectrum_weights_score(p,spectrum),p))
    l.sort()
    l.reverse()
    lb=[]
    s=l[0][0]
    c=0
    debug=0
    for i in range(len(l)):
        if i<N: lb.append(l[i][1])
        if i==N-1:
            s=l[i][0]
            i+=1
            while l[i][0]==s:
                debug+=1
                lb.append(l[i][1])
                i+=1
                if i>=len(l): break
    return lb

# from sortedcontainers import SortedSet

def trim_lp(lp,M):
    n=len(lp)
    if n<=M: return lp
    p=-1
    for i in range(n-M):
        if lp[i]!=lp[i+1]: p=i
    if p<0: return lp
    for i in range(p+1):
        lp.pop(0)
    return lp

def leaderboard_cyclopeptide_sequencing(spectrum, N):
    lb = [[]]
    lp=[]
    # lp=SortedSet()
    sc=0
    f=0
    mass=max(spectrum)
    while len(lb):
        scores=[]
        lb=expand_peptides(lb)
        for p in lb:
            score=peptide_spectrum_weights_score(p,spectrum)
            lmao = linear_spectrum_of_peptide_weights(p)
            if max(lmao)==mass:
                if score>sc:
                # if score>=sc:
                    lp=p
                    sc=score
                    # lp.append((score,p))
                    # lp.add((score,tuple(p)))
                    # trim_lp(lp,86)
            elif max(lmao)>mass:
                score=0
            scores.append(score)
        if len(lb)>N:
            cut_off=sorted(scores)[-N]
            leaders = []
            for i in range(len(lb)):
                    if scores[i] >= cut_off and scores[i] > 0:
                        leaders.append(lb[i])
            lb = leaders 
    return lp

def spectral_convolution(spectrum):
    s=spectrum
    s.sort()
    n=len(s)
    ans=[]
    for i in range(n-1):
        for j in range(i+1,n):
            x = s[j]-s[i]
            if x!=0: ans.append(x)
    return ans

def convolution_cyclopeptide_sequencing(spectrum,M,N):
    conv = spectral_convolution(spectrum)
    m={}
    for x in conv:
        if x<57 or x>200: continue
        if x not in m: m[x]=0
        m[x]+=1
    pairs=[]
    for x,f in m.items():
        pairs.append((f,x))
    pairs.sort()
    # print(pairs)
    cut_off = pairs[-M][0]
    l=len(pairs)
    global amino_weights
    temp=amino_weights.copy()
    amino_weights=set()
    for i in range(l-1,-1,-1):
        if pairs[i][0]>=cut_off: amino_weights.add(pairs[i][1])
    print(amino_weights)
    ans = leaderboard_cyclopeptide_sequencing(spectrum,N)
    amino_weights=temp.copy()
    return ans
    
