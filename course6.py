def trie_construction(patterns):
    patterns=list(set(patterns))
    n=1
    g=[{}]
    lis=[]
    for s in patterns:
        node=0
        for c in s:
            if c not in g[node]:
                g.append({})
                g[node][c]=n
                node=n
                n+=1
            else:
                node=g[node][c]
        g[node]['$']=-1
    for i in range(n):
        for c,j in g[i].items():
            if j!=-1: lis.append((i,j,c))
    return g,lis

def trie_find(g,s):
    node=0
    for c in s:
        if c not in g[node]: return False
        node=g[node][c]
    if '$' in g[node]: return True
    else: return False

def trie_matching(genome,patterns):
    g,lis=trie_construction(patterns)
    print(g)
    max=0
    m={}
    for s in patterns:
        if len(s)>max: max=len(s)
    n=len(genome)
    for i in range(n):
        for j in range(i+1,min(n,i+max)+1):
            s=genome[i:j]
            if trie_find(g,s):
                if s not in m: m[s]=[]
                m[s].append(i)
    return m
    
def suffix_tree_construction(genome):
    n=1
    g=[{}]
    p=[{}]
    nn=[{}]
    d={}
    for i in range(len(genome)):
        s=genome[i:]
        node=0
        for j in range(len(s)):
            c=s[j]
            if c not in g[node]:
                g.append({})
                g[node][c]=n
                nn.append({})
                nn[node][c]=1
                if c=='$': d[n]=i
                p.append({})
                p[node][c]=i+j
                node=n
                n+=1
            else:
                node=g[node][c]
    for i in range(n-1,-1,-1):
        if len(g[i])==0: continue
        for c,j in g[i].items():
            if len(g[j])!=1: continue
            cc=list(g[j].keys())[0]
            k=list(g[j].values())[0]
            g[i][c]=k
            nn[i][c]+=nn[j][cc]
            g[j]={}
    edges=[]
    e={}
    ans=0
    for i in range(n):
        for c,j in g[i].items():
            edges.append(genome[p[i][c]:p[i][c]+nn[i][c]])
            if i not in e: e[i]={}
            e[i][j]=edges[-1]
            if len(g[j])==0: ans+=1
    # return ans
    # return longest_repeat(0,g,e)
    return edges

def longest_repeat(v,g,e):
    max=0
    s=''
    for j in g[v].values():
        if len(g[j])==0: continue
        ss=longest_repeat(j,g,e)
        if len(e[v][j])+len(ss)>max:
            s=e[v][j]+ss
            max=len(s)
    return s

def suffix_array(genome):
    l=[]
    n=len(genome)
    for i in range(n):
        l.append((genome[i:],i))
    l.sort()
    return l

def burrows_wheeler_contruction(genome):
    l=[]
    n=len(genome)
    for i in range(n):
        l.append((genome[i:]+genome[:i],i))
    l.sort()
    s=''
    pos=[]
    for ss in l:
        s+=ss[0][-1]
        pos.append(ss[1])
    return s,pos

def inverse_burrows_wheeler_transform(bwt):
    s=bwt
    n=len(s)
    m={}
    l=[]
    for c in s:
        if c not in m: m[c]=0
        m[c]+=1
        l.append((c,m[c]))
    ll=sorted(l)
    rm={}
    for i in range(n):
        rm[l[i]]=ll[i]
    ans=''
    t=('$',1)
    f=0
    while t!=('$',1) or f==0:
        f=1
        t=rm[t]
        ans+=t[0]
    return ans

def burrows_wheeler_matching(bwt,pattern):
    s=bwt
    n=len(s)
    m={}
    l=[]
    for c in s:
        if c not in m: m[c]=0
        m[c]+=1
        l.append((c,m[c]))
    ll=sorted(l)
    rm={}
    for i in range(n):
        rm[l[i]]=ll[i]
    pn=len(pattern)
    if pattern[0] not in m: return 0
    lis=[i for i in range(1,m[pattern[0]]+1)]
    for i in range(pn-1):
        c1=pattern[i]
        c2=pattern[i+1]
        temp=[]
        for j in lis:
            t=rm[(c1,j)]
            if t[0]!=c2: continue
            temp.append(t[1])
        lis=temp.copy()
    return len(lis)

def multiple_pattern_matching(genome,patterns):
    bwt,pos=burrows_wheeler_contruction(genome)
    s=bwt
    n=len(s)
    m={}
    l=[]
    rmm={}
    index=0
    for c in s:
        if c not in m: m[c]=0
        m[c]+=1
        l.append((c,m[c]))
        rmm[(c,m[c])]=pos[index]
        index+=1
    ll=sorted(l)
    rm={}
    for i in range(n):
        rm[l[i]]=ll[i]
    ans=[]
    for pattern in patterns:
        pn=len(pattern)
        if pattern[0] not in m:
            ans.append([])
            continue
        lis=[i for i in range(1,m[pattern[0]]+1)]
        for i in range(pn-1):
            c1=pattern[i]
            c2=pattern[i+1]
            temp=[]
            for j in lis:
                t=rm[(c1,j)]
                if t[0]!=c2: continue
                if i==pn-2: temp.append((rmm[t]-pn+n)%n)
                else: temp.append(t[1])
            lis=temp.copy()
        lis.sort()
        ans.append(lis)
    return ans

def hidden_path_probability(path,transition_matrix):
    m=transition_matrix
    p=1/(len(m))
    n=len(path)
    for i in range(1,n):
        p*=m[path[i-1]][path[i]]
    return p

def outcome_given_hidden_path_probability(path,outcome,emission_matrix):
    p=1
    n=len(path)
    m=emission_matrix
    for i in range(n):
        p*=m[path[i]][outcome[i]]
    return p

def outcome_to_hidden_path(outcome,transmission_matrix,emission_matrix):
    tm=transmission_matrix
    em=emission_matrix
    n=len(outcome)
    states=list(tm.keys())
    dp=[{state:(0,'') for state in states} for i in range(n)]
    for state in states:
        dp[-1][state]=(em[state][outcome[-1]],'')
    for i in range(n-1,0,-1):
        for state1 in states:
            for state2 in states:
                x=dp[i][state1][0]*tm[state2][state1]*em[state2][outcome[i-1]]
                if x>dp[i-1][state2][0]:
                    dp[i-1][state2]=(x,state1)
    start=''
    next=''
    max=0
    for state in states:
        if max<dp[0][state][0]:
            start=state
            max=dp[0][state][0]
    ans=''
    for i in range(n-1):
        ans+=start
        next=dp[i][start][1]
        start=next
    ans+=start
    return ans

def outcome_probability_for_all_hidden_paths(outcome,transmission_matrix,emission_matrix):
    tm=transmission_matrix
    em=emission_matrix
    n=len(outcome)
    states=list(tm.keys())
    dp=[{state:0 for state in states} for i in range(n)]
    for state in states:
        dp[-1][state]=em[state][outcome[-1]]
    for i in range(n-1,0,-1):
        for state1 in states:
            for state2 in states:
                x=dp[i][state1]*tm[state2][state1]*em[state2][outcome[i-1]]
                dp[i-1][state2]+=x
    ans=0
    for state in states:
        ans+=dp[0][state]
    ans/=len(states)
    return ans

