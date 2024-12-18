import copy

def limb_length(x,ll,k=0,nodes=[]):
    n=len(ll)
    if len(nodes)==0:
        nodes=[i for i in range(n)]
    n=len(nodes)
    if n==1:
        return 0
    if n==2:
        return ll[nodes[0]][nodes[1]]
    y=0
    for i in nodes:
        if x!=i:
            y=i
            break
    m=0
    f=0
    p1=0
    for i in nodes:
        if i==x or i==y: continue
        s=(ll[i][x]+ll[y][x]-ll[i][y])/2
        if s<m or f==0:
            m=s
            f=1
            p1=i
    if k==1: return int(m),p1,y
    return int(m)

def reverse_dfs_undirected_path(x,y,g,v):
    if x not in v: v[x]=0
    if len(g[x])==0: return 0,[]
    if x==y: return 1,[y]
    for z in g[x]:
        z=z[0]
        if z in v: continue
        k,path = reverse_dfs_undirected_path(z,y,g,v)
        if k==1:
            path.append(x)
            return 1,path
    return 0,[]

def additive_phylogeny(ll,nodes,g):
    if len(nodes)==2:
        g[nodes[0]].append((nodes[1],ll[nodes[0]][nodes[1]]))
        g[nodes[1]].append((nodes[0],ll[nodes[0]][nodes[1]]))
        return
    x=nodes[-1]
    l,y,z = limb_length(x,ll,1,nodes)
    nodes.pop()
    additive_phylogeny(ll,nodes,g)
    v={}
    k,path = reverse_dfs_undirected_path(y,z,g,v)
    path.reverse()
    le = ll[x][y]-l
    if le<0: print(x,y,ll[x][y],l)
    nn=len(g)
    n=len(path)
    s=0
    # print(g)
    for i in range(1,n):
        lol=ll[path[i]][path[i-1]]
        s+=lol
        if s>le:
            g[path[i-1]].remove((path[i],lol))
            g[path[i]].remove((path[i-1],lol))
            ll[path[i-1]][path[i]]=-1
            ll[path[i]][path[i-1]]=-1
            g.append([])
            g[nn].append((path[i],s-le))
            g[nn].append((path[i-1],lol+le-s))
            g[nn].append((x,l))
            g[path[i]].append((nn,s-le))
            g[path[i-1]].append((nn,lol+le-s))
            g[x].append((nn,l))
            ll[nn][path[i]]=s-le
            ll[nn][path[i-1]]=lol+le-s
            ll[nn][x]=l
            ll[path[i]][nn]=s-le
            ll[path[i-1]][nn]=lol+le-s
            ll[x][nn]=l
            break
        if s==le:
            print(x,path[i],g,s)
            if i==n-1:
                g[path[i-1]].append((x,l))
                g[x].append((path[i-1],l))
                ll[path[i-1]][x]=l
                ll[x][path[i-1]]=l
            else:
                g[path[i]].append((x,l))
                g[x].append((path[i],l))
                ll[path[i]][x]=l
                ll[x][path[i]]=l
    # nodes.append(x)
    # nodes.append(nn)
    return

def cluster_distance_combining(c1,d1,c2,d2):
    n1=len(c1)
    n2=len(c2)
    ans = d1*n1+d2*n2
    ans/=(n1+n2)
    return ans

def upgma_cluster_tree(ll,n):
    g = [[] for i in range(2*n)]
    clusters=[[i] for i in range(n)]
    indices=set([i for i in range(n)])
    d=[0 for i in range(n)]
    while(len(indices)>1):
        nn=len(clusters)
        m,f,i1,i2=0,0,0,0
        for i in indices:
            for j in indices:
                if i>=j: continue
                if f==0 or m>ll[i][j]:
                    m=ll[i][j]
                    f=1
                    i1=i
                    i2=j
        c1=clusters[i1]
        c2=clusters[i2]
        c=c1+c2
        clusters.append(c)
        indices.remove(i1)
        indices.remove(i2)
        for i in indices:
            ll[nn][i]=cluster_distance_combining(c1,ll[i1][i],c2,ll[i2][i])
            ll[i][nn]=ll[nn][i]
        indices.add(nn)
        d.append(ll[i1][i2]/2)
        g[nn].append((i1,d[nn]-d[i1]))
        g[nn].append((i2,d[nn]-d[i2]))
        g[i1].append((nn,d[nn]-d[i1]))
        g[i2].append((nn,d[nn]-d[i2]))
    return g

def neighbor_joining(ll,n):
    g = [[] for i in range(2*n)]
    clusters=[[i] for i in range(n)]
    indices=set([i for i in range(n)])
    d=[0 for i in range(n)]
    for i in range(n):
        lol=0
        for j in range(n):
            lol+=ll[i][j]
        d[i]=lol
    while(len(indices)>1):
        if len(indices)==2:
            i1=0
            i2=0
            for i in indices:
                for j in indices:
                    if i!=j:
                        i1=i
                        i2=j
                        break
            g[i1].append((i2,ll[i1][i2]))
            g[i2].append((i1,ll[i1][i2]))
            break
        nn=len(clusters)
        m,f,i1,i2=0,0,0,0
        for i in indices:
            for j in indices:
                if i>=j: continue
                lol = (len(indices)-2)*ll[i][j]-d[i]-d[j]
                if f==0 or m>lol:
                    m=lol
                    f=1
                    i1=i
                    i2=j
        c1=clusters[i1]
        c2=clusters[i2]
        c=c1+c2
        clusters.append(c)
        indices.remove(i1)
        indices.remove(i2)
        d.append(0)
        for i in indices:
            ll[nn][i]=(1/2)*(ll[i][i1]+ll[i][i2]-ll[i1][i2])
            ll[i][nn]=ll[nn][i]
            d[i]+=(ll[nn][i]-ll[i][i1]-ll[i][i2])
            d[nn]+=ll[nn][i]
        indices.add(nn)
        delta = (d[i1]-d[i2])/(len(indices)-1)
        li = (1/2)*(ll[i1][i2]+delta)
        lj = (1/2)*(ll[i1][i2]-delta)
        g[nn].append((i1,li))
        g[nn].append((i2,lj))
        g[i1].append((nn,li))
        g[i2].append((nn,lj))
    return g

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

def nearest_neighbor_trees(s,x,y):
    s1,s2,s3,s4,s5,s6='','','','','',''
    c,d,e='','',''
    # ignoring first line
    for i in range(1,len(s)):
        a,b=s[i].split('->')
        if a==str(x) and c=='' and b!=str(y):
            c=b
        if a==str(y) and d=='' and b!=str(x):
            d=b
        elif a==str(y) and e=='' and d!='' and b!=str(x):
            e=b
    s1=str(x)+'->'+c
    rs1=c+'->'+str(x)
    s2=str(y)+'->'+d
    rs2=d+'->'+str(y)
    s5=str(y)+'->'+e
    rs5=e+'->'+str(y)
    s3=str(x)+'->'+d
    rs3=d+'->'+str(x)
    s6=str(x)+'->'+e
    rs6=e+'->'+str(x)
    s4=str(y)+'->'+c
    rs4=c+'->'+str(y)
    sa=s.copy()
    sa.remove(s1)
    sa.remove(rs1)
    sa.remove(s2)
    sa.remove(rs2)
    sa.append(s3)
    sa.append(rs3)
    sa.append(s4)
    sa.append(rs4)
    sb=s.copy()
    sb.remove(s1)
    sb.remove(rs1)
    sb.remove(s5)
    sb.remove(rs5)
    sb.append(s6)
    sb.append(rs6)
    sb.append(s4)
    sb.append(rs4)
    return sa,sb

def small_parsimony(g,v,dp,p):
    if dp[v]['A'][0]!=-1:
        if dp[v]['G'][0]==-1 or dp[v]['T'][0]==-1 or dp[v]['C'][0]==-1:
            print("eror")
            quit()
        return dp[v]
    for c in 'ATGC':
        x,y,f=0,0,0
        for i in g[v]:
            if i!=p and f==0:
                f=1
                x=i
            elif i!=p and f==1:
                y=i
                break
        dx=small_parsimony(g,x,dp,v)
        mx=0
        f=0
        ch1=''
        for cc in 'ATGC':
            if mx>dx[cc][0]+1-int(cc==c) or f==0:
                mx=dx[cc][0]+1-int(cc==c)
                ch1=cc
                f=1
        dy=small_parsimony(g,y,dp,v)
        my=0
        f=0
        ch2=''
        for cc in 'ATGC':
            if my>dy[cc][0]+1-int(cc==c) or f==0:
                my=dy[cc][0]+1-int(cc==c)
                ch2=cc
                f=1
        dp[v][c]=[mx+my,ch1,ch2]
    return dp[v]

def smol_dfs(g,dp,v,cc,ans,p):
    ans[v]+=cc
    if len(g[v])<2: return
    x,y,f=0,0,0
    for i in g[v]:
        if i!=p and f==0:
            f=1
            x=i
        elif i!=p and f==1:
            y=i
            break
    smol_dfs(g,dp,x,dp[v][cc][1],ans,v)
    smol_dfs(g,dp,y,dp[v][cc][2],ans,v)

def hamming_distance(p,q):
    h=0
    for i in range(len(p)):
        if p[i]!=q[i]: h+=1
    return h

def small_parsimony_score_path(liss):
    n=liss[0]
    n=int(n)
    lis=liss[1:]
    m=[]
    ll=[]
    ma=0
    for s in lis:
        a,b = s.split('->')
        if a[0]=='A' or a[0]=='G' or a[0]=='C' or a[0]=='T': continue
        a=int(a)
        if b[0]=='A' or b[0]=='G' or b[0]=='C' or b[0]=='T':
            ll.append((a,len(m)))
            m.append(b)
        else:
            b=int(b)
            if a>b: continue
            ll.append((a,b))
            if ma<b: ma=b
        if ma<a: ma=a
    ma+=1
    g = [[] for i in range(ma+1)]
    for t in ll:
        g[t[0]].append(t[1])
        g[t[1]].append(t[0])
    edges=set()
    f=0
    for t in ll:
        edges.add((t[0],t[1]))
        edges.add((t[1],t[0]))
        if f==0:
            f=1
            g[t[0]].remove(t[1])
            g[t[1]].remove(t[0])
            g[ma]=[t[0],t[1]]
    ls=len(m[0])
    ans=['' for i in range(ma+1)]
    erm=0
    for j in range(ls):
        dp = [{'A':[-1,'',''],'T':[-1,'',''],'G':[-1,'',''],'C':[-1,'','']} for i in range(ma+1)]
        for i in range(n):
            dp[i]['A']=[n,'','']
            dp[i]['T']=[n,'','']
            dp[i]['G']=[n,'','']
            dp[i]['C']=[n,'','']
        for i in range(n):
            dp[i][m[i][j]]=[0,'','']
        a = small_parsimony(g,ma,dp,-1)
        minn = min(dp[ma]['A'][0],dp[ma]['T'][0],dp[ma]['G'][0],dp[ma]['C'][0])
        erm+=minn
        cc=''
        for c in 'ATGC':
            if dp[ma][c][0]==minn:
                cc=c
                break
        smol_dfs(g,dp,ma,cc,ans,-1)
    output = []
    for edge in edges:
        x=edge[0]
        y=edge[1]
        output.append(f"{ans[x]}->{ans[y]}:{hamming_distance(ans[x],ans[y])}")
    return erm,output

def large_parsimony(lis):
    score,output = small_parsimony_score_path(lis)
    new_lis,new_score,new_output=lis.copy(),score,output.copy()
    f,ff=0,0
    while new_score<score or f==0:
        if f!=0:
            with open("output.txt", "a") as file:
                if ff==1: file.write("\n\n")
                file.write(f"{new_score}")
                for s in new_output: file.write(f"\n{s}")
            ff=1
        f=1
        lis=new_lis.copy()
        score=new_score
        output=new_output.copy()
        for i in range(1,len(lis)):
            a,b = lis[i].split('->')
            if a[0]=='A' or a[0]=='G' or a[0]=='C' or a[0]=='T': continue
            if b[0]=='A' or b[0]=='G' or b[0]=='C' or b[0]=='T': continue
            n1_lis,n2_lis=nearest_neighbor_trees(lis,a,b)
            n1_score,n1_output=small_parsimony_score_path(n1_lis)
            n2_score,n2_output=small_parsimony_score_path(n2_lis)
            if n1_score<new_score:
                new_lis=n1_lis.copy()
                new_score=n1_score
                new_output=n1_output.copy()
            if n2_score<new_score:
                new_lis=n2_lis.copy()
                new_score=n2_score
                new_output=n2_output.copy()
    return n2_lis


amino_acid_mass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

amino_weights = {57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186}

rm = {
    57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: ['L','I'], 114: 'N', 115: 'D', 128: ['K','Q'], 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'
}
imaginary_amino_acid_by_mass = {4:"X", 5:"Z", 57:"G", 71:"A", 87:"S", 97:"P", 99:"V", 101:"T", 103:"C", 113:"L", 114:"N", 115:"D", 128:"Q", 129:"E", 131:"M", 137:"H", 147:"F", 156:"R", 163:"Y", 186:"W"}

imaginary_mass_by_amino_acid = {"X":4, "Z":5, "G":57, "A":71, "S":87, "P":97, "V":99, "T":101, "C":103, "L":113, "I":113, "N":114, "D":115, "Q":128, "K":128, "E":129, "M":131, "H":137, "F":147, "R":156, "Y":163, "W":186, "I":113, "K":128}

def spectrum_dfs_path(g,x,y,dp):
    if dp[x]!=['lol']: return dp[x]
    dp[x]=[]
    if x not in g: return []
    for (c,a) in g[x]:
        l=spectrum_dfs_path(g,c,y,dp)
        for s in l:
            s=a+s
            dp[x].append(s)
    return dp[x]

def linear_prefix_suffix_spectrum(s):
    ans=[0]
    n=len(s)
    for i in range(n):
        p=ans[-1]
        ans.append(p+amino_acid_mass[s[i]])
    for i in range(n-1):
        p=ans[-1]
        ans.append(p-amino_acid_mass[s[i]])
    ans.sort()
    return ans

def decoding_ideal_spectrum(l):
    g={}
    s=l
    for x in s:
        for w in amino_weights:
            if x+w in s:
                if isinstance(rm[w],str):
                    if x not in g: g[x]=[]
                    g[x].append((x+w,rm[w]))
                else:
                    if x not in g: g[x]=[]
                    g[x].append((x+w,rm[w][0]))
                    # g[x].append((x+w,rm[w][1]))
    dp={}
    for k in g:
        dp[k]=['lol']
    dp[s[-1]]=['']
    ll = spectrum_dfs_path(g,0,s[-1],dp)
    ans=[]
    for ss in ll:
        if linear_prefix_suffix_spectrum(ss) == s: ans.append(ss)
    return ans

def peptide_to_peptide_vector(s):
    n=len(s)
    l=[]
    temp=0
    for x in s:
        x=imaginary_mass_by_amino_acid[x]
        temp+=x
        l.append(temp)
    ans=[0 for i in range(l[-1])]
    for x in l: ans[x-1]=1
    return ans

def peptide_vector_to_peptide(l):
    m=[]
    for i in range(len(l)):
        if l[i]==1: m.append(i+1)
    s=imaginary_amino_acid_by_mass[m[0]]
    for i in range(1,len(m)):
        s+=imaginary_amino_acid_by_mass[m[i]-m[i-1]]
    return s

def peptide_sequencing_dfs(g,x,dp,sv):
    n=len(g)
    if x>=n: return []
    if dp[x]!=['lol',0]: return dp[x]
    scm=0
    f=0
    dps=''
    for (c,a) in g[x]:
        l=peptide_sequencing_dfs(g,c,dp,sv)
        if l==[]: continue
        s=l[0]
        sc=l[1]
        if sc+sv[c]>scm or f==0:
            f=1
            dps=a+s
            scm=sc+sv[c]
    if f==1: dp[x]=[dps,scm]
    else: dp[x]=[]
    return dp[x]
    
def peptide_sequencing(sv):
    sv.insert(0,0)
    n=len(sv)
    g=[[] for i in range(n)]
    for i in range(n):
        for w,a in imaginary_amino_acid_by_mass.items():
            if i+w<n:
                g[i].append((i+w,a))
    dp=[['lol',0] for i in range(n)]
    dp[n-1]=['',0]
    ll = peptide_sequencing_dfs(g,0,dp,sv)
    return ll

def peptide_score(pv,sv):
    n=len(pv)
    ans=0
    for i in range(n):
        ans+=pv[i]*sv[i]
    return ans

def best_peptide_in_proteome(p,sv):
    n=len(p)
    w=len(sv)
    l=[]
    for i in range(n):
        temp=imaginary_mass_by_amino_acid[p[i]]
        j=i
        for jj in range(i+1,n):
            if temp>=w: break
            temp+=imaginary_mass_by_amino_acid[p[jj]]
            j=jj
        if temp!=w: continue
        l.append(p[i:j+1])
    if len(l)==0: return '',0
    ll=[]
    for i in range(len(l)):
        ll.append(peptide_score(peptide_to_peptide_vector(l[i]),sv))
    ma=max(ll)
    ans=[]
    for i in range(len(l)):
        if ll[i]==ma:
            ans.append(l[i])
    # return ans
    return ans[0],ma

def psm_search(p,l,th):
    ans=set()
    for sv in l:
        x,y=best_peptide_in_proteome(p,sv)
        if x=='': continue
        if y>=th: ans.add(x)
    return ans

def spectral_dictionary_size(v,th,ma):
    v.insert(0,0)
    shift=0
    p=0
    for x in v:
        if x<0: shift-=x
        else: p+=x
    dp=[[0 for i in range(p+shift+1)] for j in range(len(v))]
    weights = [57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]
    dp[0][shift]=1
    for i in range(1,len(v)):
        for s in range(p+shift+1):
            for w in weights:
                # if i-w<0 or s-v[i]<0 or s-v[i]>p+shift: continue
                if i-w<0 or s-v[i]<shift or s-v[i]>ma+shift: continue
                dp[i][s]+=dp[i-w][s-v[i]]
    ans=0
    for i in range(shift+th,shift+ma+1): ans+=dp[len(v)-1][i]
    return ans

def spectral_dictionary_probability(v,th,ma):
    v.insert(0,0)
    shift=0
    p=0
    for x in v:
        if x<0: shift-=x
        else: p+=x
    dp=[[0 for i in range(p+shift+1)] for j in range(len(v))]
    weights = [57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]
    dp[0][shift]=1
    for i in range(1,len(v)):
        for s in range(p+shift+1):
            for w in weights:
                if i-w<0 or s-v[i]<0 or s-v[i]>p+shift: continue
                # if i-w<0 or s-v[i]<shift or s-v[i]>ma+shift: continue
                dp[i][s]+=dp[i-w][s-v[i]]*(0.05)
    ans=0
    for i in range(shift+th,shift+ma+1): ans+=dp[len(v)-1][i]
    return ans

def spectral_alignment_with_k_mutations(s,v,k):
    p=[]
    for a in s: p.append(imaginary_mass_by_amino_acid[a])
    v.insert(0,0)
    n=len(v)
    dp=[[('lol',0,0) for i in range(k+1)] for j in range(n)]
    dp[n-1][0]=('',0,len(p))
    for ia in range(len(p)-1,-1,-1):
        w=p[ia]
        a=s[ia]
        for i in range(n):
            for j in range(k+1):
                if dp[i][j][0]=='lol' or dp[i][j][2]!=ia+1: continue
                ss=a+dp[i][j][0]
                sc=v[i]+dp[i][j][1]
                if i-w>=0:  
                    if sc>dp[i-w][j][1] or dp[i-w][j][0]=='lol':
                        dp[i-w][j]=(ss,sc,ia)
                if j==k: continue
                for q in range(n):
                    if q==i-w: continue
                    if sc>dp[q][j+1][1] or dp[q][j+1][0]=='lol':
                        x = i-w-q
                        if x<0: ss=f"{a}({x}){dp[i][j][0]}"
                        else: ss=f"{a}(+{x}){dp[i][j][0]}"
                        dp[q][j+1]=(ss,sc,ia)
    ans=''
    max=0
    f=0
    for j in range(k+1):
        if dp[0][j][0]=='lol': continue
        if max<dp[0][j][1] or f==0:
            f=1
            max=dp[0][j][1]
            ans=dp[0][j][0]
    return ans
