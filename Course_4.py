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

