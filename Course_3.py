import sys
sys.setrecursionlimit(10000)

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

def dp_change_problem(sum, coins):
    n=len(coins)
    # l=[(0,0),]
    l=[-1]*(sum+1)
    l[0]=0
    for i in range(sum+1):
        for coin in coins:
            if coin+i<=sum:
                if l[i]!=-1: 
                    if l[coin+i]>l[i]+1 or l[coin+i]==-1:
                        l[coin+i]=l[i]+1
    return l[sum]

def manhattan_tourist_problem(g,d,x,y):
    # n=4
    # m=14
    # g=[[[0,0] for i in range(m+1)] for i in range(n+1)]
    # down = ''.split()
    # for i in range(len(down)): down[i]=int(down[i])
    # x=-1
    # for i in range(len(down)):
    #     y=i%(m+1)
    #     if y==0: x+=1
    #     g[x][y][1]=down[i]
    # right=''.split()
    # for i in range(len(right)): right[i]=int(right[i])
    # x=-1
    # for i in range(len(right)):
    #     y=i%(m)
    #     if y==0: x+=1
    #     g[x][y][0]=right[i]
    # a=max(max(down),max(right))*(-1)
    # ll=[-1 for i in range(m+2)]
    # ll[m+1]=a
    # d=[[-1 for i in range(m+2)] for i in range(n+2)]
    # for i in range(n+2): d[i]=ll.copy()
    # d[n+1]=[a for i in range(m+2)]
    # d[n][m]=0
    # print(manhattan_tourist_problem(g,d,0,0))
    if d[x][y]!=-1: return d[x][y]
    d[x][y] = max(manhattan_tourist_problem(g,d,x+1,y)+g[x][y][1],manhattan_tourist_problem(g,d,x,y+1)+g[x][y][0])
    return d[x][y]

def longest_common_subsequence(s1,s2,i1,i2,dp):
    # dp = [[[-1,''] for i in range(n2+1)] for i in range(n1+1)]
    if dp[i1][i2][0]!=-1: return dp[i1][i2]
    n1=len(s1)
    n2=len(s2)
    if i1>=n1 or i2>=n2:
        dp[i1][i2][0]=0
        dp[i1][i2][1]=''
        return dp[i1][i2]
    a=0
    sa=''
    if s1[i1]==s2[i2]:
        a,sa=longest_common_subsequence(s1,s2,i1+1,i2+1,dp)[0]+1,longest_common_subsequence(s1,s2,i1+1,i2+1,dp)[1]
    b,sb=longest_common_subsequence(s1,s2,i1+1,i2,dp)[0],longest_common_subsequence(s1,s2,i1+1,i2,dp)[1]
    c,sc=longest_common_subsequence(s1,s2,i1,i2+1,dp)[0],longest_common_subsequence(s1,s2,i1,i2+1,dp)[1]
    d=max(a,b,c)
    if d==a:
        dp[i1][i2]=[a,s1[i1]+sa]
    elif d==b:
        dp[i1][i2]=[b,sb]
    else:
        dp[i1][i2]=[c,sc]
    return dp[i1][i2]

def longest_path_in_DAG(g,index,end,dp,cycle):
    if dp[index][0]!=-1: return dp[index]
    cycle[index]=1
    n=len(g)
    max=0
    l=[]
    for p in g[index]:
        if cycle[p[0]]==1:
            print("cycle",p[0])
        pp=longest_path_in_DAG(g,p[0],end,dp,cycle)
        if p[1]+pp[0]>max:
            max=p[1]+pp[0]
            l=pp[1].copy()
            l.reverse()
            l.append(index)
            l.reverse()
    dp[index]=[max,l]
    cycle[index]=0
    return dp[index]

def global_alignment_problem(s1,s2,i1,i2,match,mismatch,indel,dp):
    # bug might be present did not test with big cases
    n1=len(s1)
    n2=len(s2)
    if dp[i1][i2][0]!=0.5: return dp[i1][i2]
    print(i1,i2)
    if i1>=n1 and i2>=n2:
        print("1")
        dp[i1][i2][0]=0
        dp[i1][i2][1]=''
        dp[i1][i2][2]=''
        return dp[i1][i2]
    if i1>=n1:
        print("2")
        dp[i1][i2][0]=(n2-i2)*indel*(-1)
        dp[i1][i2][2]=s2[i2:]
        dp[i1][i2][1]='-'*(n2-i2)
        return dp[i1][i2]
    if i2>=n2:
        print("3")
        dp[i1][i2][0]=(n1-i1)*indel*(-1)
        dp[i1][i2][1]=s1[i1:]
        dp[i1][i2][2]='-'*(n1-i1)
        return dp[i1][i2]
    a=0
    sa1=''
    sa2=''
    if s1[i1]==s2[i2]:
        la=global_alignment_problem(s1,s2,i1+1,i2+1,match,mismatch,indel,dp)
        a=la[0]+match
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    else:
        la=global_alignment_problem(s1,s2,i1+1,i2+1,match,mismatch,indel,dp)
        a=la[0]-mismatch
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    lb=global_alignment_problem(s1,s2,i1+1,i2,match,mismatch,indel,dp)
    b=lb[0]-indel
    sb1=s1[i1]+lb[1]
    sb2='-'+lb[2]
    lc=global_alignment_problem(s1,s2,i1,i2+1,match,mismatch,indel,dp)
    c=lc[0]-indel
    sc1='-'+lc[1]
    sc2=s2[i2]+lc[2]
    d=max(a,b,c)
    if d==a:
        dp[i1][i2]=[a,sa1,sa2]
    elif d==b:
        dp[i1][i2]=[b,sb1,sb2]
    else:
        dp[i1][i2]=[c,sc1,sc2]
    return dp[i1][i2]

def global_alignment_problem_with_scoring_matrix(s1,s2,i1,i2,score,dp):
    # bug might be present did not test with big cases
    m = {'A':0,'G':1,'C':2,'T':3,'-':4}
    n1=len(s1)
    n2=len(s2)
    if dp[i1][i2][0]!='inf': return dp[i1][i2]
    if i1>=n1 and i2>=n2:
        dp[i1][i2][0]=0
        dp[i1][i2][1]=''
        dp[i1][i2][2]=''
        return dp[i1][i2]
    if i1>=n1:
        # dp[i1][i2][0]=(n2-i2)*indel*(-1)
        x=0
        for i in range(i2,n2):
            x+=score[m['-']][m[s2[i]]]
        dp[i1][i2][0]=x
        dp[i1][i2][2]=s2[i2:]
        dp[i1][i2][1]='-'*(n2-i2)
        return dp[i1][i2]
    if i2>=n2:
        # dp[i1][i2][0]=(n1-i1)*indel*(-1)
        x=0
        for i in range(i1,n1):
            x+=score[m[s1[i]]][m['-']]
        dp[i1][i2][0]=x
        dp[i1][i2][1]=s1[i1:]
        dp[i1][i2][2]='-'*(n1-i1)
        return dp[i1][i2]
    a=0
    sa1=''
    sa2=''
    if s1[i1]==s2[i2]:
        la=global_alignment_problem_with_scoring_matrix(s1,s2,i1+1,i2+1,score,dp)
        a=la[0]+score[m[s1[i1]]][m[s2[i2]]]
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    else:
        la=global_alignment_problem_with_scoring_matrix(s1,s2,i1+1,i2+1,score,dp)
        a=la[0]+score[m[s1[i1]]][m[s2[i2]]]
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    lb=global_alignment_problem_with_scoring_matrix(s1,s2,i1+1,i2,score,dp)
    b=lb[0]+score[m[s1[i1]]][m['-']]
    sb1=s1[i1]+lb[1]
    sb2='-'+lb[2]
    lc=global_alignment_problem_with_scoring_matrix(s1,s2,i1,i2+1,score,dp)
    c=lc[0]+score[m['-']][m[s2[i2]]]
    sc1='-'+lc[1]
    sc2=s2[i2]+lc[2]
    d=max(a,b,c)
    if d==a:
        dp[i1][i2]=[a,sa1,sa2]
    elif d==b:
        dp[i1][i2]=[b,sb1,sb2]
    else:
        dp[i1][i2]=[c,sc1,sc2]
    return dp[i1][i2]

def local_alignment_problem_with_scoring_matrix(s1,s2,i1,i2,score,indel,dp):
    n1=len(s1)
    n2=len(s2)
    if dp[i1][i2][0]!='inf': return dp[i1][i2]
    if i1>=n1 and i2>=n2:
        dp[i1][i2][0]=0
        dp[i1][i2][1]=''
        dp[i1][i2][2]=''
        return dp[i1][i2]
    if i1>=n1:
        dp[i1][i2][0]=(n2-i2)*indel*(-1)
        dp[i1][i2][2]=s2[i2:]
        dp[i1][i2][1]='-'*(n2-i2)
        return dp[i1][i2]
    if i2>=n2:
        dp[i1][i2][0]=(n1-i1)*indel*(-1)
        dp[i1][i2][1]=s1[i1:]
        dp[i1][i2][2]='-'*(n1-i1)
        return dp[i1][i2]
    a=0
    sa1=''
    sa2=''
    if s1[i1]==s2[i2]:
        la=local_alignment_problem_with_scoring_matrix(s1,s2,i1+1,i2+1,score,indel,dp)
        a=la[0]+score[s1[i1]][s2[i2]]
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    else:
        la=local_alignment_problem_with_scoring_matrix(s1,s2,i1+1,i2+1,score,indel,dp)
        a=la[0]+score[s1[i1]][s2[i2]]
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    lb=local_alignment_problem_with_scoring_matrix(s1,s2,i1+1,i2,score,indel,dp)
    b=lb[0]-indel
    sb1=s1[i1]+lb[1]
    sb2='-'+lb[2]
    lc=local_alignment_problem_with_scoring_matrix(s1,s2,i1,i2+1,score,indel,dp)
    c=lc[0]-indel
    sc1='-'+lc[1]
    sc2=s2[i2]+lc[2]
    e=0
    se1=''
    se2=''
    if i1==0 and i2==0:
        for i in range(n1):
            for j in range(n2):
                if i==0 and j==0: continue
                le=local_alignment_problem_with_scoring_matrix(s1,s2,i,j,score,indel,dp)
                ee=le[0]
                if ee>e:
                    e=ee
                    se1=le[1]
                    se2=le[2]
    else:
        e=score[s1[i1]][s2[i2]]
        se1=s1[i1]
        se2=s2[i2]
    d=max(a,b,c,e)
    if d==a:
        dp[i1][i2]=[a,sa1,sa2]
    elif d==b:
        dp[i1][i2]=[b,sb1,sb2]
    elif d==c:
        dp[i1][i2]=[c,sc1,sc2]
    else:
        dp[i1][i2]=[e,se1,se2]
    return dp[i1][i2]

def fitting_alignment_problem(s1,s2,i1,i2,score,indel,dp):
    n1=len(s1)
    n2=len(s2)
    if dp[i1][i2][0]!='inf': return dp[i1][i2]
    if i1>=n1 and i2>=n2:
        dp[i1][i2][0]=0
        dp[i1][i2][1]=''
        dp[i1][i2][2]=''
        return dp[i1][i2]
    if i1>=n1:
        dp[i1][i2][0]=(n2-i2)*indel*(-1)
        dp[i1][i2][2]=s2[i2:]
        dp[i1][i2][1]='-'*(n2-i2)
        return dp[i1][i2]
    if i2>=n2:
        dp[i1][i2][0]=(n1-i1)*indel*(-1)
        dp[i1][i2][1]=s1[i1:]
        dp[i1][i2][2]='-'*(n1-i1)
        return dp[i1][i2]
    a=0
    sa1=''
    sa2=''
    if s1[i1]==s2[i2]:
        la=fitting_alignment_problem(s1,s2,i1+1,i2+1,score,indel,dp)
        a=la[0]+score[s1[i1]][s2[i2]]
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    else:
        la=fitting_alignment_problem(s1,s2,i1+1,i2+1,score,indel,dp)
        a=la[0]+score[s1[i1]][s2[i2]]
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    lb=fitting_alignment_problem(s1,s2,i1+1,i2,score,indel,dp)
    b=lb[0]-indel
    sb1=s1[i1]+lb[1]
    sb2='-'+lb[2]
    lc=fitting_alignment_problem(s1,s2,i1,i2+1,score,indel,dp)
    c=lc[0]-indel
    sc1='-'+lc[1]
    sc2=s2[i2]+lc[2]
    e=0
    se1=''
    se2=''
    if i1==0 and i2==0:
        for i in range(n1):
            if i==0: continue
            le=fitting_alignment_problem(s1,s2,i,0,score,indel,dp)
            ee=le[0]
            if ee>e:
                e=ee
                se1=le[1]
                se2=le[2]
    elif i2==n2-1 and i1!=n1-1:
        e=score[s1[i1]][s2[i2]]
        se1=s1[i1]
        se2=s2[i2]
    else:
        d=max(a,b,c)
        if d==a:
            dp[i1][i2]=[a,sa1,sa2]
        elif d==b:
            dp[i1][i2]=[b,sb1,sb2]
        else:
            dp[i1][i2]=[c,sc1,sc2]
        return dp[i1][i2]
    d=max(a,b,c,e)
    if d==a:
        dp[i1][i2]=[a,sa1,sa2]
    elif d==b:
        dp[i1][i2]=[b,sb1,sb2]
    elif d==c:
        dp[i1][i2]=[c,sc1,sc2]
    else:
        dp[i1][i2]=[e,se1,se2]
    return dp[i1][i2]

def gap_penalties_alignemnt_problem(s1,s2,i1,i2,match,mismatch,sig,ep,k,dp):
    indel=sig
    if k==1: indel=ep
    n1=len(s1)
    n2=len(s2)
    if dp[i1][i2][0]!='inf': return dp[i1][i2]
    if i1>=n1 and i2>=n2:
        dp[i1][i2][0]=0
        dp[i1][i2][1]=''
        dp[i1][i2][2]=''
        return dp[i1][i2]
    if i1>=n1:
        dp[i1][i2][0]=(n2-i2-1)*k*(-1)-indel
        dp[i1][i2][2]=s2[i2:]
        dp[i1][i2][1]='-'*(n2-i2)
        return dp[i1][i2]
    if i2>=n2:
        dp[i1][i2][0]=(n1-i1-1)*indel*(-1)-indel
        dp[i1][i2][1]=s1[i1:]
        dp[i1][i2][2]='-'*(n1-i1)
        return dp[i1][i2]
    a=0
    sa1=''
    sa2=''
    if s1[i1]==s2[i2]:
        la=gap_penalties_alignemnt_problem(s1,s2,i1+1,i2+1,match,mismatch,sig,ep,0,dp)
        a=la[0]+match
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    else:
        la=gap_penalties_alignemnt_problem(s1,s2,i1+1,i2+1,match,mismatch,sig,ep,0,dp)
        a=la[0]-mismatch
        sa1=s1[i1]+la[1]
        sa2=s2[i2]+la[2]
    lb=gap_penalties_alignemnt_problem(s1,s2,i1+1,i2,match,mismatch,sig,ep,1,dp)
    b=lb[0]-indel
    sb1=s1[i1]+lb[1]
    sb2='-'+lb[2]
    lc=gap_penalties_alignemnt_problem(s1,s2,i1,i2+1,match,mismatch,sig,ep,1,dp)
    c=lc[0]-indel
    sc1='-'+lc[1]
    sc2=s2[i2]+lc[2]
    d=max(a,b,c)
    if d==a:
        dp[i1][i2]=[a,sa1,sa2]
    elif d==b:
        dp[i1][i2]=[b,sb1,sb2]
    else:
        dp[i1][i2]=[c,sc1,sc2]
    return dp[i1][i2]

def middle_edge_in_alignment_graph(s1,s2,x1,x2,y1,y2,match,mismatch,indel):
    # bug might be present did not test with big cases
    if x1==x2 and y1==y2: return
    if y1==y2:
        m=(x1+x2)//2
        return (m,y1),(m+1,y1)
    if x1==x2:
        m=(y1+y2)//2
        return (x1,m),(x1,m+1)
    indel*=(-1)
    mismatch*=(-1)
    m = (x1+x2)//2
    n = y2-y1+1
    f=[]
    pf=[i*indel for i in range(n)]
    if x1==m: f=pf
    for i in range(x1,m):
        f=['-' for i in range(n)]
        for j in range(n):
            if j!=n-1:
                if s1[i]==s2[y1+j]:
                    if f[j+1]=='-' or pf[j]+match>f[j+1]: f[j+1]=pf[j]+match
                else:
                    if f[j+1]=='-' or pf[j]+mismatch>f[j+1]: f[j+1]=pf[j]+mismatch
                pf[j+1]=max(pf[j]+indel,pf[j+1])
            if f[j]=='-' or pf[j]+indel>f[j]: f[j]=pf[j]+indel
        if i!=m-1: pf=f.copy()
    b=[]
    pb=[(n-i-1)*indel for i in range(n)]
    if x2==m: b=pb
    for i in range(x2,m,-1):
        b=['-' for i in range(n)]
        for j in range(n-1,-1,-1):
            if j!=0:
                if s1[i-1]==s2[y1+j-1]:
                    if b[j-1]=='-' or pb[j]+match>b[j-1]: b[j-1]=pb[j]+match
                else:
                    if b[j-1]=='-' or pb[j]+mismatch>b[j-1]: b[j-1]=pb[j]+mismatch
                pb[j-1]=max(pb[j]+indel,pb[j-1])
            if b[j]=='-' or pb[j]+indel>b[j]: b[j]=pb[j]+indel
        if i!=m+1: pb=b.copy()
    mm = []
    for i in range(n):
        try: mm.append(b[i]+f[i])
        except:
            print(i)
            quit()
    # print(f)
    # print(b)
    # print(mm)
    # print(pb)
    ma=0
    if len(mm):
        ma=max(mm)
    mx=0
    for i in range(n):
        if mm[i]==ma:
            mx=i
            break
    if mx==n-1:
        return (m,mx+y1),(m+1,mx+y1)
    if s1[m]==s2[mx+y1]:
        return (m,mx+y1),(m+1,mx+1+y1)
    a=pb[mx]+indel
    c=b[mx+1]+indel
    d=pb[mx+1]+mismatch
    if d>=a and d>=c:
        return (m,mx+y1),(m+1,mx+y1+1)
    elif a>=c and a>=d:
        return (m,mx+y1),(m+1,mx+y1)
    else:
        return (m,mx+y1),(m,mx+1+y1)

def linear_space_alignment(s1,s2,x1,x2,y1,y2,match,mismatch,indel,m):
    # bug might be present did not test with big cases
    if x1==x2 and y1==y2: return
    # print(x1,x2,y1,y2)
    p1,p2=middle_edge_in_alignment_graph(s1,s2,x1,x2,y1,y2,match,mismatch,indel)
    m[p1]=p2
    linear_space_alignment(s1,s2,x1,p1[0],y1,p1[1],match,mismatch,indel,m)
    linear_space_alignment(s1,s2,p2[0],x2,p2[1],y2,match,mismatch,indel,m)

def two_break_sorting(p,q):
    eq = set()
    vp={}
    # for l in q:
    l=q
    for i in range(len(l)):
        if i==len(l)-1:
            eq.add((l[i],-l[0]))
            eq.add((-l[0],l[i]))
        else:
            eq.add((l[i],-l[i+1]))
            eq.add((-l[i+1],l[i]))
    n=len(eq)
    n//=2
    # for l in p:
    l=p
    for i in range(len(l)):
        if i==len(l)-1:
            vp[l[i]]=-l[0]
            vp[-l[0]]=l[i]
        else:
            vp[l[i]]=-l[i+1]
            vp[-l[i+1]]=l[i]
    # ans=0
    ans = []
    s=set()
    done=set()
    for i in range(1,n+1):
        s.add(i)
        s.add(-i)
    for x in s:
        if x in done: continue
        start = x
        cycle=[]
        cycle.append(x)
        done.add(x)
        done.add(-x)
        x=vp[x]
        x*=(-1)
        ll=[]
        while x!=start:
            cycle.append(x)
            done.add(x)
            done.add(-x)
            x=vp[x]
            x*=-1
        ll.append(cycle)
        ans.append(ll)
    for e in eq:
        x=e[0]
        y=e[1]
        nx=vp[x]
        ny=vp[y]
        if nx==y and ny==x: continue
        vp[x]=y
        vp[y]=x
        vp[nx]=ny
        vp[ny]=nx
        s=set()
        done=set()
        for i in range(1,n+1):
            s.add(i)
            s.add(-i)
        ll=[]
        for x in s:
            if x in done: continue
            start = x
            cycle=[]
            cycle.append(x)
            done.add(x)
            done.add(-x)
            x=vp[x]
            x*=(-1)
            while x!=start:
                cycle.append(x)
                done.add(x)
                done.add(-x)
                x=vp[x]
                x*=-1
            ll.append(cycle)
        ans.append(ll)
        # ans+=1
    return ans
