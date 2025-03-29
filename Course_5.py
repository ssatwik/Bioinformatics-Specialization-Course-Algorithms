def eucledian_distance(point1,point2):
    m=len(point1)
    ans=0
    for i in range(m):
        ans+=(point1[i]-point2[i])**2
    ans**=(0.5)
    return ans

def farthest_first_traversal(k,points):
    n=len(points)
    m=len(points[0])
    d=[0 for i in range(n)]
    centers=[points[0]]
    while len(centers)<k:
        max=0
        cc=0
        for i in range(n):
            for c in centers:
                x=0
                for p in range(m):
                    x+=(points[i][p]-c[p])**2
                if len(centers)==1 or x<d[i]: d[i]=x
            if d[i]>max:
                max=d[i]
                cc=i
        centers.append(points[cc])
    return centers
                
def squared_error_distortion(centers,points):
    n=len(points)
    m=len(points[0])
    d=[0 for i in range(n)]
    ans=0
    for i in range(n):
        for c in centers:
            x=0
            for p in range(m):
                x+=(points[i][p]-c[p])**2
            if c==centers[0] or x<d[i]: d[i]=x
        ans+=d[i]
    ans/=n
    return ans

def lloyd_clustering(k,points):
    n=len(points)
    m=len(points[0])
    centers=[]
    for i in range(k): centers.append(points[i].copy())
    cluster=[-1 for i in range(n)]
    new_cluster=cluster.copy()
    d=[-1 for i in range(n)]
    f=1
    debug=0
    while f==1:
        debug+=1
        # print(centers)
        f=0
        d=[-1 for i in range(n)]
        new_cluster=[-1 for i in range(n)]
        for i in range(n):
            for j in range(k):
                x=0
                for p in range(m):
                    x+=(points[i][p]-centers[j][p])**2
                if d[i]==-1 or x<d[i]:
                    d[i]=x
                    new_cluster[i]=j
        if cluster!=new_cluster: f=1
        cluster=new_cluster.copy()
        cog=[[0 for i in range(m)] for j in range(k)]
        freq=[0 for i in range(k)]
        for i in range(n):
            c=cluster[i]
            freq[c]+=1
            for j in range(m): cog[c][j]+=points[i][j]
        for i in range(k):
            for j in range(m):
                # if freq[i]==0: continue
                cog[i][j]/=freq[i]
                centers[i][j]=cog[i][j]
    return centers


def cluster_distance_combining(c1,d1,c2,d2):
    n1=len(c1)
    n2=len(c2)
    ans = d1*n1+d2*n2
    ans/=(n1+n2)
    return ans

def heirarchial_clustering(ll):
    n=len(ll)
    for i in range(n):
        for j in range(n): ll[i].append(-1)
    for i in range(n):
        ll.append([])
        for j in range(2*n): ll[i+n].append(-1)
    clusters=[[i] for i in range(n)]
    indices=set([i for i in range(n)])
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
        for x in c: print(x+1,end=' ')
        print()
        indices.remove(i1)
        indices.remove(i2)
        for i in indices:
            ll[nn][i]=cluster_distance_combining(c1,ll[i1][i],c2,ll[i2][i])
            ll[i][nn]=ll[nn][i]
        indices.add(nn)
    return

