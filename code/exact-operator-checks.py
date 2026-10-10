"""Independent exact checks of the universal operator proof's normalization.

No imports from the original pilot, no writes outside this directory.
The proof is in ../paper.tex; finite checks are validation only.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import hashlib, json, math, time

BASE=Path(__file__).resolve().parent
START=time.monotonic()

def choose(n,k):
    return math.comb(n,k) if 0<=k<=n else 0

def masks_of_size(m,j):
    return [sum(1<<i for i in subset) for subset in combinations(range(m),j)]

def central_weight_layers(n):
    """E_(n-1) from the central transposition sum via exact interpolation."""
    m=n-1
    out={}
    for j in range(m+1):
        masks=masks_of_size(m,j); size=len(masks); idx={a:i for i,a in enumerate(masks)}
        ellmax=min(j,m-j)
        eig=[F(m*(m-1),2)-l*(m-l+1) for l in range(ellmax+1)]
        coeff=[F(0)]*(ellmax+1)
        for l in range(ellmax+1):
            p=[F(1)]; denominator=F(1)
            for h in range(ellmax+1):
                if h==l: continue
                temp=[F(0)]*(len(p)+1)
                for k,c in enumerate(p):
                    temp[k]-=eig[h]*c;temp[k+1]+=c
                p=temp;denominator*=eig[l]-eig[h]
            weight=F(1,n*choose(n,l))
            for k,c in enumerate(p): coeff[k]+=weight*c/denominator
        # The transposition sum has a constant diagonal and Johnson-graph edges.
        diagonal=choose(m,2)-j*(m-j)
        neighbors=[]
        for a in masks:
            ones=[i for i in range(m) if a&(1<<i)]
            zeros=[i for i in range(m) if not a&(1<<i)]
            neighbors.append([idx[a^(1<<i)^(1<<k)] for i in ones for k in zeros])
        e=[[coeff[-1] if a==b else F(0) for b in range(size)] for a in range(size)]
        for c in reversed(coeff[:-1]):
            e=[[diagonal*e[a][b]+sum((e[k][b] for k in neighbors[a]),F(0))+(c if a==b else 0)
                for b in range(size)] for a in range(size)]
        for a,aa in enumerate(masks):
            for b,bb in enumerate(masks):out[aa,bb]=e[a][b]
    return out

def delete_bit(a,pos):
    return (a&((1<<pos)-1))|((a>>(pos+1))<<pos)

def orbit_representatives(N,j,k):
    for shared in range(max(0,j+k-N),min(j,k)+1):
        a=(1<<j)-1
        b=((1<<shared)-1)|(((1<<(k-shared))-1)<<j)
        yield a,b,shared

def key_lemma_check(n,e):
    counts={'Z':0,'raising':0}
    for j in range(n+1):
        for a,b,i in orbit_representatives(n,j,j):
            lhs=F(0)
            for pos in range(n):
                aa=(a>>pos)&1;bb=(b>>pos)&1
                if aa==bb:lhs+=(1-2*aa)*e.get((delete_bit(a,pos),delete_bit(b,pos)),F(0))/n
            rhs=F(n-2*j,n*n*choose(n,j))
            assert lhs==rhs,(n,'Z',j,i,lhs,rhs)
            counts['Z']+=1
    for k in range(n):
        j=k+1
        for a,b,i in orbit_representatives(n,j,k):
            lhs=F(0)
            for pos in range(n):
                if (a>>pos)&1 and not ((b>>pos)&1):
                    lhs+=e.get((delete_bit(a,pos),delete_bit(b,pos)),F(0))/n
            rhs=F(k+1,n*n*choose(n,k))
            assert lhs==rhs,(n,'raising',j,k,i,lhs,rhs)
            counts['raising']+=1
    # E preserves subset weight. All other blocks vanish. The lowering check
    # is the transpose of raising since E is real symmetric.
    assert all(value==e.get((b,a),F(0)) for (a,b),value in e.items())
    return counts

def tensor_recurrence_check(n,e,d,x,y):
    m=n-1
    @lru_cache(None)
    def a(j,k):
        if not (0<=j<=m and 0<=k<=m):return F(0)
        return sum((choose(m,j)*choose(j,h)*choose(m-j,k-h)*d**(m-j-k+2*h)*x**(j-h)*y**(k-h)
                    for h in range(max(0,j+k-m),min(j,k)+1)),F(0))
    @lru_cache(None)
    def b(j,k,i):
        if min(j,k,i)<0 or i>min(j,k) or max(j,k)>m or j+k-i>m:return F(0)
        c=(a(j,k)-a(j-1,k-1))/(choose(n,j)*choose(n,k))
        return (c+i*b(j-1,k-1,i-1))/(n-j-k+i)
    counts=0
    for j in range(m+1):
        candidates=masks_of_size(m,j)
        for k in range(m+1):
            for aa,bb,i in orbit_representatives(m,j,k):
                direct=sum((e.get((aa,t),F(0))*d**(m-j-k+2*(t&bb).bit_count())*
                            x**(j-(t&bb).bit_count())*y**(k-(t&bb).bit_count()) for t in candidates),F(0))
                assert direct==b(j,k,i),(n,j,k,i,direct,b(j,k,i))
                counts+=1
    return counts

rows=[]
for n in range(1,9):
    e=central_weight_layers(n)
    lemma=key_lemma_check(n,e)
    nonhermitian=tensor_recurrence_check(n,e,F(3,2),F(2,3),F(5,7))
    singular=tensor_recurrence_check(n,e,F(1),F(-1),F(-1))
    rows.append({'n':n,'dimension':2**(n-1),'lemma_orbits':lemma,
                 'nonhermitian_tensor_recurrence_orbits':nonhermitian,
                 'singular_tensor_recurrence_orbits':singular,'exact_equalities':True})

# Check factorial normalization and the branch cancellation at larger n.
normalization=[]
for n in range(1,101):
    m=n-1
    for l in range(m//2+1):
        weight=F(1,n*choose(n,l))
        matching_count=F(math.factorial(m),2**l*math.factorial(l)*math.factorial(m-2*l))
        c=F(2**l*math.factorial(l)**2*math.factorial(m-2*l)*math.factorial(n-l),n*n*math.factorial(m)**2)
        assert c==weight/matching_count
    for l in range(1,(n-1)//2+1):
        twice_J=n-2*l
        first=F(choose(m,l)-choose(m,l-1),twice_J*n*choose(n,l))
        second=F(choose(m,l-1)-choose(m,l-2),(twice_J+2)*n*choose(n,l-1))
        assert first==second==F(1,n*n)
    normalization.append(n)

result={'exact_operator_and_recurrence_checks':rows,
        'normalization_and_branch_cancellation_n':normalization,
        'universal_proof':'../paper.tex','seconds':time.monotonic()-START,
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'limits':'Finite checks validate implementation and normalization. Universal reasoning is in the accompanying proof.'}
(BASE/'exact-operator-checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'all_exact_checks_passed':True,'rows':rows,'seconds':result['seconds']},indent=2))
