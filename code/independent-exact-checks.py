"""Independent bounded checks, using Casimir spectral projectors and full matrices.

This is ancillary evidence. The all-n argument is in ../paper.tex.
No scientific input is edited or imported by this checker.
"""
from pathlib import Path
import datetime, itertools, json, math, time
import sympy as S

OUT = Path(__file__).resolve().parent
START = time.monotonic()
results = []

def record(kind, **values):
    item = {'kind': kind, **values}
    results.append(item)
    (OUT/'independent-exact-results.json').write_text(json.dumps({
        'method': 'Casimir polynomial projectors, full subset matrices, exact rational-complex arithmetic',
        'elapsed_seconds': time.monotonic()-START,
        'checks': results,
    }, indent=2)+'\n')
    print(json.dumps(item), flush=True)

def without(a, s):
    return (a & ((1 << s)-1)) | ((a >> (s+1)) << s)

def indexed_matrix(rows,cols,func):
    return S.Matrix(rows,cols,lambda a,b:func(int(a),int(b)))

def assert_matrix_equal(a,b):
    assert a.shape==b.shape
    difference=(a-b).applyfunc(S.expand)
    assert difference==S.zeros(*a.shape), str(difference)

def embed_sum(Y, B, n):
    return indexed_matrix(1 << n, 1 << n, lambda a,b:
        sum(Y[(a >> s)&1, (b >> s)&1] * B[without(a,s), without(b,s)]
            for s in range(n)))

def symproj(n):
    return indexed_matrix(1 << n, 1 << n, lambda a,b:
        S.Rational(1, math.comb(n,a.bit_count())) if a.bit_count()==b.bit_count() else 0)

def weighted_casimir(m):
    size = 1 << m
    I = S.eye(size)
    # J^2 = m(4-m)/4 I + sum of all factor transpositions.
    cas = S.Rational(m*(4-m),4)*I
    for a,b in itertools.combinations(range(m),2):
        swap = S.zeros(size)
        for u in range(size):
            v = u ^ ((1 << a)|(1 << b)) if ((u >> a)&1)!=((u >> b)&1) else u
            swap[v,u] = 1
        cas += swap
    eig = [S.Rational(m-2*l,2)*(S.Rational(m-2*l,2)+1) for l in range(m//2+1)]
    projections = []
    for l,e in enumerate(eig):
        P = I
        for f in eig:
            if f != e: P = P*(cas-f*I)/(e-f)
        projections.append(P)
    E = S.zeros(size)
    for l,P in enumerate(projections): E += P/S.Integer((m+1)*math.comb(m+1,l))
    assert sum(projections, S.zeros(size)) == I
    for P in projections: assert P*P == P
    return E, projections

Z = S.diag(1,-1)
Ybasis = [Z, S.Matrix([[0,1],[0,0]]), S.Matrix([[0,0],[1,0]])]
for n in range(1,7):
    E, _ = weighted_casimir(n-1)
    ps = symproj(n)
    for Y in Ybasis:
        actual = embed_sum(Y,E,n)/n
        expected = embed_sum(Y,S.eye(1 << (n-1)),n)*ps/(n*n)
        assert actual == expected
    record('twirl_full_matrix', n=n, basis_count=3, entries_per_basis=(1 << n)**2, passed=True)

def product_tensor(mats):
    if not mats: return S.ones(1,1)
    # The bit indexing makes slot 0 the least significant bit.
    return S.kronecker_product(*mats[::-1])

def sym_tensor(mats):
    m = len(mats)
    out = S.zeros(1 << m)
    for p in itertools.permutations(range(m)):
        out += product_tensor([mats[i] for i in p])
    return out/S.factorial(m)

def coeff_compression(G,m):
    U = indexed_matrix(1 << m,m+1,lambda a,j: int(a.bit_count()==j))
    return U.T*G*U

def numerator_matrix(A,n):
    def a(j,k):
        return A[j,k] if 0<=j<n and 0<=k<n else 0
    return indexed_matrix(1 << n,1 << n,lambda x,y:
        (a(x.bit_count(),y.bit_count())-a(x.bit_count()-1,y.bit_count()-1)) /
        S.Integer(math.comb(n,x.bit_count())*math.comb(n,y.bit_count())))

for n in range(1,6):
    m=n-1
    E,_=weighted_casimir(m)
    mats=[S.Matrix([[i+2,(i+1)+S.I],[-i+2*S.I,i+2]]) for i in range(m)]
    G=sym_tensor(mats)
    A=coeff_compression(G,m)
    assert_matrix_equal(embed_sum(Z,E*G,n),numerator_matrix(A,n))
    # This compares every entry, including union-full entries, without recurrence solving.
    record('equal_diagonal_full_transfer',n=n,heterogeneous_nonhermitian=True,
           union_full_entries=sum(1 for x in range(1 << n) for y in range(1 << n) if (x|y)==(1 << n)-1),
           passed=True)

# Real assumptions make conjugation act only on coefficients; u denotes bar(w).
t,u=S.symbols('t u',real=True)
root_pool=[S.Integer(1), S.I, S.Integer(-1), S.Rational(3,5)+S.Rational(4,5)*S.I]
c=2+3*S.I

def polynomial_coeffs(p, degree):
    p=S.Poly(S.expand(p),t)
    return S.Matrix([p.nth(i) for i in range(degree+1)])

def matchings(ids, count):
    if count==0:
        yield ()
        return
    if len(ids)<2*count: return
    first,*tail=ids
    if len(tail)>=2*count:
        yield from matchings(tail,count)
    for j in tail:
        rem=[x for x in tail if x!=j]
        for rest in matchings(rem,count-1): yield ((first,j),)+rest

for n in range(1,6):
    m=n-1
    for repeated in [False,True]:
        roots=(root_pool[:m] if not repeated else [S.Integer(1)]*m)
        E,Ps=weighted_casimir(m)
        mats=[S.Matrix([1,-a])*S.Matrix([1,-a]).conjugate().T for a in roots]
        G=S.expand(c*S.conjugate(c))*sym_tensor(mats)
        B=S.simplify(E*G)
        r=c*S.prod(1-a*t for a in roots)
        rc=polynomial_coeffs(r,m)
        assert_matrix_equal(coeff_compression(G,m),rc*rc.conjugate().T)
        assert_matrix_equal(embed_sum(Z,B,n),numerator_matrix(rc*rc.conjugate().T,n))
        for ell in range(m//2+1):
            d=m-2*ell
            singlet=S.Matrix([0,1,-1,0])
            cols=[]
            for a in range(d+1):
                vec=S.Matrix([int(x.bit_count()==a) for x in range(1 << d)])
                for _ in range(ell): vec=S.kronecker_product(vec,singlet)
                cols.append(vec)
            W=S.Matrix.hstack(*cols)
            assert Ps[ell]*W==W
            H=S.simplify(W.conjugate().T*B*W)
            rhs=S.zeros(d+1)
            for pairs in matchings(list(range(m)),ell):
                used={i for pair in pairs for i in pair}
                b=polynomial_coeffs(S.prod(1-roots[i]*t for i in range(m) if i not in used),d)
                delta=S.prod((roots[i]-roots[j])*S.conjugate(roots[i]-roots[j]) for i,j in pairs)
                rhs+=delta*(b*b.conjugate().T)
            normal=S.Rational(2**ell*math.factorial(ell)**2*math.factorial(d)*math.factorial(n-ell),n*n*math.factorial(m)**2)
            assert S.simplify(H-S.expand(c*S.conjugate(c))*normal*rhs)==S.zeros(d+1)
        record('root_and_matching_full_matrix',n=n,repeated_roots=repeated,complex_scale=True,passed=True)

def reflected(q,n):
    qc=polynomial_coeffs(q,n)
    return sum(S.conjugate(qc[n-j])*t**j for j in range(n+1))

def kernel(q,n):
    star=reflected(q,n)
    num=q*S.conjugate(q).subs(t,u)-star*S.conjugate(star).subs(t,u)
    return S.cancel(num/(1-t*u))

clark_cases=[]
for n in [1,2]:
    b=(1+S.I)/4
    for k in range(n+1):
        q=S.expand(c*(1-b*t)**k)
        star=reflected(q,n)
        lam=S.cancel(star.subs(t,1)/q.subs(t,1))
        Q=S.expand(q-S.conjugate(lam)*star)
        if n==1: nodes=[S.Integer(1)]
        else:
            quotient=S.div(Q,t-1,t)[0]
            nodes=[S.Integer(1),S.cancel(-S.Poly(quotient,t).nth(0)/S.Poly(quotient,t).nth(1))]
        clark_cases.append((f'n{n}-k{k}-complex-repeated',q,n,lam,nodes))
for k in [0,4]:
    n=4
    q=c if k==0 else S.expand(c*(1-S.Rational(2,5)*t**4))
    lam=S.conjugate(c)/c
    clark_cases.append((f'n4-k{k}',q,n,lam,[1,S.I,-1,-S.I]))

for name,q,n,lam,nodes in clark_cases:
    star=reflected(q,n)
    Q=S.expand(q-S.conjugate(lam)*star)
    phi=star/q
    rhs=0
    weights=[]
    for xi in nodes:
        assert S.simplify(xi*S.conjugate(xi))==1
        assert S.simplify(phi.subs(t,xi)-lam)==0
        r=S.cancel(Q/(t-xi))
        rho=S.simplify(xi*S.conjugate(lam)*S.diff(phi,t).subs(t,xi))
        assert rho.is_positive
        weights.append(str(rho))
        rhs += r*S.conjugate(r).subs(t,u)/rho
    assert S.simplify(kernel(q,n)-rhs)==0
    record('clark_exact_phase_and_weights',case=name,n=n,weights=weights,passed=True)

record('completion',time_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),all_passed=True)
