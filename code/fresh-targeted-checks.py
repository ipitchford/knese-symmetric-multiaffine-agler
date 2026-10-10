"""Fresh exact checks targeting degree conventions, complex phases and adjoint order.

Written in this review independently of the supplied executable implementations.
These finite examples are consistency evidence, not a proof for arbitrary n.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import datetime, hashlib, json, math, platform, sys, time
import sympy as s

ROOT = Path(__file__).resolve().parent
started = time.monotonic()
t = s.Symbol('t')
c = 2+3*s.I
rows = []

def coeffs(q, n):
    p = s.Poly(s.expand(q), t)
    return [p.nth(j) for j in range(n+1)]

def A_from_difference(qc, n):
    return s.Matrix(n,n,lambda j,k: s.expand(sum(
        qc[j-h]*s.conjugate(qc[k-h])-s.conjugate(qc[n-j+h])*qc[n-k+h]
        for h in range(min(j,k)+1))))

def equal_matrices(a,b):
    assert a.shape == b.shape
    assert all(s.cancel(x)==0 for x in a-b)

def check_interpolation(label, q, n, nodes):
    qc = coeffs(q,n)
    qstar = sum(s.conjugate(qc[n-j])*t**j for j in range(n+1))
    lam = s.conjugate(c)/c
    Q = s.expand(q-s.conjugate(lam)*qstar)
    assert s.degree(Q,t)==n
    phi = qstar/q
    rhs = s.zeros(n)
    weights = []
    for xi in nodes:
        xi = s.cancel(xi)
        assert s.cancel(xi*s.conjugate(xi))==1
        assert s.cancel(phi.subs(t,xi)-lam)==0
        r,rem = s.div(Q,t-xi,t)
        assert s.cancel(rem)==0
        assert s.degree(r,t)==n-1
        rho = s.cancel(xi*s.conjugate(lam)*s.diff(phi,t).subs(t,xi))
        assert rho.is_positive is True
        assert s.cancel(r.subs(t,xi)+rho*q.subs(t,xi)/xi)==0
        v = s.Matrix(coeffs(r,n-1))
        rhs += v*v.conjugate().T/rho
        weights.append(str(rho))
    equal_matrices(A_from_difference(qc,n),rhs)
    rows.append({'kind':'fresh_exact_interpolation','label':label,'formal_degree':n,
                 'actual_univariate_degree':int(s.degree(q,t)), 'node_count':len(nodes),
                 'weights':weights,'passed':True})

etas = [s.Integer(1),s.I,s.Integer(-1),-s.I]
b = (2+s.I)/5
check_interpolation('four repeated complex reciprocal roots',c*(1-b*t)**4,4,
                    [(eta+s.conjugate(b))/(1+b*eta) for eta in etas])
check_interpolation('degree two inside formal degree four',c*(1-s.Rational(3,5)*t)**2,4,
                    [s.Integer(1),s.Integer(-1),(3+4*s.I)/5,(3-4*s.I)/5])

# A commuting but non-doubly-commuting tuple tests the hereditary ordering.
n = 3
N = s.Matrix([[0,1,0],[0,0,1],[0,0,0]])
Ts = [N/3, N/4+N*N/5, N/5-N*N/7]
Id = s.eye(3)
for T in Ts:
    # The Frobenius norm bounds the operator norm, with this exact value < 1.
    assert s.trace(T*T.conjugate().T)<1
for T in Ts:
    for U in Ts:
        assert T*U==U*T
assert Ts[0]*Ts[1].conjugate().T != Ts[1].conjugate().T*Ts[0]
q = c*(1-(1+s.I)*t/4)*(1-s.I*t/3)*(1+t/5)
qc = coeffs(q,n)
A = A_from_difference(qc,n)

def av(j,k):
    return A[j,k] if 0<=j<n and 0<=k<n else 0

@lru_cache(None)
def bv(j,k,i):
    if min(j,k,i)<0 or i>min(j,k) or max(j,k)>n-1 or j+k-i>n-1:
        return s.Integer(0)
    return s.cancel(((av(j,k)-av(j-1,k-1))/(math.comb(n,j)*math.comb(n,k))
                     +i*bv(j-1,k-1,i-1))/(n-j-k+i))

B = s.Matrix(4,4,lambda a,b: bv(int(a).bit_count(),int(b).bit_count(),(int(a)&int(b)).bit_count()))
equal_matrices(B,B.conjugate().T)
minor_count=0
for size in range(1,5):
    for ids in combinations(range(4),size):
        determinant=s.cancel(B.extract(ids,ids).det())
        assert determinant.is_nonnegative is True
        minor_count+=1

def subset_product(mask,tuple_):
    product=Id
    for j,T in enumerate(tuple_):
        if mask & (1<<j): product=product*T
    return product

pT=s.zeros(3)
ptildeT=s.zeros(3)
for mask in range(1<<n):
    k=mask.bit_count()
    monomial=subset_product(mask,Ts)
    pT+=qc[k]*monomial/math.comb(n,k)
    ptildeT+=s.conjugate(qc[n-k])*monomial/math.comb(n,k)
left=pT*pT.conjugate().T-ptildeT*ptildeT.conjugate().T
right=s.zeros(3)
for slot,T in enumerate(Ts):
    omitted=[U for j,U in enumerate(Ts) if j!=slot]
    for a in range(4):
        for bmask in range(4):
            right+=B[a,bmask]*subset_product(a,omitted)*(Id-T*T.conjugate().T)*subset_product(bmask,omitted).conjugate().T
equal_matrices(left,right)
ph=s.simplify(pT.inv()*ptildeT)
defect=s.simplify(Id-ph*ph.conjugate().T)
operator_minors=[]
for size in range(1,4):
    for ids in combinations(range(3),size):
        determinant=s.cancel(defect.extract(ids,ids).det())
        assert determinant.is_nonnegative is True
        operator_minors.append(str(determinant))
rows.append({'kind':'fresh_exact_hereditary_substitution','tuple_dimension':3,'variable_count':n,
             'tuple_commutes':True,'tuple_not_doubly_commuting':True,'B_principal_minors_checked':minor_count,
             'operator_defect_principal_minors_checked':len(operator_minors),'passed':True})

result={'python_version':sys.version,'sympy_version':s.__version__,
        'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'wall_seconds':time.monotonic()-started,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'checks':rows,'all_passed':True,'limits':'Three exact finite examples target phase, degree and adjoint-order risks; they do not establish the universal theorem.'}
(ROOT/'fresh-targeted-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
