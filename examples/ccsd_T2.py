from wick.index import Idx
from wick.expression import *
from wick.hamiltonian import one_e, two_e, E1, E2, commute
from wick.wick import apply_wick

H1 = one_e("f",["occ","vir"], norder=True)
H2 = two_e("I",["occ","vir"], norder=True)
H = H1 + H2

i = Idx(0,"occ")
a = Idx(0,"vir")
j = Idx(1,"occ")
b = Idx(1,"vir")
operators = [FOperator(i,True), FOperator(a,False), FOperator(j,True), FOperator(b,False)]
bra = Expression([Term(1.0, [], [Tensor([i,j,a,b],"")], operators, [])])
T1 = E1("t", ["occ"], ["vir"])
T2 = E2("t", ["occ"], ["vir"])
T = T1 + T2

HT = commute(H,T)
HTT = commute(HT,T)
HTTT = commute(HTT,T)
HTTTT = commute(HTTT,T)

S = bra*(H + HT + (1.0/2.0)*HTT + (1/6.0)*HTTT + (1/24.0)*HTTTT)
out = apply_wick(S)
out.resolve()
final = AExpression(Ex=out)
final.simplify()
final.sort()
print(final._print_str())
