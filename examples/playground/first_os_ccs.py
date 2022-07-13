from fractions import Fraction
from wick.expression import AExpression, Term
from wick.index import Idx
from wick.operator import Tensor
from wick.wick import apply_wick
from wick.convenience import commute, one_e_spin, two_e_spin, E1_spin, braE1_spin, E2_spin, braE2_spin, bath_projector
from sympy import sin, cos
from sympy.core.numbers import pi



def bath_projector(ospace, vspace, sin, cos, index_key=None):
    """
    Return an expression representing all pieces of a one-electron operator

    name (str): Name of the operator.
    spaces (list): List orbital subspaces
    norder (bool): Return only normal-ordered part?
    """
    terms = []
    factors = [sin*sin, cos*cos, -1*sin*cos, -1*sin*cos]
    print(factors)
    dagger_indices = [Idx(0, ospace), Idx(0, vspace), Idx(0, ospace), Idx(0, vspace)]
    indices = [Idx(0, ospace), Idx(0, vspace), Idx(0, vspace), Idx(0, ospace)]

    for p, q, f in zip(dagger_indices, indices, factors):
        operators = [FOperator(p, True), FOperator(q, False)]
        sign = 1
        t = Term(
            sign*f, [], [], operators, [], index_key=index_key)
        print(t)
        terms.append(t)

    return Expression(terms)

# sin = Tensor([], "sin(\\theta)")
# cos = Tensor([], "cos(\\theta)")
# print(Term(Float(3), [], [sin], [], [])*Term(1, [], [sin], [], []))
# print(Term(1, [], [sin], [], [])*Term(1, [], [cos], [], []))

sin_t = sin(pi/6)
cos_t = cos(pi/6)
print(f"sin: {sin_t}")
print(f"cos: {cos_t}")

o_general = "ijklmno"
v_general = "abcdefg"

o_alpha = [f"{i}_a" for i in o_general]
v_alpha = [f"{i}_a" for i in v_general]

o_beta = [f"{i}_b" for i in o_general]
v_beta = [f"{i}_b" for i in v_general]

index_dict = {"o_a": o_alpha, "v_a": v_alpha,
              "o_b": o_beta, "v_b": v_beta,
              "o": o_general, "v": v_general}

# Need several indices because we treate alpha and beta separately
O_general = "IIII"
V_general = "AAAA"
O_alpha = [f"{i}_a" for i in o_general]
V_alpha = [f"{i}_a" for i in v_general]
O_beta = [f"{i}_b" for i in o_general]
V_beta = [f"{i}_b" for i in v_general]

bath_key = {"o_a": O_alpha, "v_a": V_alpha,
            "o_b": O_beta, "v_b": V_beta,
            "o": O_general, "v": V_general}

P_a = bath_projector("o_a", "v_a", sin_t, cos_t, index_key=bath_key)
P_b = bath_projector("o_b", "v_b", sin_t, cos_t, index_key=bath_key)
P = P_b - P_b*P_a

H1 = one_e_spin("F", [("o_a", "v_a"), ("o_b", "v_b")], norder=True, index_key=index_dict)
H2 = two_e_spin("g", [("o_a", "v_a"), ("o_b", "v_b")], norder=True, index_key=index_dict)

T1 = E1_spin("t1", ["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)
print("\nT1:", T1)
singles_projection = braE1_spin(["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)

T = T1

print("\t[P T]")
PT = commute(P, T)
print("\t[[P T], T]")
PTT = commute(PT, T)

print("Project on R")
print("---------------")
print("Righthand side:")
print("\n\tApply wick")
out = apply_wick(P + PT + Fraction(1,2) * PTT)
print("\n\tResolve")
out.resolve()

print("\n\tExplicit spin:")
explicit_spin = AExpression(Ex=out)
print("\t", explicit_spin)

space_dict = {"o_a": "o", "o_b": "o", "v_a": "v", "v_b": "v"}
final = explicit_spin.update_index_spaces(space_dict)
final.simplify()

print("\n\tRHS spin-adapted expression:")
print(final)

# H = H1 + H2
# print("[H, T]")
# HT = commute(H, T)
# print("[[H, T], T]")
# HTT = commute(HT, T)
#
# Omega2 = doubles_projection*(H + HT + Fraction('1/2')*HTT)
#
# print("\nApply wick")
# out = apply_wick(Omega2)
#
# print("\nResolve")
# out.resolve()
#
# print("\nOmega2 with explicit spin contributions:")
# explicit_spin = AExpression(Ex=out)
# print(explicit_spin)
#
# space_dict = {"o_a": "o", "o_b": "o", "v_a": "v", "v_b": "v"}
#
# final = explicit_spin.update_index_spaces(space_dict)
# final.simplify()
# print("\nOmega2 spin-adapted expression:")
# print(final)
