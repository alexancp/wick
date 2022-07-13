from fractions import Fraction
from wick.expression import AExpression, Term
from wick.operator import Tensor
from wick.wick import apply_wick
from wick.convenience import commute, one_e_spin, two_e_spin, E1_spin, braE1_spin, bath_projector

o_general = "ijklmno"
v_general = "abcdefg"

o_alpha = [f"{i}_a" for i in o_general]
v_alpha = [f"{i}_a" for i in v_general]

o_beta = [f"{i}_b" for i in o_general]
v_beta = [f"{i}_b" for i in v_general]

index_dict = {"o_a": o_alpha, "v_a": v_alpha,
              "o_b": o_beta, "v_b": v_beta,
              "o": o_general, "v": v_general}

O_general = "IIII"
V_general = "AAAA"
O_alpha = [f"{i}_a" for i in O_general]
V_alpha = [f"{i}_a" for i in V_general]
O_beta = [f"{i}_b" for i in O_general]
V_beta = [f"{i}_b" for i in V_general]

bath_key = {"o_a": O_alpha, "v_a": V_alpha,
            "o_b": O_beta, "v_b": V_beta,
            "o": O_general, "v": V_general}

# Define projection operator
sin_t = Tensor([], "\\sin(\\theta)")
cos_t = Tensor([], "\\cos(\\theta)")

P_a = bath_projector("o_a", "v_a", sin_t, cos_t)#, index_key=bath_key)
P_b = bath_projector("o_b", "v_b", sin_t, cos_t)#, index_key=bath_key)
P = P_b - P_b*P_a

# Define cluster operator
T1 = E1_spin("t1", ["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)
T = T1

print("[P T]")
PT = commute(P, T)
print("[[P T], T]")
PTT = commute(PT, T)

transformed_P = P + PT + Fraction(1,2) * PTT

print("\n\nRighthand side <R| P^T |R>:")
print("---------------------------")

print("\nApply wick")
out = apply_wick(transformed_P)

print("\nResolve")
out.resolve()

print("\nRHS with explicit spin:")
explicit_spin = AExpression(Ex=out)
print(explicit_spin)

# Replace ranges with explicit range by general ones
space_dict = {"o_a": "o", "o_b": "o", "v_a": "v", "v_b": "v"}
final = explicit_spin.update_index_spaces(space_dict)
final.simplify()

print("\nRHS spin-adapted expression:")
print(final)


print("\n\nLefthand side <R| P^T H^T |R>:")
print("------------------------------")

# Define Hamiltonian
H1 = one_e_spin("F", [("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict)
H2 = two_e_spin("g", [("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict)

H = H1 + H2
print(H1)
#HT = commute(H, T)
#HTT = commute(HT, T)

transformed_H = H #+ HT + Fraction(1,2) * HTT

print("\nApply wick")
out = apply_wick(transformed_P*transformed_H)

print("\nResolve")
out.resolve()

print("\n\nLHS with explicit spin:")
explicit_spin = AExpression(Ex=out)
print(explicit_spin)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
final.simplify()

print("\nLHS spin-adapted expression:")
print(final)


singles_projection = braE1_spin(["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)
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
