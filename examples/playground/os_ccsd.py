from fractions import Fraction
from wick.expression import AExpression
from wick.operator import Tensor
from wick.wick import apply_wick
from wick.convenience import commute
from wick.convenience import one_e_spin, two_e_spin, get_g_oovo
from wick.convenience import E1_spin, braE1_spin, bath_projector, E2_spin, braE2_spin
from wick.postprocessing_tools import sort_by_deltas, simplify_sin_cos_tensors


o_general = "ijklmno"
v_general = "abcdefg"

o_alpha = [f"{i}_a" for i in o_general]
v_alpha = [f"{i}_a" for i in v_general]

o_beta = [f"{i}_b" for i in o_general]
v_beta = [f"{i}_b" for i in v_general]

index_dict = {"o_a": o_alpha, "v_a": v_alpha,
              "o_b": o_beta, "v_b": v_beta,
              "o": o_general, "v": v_general}

# Define projection operator
sin_t = Tensor([], "\\sin(\\theta)")
cos_t = Tensor([], "\\cos(\\theta)")

P_a = bath_projector("o_a", "v_a", sin_t, cos_t,)
P_b = bath_projector("o_b", "v_b", sin_t, cos_t,)
P = P_b - P_b*P_a

# Define cluster operator
T1 = E1_spin("t1", ["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)
T2 = E2_spin("t2", [("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict)
T = T1 + T2

PT = commute(P, T)
PTT = commute(PT, T)

transformed_P = P + PT + Fraction(1,2) * PTT

print("\n\nRighthand side <R| P^T |R>:")
print("---------------------------")

out = apply_wick(transformed_P)
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
space_dict = {"o_a": "o", "o_b": "o", "v_a": "v", "v_b": "v"}
final = explicit_spin.update_index_spaces(space_dict)
print(simplify_sin_cos_tensors(str(final)))


print("\n\nLefthand side <R| P^T H^T |R>:")
print("------------------------------")

# Define Hamiltonian
H1 = one_e_spin("h", [("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict, )
H2 = two_e_spin("g", [("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict, )

# Use T1 transformed H
H = H1 + H2
HT = commute(H, T)
HTT = commute(HT, T)
transformed_H = H + HT + Fraction(1, 2)*HTT

out = apply_wick(transformed_P*transformed_H)
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
L = final.introduce_Coulomb_minus_Exchange("g", (0, 3, 2, 1), "L")
F = L.introduce_Fock_matrix()
print(simplify_sin_cos_tensors(str(F)))
new = sorted(F.terms, key=lambda x: len(x.tensors))
print(new)

print("\n\nRighthand side <^a_i| P^T |R>:")
print("------------------------------")

PTTT = commute(PTT, T)
transformed_P = P + PT + Fraction(1, 2) * PTT + Fraction(1, 6) * PTTT

singles_projection = braE1_spin(["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)

out = apply_wick(singles_projection*transformed_P)
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
L = final.introduce_Coulomb_minus_Exchange("t2", (0, 3, 2, 1), "u2")
print(sort_by_deltas(simplify_sin_cos_tensors(str(L))))

print("\n\nLefthand side <^a_i| P^T H^T |R>:")
print("---------------------------------")

HTTT = commute(HTT, T)
transformed_H = H + HT + Fraction(1, 2)*HTT + Fraction(1, 6) * HTTT

out = apply_wick(singles_projection*transformed_P*transformed_H)
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
L = final.introduce_Coulomb_minus_Exchange("g", (0, 3, 2, 1), "L")
F = L.introduce_Fock_matrix()
L = F.introduce_Coulomb_minus_Exchange("t2", (0, 3, 2, 1), "u2")
print(sort_by_deltas(simplify_sin_cos_tensors(str(L))))





print("\n\nRighthand side <^{ab}_{ij}| P^T |R>:")
print("------------------------------------")

PTTT = commute(PTT, T)
transformed_P = P + PT + Fraction(1, 2) * PTT + Fraction(1, 6) * PTTT

doubles_projection = braE2_spin([("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict)

out = apply_wick(doubles_projection*transformed_P)#, occ=["o_a", "o_b"])
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
L = final.introduce_Coulomb_minus_Exchange("t2", (0, 3, 2, 1), "u2")
F = L.introduce_permutation_operators((2, 3, 0, 1))
print(sort_by_deltas(simplify_sin_cos_tensors(str(F))))

print("\n\nLefthand side <^{ab}_{ij}| P^T H^T |R>:")
print("---------------------------------------")

HTTT = commute(HTT, T)
transformed_H = H + HT + Fraction(1, 2)*HTT + Fraction(1, 6) * HTTT

out = apply_wick(doubles_projection*transformed_P*transformed_H)
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
L = final.introduce_Coulomb_minus_Exchange("g", (0, 3, 2, 1), "L")
F = L.introduce_Fock_matrix()

L = F.introduce_Coulomb_minus_Exchange("t2", (0, 3, 2, 1), "u2")
F = L.introduce_permutation_operators((2, 3, 0, 1))

new = sorted(F.terms, key=lambda x: len(x.tensors))
print(sort_by_deltas("\n".join([t._print_str() for t in new])))

