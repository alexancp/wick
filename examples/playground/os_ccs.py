from fractions import Fraction
from wick.expression import AExpression
from wick.operator import Tensor
from wick.wick import apply_wick
from wick.convenience import commute
from wick.convenience import one_e_spin, two_e_spin, get_g_oovo
from wick.convenience import E1_spin, braE1_spin, bath_projector

from re import findall

def sort_by_deltas(string: str):
    lines = string.split("\n")
    delta_dict = {}
    for line in lines:
        #n_delta = line.count("delta_{")
        deltas = "".join(sorted(findall(r"(delta_\{\w\w\})", line)))
        if deltas in delta_dict:
            delta_dict[deltas] += f"{line}\n"
        else:
            delta_dict[deltas] = f"{line}\n"
    out = ""
    for i in sorted(list(delta_dict)):
        print(i)
        out += delta_dict[i]
    return out

def introduce_exponent(string: str, tensor: str, position: int):
    exponent = string.count(tensor)
    if exponent == 1:
        return string
    temp = string.replace(tensor, "", exponent-1)
    replacement = f"{tensor[:position]}^{{{exponent}}}{tensor[position:]}"
    out = temp.replace(tensor, replacement, 1)
    return out

def make_nice_string(lines):
    out = ""
    for term in lines.split("\n"):
        temp = introduce_exponent(term, "\\sin(\\theta)_{}", 4)
        out += introduce_exponent(temp, "\\cos(\\theta)_{}", 4)
        out += "\n"
    final = out.replace("_{}", "")
    return final

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

P_a = bath_projector("o_a", "v_a", sin_t, cos_t,) # suffix="_a")
P_b = bath_projector("o_b", "v_b", sin_t, cos_t,) # suffix="_b")
P = P_b - P_b*P_a

# Define cluster operator
T1 = E1_spin("t1", ["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)
T = T1

PT = commute(P, T)
PTT = commute(PT, T)

transformed_P = P + PT + Fraction(1,2) * PTT

print("\n\nRighthand side <R| P^T |R>:")
print("---------------------------")

out = apply_wick(transformed_P)#, occ=["o_a", "o_b"])
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
space_dict = {"o_a": "o", "o_b": "o", "v_a": "v", "v_b": "v"}
final = explicit_spin.update_index_spaces(space_dict)
print(make_nice_string(str(final)))


print("\n\nLefthand side <R| P^T H^T |R>:")
print("------------------------------")

# Define Hamiltonian
H1 = one_e_spin("h", [("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict, )
H2 = two_e_spin("g", [("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict, )

# Use T1 transformed H
H = H1 + H2
transformed_H = H

out = apply_wick(transformed_P*transformed_H)#, occ=["o_a", "o_b"])
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
L = final.introduce_Coulomb_minus_Exchange("g", (0, 3, 2, 1), "L")
F = L.introduce_Fock_matrix()
print(make_nice_string(str(F)))
new = sorted(F.terms, key=lambda x: len(x.tensors))
#print("\n".join([t._print_str() for t in new]))
print("\n".join([t._print_str() for t in new]))

print("\n\nRighthand side <^a_i| P^T |R>:")
print("------------------------------")

PTTT = commute(PTT, T)
transformed_P = P + PT + Fraction(1, 2) * PTT + Fraction(1, 6) * PTTT

singles_projection = braE1_spin(["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)

out = apply_wick(singles_projection*transformed_P)#, occ=["o_a", "o_b"])
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
print(sort_by_deltas(make_nice_string(str(final))))


print("\n\nLefthand side <^a_i| P^T H^T |R>:")
print("---------------------------------")

out = apply_wick(singles_projection*transformed_P*transformed_H)#, occ=["o_a", "o_b"])
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
L = final.introduce_Coulomb_minus_Exchange("g", (0, 3, 2, 1), "L")
F = L.introduce_Fock_matrix()
print(make_nice_string(str(F)))

new = sorted(F.terms, key=lambda x: len(x.tensors))
print(sort_by_deltas("\n".join([t._print_str() for t in new])))