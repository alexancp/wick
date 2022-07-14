from fractions import Fraction
from wick.expression import AExpression, Expression, Term
from wick.operator import Tensor
from wick.wick import apply_wick
from wick.convenience import commute
from wick.convenience import one_e_spin, two_e_spin, get_g_oovo
from wick.convenience import E1_spin, braE1_spin, bath_projector

from re import finditer

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

def find_and_replace_external_indices(lines, index_dict):
    out = ""
    for i, term in enumerate(lines.split("\n")):
        is_sum = False
        if "\\sum_{" in term:
            is_sum = True

        """
        Find all indices, extract summation indices,
        replace all other indices according to the dict
        """
        matches = finditer(r"_{(\w+)}", term)
        previous_end = 0
        sum_indices = "_{}"
        new_term = ""
        for m in matches:
            # Add unmatched part of the string
            new_term += term[previous_end:m.start()]
            previous_end = m.end()
            if not is_sum:
                for char in m.group():
                    if char in sum_indices or char not in index_dict:
                        new_term += char
                    else:
                        new_term += index_dict[char]
            else:
                sum_indices = m.group()
                new_term += m.group()
            is_sum = False # Summation indices are found in first match

        out += f"{new_term}{term[previous_end:]}\n"

    return out

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
final.simplify()
final.sort()
nice = make_nice_string(str(final))
print(nice)


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
final.simplify()
final.sort()
nice = make_nice_string(str(final))
print(nice)


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
final.simplify()
final.sort()
nice = make_nice_string(str(final))
print(nice)


print("\n\nLefthand side <^a_i| P^T H^T |R>:")
print("---------------------------------")

out = apply_wick(singles_projection*transformed_P*transformed_H)#, occ=["o_a", "o_b"])
out.resolve()
explicit_spin = AExpression(Ex=out)

# Replace ranges with explicit range by general ones
final = explicit_spin.update_index_spaces(space_dict)
final.simplify()
final.sort()
nice = make_nice_string(str(final))
print(nice)