from fractions import Fraction
from wick.expression import AExpression
from wick.wick import apply_wick
from wick.convenience import commute, one_e_spin, two_e_spin, E2_spin, braE1_spin

o_general = "ijklmno"
v_general = "abcdefg"

o_alpha = [f"{i}_a" for i in o_general]
v_alpha = [f"{i}_a" for i in v_general]

o_beta = [f"{i}_b" for i in o_general]
v_beta = [f"{i}_b" for i in v_general]

index_dict = {"o_a": o_alpha, "v_a": v_alpha,
              "o_b": o_beta, "v_b": v_beta,
              "o": o_general, "v": v_general}

H1 = one_e_spin("F", [("o_a", "v_a"), ("o_b", "v_b")], norder=True, index_key=index_dict)
H2 = two_e_spin("g", [("o_a", "v_a"), ("o_b", "v_b")], norder=True, index_key=index_dict)

T2 = E2_spin("t2", [("o_a", "v_a"), ("o_b", "v_b")], index_key=index_dict)
T = T2

singles_projection = braE1_spin(["o_a", "o_b"], ["v_a", "v_b"], index_key=index_dict)

H = H1 + H2
print("[H, T]")
HT = commute(H, T)
print("[[H, T], T]")
HTT = commute(HT, T)

Omega1 = singles_projection*(H + HT + Fraction('1/2')*HTT)

print("\nApply wick")
out = apply_wick(Omega1)

print("\nResolve")
out.resolve()

print("\nOmega1 with explicit spin contributions:")
explicit_spin = AExpression(Ex=out)
print(explicit_spin)

space_dict = {"o_a": "o", "o_b": "o", "v_a": "v", "v_b": "v"}

final = explicit_spin.update_index_spaces(space_dict)
final.simplify()
print("\nOmega1 spin-adapted expression:")
print(final)
