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


def simplify_sin_cos_tensors(lines):
    out = ""
    for term in lines.split("\n"):
        temp = introduce_exponent(term, "\\sin(\\theta)_{}", 4)
        out += introduce_exponent(temp, "\\cos(\\theta)_{}", 4)
        out += "\n"
    final = out.replace("_{}", "")
    return final