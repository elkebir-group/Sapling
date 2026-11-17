"""Graphviz DOT for clone trees with events on edges (compact form of clonesim's -dot output)."""


def write_dot(tree, u, data, title="", match=None):
    """Nodes: clone index and usage. Edges: SNV cluster (with size) and copy-number events.
    `match` maps event -> True/False (placed as in the true tree, by ancestor relation) to color it."""
    def col(e):
        if match is None or e not in match:
            return "black"
        return "darkgreen" if match[e] else "red"
    out = ["digraph T {", '  node [shape=box, style=rounded, fontname="Helvetica", fontsize=10];',
           '  edge [fontname="Helvetica", fontsize=9];']
    if title:
        out.append(f'  labelloc=t; label="{title}";')
    for v in range(tree.k):
        lab = "normal" if v == 0 else f"clone {v}"
        style = ', style="rounded,dashed", color=gray' if u[v] < 1e-3 else ""
        out.append(f'  {v} [label="{lab}\\nu = {u[v]:.3f}"{style}];')
    for v in range(1, tree.k):
        parts = []
        for c in sorted(tree.clusters[v]):
            e = ("c", c)
            parts.append(f'<font color="{col(e)}"><b>C{c}</b> ({len(data.clusters[c])} SNVs)</font>')
        for s, (a, b) in sorted(tree.cn[v].items()):
            e = ("s", s, (a, b))
            parts.append(f'<font color="{col(e)}">s{s}: ({a},{b})</font>')
        out.append(f'  {tree.parent[v]} -> {v} [label=<{"<br/>".join(parts)}>];')
    out.append("}")
    return "\n".join(out)
