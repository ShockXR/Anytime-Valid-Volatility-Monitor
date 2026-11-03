"""Build ``docs/api.md``, the API reference of the ``avmonitor`` package, from its docstrings.

Only the standard library is used (``inspect``, ``importlib``, ``pkgutil``); the package itself must be
importable, so numpy, scipy and pandas must be installed. Run from anywhere:

    python scripts/build_api_docs.py              # write docs/api.md
    python scripts/build_api_docs.py --check      # exit 1 if docs/api.md differs from the docstrings
    python scripts/build_api_docs.py -o api.md    # write somewhere else

The script documents the source tree next to it (``src/avmonitor``), not an installed copy. For every
public module (name without a leading underscore) it lists the first paragraph of the module docstring,
then each public class and function defined in that module, in source order, with its signature and
the first paragraph of its docstring. Public methods and properties defined in a class are listed under
it. The top of the page lists the names that ``avmonitor`` exports (``avmonitor.__all__``).
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import inspect
import pkgutil
import re
import sys
from pathlib import Path, PurePath

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
OUT = ROOT / "docs" / "api.md"
PACKAGE = "avmonitor"

# Modules a user of the monitor needs first, in reading order. Every other public module is listed
# after them as a research module (used by the experiments and the theory notes).
CORE = ["monitor", "periodic", "multifreq", "panel", "detectors", "confseq", "robust", "simulate", "forecast",
        "data", "cli"]
WRAP = 96  # signatures longer than this are written one parameter per line


def first_paragraph(obj) -> str:
    """First paragraph of the docstring of ``obj`` on one line, or "" when it has none of its own."""
    doc = inspect.getdoc(obj)
    if not doc:
        return ""
    if inspect.isclass(obj) and dataclasses.is_dataclass(obj) and doc.startswith(obj.__name__ + "("):
        return ""  # the docstring that dataclasses generates is only the signature
    para = re.split(r"\n\s*\n", doc.strip(), maxsplit=1)[0]
    return " ".join(line.strip() for line in para.splitlines())


def fmt_default(value) -> str:
    if isinstance(value, PurePath):  # paths inside the checkout: relative, so the output does not depend on it
        try:
            return f"<repository>/{Path(value).resolve().relative_to(ROOT.resolve()).as_posix()}"
        except ValueError:
            return repr(str(value))
    text = repr(value)
    if " at 0x" in text:  # functions and objects: no memory addresses in the output
        text = getattr(value, "__name__", type(value).__name__)
    return text


def fmt_annotation(ann) -> str:
    if isinstance(ann, str):  # postponed annotations (from __future__ import annotations)
        return ann
    if isinstance(ann, type):
        return ann.__name__
    return inspect.formatannotation(ann)


def fmt_params(sig: inspect.Signature) -> tuple[list[str], str]:
    parts: list[str] = []
    star_done = False
    params = list(sig.parameters.values())
    for i, p in enumerate(params):
        if p.kind is p.KEYWORD_ONLY and not star_done:
            parts.append("*")
            star_done = True
        if p.kind is p.VAR_POSITIONAL:
            star_done = True
            text = "*" + p.name
        elif p.kind is p.VAR_KEYWORD:
            text = "**" + p.name
        else:
            text = p.name
        if p.annotation is not p.empty:
            text += ": " + fmt_annotation(p.annotation)
        if p.default is not p.empty:
            text += (" = " if p.annotation is not p.empty else "=") + fmt_default(p.default)
        parts.append(text)
        nxt = params[i + 1] if i + 1 < len(params) else None
        if p.kind is p.POSITIONAL_ONLY and (nxt is None or nxt.kind is not p.POSITIONAL_ONLY):
            parts.append("/")
    ret = "" if sig.return_annotation is sig.empty else " -> " + fmt_annotation(sig.return_annotation)
    return parts, ret


def signature(name: str, obj, drop_self: bool = False) -> str:
    """``name(params) -> ret`` (one parameter per line when it is long), or ``name(...)`` if unknown."""
    try:
        sig = inspect.signature(obj)
    except (TypeError, ValueError):
        return name + "(...)"
    if drop_self and sig.parameters:
        first = next(iter(sig.parameters))
        if first in ("self", "cls"):
            sig = sig.replace(parameters=list(sig.parameters.values())[1:])
    if inspect.isclass(obj):
        sig = sig.replace(return_annotation=inspect.Signature.empty)
    parts, ret = fmt_params(sig)
    one_line = f"{name}({', '.join(parts)}){ret}"
    if len(one_line) <= WRAP:
        return one_line
    body = "".join(f"    {p},\n" for p in parts)
    return f"{name}(\n{body}){ret}"


def source_line(obj) -> int:
    try:
        return inspect.getsourcelines(obj)[1]
    except (OSError, TypeError):
        return 10**9


def public_members(mod):
    """Public classes and functions defined in ``mod`` (not imported into it), in source order."""
    out = []
    for name, obj in vars(mod).items():
        if name.startswith("_"):
            continue
        if (inspect.isclass(obj) or inspect.isfunction(obj)) and getattr(obj, "__module__", None) == mod.__name__:
            out.append((name, obj))
    return sorted(out, key=lambda item: (source_line(item[1]), item[0]))


def public_constants(mod) -> list[tuple[str, str]]:
    """Upper-case module constants assigned in ``mod`` itself: numbers, strings and tuples.

    Values that hold a URL are left out: they are download configuration, not API.
    """
    try:
        src = inspect.getsource(mod)
    except (OSError, TypeError):
        return []
    out = []
    for name, value in vars(mod).items():
        if name.startswith("_") or not name.isupper():
            continue
        if not re.search(rf"^{re.escape(name)}\s*[:=]", src, flags=re.M):
            continue
        if type(value).__module__ == "numpy" and hasattr(value, "item"):
            value = value.item()  # numpy scalar: print the plain number
        if isinstance(value, (bool, int, float, str, tuple)) and "://" not in repr(value):
            out.append((name, repr(value)))
    return out


def class_members(cls):
    """Public methods, class methods, static methods and properties defined in ``cls`` itself."""
    out = []
    for name, raw in vars(cls).items():
        if name.startswith("_"):
            continue
        if isinstance(raw, property):
            out.append((name, raw.fget, "property"))
        elif isinstance(raw, (classmethod, staticmethod)):
            out.append((name, raw.__func__, type(raw).__name__))
        elif inspect.isfunction(raw):
            out.append((name, raw, "method"))
    return sorted(out, key=lambda item: (source_line(item[1]), item[0]))


def render_class(name: str, cls) -> list[str]:
    lines = [f"### class `{name}`", "", "```python", signature(name, cls), "```", ""]
    bases = [b for b in cls.__bases__ if b is not object]
    if bases:
        names = ", ".join(f"`{b.__module__}.{b.__qualname__}`" if b.__module__ != "builtins" else f"`{b.__name__}`"
                          for b in bases)
        lines += [f"Subclass of {names}.", ""]
    doc = first_paragraph(cls)
    if doc:
        lines += [doc, ""]
    members = class_members(cls)
    if members:
        lines += ["Members:", ""]
        for mname, fn, kind in members:
            if kind == "property":
                head = f"`{mname}` (property)"
            else:
                head = f"`{signature(mname, fn, drop_self=(kind != 'staticmethod'))}`"
                if "\n" in head:  # keep list items on one line
                    head = "`" + " ".join(head.strip("`").split()).replace("( ", "(").replace(", )", ")") + "`"
                if kind != "method":
                    head += f" ({kind})"
            mdoc = first_paragraph(fn)
            lines.append(f"- {head}" + (f": {mdoc}" if mdoc else ""))
        lines.append("")
    return lines


def render_function(name: str, fn) -> list[str]:
    lines = [f"### `{name}`", "", "```python", signature(name, fn), "```", ""]
    doc = first_paragraph(fn)
    if doc:
        lines += [doc, ""]
    return lines


def anchor(title: str) -> str:
    """GitHub heading anchor: lower case, punctuation other than '-' and '_' removed, spaces to '-'."""
    return re.sub(r"[^\w\- ]", "", title.lower()).replace(" ", "-")


def render_module(mod) -> list[str]:
    lines = [f"## {mod.__name__}", ""]
    doc = first_paragraph(mod)
    if doc:
        lines += [doc, ""]
    consts = public_constants(mod)
    if consts:
        lines += ["Constants:", ""] + [f"- `{n} = {v}`" for n, v in consts] + [""]
    for name, obj in public_members(mod):
        lines += render_class(name, obj) if inspect.isclass(obj) else render_function(name, obj)
    return lines


def where_defined(pkg, name: str, modules: dict) -> str:
    obj = getattr(pkg, name)
    mod = getattr(obj, "__module__", None)
    if mod and mod.startswith(PACKAGE + ".") and (inspect.isclass(obj) or inspect.isfunction(obj)):
        return mod
    for sub in modules.values():  # constants: the module that assigns them
        if any(n == name for n, _ in public_constants(sub)):
            return sub.__name__
    return PACKAGE


def kind_of(obj) -> str:
    if inspect.isclass(obj):
        return "class"
    if inspect.isfunction(obj):
        return "function"
    return "constant"


def build() -> str:
    sys.path.insert(0, str(SRC))
    pkg = importlib.import_module(PACKAGE)
    if not Path(pkg.__file__).resolve().is_relative_to(SRC.resolve()):
        raise SystemExit(f"imported {pkg.__file__}, expected the package under {SRC}")
    names = sorted(i.name for i in pkgutil.iter_modules(pkg.__path__) if not i.name.startswith("_"))
    modules = {n: importlib.import_module(f"{PACKAGE}.{n}") for n in names}
    documented = [n for n in names if public_members(modules[n]) or public_constants(modules[n])]
    empty = [n for n in names if n not in documented]
    core = [n for n in CORE if n in documented]
    research = [n for n in documented if n not in CORE]

    lines = [
        "<!-- Generated by scripts/build_api_docs.py from the docstrings. Edit the docstrings, then rerun"
        " the script; do not edit this file by hand. -->",
        "",
        "# API reference",
        "",
        f"Package `{PACKAGE}`, version {pkg.__version__}. "
        + first_paragraph(pkg),
        "",
        "Each entry gives the signature and the first paragraph of the docstring. The docstrings in the"
        " source, and the theory notes in `docs/` that they cite, give the full statements and conditions.",
        "",
        "## Top-level imports",
        "",
        f"`from {PACKAGE} import <name>` works for every name below. Each one is also importable from the"
        " module in the last column.",
        "",
        "| Name | Kind | Defined in |",
        "|---|---|---|",
    ]
    for name in pkg.__all__:
        obj = getattr(pkg, name)
        lines.append(f"| `{name}` | {kind_of(obj)} | `{where_defined(pkg, name, modules)}` |")
    lines += ["", "## Modules", ""]
    lines += ["Core modules:", ""] + [f"- [`{PACKAGE}.{n}`](#{anchor(f'{PACKAGE}.{n}')})" for n in core] + [""]
    lines += ["Research modules (used by the experiments and the theory notes):", ""]
    lines += [f"- [`{PACKAGE}.{n}`](#{anchor(f'{PACKAGE}.{n}')})" for n in research] + [""]
    if empty:
        lines += ["Placeholder modules with no public classes or functions: "
                  + ", ".join(f"`{PACKAGE}.{n}`" for n in empty) + ".", ""]
    for n in core + research:
        lines += render_module(modules[n])
    text = "\n".join(lines).rstrip() + "\n"
    return re.sub(r"\n{3,}", "\n\n", text)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build docs/api.md from the docstrings of the avmonitor package.")
    ap.add_argument("-o", "--output", type=Path, default=OUT, help="output file (default: docs/api.md)")
    ap.add_argument("--check", action="store_true", help="do not write; exit 1 if the output file is out of date")
    a = ap.parse_args(argv)
    text = build()
    if a.check:
        current = a.output.read_text(encoding="utf-8") if a.output.exists() else ""
        if current != text:
            print(f"{a.output} is out of date; run python scripts/build_api_docs.py")
            return 1
        print(f"{a.output} is up to date")
        return 0
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(text, encoding="utf-8")
    print(f"wrote {a.output} ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
