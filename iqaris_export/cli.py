"""Command-line interface: `python -m iqaris_export {explore|stats|list-props|plot|export} ...`."""
from __future__ import annotations
import argparse
import json

from . import (config, connect, Selection, Predicate, describe, coverage, column_stats,
               list_props, export, plot_distribution, plot_coverage, property_docs, info as _info)
from .writers import ALL_FORMATS

# Organisation banner shown atop `list-props` output and its --help.
_LISTPROPS_HEADER = (
    "Three analyses (QTAIM topology · ab initio IQA · IQA-SQM) × seven levels of theory",
    "(M06-2X, PBE0, PM7, PM6, PM6-D3H4, PM6-D3H4X, PM6-ORG). QTAIM exists at M06-2X and PBE0;",
    "ab initio IQA at PBE0 only; IQA-SQM at the five SQM levels. VIEW = the table holding the column.",
)


def _add_selection(p, root=True, level=False):
    g = p.add_argument_group("selection")
    if level:                       # count/explore: no served column, but predicates still need a level
        g.add_argument("--level", choices=["m062x", "pbe0", "dft"],
                       help="level of theory that UNSCOPED --predicate terms on shared QTAIM columns "
                            "and sanity rows bind to: m062x (M06-2X, default) or pbe0")
    if root:
        g.add_argument("--root", default=config.DB_ROOT, help="database root")
    g.add_argument("--elements", help="comma list; structures containing ALL of these (e.g. S,Cl)")
    g.add_argument("--elements-only", help="comma list; composition subset of these (e.g. C,H,O,N)")
    g.add_argument("--family", help="comma list of source families")
    g.add_argument("--natoms", help="atom-count range LO:HI (inclusive)")
    g.add_argument("--formula", help="exact Hill formula, or a LIKE pattern with %%")
    g.add_argument("--tier", choices=sorted(["sane_qtaim", "sane_mopac_all5", "sane_all",
                                             "recommended"]))
    g.add_argument("--ids", help="comma list of structure_ids")
    g.add_argument("--ids-file", help="file with one structure_id per line")
    g.add_argument("--one-per-molecule", action="store_true",
                   help="keep one conformer per molecule (inchikey_block1)")
    g.add_argument("--molecule-key", default="inchikey_block1",
                   choices=["inchikey_block1", "inchikey", "formula"])
    g.add_argument("--molecule-pick", default="min_id",
                   choices=["min_id", "max_natoms", "recommended"])
    g.add_argument("--predicate", action="append", dest="predicates", metavar="EXPR",
                   help='value filter, repeatable: "COL OP VAL [@METHOD]". COL = a sanity '
                        'residual/flag (charge_resid, L_atom_max, etot_recon_resid, '
                        'halogen_quarantine, ...), a manifest column (natoms, charge, formula, '
                        'source_family, ...), or any property column; OP = < <= > >= == != '
                        'between in is_true is_false. Examples: "charge_resid < 1e-3 @QTAIM", '
                        '"etot_recon_resid < 0.01 @PM7", "natoms between 2 30", '
                        '"halogen_quarantine is_false", "source_family in qm9,gdb". An optional '
                        '@LEVEL scopes the term: @QTAIM (= @M06-2X) or @PBE0 for sanity rows and '
                        'shared *_QTAIM columns, @PM7 / @PM6 / ... for *_SQM columns. Without a '
                        'scope the subcommand\'s --level (and --method for SQM columns) applies; '
                        'defaults M06-2X and PM7.')
    g.add_argument("--sample", type=int, help="reproducible random sample of N structures")
    g.add_argument("--seed", type=int, default=config.DEFAULT_SEED)
    g.add_argument("--limit", type=int, help="cap the number of structures")


def _coerce(x):
    """Best-effort literal: int, then float, then true/false, else the raw string."""
    try:
        return int(x)
    except ValueError:
        pass
    try:
        return float(x)
    except ValueError:
        pass
    low = x.lower()
    if low in ("true", "false"):
        return low == "true"
    return x


_OP_ALIAS = {"=": "==", "lt": "<", "le": "<=", "gt": ">", "ge": ">=", "eq": "==", "ne": "!=",
             "~": "contains", "!~": "not_contains"}


def _side_from_level(level):
    """Map a user `--level` token to the internal registry side for a shared QTAIM/QCT column:
    m062x/dft -> 'qtaim' (M06-2X), pbe0 -> 'iqadft' (PBE0).  None leaves the default (M06-2X)."""
    if not level:
        return None
    low = level.lower()
    if low in ("m062x", "m06-2x", "dft", "qtaim"):
        return "qtaim"
    if low in ("pbe0", "iqadft"):
        return "iqadft"
    return None


def _parse_predicate(expr):
    """Parse "COL OP VAL [@METHOD]" into a Predicate (raises SystemExit on a malformed spec)."""
    s = expr.strip()
    method = None
    if "@" in s:
        s, method = s.rsplit("@", 1)
        s, method = s.strip(), method.strip()
    toks = s.split()
    if len(toks) < 2:
        raise SystemExit(f"error: bad --predicate {expr!r}: expected 'COL OP VAL [@METHOD]'")
    col, op = toks[0], _OP_ALIAS.get(toks[1], toks[1])
    vals = [v for tok in toks[2:] for v in tok.split(",") if v]
    if op in ("is_true", "is_false"):
        value = None
    elif op == "between":
        if len(vals) != 2:
            raise SystemExit(f"error: bad --predicate {expr!r}: 'between' needs two values")
        value = (_coerce(vals[0]), _coerce(vals[1]))
    elif op == "in":
        if not vals:
            raise SystemExit(f"error: bad --predicate {expr!r}: 'in' needs a value list")
        value = tuple(_coerce(v) for v in vals)
    else:
        if len(vals) != 1:
            raise SystemExit(f"error: bad --predicate {expr!r}: '{op}' needs one value")
        value = _coerce(vals[0])
    try:
        return Predicate(col, op, value, method=method or None)
    except ValueError as e:
        raise SystemExit(f"error: bad --predicate {expr!r}: {e}")


def _selection(a):
    natoms = None
    if a.natoms:
        lo, hi = a.natoms.split(":")
        natoms = (int(lo), int(hi))
    ids = None
    if a.ids:
        ids = tuple(x.strip() for x in a.ids.split(",") if x.strip())
    elif a.ids_file:
        ids = tuple(l.strip() for l in open(a.ids_file) if l.strip())
    csv = lambda s: tuple(x.strip() for x in s.split(",") if x.strip())
    preds = getattr(a, "predicates", None)
    predicates = tuple(_parse_predicate(p) for p in preds) if preds else None
    # Bind unscoped predicates to the subcommand's --level (M06-2X/PBE0) and --method / SQM --level,
    # so `--level pbe0` filters on the PBE0 reading rather than silently at M06-2X.
    level = _side_from_level(getattr(a, "level", None))
    sqm = None
    for tok in (getattr(a, "method", None), getattr(a, "level", None)):
        if tok and any(str(tok).upper() == m.upper() for m in config.SQM_METHODS):
            sqm = next(m for m in config.SQM_METHODS if m.upper() == str(tok).upper())
            break
    return Selection(
        elements=csv(a.elements) if a.elements else None,
        elements_only=csv(a.elements_only) if a.elements_only else None,
        source_family=csv(a.family) if a.family else None,
        natoms=natoms, formula=a.formula, tier=a.tier, ids=ids, predicates=predicates,
        one_per_molecule=a.one_per_molecule, molecule_key=a.molecule_key,
        molecule_pick=a.molecule_pick, sample=a.sample, seed=a.seed, limit=a.limit,
        level=level, sqm_method=sqm,
    )


def main(argv=None):
    ap = argparse.ArgumentParser(prog="iqaris_export",
                                 description="Explore, visualize and export the IQARIS database.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    pc = sub.add_parser("count", help="count structures surviving a selection (fast preview)")
    _add_selection(pc, level=True)

    pe = sub.add_parser("explore", help="summarize a selection (counts, coverage)")
    _add_selection(pe, level=True)

    ps = sub.add_parser("stats", help="statistics for one column over a selection")
    _add_selection(ps)
    ps.add_argument("--column", required=True)
    ps.add_argument("--method", help="SQM method for level-aware columns (default PM7)")
    ps.add_argument("--level", choices=["m062x", "pbe0", "dft"],
                    help="for a name on both DFT sides, pick M06-2X (m062x, default) or PBE0 (pbe0)")

    pl = sub.add_parser("list-props", help="list exportable properties",
                        epilog="\n".join(_LISTPROPS_HEADER),
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    pl.add_argument("--root", default=config.DB_ROOT, help="database root")
    pl.add_argument("--prop-family", dest="prop_family", choices=list(config.FAMILIES))
    pl.add_argument("--level",
                    choices=["m062x", "pbe0", "pm7", "pm6", "pm6-d3h4", "pm6-d3h4x", "pm6-org", "all"],
                    default=None,
                    help="filter by LEVEL of theory -- one of the seven actual levels: m062x (M06-2X), "
                         "pbe0, or a semiempirical Hamiltonian pm7 / pm6 / pm6-d3h4 / pm6-d3h4x / "
                         "pm6-org. Default: all levels. (For 'all SQM levels' use --analysis iqa-sqm.)")
    pl.add_argument("--analysis", choices=["qtaim", "iqa", "iqa-sqm"], default=None,
                    help="filter by ANALYSIS type: qtaim (density topology), iqa (ab initio energy "
                         "decomposition, PBE0 only), or iqa-sqm (the semiempirical energy "
                         "decomposition, at all five Hamiltonians). Composable with --level, e.g. "
                         "`--level pbe0 --analysis iqa` = only the ab initio IQA columns.")
    pl.add_argument("--desc", action="store_true", help="also print each property's description")

    pi = sub.add_parser("info", help="describe one property (name, symbol, unit, type, equation)")
    pi.add_argument("key", help="a database column or public Key, e.g. q_QTAIM, DI_QTAIM, Einter_SQM")
    pi.add_argument("--root", default=config.DB_ROOT, help="database root (for the exact unit)")
    pi.add_argument("--no-db", action="store_true",
                    help="answer from the shipped dictionary only, without opening the database")

    pgl = sub.add_parser("glossary", help="render the whole property dictionary as a document")
    pgl.add_argument("--format", default="markdown", choices=["markdown", "html", "json"])
    pgl.add_argument("--out", help="write here instead of stdout")

    pp = sub.add_parser("plot", help="distribution of a column, or coverage")
    _add_selection(pp)
    pp.add_argument("--column", help="column to plot")
    pp.add_argument("--coverage", action="store_true", help="plot elemental coverage instead")
    pp.add_argument("--kind", default="hist", choices=["hist", "violin"])
    pp.add_argument("--by", choices=["element"])
    pp.add_argument("--method")
    pp.add_argument("--level", choices=["m062x", "pbe0", "dft"],
                    help="for a name on both DFT sides, pick M06-2X (m062x, default) or PBE0 (pbe0)")
    pp.add_argument("--out")
    pp.add_argument("--show", action="store_true",
                    help="display interactively instead of writing to a file (needs a working "
                         "GUI backend + display); combine with --out to do both")

    pj = sub.add_parser("project", help="render a per-atom property onto a structure's 3D geometry")
    pj.add_argument("--root", default=config.DB_ROOT, help="database root")
    pj.add_argument("--sid", help="structure_id to render (or the starting structure with --browse)")
    pj.add_argument("--browse", action="store_true",
                    help="interactive only: browse the filtered selection (Prev/Next + slider) "
                         "instead of a single structure; combine with the selection flags below")
    pj.add_argument("--prop", help="per-atom scalar column to color by (e.g. q_QTAIM, q_SQM)")
    pj.add_argument("--level", choices=["m062x", "pbe0", "dft"], default="m062x",
                    help="render the M06-2X (m062x, default) or PBE0 (pbe0) reading of a shared "
                         "QTAIM/QCT property")
    pj.add_argument("--method", default=config.SQM_METHODS[0],
                    help="SQM level for level-aware props / fragment engagement")
    pj.add_argument("--graph", action="store_true",
                    help="draw the QTAIM molecular graph (bond paths + BCP/RCP/CCP dots)")
    pj.add_argument("--paths", action="store_true",
                    help="draw bond paths only (no critical-point dots)")
    pj.add_argument("--ias", action="store_true",
                    help="draw interatomic (zero-flux) surfaces as area-scaled disks")
    pj.add_argument("--fragment", action="store_true",
                    help="split host+guest and color by their summed host<->guest interaction when "
                         "the structure fragments (tints the two fragment skeletons distinctly)")
    pj.add_argument("--fragment-prop", dest="fragment_prop",
                    help="interface channel aggregated in --fragment mode: Einter_SQM (default), "
                         "Eelstat_SQM, Eexchange_SQM, Eresonance_SQM, Edisp_SQM, or DI_QTAIM")
    pj.add_argument("--row-atom", type=int, dest="row_atom", metavar="A",
                    help="color every atom B by one ROW of a matrix for atom A (1-based atom_index): "
                         "the delocalization index delta(A,B) (log10 scale) or, with --row-kind iqa, "
                         "the interaction energy E_int(A,B) in kcal/mol; A is highlighted")
    pj.add_argument("--row-kind", dest="row_kind", choices=["ldm", "iqa"], default="ldm",
                    help="ldm (default): delta(A,B) at --level m062x|pbe0; iqa: E_int(A,B), ab initio "
                         "at --level pbe0, otherwise at the semiempirical --method")
    pj.add_argument("--zoom", type=float, help="Jmol camera zoom percent")
    pj.add_argument("--size", help="render/window size WxH in px, e.g. 1200x900")
    pj.add_argument("--out", help="output directory (default exports/projections)")
    pj.add_argument("--interactive", action="store_true",
                    help="open a rotatable PyVista/VTK window instead of a static Jmol PNG")
    pj.add_argument("--save", help="with --interactive: also write a screenshot here; "
                    "alone: headless PyVista screenshot to this file (PNG)")
    pj.add_argument("--points-only", action="store_true",
                    help="PyVista lightweight mode (atoms only) for very large systems")
    # --- image controls (Style; apply to BOTH the PyVista and Jmol backends) ---------
    g = pj.add_argument_group("image controls")
    g.add_argument("--cmap", help="colormap name (mako, viridis, RdBu, coolwarm, ...)")
    g.add_argument("--vmin", type=float, help="lower color-scale bound")
    g.add_argument("--vmax", type=float, help="upper color-scale bound")
    g.add_argument("--vcenter", type=float, help="center value for a diverging colormap")
    g.add_argument("--atom-scale", type=float, default=1.0, help="atom sphere size multiplier")
    g.add_argument("--bond-scale", type=float, default=1.0, help="bond thickness multiplier")
    g.add_argument("--bonds", choices=["qtaim", "di", "distance"], default="di",
                   help="bond perception: qtaim (any BCP) | di (DI>cutoff) | distance (radii)")
    g.add_argument("--bond-cutoff", type=float,
                   help="di: DI threshold (default 0.30); distance: radii tolerance (default 1.20)")
    g.add_argument("--no-h", action="store_true", help="hide hydrogens (and bonds to them)")
    g.add_argument("--labels", action="store_true", help="draw element labels at each atom")
    g.add_argument("--values", action="store_true",
                   help="annotate the projected property's value in 3D (on atoms / bonds / CPs; "
                        "interactive viewer + PyVista screenshots)")
    g.add_argument("--background", default="white", help="background color (e.g. white, black)")
    g.add_argument("--sphere-res", type=int, help="sphere facet resolution (PyVista)")
    _add_selection(pj, root=False)     # for --browse: filter the set to step through in the window

    px = sub.add_parser("export", help="export a selection to an ML format")
    _add_selection(px)
    px.add_argument("--props", required=True,
                    help="comma tokens: family shorthands / bundles / column names")
    px.add_argument("--level", default="dft",
                    help="dft | pbe0 | PM7 | PM6 | PM6-D3H4 | PM6-D3H4X | PM6-ORG | all-sqm | all")
    px.add_argument("--format", default="extxyz", choices=list(ALL_FORMATS))
    px.add_argument("--out", required=True)
    px.add_argument("--energy-target", default=config.DEFAULT_ENERGY_TARGET,
                    choices=list(config.ENERGY_TARGETS))
    px.add_argument("--with-bonus", action="store_true",
                    help="family shorthands expand to ALL columns, not just curated defaults")
    px.add_argument("--include-problematic", action="store_true", dest="include_problematic",
                    help="include the quality-flagged PBE0 records (opt-in `_all` views) and emit "
                         "pbe0_qc_status; default is clean data only")
    px.add_argument("--batch-size", type=int, default=200)
    px.add_argument("--force-stub", action="store_true",
                    help="write all-zero forces (placeholder; DB has no real forces)")

    pm = sub.add_parser(
        "matrix", help="localization-delocalization (LDM) or IQA energy matrix of structures",
        description="Atom-condensed matrices built from the stored atomic and pair columns (every "
                    "atom pair is stored, bonded or not). ldm: L[A,A]=LI_QTAIM, L[A,B]=DI_QTAIM/2, "
                    "rows sum to N_QTAIM. iqa: E[A,A]=E_intra(A), E[A,B]=E_int(A,B)/2, rows sum to "
                    "the additive atomic energy (E_IQA at PBE0, Eadd_SQM at an SQM level).")
    _add_selection(pm)
    pm.add_argument("--sid", help="one structure_id (otherwise the selection flags pick structures)")
    pm.add_argument("--kind", default="ldm", choices=["ldm", "iqa"],
                    help="ldm (QTAIM localization-delocalization matrix, default) or iqa (energy)")
    pm.add_argument("--level",
                    help="ldm: m062x (default) | pbe0;  iqa: pbe0 (default, ab initio) | PM7 | PM6 | "
                         "PM6-D3H4 | PM6-D3H4X | PM6-ORG")
    pm.add_argument("--fragments",
                    help='with --sid: condense over atom groups (1-based atom_index), e.g. "1-6;7,8" '
                         'or "ring=1-6;CO=7,8"; atoms in no group form a final "rest" group')
    pm.add_argument("--dispersion", action="store_true",
                    help="iqa at an SQM level: add the pairwise empirical dispersion (Edisp_SQM/2) "
                         "to the off-diagonal (excluded by default, as in Eadd_SQM)")
    pm.add_argument("--include-problematic", action="store_true", dest="include_problematic",
                    help="PBE0: read the quality-flagged records too (opt-in `_all` views)")
    pm.add_argument("--format", choices=["table", "json", "csv", "npy", "npz"],
                    help="table (readable, the default on screen), json (the default with --out), "
                         "csv / npy (one structure), npz (one or many)")
    pm.add_argument("--out", help="output file (required for csv, npy and npz)")
    pm.add_argument("--plot", metavar="IMAGE",
                    help="with --sid: also save a heatmap image (png, pdf or svg by extension); the "
                         "off-diagonal entries are colored and the diagonal is printed")

    a = ap.parse_args(argv)

    # Commands that do not require (or only optionally use) a database connection.
    if a.cmd == "glossary":
        return _do_glossary(a)
    if a.cmd == "info":
        return _do_info(a)

    con = connect(a.root)
    try:
        _dispatch(a, con)
    except (ValueError, RuntimeError) as e:
        raise SystemExit(f"error: {e}")


def _do_glossary(a):
    if a.format == "json":
        text = json.dumps(property_docs.as_records(), indent=2, default=str)
    else:
        text = property_docs.glossary(a.format)
    if a.out:
        with open(a.out, "w") as fh:
            fh.write(text)
        print("wrote", a.out)
    else:
        print(text)


def _do_info(a):
    con = None
    if not a.no_db:
        try:
            con = connect(a.root)
        except Exception:
            con = None                       # dictionary-only fallback
    _info(a.key, con=con)


def _style_from_args(a):
    """Build a :class:`~iqaris_export.style.Style` from the `project` subcommand's image flags."""
    from .style import Style
    kw = dict(cmap=a.cmap, vmin=a.vmin, vmax=a.vmax, vcenter=a.vcenter,
              atom_scale=a.atom_scale, bond_scale=a.bond_scale, bonds=a.bonds,
              bond_cutoff=a.bond_cutoff, show_h=not a.no_h, labels=a.labels,
              show_values=a.values, background=a.background)
    if a.sphere_res is not None:
        kw["sphere_res"] = a.sphere_res
    if a.size:
        w, h = a.size.lower().split("x")
        kw["window_size"] = (int(w), int(h))
    return Style(**kw)


def _dispatch(a, con):
    if a.cmd == "count":
        sel = _selection(a)
        print(json.dumps({"surviving": sel.count(con)}, indent=2, default=str))
    elif a.cmd == "explore":
        sel = _selection(a)
        print(json.dumps({"describe": describe(con, sel),
                          "coverage": coverage(con, sel)}, indent=2, default=str))
    elif a.cmd == "stats":
        sel = _selection(a)
        print(json.dumps(column_stats(con, a.column, sel, a.method,
                                      side=_side_from_level(getattr(a, "level", None))),
                         indent=2, default=str))
    elif a.cmd == "list-props":
        rows = list_props(con, getattr(a, "prop_family", None),
                          level=getattr(a, "level", None), analysis=getattr(a, "analysis", None))
        # Two orthogonal axes: ANALYSIS (what the column reports) and LEVEL (its level of theory).
        # PBE0 contributes both QTAIM and ab initio IQA columns; filter either with --analysis/--level.
        for line in _LISTPROPS_HEADER:
            print(line)
        print()
        last_col = "DESCRIPTION" if getattr(a, "desc", False) else "VIEW"
        header = (f"{'FAMILY':10s} {'ANALYSIS':9s} {'LEVEL':10s} {'NAME':28s} {'UNIT':16s} {last_col}")
        print(header)
        print("-" * max(len(header), 80))
        for r in rows:
            line = (f"{r['family']:10s} {r['analysis']:9s} {r['level']:10s} {r['name']:28s} "
                    f"{r['unit']:16s}")
            if getattr(a, "desc", False):
                line += f" {r.get('description', '')}"
            else:
                line += f" {r['view']}"
            print(line)
        # count distinct public property NAMES too (a *_QTAIM name listed at both M06-2X and PBE0 is
        # ONE property queryable at two levels, not two) so the tally isn't mistaken for duplication.
        print(f"\n{len(rows)} rows  ({len({r['name'] for r in rows})} distinct properties across levels)")
    elif a.cmd == "plot":
        sel = _selection(a)
        if a.coverage:
            path = plot_coverage(con, sel, out=a.out, show=a.show)
        elif a.column:
            path = plot_distribution(con, a.column, sel, kind=a.kind, by=a.by,
                                     method=a.method, out=a.out, show=a.show,
                                     side=_side_from_level(getattr(a, "level", None)))
        else:
            raise SystemExit("plot: pass --column COL or --coverage")
        print("wrote", path) if path else print("displayed (not saved)")
    elif a.cmd == "project":
        style = _style_from_args(a)
        pbe0 = _side_from_level(getattr(a, "level", None)) == "iqadft"
        row = {}
        if a.row_atom is not None:
            if a.browse:
                raise SystemExit("project: --row-atom picks one atom of one structure; it cannot be "
                                 "combined with --browse (double-click an atom in the window instead)")
            row = dict(row_atom=a.row_atom, row_kind=a.row_kind,
                       row_level=("pbe0" if pbe0 else "m062x") if a.row_kind == "ldm"
                       else ("pbe0" if pbe0 else a.method))
        if a.interactive or a.save:                     # interactive PyVista / VTK path
            from .viz3d import view as _view
            sids = None
            if a.browse:
                if not a.interactive:                    # browsing needs the live window
                    raise SystemExit("project --browse requires --interactive (Prev/Next has no "
                                     "meaning in a headless --save screenshot)")
                sids = _selection(a).resolve(con)
                if not sids:
                    raise SystemExit("browse: the selection matched 0 structures")
            elif not a.sid:
                raise SystemExit("project: pass --sid, or --browse with selection filters")
            img = _view(a.sid, sids=sids, prop=a.prop, method=a.method, graph=a.graph,
                        paths=a.paths, ias=a.ias, fragment=a.fragment, fragment_prop=a.fragment_prop,
                        con=con, style=style, interactive=a.interactive, save=a.save,
                        points_only=a.points_only, pbe0=pbe0, **row)
            print(json.dumps({"backend": "pyvista", "sid": a.sid or (sids[0] if sids else None),
                              "browse_n": (len(sids) if sids else 1),
                              "interactive": a.interactive, "saved": img}, indent=2, default=str))
        else:                                           # static publication Jmol path (default)
            if not a.sid:
                raise SystemExit("project: --sid is required for the static Jmol render")
            from .project import project as _project, DEFAULT_VIEW
            view = dict(DEFAULT_VIEW, zoom=a.zoom) if a.zoom is not None else None
            size = None
            if a.size:
                w, h = a.size.lower().split("x")
                size = (int(w), int(h))
            png = _project(a.sid, prop=a.prop, method=a.method, graph=a.graph, paths=a.paths,
                           ias=a.ias, fragment=a.fragment, fragment_prop=a.fragment_prop,
                           out=a.out, con=con, view=view, size=size, style=style, pbe0=pbe0, **row)
            print(json.dumps({"backend": "jmol", "rendered": png, "sid": a.sid}, indent=2, default=str))
    elif a.cmd == "export":
        sel = _selection(a)
        stats = export(con, sel, props=a.props, level=a.level, fmt=a.format, out=a.out,
                       energy_target=a.energy_target, with_bonus=a.with_bonus,
                       batch_size=a.batch_size, force_stub=a.force_stub,
                       include_problematic=getattr(a, "include_problematic", False))
        print(json.dumps({"done": stats, "out": a.out}, indent=2, default=str))
    elif a.cmd == "matrix":
        _do_matrix(a, con)


_SELECTION_FLAGS = ("elements", "elements_only", "family", "natoms", "formula", "tier", "ids",
                    "ids_file", "predicates", "sample", "limit", "one_per_molecule")


def _do_matrix(a, con):
    """`matrix` subcommand: one structure (--sid) or a selection, written as json/csv/npy/npz."""
    import numpy as np
    from .matrices import atom_matrices, condense, parse_groups, _norm_level
    level = a.level or ("m062x" if a.kind == "ldm" else "pbe0")
    if a.format is None:
        a.format = "json" if a.out else "table"
    if a.sid:
        ids = [a.sid]
    else:
        if not any(getattr(a, f, None) for f in _SELECTION_FLAGS):
            raise SystemExit("matrix: pass --sid ID, or selection flags (e.g. --ids, --family, "
                             "--natoms, --limit) to pick the structures")
        ids = _selection(a).resolve(con)
    if a.fragments and not a.sid:
        raise SystemExit("matrix: --fragments needs a single --sid (atom indices are per structure)")
    if a.plot and not a.sid:
        raise SystemExit("matrix: --plot needs a single --sid (one image per structure)")
    if a.format in ("csv", "npy", "npz") and not a.out:
        raise SystemExit(f"matrix: --format {a.format} needs --out FILE")
    if a.format in ("csv", "npy") and len(ids) != 1:
        raise SystemExit(f"matrix: --format {a.format} writes one structure; use npz or json for "
                         f"{len(ids)} structures")
    mats, missing = atom_matrices(con, ids, kind=a.kind, level=level,
                                  include_problematic=a.include_problematic,
                                  dispersion=a.dispersion)
    if a.sid and missing:
        if not con.execute("SELECT count(*) FROM manifest WHERE structure_id = ?",
                           [a.sid]).fetchone()[0]:
            raise SystemExit(f"error: structure_id {a.sid!r} not found")
        hint = ""
        if _norm_level(level) == "PBE0" and not a.include_problematic:
            hint = " (no clean PBE0 record; --include-problematic reads the quality-flagged ones)"
        raise SystemExit(f"error: structure_id {a.sid!r} has no {_norm_level(level)} record{hint}")
    if not mats and a.format not in ("json", "table"):
        raise SystemExit(f"matrix: none of the {len(ids)} selected structures has a "
                         f"{_norm_level(level)} record; nothing written")
    if a.fragments:
        groups, names = parse_groups(a.fragments)
        m = mats[a.sid]
        mats = {a.sid: condense(m, groups, names)}
    if a.plot:
        import sys
        from .matrices import plot_matrix
        try:
            plot_matrix(mats[a.sid], out=a.plot)
        except ImportError as e:
            raise SystemExit(f"matrix: --plot needs matplotlib (pip install 'iqaris-explorer[plots]'): {e}")
        print(f"image: {a.plot}", file=sys.stderr)
    devs = [m.max_rowsum_deviation() for m in mats.values()]
    devs = [d for d in devs if np.isfinite(d)]
    summary = {"kind": a.kind, "level": _norm_level(level), "n_structures": len(mats),
               "missing": missing, "max_rowsum_deviation": max(devs) if devs else None}
    if a.format == "table":
        blocks = [m.to_text() for m in mats.values()]
        if missing:
            blocks.append(f"no {_norm_level(level)} record for {len(missing)} structure(s): "
                          + ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else ""))
        text = "\n\n".join(blocks)
        if a.out:
            with open(a.out, "w") as fh:
                fh.write(text + "\n")
            print(json.dumps(dict(summary, out=a.out), indent=2, default=str))
        else:
            print(text)
        return
    if a.format == "json":
        payload = (next(iter(mats.values())).as_dict() if (a.sid and mats)
                   else dict(summary, matrices=[m.as_dict() for m in mats.values()]))
        text = json.dumps(payload, indent=2, default=str)
        if a.out:
            with open(a.out, "w") as fh:
                fh.write(text)
            print(json.dumps(dict(summary, out=a.out), indent=2, default=str))
        else:
            print(text)
        return
    m0 = next(iter(mats.values())) if mats else None
    if a.format == "csv":
        m0.to_dataframe().to_csv(a.out)
    elif a.format == "npy":
        np.save(a.out, m0.values)
    else:                                                   # npz: one array + labels per structure
        arrays = {}
        for sid, m in mats.items():
            arrays[sid] = m.values
            arrays[f"{sid}__labels"] = np.array(m.labels)
        arrays["__meta__"] = np.array(json.dumps(
            {"kind": a.kind, "level": _norm_level(level),
             "unit": m0.unit if m0 else None, "structures": list(mats), "missing": missing}))
        np.savez_compressed(a.out, **arrays)
    print(json.dumps(dict(summary, out=a.out), indent=2, default=str))


if __name__ == "__main__":
    main()
