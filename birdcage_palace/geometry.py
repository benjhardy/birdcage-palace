"""Gmsh (OpenCASCADE) model of a birdcage coil for Palace.

Conductors are zero-thickness PEC sheets on the coil cylinder (rungs and end rings as curved strips). Every capacitor is a
small curved patch filling a gap in a rung (low-pass, band-pass) or an end-ring segment (high-pass, band-pass); Palace
treats each patch as a lumped element (LumpedPort with C). The domain is the inside of the RF shield (a closed PEC can by
default) or, unshielded, an air cylinder with a first-order absorbing outer boundary. An optional phantom is a
cylinder of lossy dielectric.

Mesh length unit is the metre (Palace Model.L0 = 1).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import gmsh

from .design import CoilDesign

ATTR_AIR, ATTR_PHANTOM = 1, 2
ATTR_COIL, ATTR_SHIELD, ATTR_OUTER = 10, 11, 12
ATTR_SYM = {"x": 13, "y": 14, "z": 15}   # symmetry planes x=0, y=0, z=0
ATTR_LEG_CAP0, ATTR_RING_TOP0, ATTR_RING_BOT0 = 100, 200, 300


def leg_angles(N: int, phi0: float = 0.0) -> list[float]:
    return [phi0 + 2 * math.pi * n / N for n in range(N)]


def _strip(occ, R, a0, a1, z0, z1):
    """Curved strip on the cylinder of radius R, angles a0..a1 (rad), z0..z1."""
    if a1 - a0 >= 2 * math.pi - 1e-12:
        c = occ.addCircle(0, 0, z0, R)
    else:
        c = occ.addCircle(0, 0, z0, R, -1, a0, a1)
    ext = occ.extrude([(1, c)], 0, 0, z1 - z0)
    return [t for d, t in ext if d == 2][0]


def _flat_gap(occ, R, a0, a1, z0, z1):
    """Planar rectangle between the vertical lines at angles a0 and a1 on the cylinder of radius R."""
    p = [occ.addPoint(R * math.cos(a), R * math.sin(a), z) for a, z in ((a0, z0), (a1, z0), (a1, z1), (a0, z1))]
    ln = [occ.addLine(p[i], p[(i + 1) % 4]) for i in range(4)]
    return occ.addPlaneSurface([occ.addCurveLoop(ln)])


def build(des: CoilDesign, out_dir: Path, name: str = "birdcage", verbosity: int = 2) -> dict:
    des.validate()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    N, R, l = des.legs, des.radius, des.length
    wl, wr, g = des.leg_width, des.ring_width, des.gap
    zr = l / 2  # ring centrelines at +-zr
    has_leg_caps = des.coil_type in ("lowpass", "bandpass")
    has_ring_caps = des.coil_type in ("highpass", "bandpass")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("General.Verbosity", verbosity)
    gmsh.model.add(name)
    occ = gmsh.model.occ

    phi0 = des.phi0
    sheets = []  # (tag, kind, index)
    for n, ph in enumerate(leg_angles(N, phi0)):
        da = wl / (2 * R)
        z0, z1 = -zr + wr / 2, zr - wr / 2
        if has_leg_caps:
            sheets.append((_strip(occ, R, ph - da, ph + da, z0, -g / 2), "coil", n))
            sheets.append((_strip(occ, R, ph - da, ph + da, -g / 2, g / 2), "legcap", n))
            sheets.append((_strip(occ, R, ph - da, ph + da, g / 2, z1), "coil", n))
        else:
            sheets.append((_strip(occ, R, ph - da, ph + da, z0, z1), "coil", n))
    for side, zc in (("top", zr), ("bot", -zr)):
        if not has_ring_caps:
            sheets.append((_strip(occ, R, 0, 2 * math.pi, zc - wr / 2, zc + wr / 2), "coil", -1))
            continue
        dg = g / (2 * R)
        for n, ph in enumerate(leg_angles(N, phi0)):  # gap n sits between rung n and n+1
            mid = ph + math.pi / N
            nxt = ph + 2 * math.pi / N + math.pi / N
            # flat patch across the gap (its two z-parallel edges are exactly the segment ends), so the lumped
            # element is planar and its Palace direction is exact
            sheets.append((_flat_gap(occ, R, mid - dg, mid + dg, zc - wr / 2, zc + wr / 2), f"ringcap_{side}", n))
            sheets.append((_strip(occ, R, mid + dg, nxt - dg, zc - wr / 2, zc + wr / 2), "coil", n))

    if des.shield_diameter:
        Rd, Ld = des.shield_diameter / 2, des.shield_len
    else:
        Rd, Ld = max(2.0 * R, R + 0.15), l + 0.6
    air = occ.addCylinder(0, 0, -Ld / 2, 0, 0, Ld, Rd)
    ph_tag = None
    if des.phantom:
        p = des.phantom
        ph_tag = occ.addCylinder(0, 0, -p.length / 2, 0, 0, p.length, p.radius)
    planes = des.planes
    if planes:
        # Keep x >= 0 (half), x, y >= 0 (quarter), x, y, z >= 0 (eighth): intersect everything with a box.
        big = 2 * max(Rd, Ld)
        lo = [0.0 if a in planes else -big for a in "xyz"]
        def _box():
            return occ.addBox(lo[0], lo[1], lo[2], big - lo[0], big - lo[1], big - lo[2])
        kept = []
        for t, kind, n in sheets:
            out, _ = occ.intersect([(2, t)], [(3, _box())], removeObject=True, removeTool=True)
            kept += [(tt, kind, n) for d, tt in out if d == 2]
        sheets = kept
        out, _ = occ.intersect([(3, air)], [(3, _box())], removeObject=True, removeTool=True)
        air = [t for d, t in out if d == 3][0]
        if ph_tag is not None:
            out, _ = occ.intersect([(3, ph_tag)], [(3, _box())], removeObject=True, removeTool=True)
            ph_tag = [t for d, t in out if d == 3][0]
    tools = [(2, t) for t, _, _ in sheets]
    if ph_tag is not None:
        tools.append((3, ph_tag))
    _, omap = occ.fragment([(3, air)], tools)
    occ.synchronize()

    vol_air = sorted({t for d, t in omap[0] if d == 3})
    vol_ph = []
    if ph_tag is not None:
        vol_ph = sorted({t for d, t in omap[-1] if d == 3})
        vol_air = [v for v in vol_air if v not in vol_ph]
    coil, legcap, ring_top, ring_bot = [], {}, {}, {}
    for i, (_, kind, n) in enumerate(sheets):
        pieces = [t for d, t in omap[1 + i] if d == 2]
        if kind == "coil":
            coil += pieces
        elif kind == "legcap":
            legcap.setdefault(n, []).extend(pieces)
        elif kind == "ringcap_top":
            ring_top.setdefault(n, []).extend(pieces)
        else:
            ring_bot.setdefault(n, []).extend(pieces)
    inner = set(coil) | {s for v in (legcap, ring_top, ring_bot) for ss in v.values() for s in ss}
    outer = []
    for v in vol_air + vol_ph:
        for _, s in gmsh.model.getBoundary([(3, v)], oriented=False):
            s = abs(s)
            up, _ = gmsh.model.getAdjacencies(2, s)
            if len(up) == 1 and s not in inner:
                outer.append(s)
    outer = sorted(set(outer))
    sym_faces = {a: [] for a in planes}
    if planes:
        rest = []
        for s_ in outer:
            bb = gmsh.model.getBoundingBox(2, s_)
            for i, a in enumerate("xyz"):
                if a in planes and abs(bb[i]) < 1e-6 and abs(bb[3 + i]) < 1e-6:
                    sym_faces[a].append(s_)
                    break
            else:
                rest.append(s_)
        outer = rest

    gmsh.model.addPhysicalGroup(3, vol_air, ATTR_AIR, name="air")
    if vol_ph:
        gmsh.model.addPhysicalGroup(3, vol_ph, ATTR_PHANTOM, name="phantom")
    gmsh.model.addPhysicalGroup(2, sorted(set(coil)), ATTR_COIL, name="coil_pec")
    if des.shield_diameter:
        lateral, caps = [], []
        for s in outer:
            bb = gmsh.model.getBoundingBox(2, s)
            (caps if abs(bb[5] - bb[2]) < 1e-9 else lateral).append(s)
        if des.shield_ends == "pec":
            gmsh.model.addPhysicalGroup(2, outer, ATTR_SHIELD, name="shield_pec")
        else:
            gmsh.model.addPhysicalGroup(2, lateral, ATTR_SHIELD, name="shield_pec")
            gmsh.model.addPhysicalGroup(2, caps, ATTR_OUTER, name="ends_absorbing")
    else:
        gmsh.model.addPhysicalGroup(2, outer, ATTR_OUTER, name="outer_absorbing")
    for a, ff in sym_faces.items():
        gmsh.model.addPhysicalGroup(2, sorted(ff), ATTR_SYM[a], name=f"sym_{a}0")
    ports = []
    for n, ss in sorted(legcap.items()):
        gmsh.model.addPhysicalGroup(2, ss, ATTR_LEG_CAP0 + n, name=f"cap_leg{n}")
        # A capacitor patch cut by a symmetry plane: lengthwise (through-rung cut, half the width) each half carries
        # half the current at the same voltage -> C/2; across (z=0 plane, half the gap length) each half sees half the
        # voltage at the same current -> 2C. Both keep Palace's per-square value Cs = C * length / width unchanged.
        area = sum(gmsh.model.occ.getMass(2, t) for t in ss)
        bb = [gmsh.model.getBoundingBox(2, t) for t in ss]
        zlen = max(b[5] for b in bb) - min(b[2] for b in bb)
        w_frac = area / zlen / wl
        l_frac = zlen / g
        ports.append({"attr": ATTR_LEG_CAP0 + n, "kind": "leg", "n": n, "phi": leg_angles(N, phi0)[n],
                      "dir": [0.0, 0.0, 1.0], "w_frac": round(w_frac, 4), "l_frac": round(l_frac, 4),
                      "c_scale": round(w_frac, 3) / round(l_frac, 3)})
    for base, side, d in ((ATTR_RING_TOP0, "ring_top", ring_top), (ATTR_RING_BOT0, "ring_bot", ring_bot)):
        for n, ss in sorted(d.items()):
            mid = leg_angles(N, phi0)[n] + math.pi / N
            # +phi direction on the top ring, -phi on the bottom ring: the same sense as a mesh current.
            sgn = 1.0 if side == "ring_top" else -1.0
            gmsh.model.addPhysicalGroup(2, ss, base + n, name=f"cap_{side}{n}")
            ports.append({"attr": base + n, "kind": side, "n": n, "phi": mid,
                          "dir": [-sgn * math.sin(mid), sgn * math.cos(mid), 0.0]})

    # Mesh size: fine along the strip edges (current crowding) and finer at the capacitor gaps.
    fld = gmsh.model.mesh.field
    coil_curves = sorted({abs(c) for s in inner for _, c in gmsh.model.getBoundary([(2, s)], oriented=False)})
    cap_surfs = sorted({s for v in (legcap, ring_top, ring_bot) for ss in v.values() for s in ss})
    d1 = fld.add("Distance")
    fld.setNumbers(d1, "CurvesList", coil_curves)
    fld.setNumber(d1, "Sampling", 40)
    t1 = fld.add("Threshold")
    fld.setNumber(t1, "InField", d1)
    for k, v in (("SizeMin", des.h_coil), ("SizeMax", des.h_max), ("DistMin", des.h_coil), ("DistMax", 0.12)):
        fld.setNumber(t1, k, v)
    fields = [t1]
    if cap_surfs:
        d2 = fld.add("Distance")
        fld.setNumbers(d2, "SurfacesList", cap_surfs)
        fld.setNumber(d2, "Sampling", 10)
        t2 = fld.add("Threshold")
        fld.setNumber(t2, "InField", d2)
        for k, v in (("SizeMin", des.h_gap), ("SizeMax", des.h_max), ("DistMin", des.h_gap), ("DistMax", 0.05)):
            fld.setNumber(t2, k, v)
        fields.append(t2)
    if vol_ph:  # resolve the wavelength in the phantom (lambda/sqrt(eps_r) ~ 0.26 m at 128 MHz)
        bx = fld.add("Constant")
        fld.setNumbers(bx, "VolumesList", vol_ph)
        fld.setNumber(bx, "VIn", min(des.h_max, 0.02))
        fld.setNumber(bx, "VOut", des.h_max)
        fields.append(bx)
    mn = fld.add("Min")
    fld.setNumbers(mn, "FieldsList", fields)
    fld.setAsBackgroundMesh(mn)
    for k, v in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                 ("Mesh.MeshSizeFromCurvature", 24), ("Mesh.MeshSizeMax", des.h_max), ("Mesh.Algorithm", 6),
                 ("Mesh.Algorithm3D", 10), ("Mesh.MshFileVersion", 2.2), ("Mesh.SaveAll", 0)):
        gmsh.option.setNumber(k, v)
    gmsh.model.mesh.generate(3)
    n_tet = len(gmsh.model.mesh.getElementsByType(4)[0])
    msh = out_dir / f"{name}.msh"
    gmsh.write(str(msh))
    # Surface triangles of the coil, shield and caps for rendering (small, numpy friendly).
    surf = {}
    for label, tags in [("coil", sorted(set(coil))), ("caps", cap_surfs), ("outer", outer)] +             [(f"sym_{a}", ff) for a, ff in sym_faces.items()]:
        tris, nodes = [], {}
        for s in tags:
            et, _, en = gmsh.model.mesh.getElements(2, s)
            for t, conn in zip(et, en):
                if t == 2:
                    tris += list(map(int, conn))
        ids = sorted(set(tris))
        for i in ids:
            nodes[i] = [round(x, 6) for x in gmsh.model.mesh.getNode(i)[0].tolist()]
        idx = {i: j for j, i in enumerate(ids)}
        surf[label] = {"nodes": [nodes[i] for i in ids], "tris": [idx[i] for i in tris]}
    gmsh.finalize()

    info = {"name": name, "mesh": msh.name, "tets": int(n_tet), "ports": ports, "domain_radius": Rd, "domain_length": Ld,
            "attrs": {"air": ATTR_AIR, "phantom": ATTR_PHANTOM if vol_ph else None, "coil": ATTR_COIL,
                      "shield": ATTR_SHIELD if des.shield_diameter else None,
                      "absorbing": ATTR_OUTER if (not des.shield_diameter or des.shield_ends != "pec") else None,
                      "sym": {a: ATTR_SYM[a] for a in planes}},
            "symmetry": des.symmetry, "planes": planes, "phi0": phi0}
    (out_dir / "build.json").write_text(json.dumps(info, indent=1))
    (out_dir / "surfaces.json").write_text(json.dumps(surf))
    des.to_json(out_dir / "design.json")
    return info
