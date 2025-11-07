from .simulation import Simulation
from .triangle import FEATriangle
from .node import Node
from .vector3 import Vector3
import plotly.graph_objects as go
import numpy as np
import scipy




sim = Simulation(36, 120)




bad = 0
for tri in sim.featriangles:
    for (zeta,eta) in [(1/6,1/6),(2/3,1/6),(1/6,2/3)]:
        detJ = np.linalg.det(tri.get_jacobian(zeta,eta))
        if detJ <= 1e-8:
            print("Small detJ", tri, detJ)
            bad += 1
print("bad detJ count:", bad)

print()

d = sim.calculate_prestress_deflections()

jumps = []
for tri in sim.featriangles:
    nodes = tri.get_node_indexes()
    # extract rotation vector (global) at each node from solution d
    rvecs = []
    for nid in nodes:
        r_global = np.array([d[nid*6+3], d[nid*6+4], d[nid*6+5]])
        # transform into element-local for comparison (optional)
        r_local = np.transpose(tri.Re) @ r_global
        rvecs.append(r_local)
    # edges (0-1),(1-2),(2-0),(0-3),(1-4),(2-5) — compare adjacent nodes on element
    edge_pairs = [(0,1),(1,2),(2,0),(0,3),(1,4),(2,5)]
    for a,b in edge_pairs:
        jump = np.linalg.norm(rvecs[a]-rvecs[b])
        jumps.append(jump)
# Sort and print largest jumps
jumps_sorted = sorted(jumps, reverse=True)
print("top rotation jumps:", jumps_sorted[:20])

print()

bending_energies = []
for tri in sim.featriangles:
    e = 0.0
    for (zeta,eta,w) in [(1/6,1/6,1/6),(2/3,1/6,1/6),(1/6,2/3,1/6)]:
        B_b = tri.get_bending_B(zeta,eta)
        D_b = tri.get_bending_D(zeta,eta)
        # make local displacement vector for tri
        local_d = []
        for nid in tri.get_node_indexes():
            ug = np.array([d[nid*6+0], d[nid*6+1], d[nid*6+2]])
            rg = np.array([d[nid*6+3], d[nid*6+4], d[nid*6+5]])
            u_local = np.transpose(tri.Re) @ ug
            r_local = np.transpose(tri.Re) @ rg
            local_d.extend([u_local[0], u_local[1], u_local[2], r_local[0], r_local[1]])
        local_d = np.array(local_d)
        curv = B_b @ local_d   # curvature vector
        e += 0.5 * (curv.T @ D_b @ curv) * np.linalg.det(tri.get_jacobian(zeta,eta)) * w
    bending_energies.append(e)
# inspect distribution
print("bending energy min/median/max:", min(bending_energies), np.median(bending_energies), max(bending_energies))
