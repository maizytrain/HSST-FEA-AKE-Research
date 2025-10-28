from .simulation import Simulation
from .triangle import FEATriangle
from .node import Node
from .vector3 import Vector3
import plotly.graph_objects as go
import numpy as np
import scipy




sim = Simulation(36, 120)

# Ks = sim.get_structural_stiffness()
# print("K conditional number:", np.linalg.cond(Ks))
# F = np.zeros(Ks.shape[0])
# F[20] = 100
# d = np.linalg.solve(Ks, F)
# print("Max Deflection:", np.max(np.abs(d)))



# Ke = sim.featriangles[0].get_Ke_global()
# print("Ke symmetry error:", np.max(np.abs(Ke - Ke.T)))
# print("Ke condition number:", np.linalg.cond(Ke))


# Ke_local = sim.featriangles[0].get_Ke()
# Ke_global = sim.featriangles[0].get_Ke_global()
# print("Local cond:", np.linalg.cond(Ke_local))
# print("Global cond:", np.linalg.cond(Ke_global))



# print("Re:\n", sim.featriangles[0].Re)
# print("Re * Re.T:\n", sim.featriangles[0].Re @ sim.featriangles[0].Re.T)


# print()
# print("Membrane D", sim.featriangles[0].get_membrane_D(0,0))
# print("Bending D", sim.featriangles[0].get_bending_D(0,0))
# print("Shear D", sim.featriangles[0].get_shear_D(0,0))

# locked_count = 0
# for i in range(len(sim.nodes)):
#     if (sim.nodes[i].fixed == [True, True, True, True, True, True]):
#         locked_count += 1
# print(locked_count)


# print()
# A = sim.featriangles[0].approximate_area()
# zetas = [1/6, 2/3, 1/6]
# etas = [1/6, 1/6, 2/3]
# sum_detJ = 0
# for i in range(3):
#     detJ = np.linalg.det(sim.featriangles[0].get_jacobian(zetas[i], etas[i]))
#     sum_detJ += detJ
# print("Area from detJ:", sum_detJ/2, "Analytical area:", A)



# Ks = sim.get_structural_stiffness()
# Ks_r = sim.remove_fixities(Ks)

# def power_iteration(A, num_iter=50):
#     b = np.random.rand(A.shape[0])
#     b /= np.linalg.norm(b)
#     for _ in range(num_iter):
#         b = A @ b
#         b /= np.linalg.norm(b)
#     return b

# def smallest_mode_estimate(A, num_iter=100):
#     # inverse iteration estimate for smallest mode using Rayleigh quotient
#     x = np.random.rand(A.shape[0])
#     for _ in range(num_iter):
#         Ax = A @ x
#         lam = x.dot(Ax) / x.dot(x)
#         x = Ax / np.linalg.norm(Ax)
#     return lam

# lam_est = smallest_mode_estimate(Ks_r)
# print("Approx smallest eigenvalue:", lam_est)


# bad = []
# for i, tri in enumerate(sim.featriangles):
#     dets = []
#     zetas = [1/6, 2/3, 1/6]
#     etas = [1/6, 1/6, 2/3]
#     for z, e in zip(zetas, etas):
#         detJ = np.linalg.det(tri.get_jacobian(z, e))
#         dets.append(detJ)
#     if any(abs(d) < 1e-6 for d in dets):
#         bad.append((i, dets, tri.p1, tri.p2, tri.p3))
# print("Num degenerate elements:", len(bad))
# if bad:
#     print("First few bad elements:", bad[:3])


# mid_nodes = {}
# for tri in sim.featriangles:
#     mids = [tri.n4, tri.n5, tri.n6]
#     for m in mids:
#         mid_nodes[m] = mid_nodes.get(m, 0) + 1

# unshared = [n for n, count in mid_nodes.items() if count == 1]
# print("Total midside nodes:", len(mid_nodes))
# print("Unshared midside nodes (should be 0 ideally):", len(unshared))
# if len(unshared) < 20:
#     print("Example unshared midside node indexes:", unshared[:10])

# magnitudes = [np.max(np.abs(tri.get_Ke_global())) for tri in sim.featriangles]
# print("Ke magnitude stats → min:", np.min(magnitudes),
#       "median:", np.median(magnitudes),
#       "max:", np.max(magnitudes))

# for i in range(3):
#     Te = sim.featriangles[i].get_Te()
#     print(f"Element {i}: Te^T*Te symmetry error =", np.max(np.abs(Te.T @ Te - np.eye(Te.shape[1]))))

# print()
# Re = sim.featriangles[0].Re
# print("Re @ Re.T:", Re @ Re.T)
# Te = sim.featriangles[0].get_Te()
# print("Te^T @ Te symmetry error:", np.max(np.abs(Te.T @ Te - np.eye(Te.shape[1]))))
# print()
# print(sim.featriangles[0].get_normal())
# # parr = Te.T @ Te
# # strs = []
# # for i in range(parr.shape[0]):
# #     strs.append("")
# #     for j in range(parr.shape[1]):
# #         strs[i] += str(parr[i,j]) + ", "
# #     print(strs[i])


# from numpy.linalg import norm
# Ks = sim.get_structural_stiffness()
# Ks_r = sim.remove_fixities(Ks)

# # power iteration to estimate largest eigenvalue
# def power_iter(A, steps=200):
#     b = np.random.rand(A.shape[0])
#     b /= norm(b)
#     for _ in range(steps):
#         b = A @ b
#         b /= norm(b)
#     Ab = A @ b
#     return b.dot(Ab)  # Rayleigh quotient ~ largest eigenvalue

# lam_max_est = power_iter(Ks_r, steps=200)
# print("Estimated largest eigenvalue:", lam_max_est)

# print()

# Ks_r = sim.remove_fixities(Ks)
# absK = np.abs(Ks_r)
# maxval = absK.max()
# print("Ks_reduced max absolute entry:", maxval)

# # list top 10 largest entries
# flat = absK.flatten()
# idxs = np.argsort(flat)[-10:][::-1]
# n = Ks_r.shape[0]
# for idx in idxs:
#     r = idx // n
#     c = idx % n
#     print(f"Entry ({r},{c}) = {Ks_r[r,c]}")


# print()

# def dof_to_node_dof(idx):
#     node = idx // 6
#     dof = idx % 6
#     return node, dof

# # Example: map the top offending entries from step 2:
# offending = [ 6062, 5588, 3481, 5048, 4592 ]  # copy the r values from step 2's printout
# for r in offending:
#     node, dof = dof_to_node_dof(r)
#     print("DOF", r, "=> node", node, "dof_index", dof, "position:", sim.nodes[node].position)


# print()


# suspect_row = r  # one of the row indices from step 2
# # assemble contributions per element to that global row (fast)
# contribs = []
# for ei, tri in enumerate(sim.featriangles):
#     Ke_g = tri.get_Ke_global()
#     nodes = tri.get_node_indexes()
#     # map local entries for each local DOF j to global index
#     for local_i in range(6):
#         for local_d in range(6):
#             global_idx = nodes[local_i] * 6 + local_d
#             val = Ke_g[local_i*6 + local_d, :].copy()  # whole row for this local DOF
#             # check the absolute contribution to suspect_row (if matches)
#             contrib = Ke_g[local_i*6 + local_d, :]  # full vector
#             # contribution to suspect_row entry:
#             contrib_val = Ke_g[local_i*6 + local_d, (suspect_row)]
#             if abs(contrib_val) > 1e-6 * maxval:  # sizable fraction of max
#                 contribs.append((ei, local_i, local_d, contrib_val))
# # print findings
# print("Element contributions to suspect row (>1e-6*max):")
# for c in contribs[:20]:
#     print(c)

# print()

# diag = np.diag(Ks_r)
# min_diag = diag.min()
# zero_diag_indices = np.where(diag < 1e-8)[0]
# print("min diag:", min_diag, "num tiny diag <1e-8:", len(zero_diag_indices))

# print()

# print("Symmetry error (reduced):", np.max(np.abs(Ks_r - Ks_r.T)))
# print("Any NaNs? ", np.isnan(Ks_r).any(), "Any Infs? ", np.isinf(Ks_r).any())


# Ks = sim.get_structural_stiffness()
# Ks_r = sim.remove_fixities(Ks)
# diag = np.diag(Ks_r)
# tol = 1e-8 * np.median(np.abs(diag))  # tiny threshold relative to typical scale
# tiny_idxs = np.where(np.abs(diag) < tol)[0]
# print("Total reduced DOFs:", Ks_r.shape[0])
# print("Number of tiny-diagonal DOFs:", len(tiny_idxs))
# # map a few to node,dof
# def map_idx(idx):
#     node = idx // 6
#     dof = idx % 6
#     return node, dof
# print("Examples (global_index -> node, dof):", [ (i, map_idx(i)) for i in tiny_idxs[:20] ])


# print()
# Ks = add_drilling_stiffness(sim.get_structural_stiffness(), sim, 1e-6)
# Ks_r = sim.remove_fixities(Ks)
# vals = np.linalg.eigvals(Ks_r)
# vals = np.real_if_close(vals)
# vals_sorted = np.sort(vals)
# print("Smallest 20 eigenvalues:", vals_sorted[:20])


# Ks = sim.get_structural_stiffness()
# Ks_r = sim.remove_fixities(Ks)
# absK = np.abs(Ks_r)
# maxval = absK.max()
# print("Global max abs entry:", maxval)

# # list the top global diag indices like you saw earlier
# diag = np.diag(Ks_r)
# # get top 20 diagonal indices by value
# top_diags_idx = np.argsort(diag)[-20:][::-1]
# print("Top diag indices and values (global_index: value):")
# for idx in top_diags_idx[:20]:
#     print(idx, diag[idx])

# # build per-element contributions to global diagonals
# # for speed, precompute each element's local->global mapping and its Ke_global
# element_maps = []
# for ei, tri in enumerate(sim.featriangles):
#     nodes = tri.get_node_indexes()  # length 6
#     # if node index -1 present, skip
#     if any(n < 0 for n in nodes):
#         element_maps.append(None)
#         continue
#     local_to_global = []
#     for local_node in nodes:
#         for d in range(6):
#             local_to_global.append(local_node * 6 + d)  # global index
#     # local_to_global is length 36
#     Ke_g = tri.get_Ke_global()  # should be 36x36
#     element_maps.append((local_to_global, Ke_g))

# # function to get element contributions to a particular global diagonal
# def element_diag_contribs(global_idx, threshold_fraction=1e-6):
#     contribs = []
#     for ei, entry in enumerate(element_maps):
#         if entry is None: continue
#         local_to_global, Ke_g = entry
#         # find local index (or indices) mapping to this global_idx
#         matches = [li for li, g in enumerate(local_to_global) if g == global_idx]
#         if not matches:
#             continue
#         # sum diagonal contributions from those local indices (should be one per match)
#         val = 0.0
#         for li in matches:
#             val += Ke_g[li, li]
#         if abs(val) > threshold_fraction * maxval:
#             contribs.append((ei, matches, val))
#     # sort by magnitude descending
#     contribs.sort(key=lambda x: abs(x[2]), reverse=True)
#     return contribs

# # inspect the top diagonal indices you already saw (or the top 10)
# interesting = top_diags_idx[:10]
# for gidx in interesting:
#     c = element_diag_contribs(gidx, threshold_fraction=1e-8)
#     print("\nGlobal DOF", gidx, "-> contributions (element_index, local_indices, value):")
#     for item in c[:20]:
#         print(item)


# print()

# coords = {}
# dups = []
# for n in sim.nodes:
#     key = (round(n.position.x,9), round(n.position.y,9), round(n.position.z,9))
#     if key in coords:
#         dups.append((coords[key], n.index, key))
#     else:
#         coords[key] = n.index
# print("Duplicate node coordinate count:", len(dups))
# if len(dups) > 0:
#     print("First duplicates:", dups[:10])

#     print()

# g = 4393
# n_global = Ks.shape[0]
# # compute row g by summing element contributions
# Ks_row_from_elems = np.zeros(n_global)
# for ei, tri in enumerate(sim.featriangles):
#     nodes = tri.get_node_indexes()
#     if any(n < 0 for n in nodes):
#         continue
#     # local->global mapping
#     local_to_global = []
#     for local_node in nodes:
#         for d in range(6):
#             local_to_global.append(local_node * 6 + d)
#     Ke_g = tri.get_Ke_global()  # should be 36x36
#     # add contributions
#     for li, g_i in enumerate(local_to_global):
#         for lj, g_j in enumerate(local_to_global):
#             Ks_row_from_elems[g_j] += Ke_g[li, lj] if g_i == g else 0.0

# print("Assembled Ks[g,g]:", Ks[g, g])
# print("Recomputed Ks[g,g]:", Ks_row_from_elems[g])
# # show a small mismatch metric
# print("Max abs diff between assembled Ks row and recomputed row:", np.max(np.abs(Ks[g, :] - Ks_row_from_elems)))

# print()

# dof_counts = {}
# for idx in top_diags_idx[:50]:
#     node = idx // 6
#     dof = idx % 6
#     dof_counts[dof] = dof_counts.get(dof, 0) + 1
# print("Top diag counts per DOF index:", dof_counts)

# print()

# from scipy.sparse.linalg import eigsh
# from numpy.linalg import norm

# # full Ks and then reduced (remove fixities)
# Ks = sim.get_structural_stiffness()
# Ks_reduced = sim.remove_fixities(Ks)   # size ~8856

# print("Ks shape:", Ks.shape, "Ks_reduced shape:", Ks_reduced.shape)

# # compute k smallest eigenpairs (k=6 or 8)
# k = 8
# # Use eigsh with sigma=0 (shift-invert) to get smallest magnitude eigenvalues quickly.
# vals, vecs = eigsh(Ks_reduced, k=k, sigma=0.0, which='LM', tol=1e-6, maxiter=5000)
# # eigsh may return eigenvalues out of order; sort them
# idx_sort = np.argsort(vals)
# vals = vals[idx_sort]
# vecs = vecs[:, idx_sort]

# print("Smallest eigenvalues (Ks_reduced):", vals)

# # For each eigenvector, find DOFs with largest magnitude
# def map_idx_to_node_dof(idx):
#     node = idx // 6
#     dof = idx % 6
#     return node, dof

# for ev_i in range(k):
#     v = vecs[:, ev_i]
#     # normalize so we can compare magnitudes
#     v = v / norm(v)
#     abs_v = np.abs(v)
#     top_idxs = np.argsort(abs_v)[-12:][::-1]  # top 12 contributing DOFs
#     print(f"\nMode {ev_i}, eigenvalue={vals[ev_i]:.3e}, top DOFs (global_index, node, dof, value):")
#     for gi in top_idxs:
#         # gi is index in reduced matrix; map back to full global index by using replace_fixities
#         # replace_fixities will place zeros at fixed DOFs — but we need mapping from reduced index -> full index
#         # Build a mapping once:
#         pass


# full_n = Ks.shape[0]
# # determine which indices were removed
# fixes = []
# for node in sim.nodes:
#     base = node.index * 6
#     for i in range(6):
#         if node.fixed[i]:
#             fixes.append(base + i)
# fixes = sorted(fixes)
# # Build list of kept indices
# all_indices = list(range(full_n))
# kept = [i for i in all_indices if i not in fixes]
# # kept[j] is the full global index that corresponds to reduced index j

# # Now finish printing top DOFs for each mode:
# for ev_i in range(k):
#     v = vecs[:, ev_i]
#     v = v / norm(v)
#     abs_v = np.abs(v)
#     top_idxs = np.argsort(abs_v)[-12:][::-1]  # top 12
#     print(f"\nMode {ev_i}, eigenvalue={vals[ev_i]:.3e}")
#     for ri in top_idxs:
#         full_idx = kept[ri]
#         node, dof = map_idx_to_node_dof(full_idx)
#         print(f" reduced_idx={ri} -> global={full_idx}, node={node}, dof={dof}, value={v[ri]:+.4f}")


# tri = sim.featriangles[0]
# ke = tri.get_Ke_global()
# print("Ke shape:", ke.shape)
# print("Top-left 6x6 block:\n", ke[:6, :6])

# nodeIndexes = tri.get_node_indexes()
# for i, n in enumerate(nodeIndexes):
#     for dof in range(6):
#         global_index = n*6 + dof
#         print("Node", n, "local DOF", dof, "-> global DOF", global_index)




positions = np.array([[n.position.x, n.position.y, n.position.z] for n in sim.nodes])
for n1 in sim.nodes:
    p1 = n1.position
    for n2 in sim.nodes:
        p2 = n2.position
        dist = (p1 - p2).magnitude()
        if dist < .001 and dist != 0:
            print("Close dist! Dist:", dist, "p1:", p1, "p2:", p2)
print("Total nodes:", len(sim.nodes))
print("Unique nodes:", len(np.unique(positions, axis=0)))