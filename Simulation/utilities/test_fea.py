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



Ks = sim.get_structural_stiffness()
Ks_r = sim.remove_fixities(Ks)

def power_iteration(A, num_iter=50):
    b = np.random.rand(A.shape[0])
    b /= np.linalg.norm(b)
    for _ in range(num_iter):
        b = A @ b
        b /= np.linalg.norm(b)
    return b

def smallest_mode_estimate(A, num_iter=100):
    # inverse iteration estimate for smallest mode using Rayleigh quotient
    x = np.random.rand(A.shape[0])
    for _ in range(num_iter):
        Ax = A @ x
        lam = x.dot(Ax) / x.dot(x)
        x = Ax / np.linalg.norm(Ax)
    return lam

lam_est = smallest_mode_estimate(Ks_r)
print("Approx smallest eigenvalue:", lam_est)


bad = []
for i, tri in enumerate(sim.featriangles):
    dets = []
    zetas = [1/6, 2/3, 1/6]
    etas = [1/6, 1/6, 2/3]
    for z, e in zip(zetas, etas):
        detJ = np.linalg.det(tri.get_jacobian(z, e))
        dets.append(detJ)
    if any(abs(d) < 1e-6 for d in dets):
        bad.append((i, dets, tri.p1, tri.p2, tri.p3))
print("Num degenerate elements:", len(bad))
if bad:
    print("First few bad elements:", bad[:3])


mid_nodes = {}
for tri in sim.featriangles:
    mids = [tri.n4, tri.n5, tri.n6]
    for m in mids:
        mid_nodes[m] = mid_nodes.get(m, 0) + 1

unshared = [n for n, count in mid_nodes.items() if count == 1]
print("Total midside nodes:", len(mid_nodes))
print("Unshared midside nodes (should be 0 ideally):", len(unshared))
if len(unshared) < 20:
    print("Example unshared midside node indexes:", unshared[:10])

magnitudes = [np.max(np.abs(tri.get_Ke_global())) for tri in sim.featriangles]
print("Ke magnitude stats → min:", np.min(magnitudes),
      "median:", np.median(magnitudes),
      "max:", np.max(magnitudes))

for i in range(3):
    Te = sim.featriangles[i].get_Te()
    print(f"Element {i}: Te^T*Te symmetry error =", np.max(np.abs(Te.T @ Te - np.eye(Te.shape[1]))))

print()
Re = sim.featriangles[0].Re
print("Re @ Re.T:", Re @ Re.T)
Te = sim.featriangles[0].get_Te()
print("Te^T @ Te symmetry error:", np.max(np.abs(Te.T @ Te - np.eye(Te.shape[1]))))