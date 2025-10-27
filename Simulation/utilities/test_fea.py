from .simulation import Simulation
from .triangle import FEATriangle
from .node import Node
from .vector3 import Vector3
import plotly.graph_objects as go




sim = Simulation(36, 120)

def replace_with_test_FEA(simulation):
    p1 = Node(Vector3(0,0,0))
    p2 = Node(Vector3(120,0,0))
    p3 = Node(Vector3(60,36,0))
    p4 = Node(Vector3(60,18,36))
    sim.triangles = [FEATriangle(p1, p2, p3), 
                     FEATriangle(p1, p2, p4), 
                     FEATriangle(p2, p3, p4), 
                     FEATriangle(p1, p3, p4)]
    sim.featriangles = sim.triangles_to_featriangles()
    sim.get_nodes()
    sim.nodes[0].fix_all()
    sim.nodes[2].fix_all()
    sim.nodes[4].fix_all()

replace_with_test_FEA(sim)

fig = go.Figure()
sim.draw(fig)
sim.draw_fixities(fig)
sim.draw_prestress_deflected(fig, scale=.05)
fig.update_layout(scene=dict(aspectmode='data'))
fig.show()

def test_simulation_creation():
    assert isinstance(sim, Simulation)

Ks = sim.get_structural_stiffness()
def test_structural_stiffness():
    assert Ks.shape == (60,60)

Ms = sim.get_mass_matrix()
def test_mass_matrix():
    assert Ms.shape == (60,60)


f = sim.get_f()
def test_get_f():
    assert len(f) == 60