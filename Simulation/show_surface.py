from utilities.simulation import Simulation
import matplotlib.pyplot as plt
import plotly.graph_objects as go



# ax = plt.figure().add_subplot(projection='3d')
# sim = Simulation(2, 3)
# sim.get_points()
# sim.draw(ax)
# print(sim.triangles)


# # ax = plt.figure().add_subplot(projection='3d')
# # ax.plot_trisurf([0,1,2],[1,2,0],[2,0,1])

# plt.show()


fig = go.Figure()
sim = Simulation(36, 120, saddles = 3)
sim.get_tris()
sim.draw(fig, .3)
sim.draw_prestress_deflected(fig, .7, .1)
sim.draw_fixities(fig)

fig.update_layout(scene=dict(aspectmode='data'))

fig.show()