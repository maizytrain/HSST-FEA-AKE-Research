# import matplotlib.pyplot as plt
import numpy as np
from .vector3 import Vector3
from .triangle import Triangle
from .triangle import FEATriangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go
from .node import Node








class Square:
    def __init__(self, p1:Vector3, p2:Vector3, p3:Vector3, p4:Vector3):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4

    def subdivide(self):
        middle = (self.p1 + self.p2 + self.p3 + self.p4) * .25
        p12 = (self.p1 + self.p2) * .5 #+ Vector3(.01, .01, .01)
        p23 = (self.p2 + self.p3) * .5 #+ Vector3(.01, .01, .01)
        p34 = (self.p3 + self.p4) * .5
        p41 = (self.p4 + self.p1) * .5

        return [Square(self.p1, p12, middle, p41), Square(self.p2, p23, middle, p12), Square(self.p3, p34, middle, p23), Square(self.p4, p41, middle, p34)]
    

    def give_triangles(self):
        return [Triangle(self.p1,self.p2,self.p3), Triangle(self.p1,self.p3,self.p4)]





class Simulation:
    def __init__(self, cylinderDiameter, cylinderLength, initialPressure = 100, resolution=3, saddles=2, saddle_inset=0, saddle_up=None):
        self.cylinderDiameter = cylinderDiameter
        self.cylinderLength = cylinderLength
        self.resolution = resolution
        self.saddles = saddles
        self.saddle_inset = saddle_inset
        if saddle_up == None:
            self.saddle_up = cylinderDiameter * .5
        else:
            self.saddle_up = saddle_up
        self.nodes = []
        self.triangles = []
        self.get_tris()
        self.featriangles = self.triangles_to_featriangles()
        self.get_nodes()
        self.define_restraints()
        self.original_volume = 4/3 * np.pi * (cylinderDiameter * .5)**3 + np.pi * (cylinderDiameter * .5)**2 * cylinderLength
        self.now_volume = self.original_volume
        self.original_pressure = initialPressure
        self.now_pressure = self.original_pressure
        self.debug_lines = []


    def get_tris(self):
        points = []
        for z in [0, self.cylinderDiameter]:
            for y in [0, self.cylinderDiameter]:
                for x in [0, self.cylinderLength]:
                    points.append(Vector3(x, y, z))

        front_squares = [Square(points[0], points[4], points[6], points[2])]
        length_squares = [Square(points[0], points[2], points[3], points[1]), 
                          Square(points[2], points[6], points[7], points[3]), 
                          Square(points[4], points[5], points[7], points[6]), 
                          Square(points[0], points[1], points[5], points[4])]
        back_squares = [Square(points[1], points[3], points[7], points[5])]


        for time in range(self.resolution):
            new_squares = []
            for i in range(len(front_squares)):
                new = front_squares[i].subdivide()
                for j in range(len(new)):
                    new_squares.append(new[j])
            front_squares = new_squares

        for time in range(self.resolution):
            new_squares = []
            for i in range(len(length_squares)):
                new = length_squares[i].subdivide()
                for j in range(len(new)):
                    new_squares.append(new[j])
            length_squares = new_squares

        for time in range(self.resolution):
            new_squares = []
            for i in range(len(back_squares)):
                new = back_squares[i].subdivide()
                for j in range(len(new)):
                    new_squares.append(new[j])
            back_squares = new_squares


        
        for i in range(len(length_squares)):
            ps = [length_squares[i].p1, length_squares[i].p2, length_squares[i].p3, length_squares[i].p4]
            for j in range(len(ps)):
                p_cent = Vector3(ps[j].x, self.cylinderDiameter*.5, self.cylinderDiameter*.5)
                psminp = ps[j] - p_cent
                ang = np.pi * .5
                if not (psminp.y == 0):
                    ang = abs(np.arctan(psminp.z / psminp.y))
                ls = self.cylinderDiameter * .5
                if (ang == 0 or ang == np.pi/2):
                    ls *= 1
                if (ang == np.pi/4):
                    ls *= np.sqrt(2)
                elif (ang < np.pi/4):
                    ls *= 1 / np.cos(ang)
                elif (ang > np.pi/4):
                    ls *= 1 / np.sin(ang)
                ps[j] = psminp * (self.cylinderDiameter * .5 / ls) + p_cent
            length_squares[i].p1 = ps[0]
            length_squares[i].p2 = ps[1]
            length_squares[i].p3 = ps[2]
            length_squares[i].p4 = ps[3]

        for i in range(len(front_squares)):
            ps = [front_squares[i].p1, front_squares[i].p2, front_squares[i].p3, front_squares[i].p4]
            for j in range(len(ps)):
                p_cent = Vector3(ps[j].x, self.cylinderDiameter*.5, self.cylinderDiameter*.5)
                psminp = ps[j] - p_cent
                ang = np.pi * .5
                if not (psminp.y == 0):
                    ang = abs(np.arctan(psminp.z / psminp.y))
                ls = self.cylinderDiameter * .5
                if (ang == 0 or ang == np.pi/2):
                    ls *= 1
                if (ang == np.pi/4):
                    ls *= np.sqrt(2)
                elif (ang < np.pi/4):
                    ls *= 1 / np.cos(ang)
                elif (ang > np.pi/4):
                    ls *= 1 / np.sin(ang)
                ps[j] = psminp * (self.cylinderDiameter * .5 / ls) + p_cent
                vs = np.sqrt(max(self.cylinderDiameter**2 * .25 - (ps[j] - p_cent).magnitude()**2, 0))
                ps[j].x -= vs
            front_squares[i].p1 = ps[0]
            front_squares[i].p2 = ps[1]
            front_squares[i].p3 = ps[2]
            front_squares[i].p4 = ps[3]

        for i in range(len(back_squares)):
            ps = [back_squares[i].p1, back_squares[i].p2, back_squares[i].p3, back_squares[i].p4]
            for j in range(len(ps)):
                p_cent = Vector3(ps[j].x, self.cylinderDiameter*.5, self.cylinderDiameter*.5)
                psminp = ps[j] - p_cent
                ang = np.pi * .5
                if not (psminp.y == 0):
                    ang = abs(np.arctan(psminp.z / psminp.y))
                ls = self.cylinderDiameter * .5
                if (ang == 0 or ang == np.pi/2):
                    ls *= 1
                if (ang == np.pi/4):
                    ls *= np.sqrt(2)
                elif (ang < np.pi/4):
                    ls *= 1 / np.cos(ang)
                elif (ang > np.pi/4):
                    ls *= 1 / np.sin(ang)
                ps[j] = psminp * (self.cylinderDiameter * .5 / ls) + p_cent
                vs = np.sqrt(max(self.cylinderDiameter**2 * .25 - (ps[j] - p_cent).magnitude()**2,0))
                ps[j].x += vs
            back_squares[i].p1 = ps[0]
            back_squares[i].p2 = ps[1]
            back_squares[i].p3 = ps[2]
            back_squares[i].p4 = ps[3]




        squares = []
        for sq in front_squares:
            squares.append(sq)
        for sq in length_squares:
            squares.append(sq)
        for sq in back_squares:
            squares.append(sq)
        
        tris = []
        for i in range(len(squares)):
            temp = squares[i].give_triangles()
            tris.append(temp[0])
            tris.append(temp[1])

        self.triangles = tris


    def define_restraints(self, tolerance = .001):
        xs = np.linspace(0 + self.saddle_inset, self.cylinderLength - self.saddle_inset, self.saddles)
        for x in xs:
            closest_x = -100000.0
            closest_dist = 100000
            for node in self.nodes:
                dist = abs(node.position.x - x)
                if dist < closest_dist:
                    closest_dist = dist
                    closest_x = node.position.x
            for node in self.nodes:
                if abs(node.position.x - closest_x) < tolerance:
                    if node.position.z < self.saddle_up:
                        node.fix_all()


    def triangles_to_featriangles(self, E=29000, v=.3, h=.25, rho=.28, shear_correction = 5):
        tris = []
        for tri in self.triangles:
            tris.append(FEATriangle(Node(tri.p1), Node(tri.p2), Node(tri.p3), E, v, h, rho, shear_correction))
        return tris
    

    def get_nodes(self):
        nodes = []
        for tri in self.featriangles:
            pos = [tri.p1, tri.p2, tri.p3, tri.p4, tri.p5, tri.p6]
            idxs = []
            for p in pos:
                found = False
                for n in nodes:
                    if p == n.position:
                        found = True
                        idxs.append(n.index)
                        break
                if not found:
                    nodes.append(Node(p, len(nodes)))
                    idxs.append(len(nodes) - 1)
            tri.set_node_indexes(idxs[0], idxs[1], idxs[2], idxs[3], idxs[4], idxs[5])
        self.nodes = nodes
                    



    def get_modal_frequencies(self):
        Ks = self.get_structural_stiffness()



    def get_structural_stiffness(self):
        Ks = np.zeros((6 * len(self.nodes), 6 * len(self.nodes)))
        for tri in self.featriangles:
            ks = tri.get_Ke_global()
            nodeIndexes = tri.get_node_indexes()
            for i in range(len(nodeIndexes)):
                for j in range(len(nodeIndexes)):
                    for k in range(6):
                        for l in range(6):
                            Ks[nodeIndexes[i] * 6 + k, nodeIndexes[j] * 6 + l] += ks[i * 6 + k, j * 6 + l]
        return Ks


    def get_mass_matrix(self):
        Ms = np.zeros((6 * len(self.nodes), 6 * len(self.nodes)))
        for tri in self.featriangles:
            m = tri.get_Me_global()
            nodeIndexes = tri.get_node_indexes()
            for i in range(len(nodeIndexes)):
                for j in range(len(nodeIndexes)):
                    for k in range(6):
                        for l in range(6):
                            Ms[nodeIndexes[i] * 6 + k, nodeIndexes[j] * 6 + l] += m[i * 6 + k, j * 6 + l]
        return Ms
    
    def get_center(self):
        return Vector3((self.cylinderLength + self.cylinderDiameter) * .5, self.cylinderDiameter * .5, self.cylinderDiameter * .5)
    

    def get_f(self):
        f = np.zeros(len(self.nodes) * 6)
        cent = self.get_center()
        for tri in self.featriangles:
            a, n = tri.get_f_helpers(cent)
            nodes = tri.get_node_indexes()
            for node in nodes:
                f[node * 6 + 0] += a * n.x
                f[node * 6 + 1] += a * n.y
                f[node * 6 + 2] += a * n.z
        return f

    def get_fluid_stiffness(self):
        kf = self.original_pressure / self.original_volume
        f = self.get_f()
        Kf = kf * np.outer(f,f)
        return Kf
    

    def calculate_prestress_forces(self):
        F = np.zeros(len(self.nodes) * 6)
        cent = self.get_center()
        for tri in self.featriangles:
            fs = tri.get_pressure_force(self.original_pressure, cent)
            nodes = tri.get_node_indexes()
            for i in range(len(nodes)):
                F[nodes[i] * 6 + 0] += fs[i].x
                F[nodes[i] * 6 + 1] += fs[i].y
                F[nodes[i] * 6 + 2] += fs[i].z
        return F
    
    def remove_fixities(self, mat):
        fixes = []
        new_mat = mat
        for node in self.nodes:
            baseval = node.index * 6
            for i in range(6):
                if node.fixed[i]:
                    fixes.append(baseval + i)
        fixes.sort(reverse=True)
        # print("Fixities Removed:", len(fixes))
        # if len(fixes) != len(np.unique(fixes)):
            # raise Exception("Multiple fixes of same type")
        new_mat = np.delete(new_mat, fixes, axis=0)
        if new_mat.ndim != 1:
            new_mat = np.delete(new_mat, fixes, axis=1)
        return new_mat
    
    def replace_fixities(self, rmat):
        fixes = []
        new_mat = rmat
        for node in self.nodes:
            baseval = node.index * 6
            for i in range(6):
                if node.fixed[i]:
                    fixes.append(baseval + i)
        fixes.sort(reverse=True)
        # print("Fixities Replaced:", len(fixes))
        zeros = np.zeros(new_mat.shape[0])
        if new_mat.ndim == 1:
            zeros = 0
        for i in fixes:
            new_mat = np.insert(new_mat, i, zeros, axis=0)
        if new_mat.ndim == 1:
            return new_mat
        zeros = np.zeros(new_mat.shape[1])
        for i in fixes:
            new_mat = np.insert(new_mat, i, zeros, axis=1)
        return new_mat



    
    def calculate_prestress_deflections(self):
        F = self.calculate_prestress_forces()
        F = self.remove_fixities(F)
        Ks = self.get_structural_stiffness()
        Ks = self.remove_fixities(Ks)
        d = np.linalg.solve(Ks, np.transpose(F))
        print("Before:", len(d))
        d = self.replace_fixities(d)
        print("After:", len(d))
        return d
    


    def get_prestress_stiffness(self):
        d = self.calculate_prestress_deflections()
        Ksigma = np.zeros((len(self.nodes) * 6, len(self.nodes) * 6))
        for tri in self.featriangles:
            ksigma = tri.get_Ks_global(d)
            nodes = tri.get_node_indexes()
            for n1 in nodes:
                for n2 in nodes:
                    for i in range(6):
                        for j in range(6):
                            Ksigma[n1 * 6 + i, n2 * 6 + j] += ksigma[i, j]
        return Ksigma
    
    




    def draw(self, fig, opacity=1.0, color="blue"):
        x_coords = []
        y_coords = []
        z_coords = []
        ivals = []
        jvals = []
        kvals = []

        for i in range(len(self.triangles)):
            tri = self.triangles[i]
            x_coords.append(tri.p1.x)
            x_coords.append(tri.p2.x)
            x_coords.append(tri.p3.x)
            y_coords.append(tri.p1.y)
            y_coords.append(tri.p2.y)
            y_coords.append(tri.p3.y)
            z_coords.append(tri.p1.z)
            z_coords.append(tri.p2.z)
            z_coords.append(tri.p3.z)

            ivals.append(3*i)
            jvals.append(3*i+1)
            kvals.append(3*i+2)

        fig.add_trace(go.Mesh3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            i=ivals,
            j=jvals,
            k=kvals,
            opacity=opacity,
            color=color,
            flatshading=False
        ))


    def draw_prestress_deflected(self, fig, opacity=1.0, scale=10.0, color="green"):
        x_coords = []
        y_coords = []
        z_coords = []
        ivals = []
        jvals = []
        kvals = []

        d = self.calculate_prestress_deflections()
        # print(d)
        print("Max:", max(d), "Min:", min(d))

        # F = self.calculate_prestress_forces()
        # F = self.remove_fixities(F)

        # # print(F)

        for i in range(len(self.triangles)):
            tri = self.featriangles[i]
            x_coords.append(tri.p1.x + d[tri.n1 * 6] * scale)
            x_coords.append(tri.p2.x + d[tri.n2 * 6] * scale)
            x_coords.append(tri.p3.x + d[tri.n3 * 6] * scale)
            y_coords.append(tri.p1.y + d[tri.n1 * 6 + 1] * scale)
            y_coords.append(tri.p2.y + d[tri.n2 * 6 + 1] * scale)
            y_coords.append(tri.p3.y + d[tri.n3 * 6 + 1] * scale)
            z_coords.append(tri.p1.z + d[tri.n1 * 6 + 2] * scale)
            z_coords.append(tri.p2.z + d[tri.n2 * 6 + 2] * scale)
            z_coords.append(tri.p3.z + d[tri.n3 * 6 + 2] * scale)

            # if i % 30 == 0 or i == len(self.triangles) - 1:
            #     print("n1:", tri.n1, "n2:", tri.n2, "n3:", tri.n3)
            #     print("n1+0:", d[tri.n1 * 6], "n2+0:", d[tri.n2 * 6], "n3+0:", d[tri.n3 * 6])
            #     print("n1+1:", d[tri.n1 * 6], "n2+1:", d[tri.n2 * 6], "n3+1:", d[tri.n3 * 6])
            #     print("n1+2:", d[tri.n1 * 6], "n2+2:", d[tri.n2 * 6], "n3+2:", d[tri.n3 * 6])

            ivals.append(3*i)
            jvals.append(3*i+1)
            kvals.append(3*i+2)

        # poi = 1476
        # print("1476:", d[poi * 6], d[poi * 6 + 1], d[poi * 6 + 2], d[poi * 6 + 3], d[poi * 6 + 4], d[poi * 6 + 5])

        # print(d[-6])
        # print(d[-5])
        # print(d[-4])
        # print(d[-3])
        # print(d[-2])
        # print(d[-1])

        fig.add_trace(go.Mesh3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            i=ivals,
            j=jvals,
            k=kvals,
            opacity=opacity,
            color=color,
            flatshading=False
        ))


    def draw_fixities(self, fig):
        scatter_x = []
        scatter_y = []
        scatter_z = []

        for node in self.nodes:
            if node.fixed[0]:
                pos = node.position
                scatter_x.append(pos.x)
                scatter_y.append(pos.y)
                scatter_z.append(pos.z)

        fig.add_trace(go.Scatter3d(
            x=scatter_x,
            y=scatter_y,
            z=scatter_z,
            mode='markers',
            name='Fixities',
            marker=dict(
                size=5,
                color='red',
                opacity=0.8
            )
        ))


    def debug_draw_line(self, start:Vector3, end:Vector3):
        self.debug_lines.append([start, end])

    def draw_debug_lines(self, fig, color="black"):
        xs = []
        ys = []
        zs = []
        for i in range(len(self.debug_lines)):
            xs.extend([self.debug_lines[i][0].x, self.debug_lines[i][1].x, None])
            ys.extend([self.debug_lines[i][0].y, self.debug_lines[i][1].y, None])
            zs.extend([self.debug_lines[i][0].z, self.debug_lines[i][1].z, None])

        fig.add_trace(go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode='lines',
            line=dict(color=color)
        ))