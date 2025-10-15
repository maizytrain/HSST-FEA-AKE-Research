import matplotlib.pyplot as plt
import numpy as np
from .vector3 import Vector3
from .triangle import Triangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import plotly.graph_objects as go








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
    def __init__(self, cylinderRadius, cylinderLength, resolution=4):
        self.cylinderDiameter = cylinderRadius
        self.cylinderLength = cylinderLength
        self.resolution = resolution
        self.nodes = []
        self.triangles = []
        self.get_points()

    def get_points(self):
        points = []
        for z in [0, self.cylinderDiameter]:
            for y in [0, self.cylinderDiameter]:
                for x in [0, self.cylinderLength]:
                    points.append(Vector3(x, y, z))
        # squares = [Square(points[0], points[2], points[3], points[1]), 
        #            Square(points[0], points[4], points[6], points[2]),
        #            Square(points[2], points[6], points[7], points[3]),
        #            Square(points[4], points[5], points[7], points[6]),
        #            Square(points[0], points[1], points[5], points[4]),
        #            Square(points[1], points[3], points[7], points[5])]

        # # squares = [Square(points[0], points[2], points[3], points[1])]

        # for time in range(self.resolution):
        #     new_squares = []
        #     for i in range(len(squares)):
        #         new = squares[i].subdivide()
        #         for j in range(len(new)):
        #             new_squares.append(new[j])
        #     squares = new_squares

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
                vs = np.sqrt(self.cylinderDiameter**2 * .25 - (ps[j] - p_cent).magnitude()**2)
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
                vs = np.sqrt(self.cylinderDiameter**2 * .25 - (ps[j] - p_cent).magnitude()**2)
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

    def draw(self, fig):
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
            flatshading=False
        ))