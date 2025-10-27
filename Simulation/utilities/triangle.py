from .vector3 import Vector3, create_Vector3_from_numpy_array
import numpy as np
from .node import Node



class Triangle:
    def __init__(self, p1, p2, p3):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    def __str__(self):
        return f"[{self.p1}], [{self.p2}], [{self.p3}]"
    def __repr__(self):
        return f"[{self.p1}], [{self.p2}], [{self.p3}]"


class FEATriangle:
    def __init__(self, node1:Node, node2:Node, node3:Node, E=29000000, v=.3, h=.25, rho=.28, shear_correction = 5):
        self.n1 = node1.index
        self.n2 = node2.index
        self.n3 = node3.index
        self.n4 = -1
        self.n5 = -1
        self.n6 = -1
        p1 = node1.position
        p2 = node2.position
        p3 = node3.position
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p1.average(p2)
        self.p5 = p2.average(p3)
        self.p6 = p3.average(p1)
        self.node4 = Node(self.p4)
        self.node5 = Node(self.p5)
        self.node6 = Node(self.p6)
        self.E = E
        self.v = v
        self.h = h
        self.rho = rho
        self.shear_correction = shear_correction
        self.set_basis_vectors()
        self.set_rotation_matrix()
        self.translate_points_to_local()

    def __str__(self):
        return f"[{self.p1}], [{self.p2}], [{self.p3}]"
    def __repr__(self):
        return f"[{self.p1}], [{self.p2}], [{self.p3}]"
    


    def set_node_indexes(self, n1, n2, n3, n4, n5, n6):
        self.n1 = n1
        self.n2 = n2
        self.n3 = n3
        self.n4 = n4
        self.n5 = n5
        self.n6 = n6
    

    def approximate_area(self):
        return abs(.5 * (self.p2 - self.p1).cross(self.p3 - self.p1).magnitude())


    def set_basis_vectors(self):
        v1 = self.p2 - self.p1
        v2 = self.p3 - self.p1
        n = v1.cross(v2)
        self.n = n.unit_vector()
        self.t1 = v1.unit_vector()
        self.t2 = n.cross(self.t1)

    def get_normal(self):
        v1 = self.p2 - self.p1
        v2 = self.p3 - self.p1
        n = v1.cross(v2)
        n = n.unit_vector()
        return n

    def set_rotation_matrix(self):
        self.Re = np.zeros((3,3))

        self.Re[0,0] = self.t1.x
        self.Re[1,0] = self.t1.y
        self.Re[2,0] = self.t1.z

        self.Re[0,1] = self.t2.x
        self.Re[1,1] = self.t2.y
        self.Re[2,1] = self.t2.z

        self.Re[0,2] = self.n.x
        self.Re[1,2] = self.n.y
        self.Re[2,2] = self.n.z

    def shape_functions(self, zeta, eta):
        shapes = []
        shapes.append((1 - zeta - eta) * (1 - 2 * zeta - 2 * eta))
        shapes.append(zeta * (2 * zeta - 1))
        shapes.append(eta * (2 * eta - 1))
        shapes.append(4 * zeta * (1 - zeta - eta))
        shapes.append(4 * zeta * eta)
        shapes.append(4 * eta * (1 - zeta - eta))
        return shapes

    def shape_functions_dzeta(self, zeta, eta):
        shapes = []
        shapes.append(4 * zeta + 4 * eta - 3)
        shapes.append(4 * zeta - 1)
        shapes.append(0)
        shapes.append(-8 * zeta - 4 * eta + 4)
        shapes.append(4 * eta)
        shapes.append(-4 * eta)
        return shapes

    def shape_functions_deta(self, zeta, eta):
        shapes = []
        shapes.append(4 * zeta + 4 * eta - 3)
        shapes.append(0)
        shapes.append(4 * eta - 1)
        shapes.append(-4 * zeta)
        shapes.append(4 * zeta)
        shapes.append(-4 * zeta - 8 * eta + 4)
        return shapes
    
    def shape_functions_dx_dy(self, zeta, eta):
        shapes_dzeta = self.shape_functions_dzeta(zeta, eta)
        shapes_deta = self.shape_functions_deta(zeta, eta)
        shapes_dx = []
        shapes_dy = []
        J = self.get_jacobian(zeta, eta)
        for i in range(len(shapes_dzeta)):
            temp = np.linalg.solve(J, np.array([shapes_dzeta[i], shapes_deta[i]]))
            shapes_dx.append(temp[0])
            shapes_dy.append(temp[1])
        return shapes_dx, shapes_dy

    def get_points_as_list(self):
        return [self.p1, self.p2, self.p3, self.p4, self.p5, self.p6]

    def translate_points_to_local(self):
        self.locals = []
        points = self.get_points_as_list()
        for i in range(6):
            vec = (points[i] - self.p1).give_np_array()
            local_coords_np = self.Re.T @ vec
            self.locals.append(create_Vector3_from_numpy_array(local_coords_np))

    def get_jacobian(self, zeta, eta):
        J = np.zeros((2,2))

        shapes_dzeta = self.shape_functions_dzeta(zeta, eta)
        shapes_deta = self.shape_functions_deta(zeta, eta)

        for i in range(6):
            J[0,0] += shapes_dzeta[i] * self.locals[i].x
            J[0,1] += shapes_dzeta[i] * self.locals[i].y
            J[1,0] += shapes_deta[i] * self.locals[i].x
            J[1,1] += shapes_deta[i] * self.locals[i].y

        return J
    
    def get_membrane_B(self, zeta, eta):
        dxs, dys = self.shape_functions_dx_dy(zeta, eta)
        B = np.zeros((3, 30))
        for i in range(6):
            B[0, 5 * i] = dxs[i]
            B[1, 5 * i + 1] = dys[i]
            B[2, 5 * i] = dys[i]
            B[2, 5 * i + 1] = dxs[i]
        return B
    
    def get_bending_B(self, zeta, eta):
        dxs, dys = self.shape_functions_dx_dy(zeta, eta)
        B = np.zeros((3, 30))
        for i in range(6):
            B[0, 5 * i + 3] = -dxs[i]
            B[1, 5 * i + 4] = -dys[i]
            B[2, 5 * i + 3] = -dys[i]
            B[2, 5 * i + 4] = -dxs[i]
        return B

    def get_shear_B(self, zeta, eta):
        Ns = self.shape_functions(zeta, eta)
        dxs, dys = self.shape_functions_dx_dy(zeta, eta)
        B = np.zeros((2, 30))
        for i in range(6):
            B[0, 5 * i + 2] = dxs[i]
            B[1, 5 * i + 2] = dys[i]
            B[0, 5 * i + 3] = Ns[i]
            B[1, 5 * i + 4] = Ns[i]
        return B
    
    def get_membrane_D(self, zeta, eta):
        D = np.zeros((3,3))
        D[0,0] = 1
        D[0,1] = self.v
        D[1,0] = self.v
        D[1,1] = 1
        D[2,2] = (1 - self.v) * .5

        D = D * (self.E * self.h / (1 - self.v**2))
        return D
    
    def get_bending_D(self, zeta, eta):
        D = np.zeros((3,3))
        D[0,0] = 1
        D[0,1] = self.v
        D[1,0] = self.v
        D[1,1] = 1
        D[2,2] = (1 - self.v) * .5

        D = D * (self.E * self.h**3 / (12 * (1 - self.v**2)))
        return D
    
    def get_shear_D(self, zeta, eta):
        D = np.zeros((2,2))
        D[0,0] = 1
        D[1,1] = 1

        D = D * (self.E * self.h / (2 * (1 - self.v)) * self.shear_correction)
        return D
    

    def get_Ke(self):
        zetas = [1/6, 2/3, 1/6]
        etas = [1/6, 1/6, 2/3]
        Ke = np.zeros((30,30))
        for i in range(len(zetas)):
            B_m = self.get_membrane_B(zetas[i], etas[i])
            D_m = self.get_membrane_D(zetas[i], etas[i])
            B_b = self.get_bending_B(zetas[i], etas[i])
            D_b = self.get_bending_D(zetas[i], etas[i])
            B_s = self.get_shear_B(zetas[i], etas[i])
            D_s = self.get_shear_D(zetas[i], etas[i])
            detJ = np.linalg.det(self.get_jacobian(zetas[i], etas[i]))
            if detJ <= 0:
                raise Exception("Dejenerate Plate" + '\n' + "DetJ:" + str(detJ) + '\n' +
                                                            "POINTS: " + str(self.p1) + ", " + str(self.p2) + ", " + 
                                                                        str(self.p3) + ", " + str(self.p4) + ", " + 
                                                                        str(self.p5) + ", " + str(self.p6) + '\n' + 
                                                                        "Normal: " + str(self.get_normal()))
            Km = np.transpose(B_m) @ D_m @ B_m
            Kb = np.transpose(B_b) @ D_b @ B_b
            Ks = np.transpose(B_s) @ D_s @ B_s
            Ke += (Km + Kb + Ks) * detJ / 6
        return Ke

    def get_Te(self):
        Re = self.Re
        Te = np.zeros((30,36))
        for i in range(6):
            Ti = np.zeros((5,6))
            for j in range(3):
                for k in range(3):
                    Ti[j,k] = Re[k,j] #Something is wrong here I think? It may be time to implement more rigorous tests.
                    if not j==2:
                        Ti[j+3,k+3] = Re[k,j]
            for j in range(5):
                for k in range(6):
                    Te[j + 5 * i, k + 6 * i] = Ti[j, k]
        return Te

    def get_Ke_global(self):
        Te = self.get_Te()
        Ke = self.get_Ke()
        Keglobal = np.transpose(Te) @ Ke @ Te
        return Keglobal

    def get_mass_matrix(self):
        M = np.zeros((30,30))
        zetas = [1/6, 2/3, 1/6]
        etas = [1/6, 1/6, 2/3]
        for i in range(len(zetas)):
            detJ = np.linalg.det(self.get_jacobian(zetas[i], etas[i]))
            shapes = self.shape_functions(zetas[i], etas[i])
            Nt = np.zeros((3,18))
            Nr = np.zeros((2,12))
            for i in range(6):
                Nt[0, 3 * i] = shapes[i]
                Nt[1, 3 * i + 1] = shapes[i]
                Nt[2, 3 * i + 2] = shapes[i]
                Nr[0, 2 * i] = shapes[i]
                Nr[1, 2 * i + 1] = shapes[i]
            Mtt = self.rho * self.h * np.transpose(Nt) @ Nt
            Mrr = self.rho * self.h**3 / 12 * np.transpose(Nr) @ Nr
            for i in range(Mtt.shape[0]):
                for j in range(Mtt.shape[1]):
                    M[i,j] += detJ * Mtt[i,j]
            for i in range(Mrr.shape[0]):
                for j in range(Mrr.shape[1]):
                    M[i + Mtt.shape[0], j + Mtt.shape[1]] += detJ * Mrr[i,j]

        M *= 1/6
        return M

    def get_Me_global(self):
        Te = self.get_Te()
        Me = self.get_mass_matrix()
        Meglobal = np.transpose(Te) @ Me @ Te
        return Meglobal
    
    def get_node_indexes(self):
        return [self.n1, self.n2, self.n3, self.n4, self.n5, self.n6]
    
    def get_normal_towards_center(self, center:Vector3) -> Vector3:
        n = self.get_normal()
        c = (self.p1 + self.p2 + self.p3) * (1/3)
        d = center - c
        val = n.dot(d)
        if val < 0:
            n = n * -1
        return n
    
    def get_f_helpers(self, center:Vector3):
        A = self.approximate_area()
        a = A / 6
        n = self.get_normal_towards_center(center)
        return a, n
    

    def get_pressure_force(self, pressure, center):
        zetas = [1/6, 2/3, 1/6]
        etas = [1/6, 1/6, 2/3]
        n = self.get_normal_towards_center(center)
        F = [Vector3(0,0,0), Vector3(0,0,0), Vector3(0,0,0), Vector3(0,0,0), Vector3(0,0,0), Vector3(0,0,0)]
        for i in range(len(zetas)):
            detJ = np.linalg.det(self.get_jacobian(zetas[i], etas[i]))
            shapes = self.shape_functions(zetas[i], etas[i])
            for j in range(len(shapes)):
                F[j] += n * (shapes[j] * pressure * detJ / 6)
        return F
    

    def get_membrane_strain(self, zeta, eta, d, B=None):
        if B == None:
            B = self.get_membrane_B(zeta, eta)
        return B @ d
    
    def get_geometric_S(self, zeta, eta, strains, D=None):
        if D == None:
            D = self.get_membrane_D(zeta, eta)
        N = D @ strains
        S = np.array([[N[0], N[2], 0], [N[2], N[1], 0], [0, 0, 0]])
        return S
    

    def get_geometric_stiffness(self, d):
        myd = []
        nodes = self.get_node_indexes()
        for n in nodes:
            for i in range(6):
                myd.append(d[n * 6 + i])

        zetas = [1/6, 2/3, 1/6]
        etas = [1/6, 1/6, 2/3]
        Ksigma = np.zeros((30,30))
        for i in range(len(zetas)):
            Bm = self.get_membrane_B(zetas[i], etas[i])
            detJ = np.linalg.det(self.get_jacobian(zetas[i], etas[i]))
            strains = self.get_membrane_strain(zetas[i], etas[i], myd, Bm)
            S = self.get_geometric_S(zetas[i], etas[i], strains)
            ks = np.transpose(Bm) @ S @ Bm * detJ / 6
            Ksigma += ks
        return Ksigma
    
    def get_Ks_global(self, d):
        Te = self.get_Te()
        Ksigma = self.get_geometric_stiffness(d)
        Ksigmaglobal = np.transpose(Te) @ Ksigma @ Te
        return Ksigmaglobal

            