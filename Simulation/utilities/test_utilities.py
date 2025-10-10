from .vector3 import Vector3
from .triangle import Triangle






def test_vector3_creation():
    p1 = Vector3(1,2,3)

    assert isinstance(p1, Vector3)
    assert p1.x == 1
    assert p1.y == 2
    assert p1.z == 3

def test_vector3_add():
    p1 = Vector3(1,2,3)
    p2 = Vector3(4,-1,3)

    calc = p1 + p2
    ans = Vector3(5, 1, 6)

    assert calc == ans

def test_vector3_subtract():
    p1 = Vector3(1,2,3)
    p2 = Vector3(4,-1,3)

    calc = p1 - p2
    ans = Vector3(-3, 3, 0)

    assert calc == ans

def test_vector3_dot():
    p1 = Vector3(1,2,3)
    p2 = Vector3(4,-1,3)

    calc = p1.dot(p2)
    ans = 11

    assert calc == ans

def test_vector3_cross():
    p1 = Vector3(1,2,3)
    p2 = Vector3(4,-1,3)

    calc = p1.cross(p2)
    ans = Vector3(9,9,-9)

    assert calc == ans

def test_vector3_average():
    p1 = Vector3(1,2,3)
    p2 = Vector3(4,-1,3)

    calc = p1.average(p2)
    ans = Vector3(2.5,0.5,3)

    assert calc == ans

def test_vector3_unit_vector():
    p1 = Vector3(1,2,3)

    calc = p1.unit_vector()
    ans = Vector3(0.26726,0.53452,0.80178)

    tolerance = .001
    assert abs(calc.x - ans.x) < tolerance
    assert abs(calc.y - ans.y) < tolerance
    assert abs(calc.z - ans.z) < tolerance







def test_triangle_creation():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    calc = [tri.p4, tri.p5, tri.p6]
    ans = [Vector3(0.65,.6,.85), Vector3(0.45,0.6,1), Vector3(.4,.3,.95)]

    tolerance = .001
    
    for i in range(len(calc)):
        assert abs(calc[i].x - ans[i].x) < tolerance
        assert abs(calc[i].y - ans[i].y) < tolerance
        assert abs(calc[i].z - ans[i].z) < tolerance


def test_dx_dy_shapes():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    dxs, dys = tri.shape_functions_dx_dy(0,0)

    assert len(dxs) == 6
    assert len(dys) == 6

def test_B_membrane():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    B_m = tri.get_membrane_B(0,0)

    assert B_m.shape == (3, 30)

def test_B_bending():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    B_b = tri.get_bending_B(0,0)

    assert B_b.shape == (3, 30)

def test_B_shear():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    B_s = tri.get_shear_B(0,0)

    assert B_s.shape == (2, 30)

def test_D_membrane():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    D_m = tri.get_membrane_D(0,0)

    assert D_m.shape == (3,3)

def test_D_bending():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    D_b = tri.get_bending_D(0,0)

    assert D_b.shape == (3,3)

def test_D_shear():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    D_s = tri.get_shear_D(0,0)

    assert D_s.shape == (2,2)

def test_get_Ke():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    Ke = tri.get_Ke()

    assert Ke.shape == (30,30)

def test_get_global_Ke():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    Keglobal = tri.get_Ke_global()
    # print(Keglobal)

    assert Keglobal.shape == (36,36)
    # assert False

def test_get_mass_matrix_M():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    M= tri.get_mass_matrix()
    # print(M)

    assert M.shape == (30,30)
    # assert False

def test_get_global_Me():
    p1 = Vector3(.6, .3, .8)
    p2 = Vector3(.7, .9, .9)
    p3 = Vector3(.2, .3, 1.1)

    tri = Triangle(p1, p2, p3)

    Meglobal = tri.get_Me_global()
    # print(Keglobal)

    assert Meglobal.shape == (36,36)
    # assert False