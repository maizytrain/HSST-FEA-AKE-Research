import numpy as np



class Vector3:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    
    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"
    
    def __repr__(self):
        return f"({self.x}, {self.y}, {self.z})"
    
    
    def __eq__(self, other):
        if isinstance(other, Vector3):
            return (self.x == other.x and self.y == other.y and self.z == other.z)
        else:
            raise TypeError("Unsupported operand type for ==")

    def __add__(self, other):
        if isinstance(other, Vector3):
            return Vector3(self.x + other.x, self.y + other.y, self.z + other.z)
        if isinstance(other, (int, float)):
            return Vector3(self.x + other, self.y + other, self.z + other)
        else:
            raise TypeError("Unsupported operand type for +")
        
    def __sub__(self, other):
        if isinstance(other, Vector3):
            return Vector3(self.x - other.x, self.y - other.y, self.z - other.z)
        if isinstance(other, (int, float)):
            return Vector3(self.x - other, self.y - other, self.z - other)
        else:
            raise TypeError("Unsupported operand type for -")
        
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Vector3(self.x * other, self.y * other, self.z * other)
        else:
            raise TypeError("Unsupported operand type for *")
        
    def dot(self, other):
        if isinstance(other, Vector3):
            return self.x * other.x + self.y * other.y + self.z * other.z
        else:
            raise TypeError("Unsupported operand type for dot product")
        
    def cross(self, other):
        if isinstance(other, Vector3):
            return Vector3(self.y * other.z - self.z * other.y, self.z * other.x - self.x * other.z, self.x * other.y - self.y * other.x)
        else:
            raise TypeError("Unsupported operand type for cross product")
        
    def average(self, other):
        if isinstance(other, Vector3):
            return Vector3(self.x + other.x, self.y + other.y, self.z + other.z) * .5
        else:
            raise TypeError("Unsupported operand type for average")
        
    def magnitude(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)
        
    def unit_vector(self):
        return self * (1 / self.magnitude())
    
    def give_np_array(self):
        return np.array([self.x, self.y, self.z])
    

def create_Vector3_from_numpy_array(array):
    return Vector3(array[0], array[1], array[2])