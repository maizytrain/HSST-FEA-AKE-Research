from .vector3 import Vector3




class Node:
    def __init__(self, position:Vector3, index:int=-1, simulation=None):
        self.position = position
        self.index = index
        self.fixed = [False, False, False, False, False, False]

    def fix_degree_of_freedom(self, index):
        self.fixed[index] = True
    
    def release_degree_of_freedom(self, index):
        self.fixed[index] = False

    def fix_all(self):
        for i in range(len(self.fixed)):
            self.fixed[i] = True

    def release_all(self):
        for i in range(len(self.fixed)):
            self.fixed[i] = False

    