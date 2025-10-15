from .vector3 import Vector3




class Node:
    def __init__(self, position:Vector3, key:str="NoneAssigned", simulation=None):
        self.position = position
        self.key = key

    