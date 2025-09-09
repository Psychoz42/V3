from V3_classesSpace import *

class Cube:
    def __init__(self, cenPos:Position, size:float, angleVec:Vector, color:list, reflectionK:float):
        self.cenPos = cenPos
        self.size = size
        self.angleVec = angleVec
        self.color = color

        v1 = Vertex(Position(round(self.cenPos.x - 0.5 * self.size, 5), round(self.cenPos.y - 0.5 * self.size, 5), round(self.cenPos.z - 0.5 * self.size, 5)))
        v2 = Vertex(Position(round(self.cenPos.x + 0.5 * self.size, 5), round(self.cenPos.y - 0.5 * self.size, 5), round(self.cenPos.z - 0.5 * self.size, 5)))
        v3 = Vertex(Position(round(self.cenPos.x + 0.5 * self.size, 5), round(self.cenPos.y + 0.5 * self.size, 5), round(self.cenPos.z - 0.5 * self.size, 5)))
        v4 = Vertex(Position(round(self.cenPos.x - 0.5 * self.size, 5), round(self.cenPos.y + 0.5 * self.size, 5), round(self.cenPos.z - 0.5 * self.size, 5)))
        v5 = Vertex(Position(round(self.cenPos.x - 0.5 * self.size, 5), round(self.cenPos.y - 0.5 * self.size, 5), round(self.cenPos.z + 0.5 * self.size, 5)))
        v6 = Vertex(Position(round(self.cenPos.x + 0.5 * self.size, 5), round(self.cenPos.y - 0.5 * self.size, 5), round(self.cenPos.z + 0.5 * self.size, 5)))
        v7 = Vertex(Position(round(self.cenPos.x + 0.5 * self.size, 5), round(self.cenPos.y + 0.5 * self.size, 5), round(self.cenPos.z + 0.5 * self.size, 5)))
        v8 = Vertex(Position(round(self.cenPos.x - 0.5 * self.size, 5), round(self.cenPos.y + 0.5 * self.size, 5), round(self.cenPos.z + 0.5 * self.size, 5)))

        Polygon(v1, v2, v4, color, reflectionK)
        Polygon(v2, v3, v4, color, reflectionK)
        Polygon(v1, v5, v6, color, reflectionK)
        Polygon(v6, v2, v1, color, reflectionK)
        Polygon(v5, v6, v7, color, reflectionK)
        Polygon(v7, v8, v5, color, reflectionK)
        Polygon(v1, v4, v8, color, reflectionK)
        Polygon(v8, v5, v1, color, reflectionK)
        Polygon(v2, v6, v7, color, reflectionK) #-
        Polygon(v7, v3, v2, color, reflectionK) #-
        Polygon(v3, v4, v8, color, reflectionK)
        Polygon(v8, v7, v3, color, reflectionK)

class Cuboid:
    def __init__(self, cenPos:Position, Xsize:float, Ysize:float, Zsize:float, angleVec:Vector, color:list, reflectionK:float):
        self.cenPos = cenPos
        self.Xsize = Xsize
        self.Ysize = Ysize
        self.Zsize = Zsize
        self.angleVec = angleVec
        self.color = color

        v1 = Vertex(Position(round(self.cenPos.x - 0.5 * self.Xsize, 5), round(self.cenPos.y - 0.5 * self.Ysize, 5), round(self.cenPos.z - 0.5 * self.Zsize, 5)))
        v2 = Vertex(Position(round(self.cenPos.x + 0.5 * self.Xsize, 5), round(self.cenPos.y - 0.5 * self.Ysize, 5), round(self.cenPos.z - 0.5 * self.Zsize, 5)))
        v3 = Vertex(Position(round(self.cenPos.x + 0.5 * self.Xsize, 5), round(self.cenPos.y + 0.5 * self.Ysize, 5), round(self.cenPos.z - 0.5 * self.Zsize, 5)))
        v4 = Vertex(Position(round(self.cenPos.x - 0.5 * self.Xsize, 5), round(self.cenPos.y + 0.5 * self.Ysize, 5), round(self.cenPos.z - 0.5 * self.Zsize, 5)))
        v5 = Vertex(Position(round(self.cenPos.x - 0.5 * self.Xsize, 5), round(self.cenPos.y - 0.5 * self.Ysize, 5), round(self.cenPos.z + 0.5 * self.Zsize, 5)))
        v6 = Vertex(Position(round(self.cenPos.x + 0.5 * self.Xsize, 5), round(self.cenPos.y - 0.5 * self.Ysize, 5), round(self.cenPos.z + 0.5 * self.Zsize, 5)))
        v7 = Vertex(Position(round(self.cenPos.x + 0.5 * self.Xsize, 5), round(self.cenPos.y + 0.5 * self.Ysize, 5), round(self.cenPos.z + 0.5 * self.Zsize, 5)))
        v8 = Vertex(Position(round(self.cenPos.x - 0.5 * self.Xsize, 5), round(self.cenPos.y + 0.5 * self.Ysize, 5), round(self.cenPos.z + 0.5 * self.Zsize, 5)))

        Polygon(v1, v2, v4, color, reflectionK)
        Polygon(v2, v3, v4, color, reflectionK)
        Polygon(v1, v5, v6, color, reflectionK)
        Polygon(v6, v2, v1, color, reflectionK)
        Polygon(v5, v6, v7, color, reflectionK)
        Polygon(v7, v8, v5, color, reflectionK)
        Polygon(v1, v4, v8, color, reflectionK)
        Polygon(v8, v5, v1, color, reflectionK)
        Polygon(v2, v6, v7, color, reflectionK)
        Polygon(v7, v3, v2, color, reflectionK)
        Polygon(v3, v4, v8, color, reflectionK)
        Polygon(v8, v7, v3, color, reflectionK)