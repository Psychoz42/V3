import math

import numpy as np
from PIL import Image
import subprocess
import os
import time

polygons = []
spheres = []
lights = []
polTextures = []

class Position:
    def __init__(self, x:float, y:float, z:float):
        self.x = x
        self.y = y
        self.z = z

class Vector:
    def __init__(self, vx:float, vy:float, vz:float):
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.npVec = np.array([vx, vy, vz])


class Object:
    def __init__(self, pos:Position, vec:Vector):
        self.pos = pos
        self.vec = vec

class Ray:
    def __init__(self,pos:Position, vec:Vector):
        self.pos = pos
        self.vec = vec

class Vertex:
    def __init__(self, pos:Position):
        self.pos = pos
        #self.connections = []

class Plane:
    def __init__(self, v1:Vertex, v2:Vertex, v3:Vertex):
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3

        self.normalV = Vector(0, 0, 0)
        self.normalV.vx = np.linalg.det(np.array([[self.v2.pos.y - self.v1.pos.y, self.v3.pos.y - self.v1.pos.y],
                                                  [self.v2.pos.z - self.v1.pos.z, self.v3.pos.z - self.v1.pos.z]], dtype=float))
        self.normalV.vy = -1*np.linalg.det(np.array([[self.v2.pos.x - self.v1.pos.x, self.v3.pos.x - self.v1.pos.x],
                                                  [self.v2.pos.z - self.v1.pos.z, self.v3.pos.z - self.v1.pos.z]], dtype=float))
        self.normalV.vz = np.linalg.det(np.array([[self.v2.pos.x - self.v1.pos.x, self.v3.pos.x - self.v1.pos.x],
                                                  [self.v2.pos.y - self.v1.pos.y, self.v3.pos.y - self.v1.pos.y]], dtype=float))

        self.free = 0
        self.free += self.normalV.vx * -1*self.v1.pos.x
        self.free += self.normalV.vy * -1*self.v1.pos.y
        self.free += self.normalV.vz * -1*self.v1.pos.z


class Camera:
    def __init__(self, pos:Position, vec:Vector, mW, mH, FOV = 90):
        self.pos = pos
        self.vec = vecToUvec(vec)
        self.matrixW = mW
        self.matrixH = mH
        self.hFOVK = mH/mW
        self.FOV = FOV
        self.FOVk = math.tan(FOV * math.pi / 360)
        self.moveVec = vecToUvec(Vector(self.vec.vx, self.vec.vy, 0))

        self.xyPerPlane = Plane(Vertex(pos), Vertex(Position(pos.x, pos.y, pos.z + 1.0)), Vertex(Position(pos.x + vec.vx, pos.y + vec.vy, pos.z + vec.vz)))
        self.xyPerV = vecToUvec(self.xyPerPlane.normalV)
        self.zPerPlane = Plane(Vertex(pos), Vertex(Position(pos.x + vec.vx, pos.y + vec.vy, pos.z + vec.vz)), Vertex(Position(pos.x + self.xyPerV.vx, pos.y + self.xyPerV.vy, pos.z)))
        self.zPerV = vecToUvec(self.zPerPlane.normalV)

    def recalculateVecs(self):
        self.xyPerPlane = Plane(Vertex(self.pos), Vertex(Position(self.pos.x, self.pos.y, self.pos.z + 1.0)), Vertex(Position(self.pos.x + self.vec.vx, self.pos.y + self.vec.vy, self.pos.z + self.vec.vz)))
        self.xyPerV = vecToUvec(self.xyPerPlane.normalV)
        self.zPerPlane = Plane(Vertex(self.pos), Vertex(Position(self.pos.x + self.vec.vx, self.pos.y + self.vec.vy, self.pos.z + self.vec.vz)), Vertex(Position(self.pos.x + self.xyPerV.vx, self.pos.y + self.xyPerV.vy, self.pos.z)))
        self.zPerV = vecToUvec(self.zPerPlane.normalV)
        self.moveVec = vecToUvec(Vector(self.vec.vx, self.vec.vy, 0))
        #print(self.xyPerV.vx, self.xyPerV.vy, self.xyPerV.vz, '|', self.zPerV.vx, self.zPerV.vy, self.zPerV.vz)


    def render(self, fileName:str, printData = False, saveTxt = True, straightMode = False):
        if printData:
            print("Start rendering...")
        #startTime = time.time()

        input_data = load(self, saveTxt)

        process = subprocess.Popen(
            ["main.exe"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, shell=True
        )

        output, _ = process.communicate(input=input_data)
        if printData:
            print("Image was rendered, start filling png file...")

        img = Image.new("RGB", (self.matrixW, self.matrixH), (0, 0, 0))
        newPix = img.load()
        #print(output)

        x = 0
        y = 0
        outLines = output.split('\n')
        for i in outLines[:-1]:
            if i == '-':
                y+=1
                x = 0
                continue

            col = i.split('|')
            newPix[x, y] = (int(col[0]), int(col[1]), int(col[2]))
            x+=1

        if straightMode:
            if printData:
                print("Finished")
            return img
        else:
            img.save(fileName)
            if printData:
                print("Finished")



        #endTime = time.time()
        #print("Time:", round(endTime - startTime, 2), 'seconds')
        #print(output)



class Light:
    def __init__(self, pos:Position, brightness:float = 1.0, isDir:bool = False, dirVec:Vector = Vector(0, 0, -1)):
        self.pos = pos
        self. brightness = brightness
        self.isDir = isDir
        self.dirVec = dirVec

        lights.append(self)

class Polygon:
    def __init__(self, v1:Vertex, v2:Vertex, v3:Vertex, color:list, reflectionK:float): #, haveTexture:bool = False
        self.color = color
        self.pl = Plane(v1, v2, v3)
        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.reflectionK = reflectionK
        #self.haveTexture = haveTexture

        polygons.append(self)


#class polTexture:
#    def __init__(self, polIndex:int, texturePath:str, anchorDot:Position, anchorVec:Vector, size:float):
#        self.polIndex = polIndex
#        self.texturePath = texturePath
#        self.anchorDot = anchorDot
#        self.anchorVec = anchorVec
#        self.size = size
#        polygons[polIndex].haveTexture = True
#        polTextures.append((self))

class Sphere:
    def __init__(self, pos:Position, size:float, color:list):
        self.pos = pos
        self.size = size
        self.color = color

        spheres.append(self)


def plRes(v1:Vertex, v2:Vertex, v3:Vertex, pos:Position):
    pl = np.array([
        [pos.x - v1.pos.x, v2.pos.x - v1.pos.x, v3.pos.x - v1.pos.x],
        [pos.y - v1.pos.y, v2.pos.y - v1.pos.y, v3.pos.y - v1.pos.y],
        [pos.z - v1.pos.z, v2.pos.z - v1.pos.z, v3.pos.z - v1.pos.z]
    ], dtype=float)
    return np.linalg.det(pl)

def vecsDeg(vec1:Vector, vec2:Vector):
    return round(np.arccos(np.dot(vec1.npVec, vec2.npVec) / (np.linalg.norm(vec1.npVec) * np.linalg.norm(vec2.npVec))), 6)

def makeVec(pos1:Position, pos2:Position):
    return Vector(pos2.x - pos1.x, pos2.y - pos1.y, pos2.z - pos1.z)

def dist(pos1:Position, pos2:Position):
    return np.linalg.norm(makeVec(pos1, pos2).npVec)

def lenV(vec:Vector):
    return (vec.vx * vec.vx + vec.vy * vec.vy + vec.vz * vec.vz) ** 0.5

def vecToUvec(vec:Vector):
    vecL = lenV(vec)
    return Vector(round(vec.vx / vecL, 6), round(vec.vy / vecL, 6), round(vec.vz / vecL, 6))

def load(cam:Camera, saveTxt):
    data = str()
    if len(polygons) > 0:
        data += "pols\n"
        for i in polygons:
            data += f"{i.v1.pos.x} {i.v1.pos.y} {i.v1.pos.z} {i.v2.pos.x} {i.v2.pos.y} {i.v2.pos.z} {i.v3.pos.x} {i.v3.pos.y} {i.v3.pos.z} {i.color[0]} {i.color[1]} {i.color[2]} {i.reflectionK}\n" # {i.haveTexture}
    data += "lights\n"
    for i in lights:
        data += f"{i.pos.x} {i.pos.y} {i.pos.z} {i.brightness} {int(i.isDir)} {i.dirVec.vx} {i.dirVec.vy} {i.dirVec.vz}\n"
    if len(spheres) > 0:
        data += "sphs\n"
        for i in spheres:
            data += f"{i.pos.x} {i.pos.y} {i.pos.z} {i.size} {i.color[0]} {i.color[1]} {i.color[2]}\n"
    #if len(polTextures) > 0:
    #    data += "polT\n"
    #    for i in polTextures:
    #        data += f"{i.polIndex} {i.anchorDot.x} {i.anchorDot.y} {i.anchorDot.z} {i.anchorVec.vx} {i.anchorVec.vy} {i.anchorVec.vz} {i.size}\n"
    #        with Image.open(f"textures/{i.texturePath}") as img:
    #            imgArr = img.load()
    #            xL = img.size[0]
    #            yL = img.size[1]
    #            for y in range(yL):
    #                for x in range(xL):
    #                    data += f"{imgArr[x, y][0]} {imgArr[x, y][1]} {imgArr[x, y][2]} "
    #                data += "| "
    #            data += '\n'

    data += "cam\n"
    data += f"{cam.pos.x} {cam.pos.y} {cam.pos.z} {cam.vec.vx} {cam.vec.vy} {cam.vec.vz} {cam.matrixW} {cam.matrixH} {cam.FOV}\n"

    if saveTxt:
        with open("scene.txt", "w") as file:
            file.write(data)

    return data

def pngDecode(imgPath:str):
    with Image.open(imgPath) as img:
        pixels = img.load()
        with open(imgPath[:-4]+"_dpng"+".txt", 'w') as file:
            for y in range(img.size[1]):
                for x in range(img.size[0]):
                    file.write(f"{pixels[x, y][0]} {pixels[x, y][1]} {pixels[x, y][2]}\n")
                file.write("-\n")

def pngDecodeDir(directory:str):
    for f in os.listdir(directory):
        pngDecode(directory + '/' + f)