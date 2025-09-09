from V3_classesSpace import *
from tkinter import *
from tkinter import ttk
import os

class Scene:
    def __init__(self, pols:list, lights:list, sphs:list, cam:Camera, name:str):
        self.pols = pols
        self.lights = lights
        self.sphs = sphs
        self.cam = cam
        self.name = name

scenes = []

if not os.path.isdir('scenes'):
    os.mkdir('scenes')
else:
    files = os.listdir('scenes')
    for i in files:
        pols_ = []
        lights_ = []
        sphs_ = []
        cam_ = 0
        with open(f'scenes/{i}', 'r') as f:
            fLines = f.readlines()
            name_ = i[:-4]
            rType = 0 # 0 - pols, 1 - lights, 2 - cam
            for j in fLines:
                if j == 'pols\n':
                    rType = 0
                    continue
                if j == 'lights\n':
                    rType = 1
                    continue
                if j == 'sphs\n':
                    rType = 2
                    continue
                if j == 'cam\n':
                    rType = 3
                    continue

                match rType:
                    case 0:
                        v1 = Vertex(Position(float(j[0:j.find('a')]), float(j[j.find('a')+1:j.find('b')]), float(j[j.find('b')+1:j.find('c')])))
                        v2 = Vertex(Position(float(j[j.find('c')+1:j.find('d')]), float(j[j.find('d')+1:j.find('e')]), float(j[j.find('e')+1:j.find('f')])))
                        v3 = Vertex(Position(float(j[j.find('f')+1:j.find('g')]), float(j[j.find('g')+1:j.find('h')]), float(j[j.find('h')+1:j.find('i')])))
                        color = [int(j[j.find('i')+1:j.find('j')]), int(j[j.find('j')+1:j.find('k')]), int(j[j.find('k')+1:-1])]
                        pols_.append(Polygon(v1, v2, v3, color))
                    case 1:
                        pos = Position(float(j[0:j.find('a')]), float(j[j.find('a') + 1:j.find('b')]), float(j[j.find('b') + 1:j.find('c')]))
                        brightness = float(j[j.find('c') + 1:j.find('d')])
                        isDir = bool(j[j.find('d') + 1:j.find('e')])
                        dirVec = Vector(float(j[j.find('e')+1:j.find('f')]), float(j[j.find('f')+1:j.find('g')]), float(j[j.find('g')+1:-1]))
                        lights_.append(Light(pos, brightness, isDir, dirVec))
                    case 2:
                        pos = Position(float(j[0:j.find('a')]), float(j[j.find('a') + 1:j.find('b')]), float(j[j.find('b') + 1:j.find('c')]))
                        size = float(j[j.find('c') + 1:j.find('d')])
                        color = [float(j[j.find('d') + 1:j.find('e')]), float(j[j.find('e') + 1:j.find('f')]), float(j[j.find('f') + 1:-1])]
                        sphs_.append(Sphere(pos, size, color))
                    case 3:
                        pos = Position(float(j[0:j.find('a')]), float(j[j.find('a') + 1:j.find('b')]), float(j[j.find('b') + 1:j.find('c')]))
                        vec = Vector(float(j[j.find('c')+1:j.find('d')]), float(j[j.find('d')+1:j.find('e')]), float(j[j.find('e')+1:j.find('f')]))
                        mW = int(j[j.find('f')+1:j.find('g')])
                        mH = int(j[j.find('g')+1:j.find('h')])
                        FOV = int(j[j.find('h')+1:j.find('i')])
                        stFOV = int(j[j.find('i')+1:-1])
                        cam_ = Camera(pos, vec, mW, mH, FOV, stFOV)

        scenes.append(Scene(pols_, lights_, sphs_, cam_, name_))


root = Tk()
root.title("V3")
root.geometry("700x700+400+200")
root.minsize(700, 700)

polsVar = Variable()
lightsVar = Variable()

namesList = [i.name for i in scenes]
def sceneChange(event):
    global polsVar
    global lightsVar
    global  sphsVar
    #polsVar.set([f"Pol {n}, color {i.color}" for n, i in enumerate(scenes[namesList.index(nameVar.get())].pols)])
    pl = []
    lg = []
    sp = []
    for n, i in enumerate(scenes[namesList.index(nameVar.get())].pols):
        pl.append(f"Pol {n}, color {i.color}")
        pl.append(f'    v1: x={i.v1.pos.x}, y={i.v1.pos.y}, z={i.v1.pos.z}')
        pl.append(f'    v2: x={i.v2.pos.x}, y={i.v2.pos.y}, z={i.v2.pos.z}')
        pl.append(f'    v3: x={i.v3.pos.x}, y={i.v3.pos.y}, z={i.v3.pos.z}')
    polsVar.set(pl)

    for n, i in enumerate(scenes[namesList.index(nameVar.get())].lights):
        lg.append(f"Light {n}")
        lg.append(f'    pos: x={i.pos.x}, y={i.pos.y}, z={i.pos.z}')
        lg.append(f'    brightness = {i.brightness}')
        lg.append(f'    is direct = {i.isDir}')
        lg.append(f'    direction vector: vx={i.dirVec.vx}, vy={i.dirVec.vy}, vz={i.dirVec.vz}')
    lightsVar.set(lg)

    for n, i in enumerate(scenes[namesList.index(nameVar.get())].sphs):
        sp.append(f"Sphere {n}, color {i.color}")
        sp.append(f'    pos: x={i.pos.x}, y={i.pos.y}, z={i.pos.z}')
        sp.append(f'    size: {i.size}')
    sphsVar.set(sp)

nameVar = StringVar()

#nameLabel = ttk.Label(textvariable=namesVar)
#nameLabel.place(height=25, width=100, x=0, y=0)

scenesCom = ttk.Combobox(values=namesList, state="readonly", textvariable=nameVar)

ttk.Label(text="Scene:").place(height=25, width=50, x=2, y=2)
scenesCom.place(height=20, width=100, x=52, y=2)
scenesCom.bind("<<ComboboxSelected>>", sceneChange)

polsListb = Listbox(listvariable=polsVar)
polsListb.place(relheight=0.31, width=200, x=2, y=25)
polsListb_scrollbar = ttk.Scrollbar(orient="vertical", command=polsListb.yview)
polsListb_scrollbar.place(relheight=0.31, width=20, x=200, y=25)
polsListb["yscrollcommand"]=polsListb_scrollbar.set
polsListb_scrollbarx = ttk.Scrollbar(orient="horizontal", command=polsListb.xview)
polsListb_scrollbarx.place(height=20, width=200, x=2, rely=0.33)
polsListb["xscrollcommand"]=polsListb_scrollbarx.set

lightsListb = Listbox(listvariable=lightsVar)
lightsListb.place(relheight=0.3, width=200, x=2, rely=0.355)
lightsListb_scrollbary = ttk.Scrollbar(orient="vertical", command=lightsListb.yview)
lightsListb_scrollbary.place(relheight=0.3, width=20, x=200, rely=0.355)
lightsListb["yscrollcommand"]=lightsListb_scrollbary.set
lightsListb_scrollbarx = ttk.Scrollbar(orient="horizontal", command=lightsListb.xview)
lightsListb_scrollbarx.place(height=20, width=200, x=2, rely=0.655)
lightsListb["xscrollcommand"]=lightsListb_scrollbarx.set

root.mainloop()

