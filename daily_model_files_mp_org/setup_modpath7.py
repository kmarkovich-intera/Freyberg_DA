import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import platform
import flopy
import pyemu

model_ws=os.path.join(".")

#load model
sim = flopy.mf6.MFSimulation.load(sim_ws=model_ws)
m = sim.get_model("freyberg6")


#set up node information
def get_nodes(locs):
    nodes = []
    for k, i, j in locs:
        nodes.append(k * m.dis.nrow.data * m.dis.ncol.data + i * m.dis.ncol.data + j)
    return nodes

df = pd.read_csv("particle_cells.csv")
kij = df.loc[:,["row","column"]].apply(lambda x: [2,x.row-1,x.column-1],axis=1).values

nodes = get_nodes(kij)
print(nodes[0])

# create modpath files
mpnamf = "freyberg6" + '_mp_forward'

# create basic forward tracking modpath simulation
mp = flopy.modpath.Modpath7.create_mp7(modelname=mpnamf, trackdir='forward', flowmodel=m, model_ws=model_ws,nodes=nodes[0], 
                                       rowcelldivisions=1, columncelldivisions=1, layercelldivisions=1)

# write modpath datasets
mp.write_input()