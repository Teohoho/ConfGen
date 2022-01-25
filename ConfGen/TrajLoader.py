import mdtraj as md
import numpy as np

##  Function that handles the input of the files
def TrajectoryLoader(topology, coordinates):
    MDTrajObj = md.load(coordinates,top=topology)

    return (MDTrajObj)
