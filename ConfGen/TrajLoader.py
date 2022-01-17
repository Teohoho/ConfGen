import mdtraj as md
import numpy as np

##  Function that handles the input of the files
def TrajectoryLoader(topology, coordinates):
    MDTrajObj = md.load(coordinates,top=topology)

# Since we'll be using it later, we also want
# to print the largest distance between CA and
# any sidechain atom of the same residue
    #print (MDTrajObj.n_residues)
    atoms_CA = [x.index for x in MDTrajObj.topology.atoms_by_name("CA")]
    max_distance = 0
    for i in range(MDTrajObj.n_residues):
        atoms = [x.index for x in MDTrajObj.topology.residue(i).atoms]
        distancePairs = [[atoms_CA[i],x] for x in atoms]
        distances = md.compute_distances(MDTrajObj, distancePairs)
        if (np.max(distances) > max_distance):
            max_distance = np.max(distances)

    return (MDTrajObj,max_distance)
