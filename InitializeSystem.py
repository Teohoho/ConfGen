import numpy as np
from scipy.spatial.transform import Rotation as R
import mdtraj as md

import sys

def InitializeSys (Sys1, Sys2):
    
    """
    Function that generates the starting conformation of the two helices.

    Parameters
    ----------

    Sys1,Sys2: MDTrajTrajectory Object
                Trajectory objects of the two systems

    """

    ##  Get positions
    Sys1Pos= Sys1.xyz
    Sys2Pos= Sys2.xyz
    
#    print (Sys1Pos)
    
    ## Define the unit vector that defines the 
    ## Y axis
    YAxis = np.array([0,1,0])
    
    ##  Get central axis
    ##  Define end points, generate vectors for the rotation
    ##  axes, normalize them
    Sys1Point1 = md.compute_center_of_mass(Sys1, "(resid 0 to 3) and name CA")
    Sys1Point2 = md.compute_center_of_mass(Sys1, "(resid {} to {}) and name CA".format(Sys1.n_residues-4,Sys1.n_residues-1))
    Sys1CentralAxis = np.array(Sys1Point1) - np.array(Sys1Point2)
    Sys1CentralAxisNormalized = Sys1CentralAxis/np.linalg.norm(Sys1CentralAxis)
    
    Sys2Point1 = md.compute_center_of_mass(Sys2, "(resid 0 to 3) and name CA")
    Sys2Point2 = md.compute_center_of_mass(Sys2, "(resid {} to {}) and name CA".format(Sys1.n_residues-4,Sys1.n_residues-1))
    Sys2CentralAxis = np.array(Sys2Point1) - np.array(Sys2Point2)
    Sys2CentralAxisNormalized = Sys2CentralAxis/np.linalg.norm(Sys2CentralAxis)
    
    #print (Sys1CentralAxisNormalized)
    #print (Sys2CentralAxisNormalized)

    #print (Sys1CentralAxisNormalized)
    #print (np.linalg.norm(Sys1CentralAxisNormalized))
    #print (Sys1CentralAxisNormalized*2)
    #print (np.linalg.norm(Sys1CentralAxisNormalized*2))
    
    ## Compute the angle between the two vectors
   
    #print(YAxis)
    #print(Sys1CentralAxisNormalized)

    RotAng1 = np.arccos(np.dot(Sys1CentralAxisNormalized, YAxis))
    RotAng2 = np.arccos(np.dot(Sys2CentralAxisNormalized, YAxis))

    ## Generate a orthogonal vector to the YAxis-Sys1/2 plane, 
    ## in order to compute the rotation of Sys1/2 to align it to the Y Axis
    Sys1Orthogonal = np.cross(Sys1CentralAxisNormalized,YAxis)
    Sys2Orthogonal = np.cross(Sys2CentralAxisNormalized,YAxis)

    ## Normalize the above vector, so it's easier to work with
    Sys1Orthogonal = Sys1Orthogonal/np.linalg.norm(Sys1Orthogonal)
    Sys2Orthogonal = Sys2Orthogonal/np.linalg.norm(Sys2Orthogonal)

#    print (np.linalg.norm(Sys1Orthogonal), np.linalg.norm(Sys2Orthogonal))
#    print (np.linalg.norm(Sys1Orthogonal), np.linalg.norm(Sys2Orthogonal))

#    print (Sys1CentralAxisNormalized, Sys2CentralAxisNormalized) 
#    print (Sys1Orthogonal, Sys2Orthogonal) 
#    print ((RotAng1), (RotAng2))
#    print (np.degrees(RotAng1), np.degrees(RotAng2))

    ## Compute the rotation vectors
    rotVec1 = R.from_rotvec(Sys1Orthogonal*RotAng1)
    rotVec2 = R.from_rotvec(Sys2Orthogonal*RotAng2)

    ## Compute the rotation matrix, so it's easier to work with

    ## Apply rotation
    AlignedVector1 = rotVec1.apply(Sys1CentralAxisNormalized)
    AlignedVector2 = rotVec2.apply(Sys2CentralAxisNormalized)

    #print (AlignedVector1)
    #print (AlignedVector2)
