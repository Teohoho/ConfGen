import mdtraj as md
import mdtraj.formats
import os

def writeTraj(positions, out_FN, topology=None, verbosity=True):
    """
    Function that takes an array of positions and generates
    an INPCRD file or a DCD trajectory file

    Parameters
    ----------
    positions:      numpy.ndarray

            if Array of shape (N,3), write rst7
            if Array of shape (M,N,3), write DCD with M frames

    out_FN:         str
            RST7/DCD output file name

    topology:       MDTrajTrajectory.Topology Object
            Topology of system. Needed for writing PDBs

    verbosity:      bool
            Print succesful writing message.

    Returns
    -------
    None

    Notes
    -----
    The INPCRD format states that each line has 6 fields, each one containing
    a float with 12 characters and 7 decimal places (12.7f)
    """

    # Open the out_FN file, after checking if it already exists
    #if (os.path.exists(out_FN)):
    #    raise ValueError("File {} already exists! Please choose a different name "
    #                    "or move to a new directory".format(out_FN))
    positions = positions * 10  # nm to Angstrom

    if (out_FN.split(".")[-1]) in ["rst7", "inpcrd"]:
        with mdtraj.formats.AmberRestartFile(out_FN, "w") as f:
            f.write(positions)
        if (verbosity):
            print ("AMBER RST7 file {} successfully written!".format(out_FN))

    elif (out_FN.split(".")[-1]) == "pdb":
        positions = positions.reshape((positions.shape[-2],positions.shape[-1]))
        currentFrame = md.formats.PDBTrajectoryFile(out_FN, "w")
        currentFrame.write(positions, topology=topology)
        if (verbosity):
            print("PDB file {} successfully written!".format(out_FN))

    elif (out_FN.split(".")[-1]) == "dcd":
        with md.formats.DCDTrajectoryFile(out_FN, "w") as f:
            f.write(positions)
        if (verbosity):
            print("CHARMM DCD file {} successfully written!".format(out_FN))




