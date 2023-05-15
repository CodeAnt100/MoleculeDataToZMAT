# Created by: Anant Gupta
# N.B. This code is still buggy, but contains the backbone for the process of constructing a Z-matrix in this way.

import math
import numpy as np
from scipy.optimize import fsolve


def calc_bond_angle(p0, p1, p2):
    """
    Calculates the angle between a set of three cartesian points, where p1 is the point in common.

    :param p0: A numpy array for the first point
    :param p1: A numpy array for the second point (shared)
    :param p2: A numpy array for the third point
    :return: The bond angle in degrees
    """

    # Defines the normalised vectors for the two bond directions
    u = p0 - p1
    u /= np.linalg.norm(u)

    v = p2 - p1
    v /= np.linalg.norm(v)

    # Returns the angle between the vector directions
    return np.degrees(np.arccos(np.dot(u, v)))


def dihedral_angle(p0, p1, p2, p3):
    """
    Calculates the dihedral angle for a set of four cartesian points, where p1 and p2 are the points in common.

    :param p0: A numpy array for the first point
    :param p1: A numpy array for the second point (shared)
    :param p2: A numpy array for the third point (shared)
    :param p3: A numpy array for the fourth point
    :return: The dihedral angle in degrees
    """

    # Defines the normalised vectors for the intersection line and the perpendicular half-plane directions
    n = p2 - p1
    b0 = p0 - p1
    b1 = p3 - p2

    u = b0 - n * np.dot(b0, n)
    u /= np.linalg.norm(u)

    v = b1 - n * np.dot(b1, n)
    v /= np.linalg.norm(v)

    # Returns the angle between the half-plane directions
    return np.degrees(np.arccos(np.dot(u, v)))


def z_matrix(atoms, bond_lengths, bond_angles):
    """
    Constructs a z-matrix from a given collection of atoms, bond lengths and bond angles.
    The function sequentially goes through each of the atoms and adds their corresponding row to the Z-matrix.

    :param atoms: A list of strings for the atoms of the molecule
    :param bond_lengths: A dictionary of bond lengths for the molecule, from string pairs to floats
    :param bond_angles: A dictionary of bond angles for the molecule, from string triples to floats
    :return: A string containing the molecule data in Z-matrix form
    """

    # EDIT: Modify to deal with creating bond angles that don't already exist
    # EDIT: Be careful of string identifiers for atoms and the actual coordinates for them
    # Defines an inner function specific to bond angles
    def get_bond_angle(p0, p1, p2):
        """
        Defines a function to retrieve the bond angle of a set of 3 points, or to calculate it if not known.

        :param p0: A numpy array for the first point
        :param p1: A numpy array for the second point (shared)
        :param p2: A numpy array for the third point
        :return: A tuple of the bond angle array and the current calculated bond angle
        """

        # Checks if the point is in the bond angle array
        if (p0, p1, p2) in bond_angles:
            return bond_angles[(p0, p1, p2)]

        # Ensures there are bonds connecting p0 to p1 and connecting p1 to p2
        elif tuple(sorted((p0, p1))) in bond_lengths and tuple(sorted((p1, p2))) in bond_lengths:
            bondAngle = calc_bond_angle(p0, p1, p2)

            # Adds the bond angle to the bond angle array and returns the value
            bond_angles[(p0, p1, p2)] = bondAngle
            bond_angles[(p2, p1, p0)] = bondAngle
            return bondAngle

        # Throws an error if the required bonds do not exist
        else:
            print(f"Attempted to find bond angle for {p0}, {p1}, {p2} but angle does not exist.")
            return 0

    # Formats the inputs
    bond_lengths = dict(map(lambda kv: (tuple(sorted(kv[0])), kv[1]), bond_lengths.items()))
    bond_angles.update(dict(map(lambda kv: (tuple(reversed(kv[0])), kv[1]), bond_angles.items())))

    # Formats the inputs so that the first atom is connected to both the second and third atoms
    if tuple(sorted((atoms[0], atoms[2]))) not in bond_lengths:
        tempAtom = atoms[0]
        atoms[0] = atoms[2]
        atoms[2] = tempAtom

    # Initialises the XYZ-coordinate dictionary for the first 3 coordinates
    atom_coordinates = {atoms[0]: np.array([0, 0, 0]),
                        atoms[1]: np.array([bond_lengths[tuple(sorted((atoms[0], atoms[1])))], 0, 0])}
    tempAngle = bond_angles[atoms[1], atoms[0], atoms[2]]
    atom_coordinates[atoms[2]] = np.array([math.cos(np.radians(tempAngle)), math.sin(np.radians(tempAngle)), 0])
    atom_coordinates[atoms[2]] *= bond_lengths[tuple(sorted((atoms[0], atoms[2])))]

    # Initialises the z-matrix output for the first 3 coordinates
    zMatrix = f"{atoms[0]}\n"
    zMatrix += f"{atoms[1]}    1    {bond_lengths[tuple(sorted((atoms[0], atoms[1])))]}\n"
    zMatrix += f"{atoms[2]}    1    {bond_lengths[tuple(sorted((atoms[0], atoms[2])))]}    2    {tempAngle}\n"

    # Loops through each of the following atoms for which a dihedral angle entry is required
    for currAtom in atoms[3:]:

        # 1. Find a set of 3 atoms for which a dihedral angle can be defined
        #    We need to find an atom bonded to part of a bonded pair for which there exists another atom bonded
        #    to an atom for which their coordinates are already known.
        #    EDIT: This doesn't work for 180 degree bond angles (only one possibility)

        # Initialises the atom names for easier reference for calculating the dihedral angle
        currSharedAtom = None
        farSharedAtom = None
        farAtom = None

        for loopAtom1, loopAtom2 in bond_lengths:
            if loopAtom1 in atom_coordinates and loopAtom2 in atom_coordinates:

                for sharedAtom1, sharedAtom2 in [(loopAtom1, loopAtom2), (loopAtom2, loopAtom1)]:
                    if tuple(sorted((currAtom, sharedAtom1))) in bond_lengths:

                        for branchAtom1, branchAtom2 in bond_lengths:
                            if branchAtom1 != currAtom and branchAtom2 != currAtom:

                                if sharedAtom1 == branchAtom1 and branchAtom2 in atom_coordinates:
                                    currSharedAtom = sharedAtom1
                                    farSharedAtom = sharedAtom2
                                    farAtom = branchAtom2
                                    break

                                elif sharedAtom1 == branchAtom2 and branchAtom1 in atom_coordinates:
                                    currSharedAtom = sharedAtom1
                                    farSharedAtom = sharedAtom2
                                    farAtom = branchAtom1
                                    break

        # Moves the current atom to end of the atom array if a suitable set of 3 atoms is not found
        # EDIT: Iterated structure is spliced and hence is not dynamically updated
        if currSharedAtom is None:
            atoms.append(currAtom)
            continue

        # 2. Use simultaneous equations with bond angles to find XYZ coordinates for atom
        #    For this, we will use SciPy's fsolve() for non-linear simultaneous equations

        # Gets the coordinates for each of the chosen set of 3 atoms
        currSharedAtomCoord = atom_coordinates[currSharedAtom]
        farSharedAtomCoord = atom_coordinates[farSharedAtom]
        farAtomCoord = atom_coordinates[farAtom]

        # Defines the equations which constrain the position of the current atom
        # EDIT: potential problem with assuming farAtom will always attach to currSharedAtom (and not farSharedAtom)
        def equations(p):

            # TESTING
            # print("u:", p - currSharedAtomCoord)
            # print("u len:", np.linalg.norm(p - currSharedAtomCoord))
            # print("v:", farSharedAtomCoord - currSharedAtomCoord)
            # print("v len:", np.linalg.norm(farSharedAtomCoord - currSharedAtomCoord))
            # print("w:", farAtomCoord - currSharedAtomCoord)
            # print("w len:", np.linalg.norm(farAtomCoord - currSharedAtomCoord), "\n")

            eq1 = np.linalg.norm(currSharedAtomCoord - p) - bond_lengths[tuple(sorted((currAtom, currSharedAtom)))]

            u = p - currSharedAtomCoord

            # Casting error occurs if syntactic sugar version used
            u = u / np.linalg.norm(u)

            v = farSharedAtomCoord - currSharedAtomCoord
            v = v / np.linalg.norm(v)

            w = farAtomCoord - currSharedAtomCoord
            w = w / np.linalg.norm(w)

            eq2 = np.degrees(np.arccos(np.dot(u, v))) - bond_angles[(currAtom, currSharedAtom, farSharedAtom)]
            eq3 = np.degrees(np.arccos(np.dot(u, w))) - bond_angles[(currAtom, currSharedAtom, farAtom)]

            return np.array([eq1, eq2, eq3])

        # Solves the equations numerically and adds the coordinates to the coordinate dictionary
        currAtomCoord = fsolve(equations, np.array([1, 1, 1]))
        atom_coordinates[currAtom] = currAtomCoord

        # 3. Use dihedral angle formula to construct and add the corresponding Z-matrix row
        currentDihedralAngle = dihedral_angle(atom_coordinates[currAtom],
                                              atom_coordinates[currSharedAtom],
                                              atom_coordinates[farSharedAtom],
                                              atom_coordinates[farAtom])

        zMatrix += f"{currAtom}    {atoms.index(currSharedAtom) + 1}    {bond_lengths[tuple(sorted((currAtom, currSharedAtom)))]}"
        zMatrix += f"    {atoms.index(farSharedAtom) + 1}    {bond_angles[(currAtom, currSharedAtom, farSharedAtom)]}"
        zMatrix += f"    {atoms.index(farAtom) + 1}    {currentDihedralAngle}\n"

    # TESTING
    print(atom_coordinates)

    return zMatrix


# TESTING
# AsH2O
atomArray = ["As", "O", "H1", "H2"]

AsHLength = 1.5134
AsOLength = 1.672011

HAsHAngle = 101.84
HAsOAngle = 106.65

bondLengthDict = {
    ("As", "O"): AsOLength,
    ("As", "H1"): AsHLength,
    ("As", "H2"): AsHLength,
}

bondAngleDict = {
    ("H1", "As", "O"): HAsOAngle,
    ("H2", "As", "O"): HAsOAngle,
    ("H1", "As", "H2"): HAsHAngle,
}

zMatrixString = z_matrix(atomArray, bondLengthDict, bondAngleDict)
print(zMatrixString)

