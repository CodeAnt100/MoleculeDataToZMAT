# MoleculeDataToZMAT
This contains the start of a program to convert molecule data to a ZMAT format.

## Method
The program works by taking data about the bond angles and bond lengths for a molecule and determining a possible XYZ coordinate configuration for the atoms in the molecule.
The program tries to do this alongside iteratively constructing the Z-matrix, although it would also be possible to create an XYZ coordinate configuration as a standalone output, which could then be fed into existing XYZ-to-ZMAT converters.
