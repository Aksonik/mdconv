# mdconv: CG-contacts

#### mdconv is extended by an option for determination of contacts between coarse-graned molecules.

#### Command
```
mdconv -contacts A D file.dat file.pdb file.dcd
```
*A* and *D* -- parameters in the contact criterion (see below). 

*file.dat* -- residue names of molecules and their sizes (*R<sub>c</sub>*) in [A].

#### Output

*con.dat* -- total number of contacts for each trajectory frame.

*pairs.dat* -- frame number, pair of molecules that form a contact.

#### Contact criterion

*d* = *A<sub>c</sub>* ( *R<sub>c</sub><sup>A</sup>* + *R<sub>c</sub><sup>B</sup>* ) 0.5 + *D<sub>c</sub>*

#### Remarks

See the sample.

The work was done at the Michigan State University.
