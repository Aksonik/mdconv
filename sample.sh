#!/bin/bash

./compile.sh
time ./mdconv -contacts 1 7 sample/sigmas.dat sample/sys.pdb sample/sys.dcd  

echo "consistent with processdcd if no different:"
diff pairs.dat sample/calc_with_processdcd/pairs.res
diff con.dat sample/calc_with_processdcd/con.res

