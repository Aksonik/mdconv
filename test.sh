#!/bin/bash

./compile.sh
time ./mdconv -contacts 1 7 test/sigmas.dat test/sys.pdb test/sys.dcd  

echo "consistent with processdcd if no different:"
diff pairs.dat test/calc_with_processdcd/pairs.res
diff con.dat test/calc_with_processdcd/con.res

