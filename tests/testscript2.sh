#!/bin/bash
for i in `seq -w 1.0 0.1 2.2`
do

cat > tmp.inp <<EOF
Atoms
    Cl   0.0000  0.0000  0.0000
    H  0.0000  0.0000  $i
End

Basis
    STO-3G
End

EOF

cd ..
go run . tests/tmp.inp >> tests/log2

cd tests

done

