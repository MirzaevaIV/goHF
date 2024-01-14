#!/bin/bash
for i in `seq -w 1.0 0.1 2.2`
do

cat > tmp.inp <<EOF
Atoms
    Cl   0.0000  0.0000  0.0000
    H  0.0000  0.0000  $i
End

Basis
    def2-svp
End

EOF

cd ..
go run . tests/tmp.inp >> tests/log1

cd tests

done

