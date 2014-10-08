#!/bin/bash

dir=`pwd`
ldir=$dir/joblogs
outfile=$2-$3-$4-output-ffrom$5-fto$6-zvtxmin$7-zvtxmax$8-nmin$9-nmax${10}-centmin${11}-centmax${12}
subfile=$dir/subfiles/subfile_$outfile

cat > $subfile <<EOF

Universe     = vanilla

Notification = Error
Requirements = (Machine != "hibat0009.cmsaf.mit.edu") && (Machine !="t2bat0236.cmsaf.mit.edu") && (Machine !="lnbat0002.cmsaf.mit.edu")
Executable   = $dir/$1
Arguments    = $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12}
Rank         = Mips
GetEnv       = True
priority     = 2000

#initialdir   = $dir
input        = /dev/null
output       = $ldir/$outfile.out
error        = $ldir/$outfile.err
log          = $ldir/$outfile.log

+AccountingGroup = "group_cmshi.$USER"

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT

Queue
EOF
