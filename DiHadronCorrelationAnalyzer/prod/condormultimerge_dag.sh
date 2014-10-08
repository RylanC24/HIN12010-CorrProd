#!/bin/bash

dir=`pwd`
ldir=$dir/joblogs
outfile=$2-$3-$4-final
subfile=./subfiles/subfile_$outfile

echo $subfile

cat > $subfile <<EOF

Universe     = vanilla

Notification = Error
Executable   = $dir/$1
Arguments    = $2 $3 $4
Rank         = Mips
GetEnv       = True
priority     = 20000

initialdir   = $dir
input        = /dev/null
output       = $ldir/$outfile.out
error        = $ldir/$outfile.err
log          = $ldir/$outfile.log

+AccountingGroup = "group_cmshi.$USER"

#transfer_input_files = 
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
notification = Complete
notify_user = davidlw@rice.edu

Queue
EOF
