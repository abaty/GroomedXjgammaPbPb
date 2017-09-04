#!/bin/sh

#useful debug info
echo echo $HOSTNAME

#extract inputs and print them
tar -zxvf inputs.tar.gz
echo $(ls -la .)

#grid certificate
proxyFile=$(ls | grep x509)
export X509_USER_PROXY=$PWD/$proxyFile
echo $X509_USER_PROXY

#DON'T CHANGE ANYTHING ABOVE unless you know what you are doing

#your code here, and change output file
echo | awk -v j=$1 -v k=$2 -v l=$3 -v m=$4 '{print "./run.exe "j" "k" "l" "m}'
echo | awk -v j=$1 -v k=$2 -v l=$3 -v m=$4 '{print "./run.exe "j" "k" "l" "m}' | bash


echo "copying output..."
gfal-copy -f output_$1.root gsiftp://se01.cmsaf.mit.edu:2811//cms/store/user/$USER/condor_tests/
echo "Cleaning up..."
rm output_$1.root
echo "done!"

