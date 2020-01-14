#!/bin/bash

# Xin Qiao, Jan 14 2020

path=`pwd`
echo -e "The current directory will be added to your \$PATH environmental variable\n\n>>> $path\n"
sleep 1

echo Setting the \$PATH variable in the \~\/.bashrc file
echo -e "...\n"
echo -e "\n# added by DupGen_finder" >> ~/.bashrc
echo PATH=$path:\$PATH >> ~/.bashrc
#sed -i '$a PATH='$path':\$PATH' ~/.bashrc
sleep 1

echo Loading the new \$PATH into the current shell session
echo -e "...\n"
source ~/.bashrc
sleep 1

echo Done!
