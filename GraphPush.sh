#!/bin/bash

find . -name "*~" -type f -delete

RUTA=`pwd`
#echo $RUTA
cd $RUTA
git-cola 
