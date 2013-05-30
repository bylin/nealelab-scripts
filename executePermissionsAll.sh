#!/bin/bash
#gives user execute permission to python,perl and shell scripts
for py in *.py
do
	if [ -f $py ];  then
		chmod u+x $py
	fi	
done

for pl in *.pl
do
	if [ -f $pl ];  then
		chmod u+x $pl
	fi	
done

for sh in *.sh
do
	if [ -f $sh ];  then
		chmod u+x $sh
	fi	
done
