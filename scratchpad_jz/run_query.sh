#!/bin/sh
#$ -S /bin/bash
#$ -cwd
#$ -N run_query
#$ -q bigmem1.q
#$ -m abe
#$ -M jjzieve@ucdavis.edu

#log starting time
date


module load bigmem

#echo "select sum(len) from Vvinifera_censor c, repbase r where c.name=r.name and r.family='Copia';" | mysql --password="foo" loblolly_repeat_project > vv_copia
#echo "select sum(len) from Vvinifera_censor c, repbase r where c.name=r.name and r.family='Gypsy';" | mysql --password="foo" loblolly_repeat_project > vv_gypsy
#echo "select sum(len) from Vvinifera_censor c, repbase r where c.name=r.name and r.family='LINE';" | mysql --password="foo" loblolly_repeat_project > vv_line
echo "select seq_name,count(*) from Ptaeda_trf where sequence = 'TTTAGGG' group by seq_name;"| mysql --password="foo" loblolly_repeat_project > tel_scaffs

date

