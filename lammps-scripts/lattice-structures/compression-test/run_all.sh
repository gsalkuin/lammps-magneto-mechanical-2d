#!/usr/bin/env bash

Bs=(0 20 40 60 80 100)
for B in ${Bs[@]}
do
	rsync -av --delete template/ B${B}mT
	fname=square-achiral-compact-3x3-B${B}mT-relaxed.lam
	cp lam-files/$fname B${B}mT/ 
	cd B${B}mT/
	sed -i "s/BMT/$B/g" in.compress
	sed -i "s/BMT/$B/g" lammps-slurm.sb
	sbatch lammps-slurm.sb	
	cd ..			
done

