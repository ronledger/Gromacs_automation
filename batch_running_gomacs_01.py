#!/usr/bin/env python3 
#................................................................. GROMACS SIMULATION (2019).............................................................................

import os
import numpy as np
import matplotlib.pyplot as plt
import statistics

print("My code DOESN'T work, I have no idea why. My code WORKS, I have no idea why.")
rmsd_file= open ("Average_RMSD.txt",'w')
mindist_file = open("AVERAGE_MINDIST.txt",'w') 

# Before running Change the Directory according to your need......

directory = '/media/quinn/Joker/Gromacs_using_python/Delete/KCATBH00'

# For Moving in Directory with range 
for x in range (8, 11):
    os.chdir(directory + str(x) +'/')
    print(os.getcwd())

# Functions for GROMACS specified tasks:
    def run_gromacs_simulation():
        pdb_2gmx()
        add_ions()
        making_index()
        traj_cavity()
        rmsd_calculate()
        mindist_calculate()
        
# Function for creating topology    
    def pdb_2gmx():
       pdb2gmx="gmx pdb2gmx -f *.pdb -o conf.gro -ignh <<END \n"
       pdb2gmx+="1\n"
       pdb2gmx+="1\n"
       os.system(pdb2gmx)
       os.system("gmx editconf -f conf.gro -o pbc.gro -c -d 0.8")
       os.system("gmx solvate -cp pbc.gro -cs spc216.gro -o solvated.gro -p topol.top")
       os.system("gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr")

#   Function for Adding ions replacing solvent    
    def add_ions():
        genion="gmx genion -s ions.tpr -o protein_solv_ions.gro -p topol.top -pname K -nname CL -neutral <<END \n"
        genion+="15"
        os.system(genion)
    
        os.system("gmx grompp -f em.mdp -c protein_solv_ions.gro -p topol.top -o em.tpr -maxwarn 1")
        os.system("gmx mdrun -v -deffnm em")

# Function for making index of Protein and Ligand Molecule (BSS)    
    def making_index():
        make_index="gmx make_ndx -f em.tpr <<END \n"
        make_index+='"PROTEIN" | r BSS\n'
        make_index+="r CHO & a O1\n"
        make_index+="r BSS & a C22\n"
        make_index+="q\n"
        os.system(make_index)

# Checking Ion in system for removal of tc-grps coupling error        
        with open ('index.ndx') as ions_search:
            if 'Ion' in ions_search.read():
                print("Ions are added in the System")
                ions_search.close()
            else:
                print("Ions are NOT added in the System")
                ions_search.close()
                with open ('npt.mdp','r+') as npt:
                    p=npt.read()
                    if 'Protein_BSS Water_and_ions' in p:
                        p=p.replace('Protein_BSS Water_and_ions', 'Protein_BSS Water')
                        npt.close()
                        with open ('npt.mdp','w') as f1:
                            f1.write (p)
                            f1.close()
                with open ('nvt.mdp','r+') as nvt:
                    v=nvt.read()
                    if 'Protein_BSS Water_and_ions' in v:
                        v=v.replace('Protein_BSS Water_and_ions','Protein_BSS Water')
                        nvt.close()
                        with open ('nvt.mdp','w') as f2:
                            f2.write(v)
                            f2.close()
                with open ('md.mdp', 'r+') as md:
                    m=md.read()
                    if 'Protein_BSS Water_and_ions' in m:
                        m=m.replace('Protein_BSS Water_and_ions', 'Protein_BSS Water')
                        md.close()
                        with open ('md.mdp', 'w') as f3:
                            f3.write(m)
                            f3.close()
    
# NVT NPT Equilibration followed by the MD RUN with four threads *** Change threads according to your need        
        os.system("gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1")
        os.system("gmx mdrun -v -deffnm nvt -nt 4")
        os.system("gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 1")
        os.system("gmx mdrun -v -deffnm npt -nt 4")
        os.system("gmx grompp -f md.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o md.tpr -maxwarn 1")
        os.system("gmx mdrun -v -deffnm md -nt 4")
        
# Removing PBC boundary condition with nojump option      
    def traj_cavity():
        trajconv ="gmx trjconv -f md.xtc -s md.tpr -o md_no_pbc.xtc -n index.ndx -pbc nojump -center <<END \n"
        trajconv += "1\n"
        trajconv +="1\n"
        os.system(trajconv)

# RMSD Calculation     
    def rmsd_calculate():
        rmsd= "gmx rms -n index.ndx -s md.tpr -f md.xtc -o rmsd_protein.xvg -tu ns -fit rot+trans <<END \n"
        rmsd+="1\n"
        rmsd+="1\n"
        os.system(rmsd)
        
        col1, col2=[],[]
        with open ('rmsd_protein.xvg', 'r') as rmsd:
            line_19_to_end = rmsd.readlines()[18:]
            for line in line_19_to_end:
                line=list(map(float,line.split()))
                col1.append(line[0])
                col2.append(line[1])
            me=statistics.mean(col2)
            sd=statistics.stdev(col2)
            plt.plot(col1,col2)          #RMSD Plotting using matplotlib
            plt.title("RMSD")
            plt.xlabel("Time (ns)")
            plt.ylabel("RMSD (nm)")
            plt.savefig("RMSD.png")
            plt.close()
            rmsd_file.write("KCATBH00")
            rmsd_file.write(str(x))
            rmsd_file.write("        ")
            rmsd_file.write(str(me))
            rmsd_file.write("        ")
            rmsd_file.write(str(sd))
            rmsd_file.write(str("\n"))
            rmsd.close()
            print("Always code as if the guy who ends up maintaining your code will be a violent psychopath who knows where you live. - Martin Golding")
#   Calculation of mindist 
    
    def mindist_calculate():
        ions_search= open ('index.ndx')
        if 'Ion'in ions_search.read():
            print ("Ions Are PRESENT in the System Group 23 & 24 is Chosen")
            mindist="gmx mindist -f md.xtc -s md.tpr -od mindist.xvg -n index.ndx <<END \n"
            mindist+="23\n"
            mindist+="24\n"
            os.system(mindist)
            
        else:
            print("Ions are NOT Present in the System Group 17 & 18 is Chosen")
        
            def mindist_calculate():
                mindist="gmx mindist -f md.xtc -s md.tpr -od mindist.xvg -n index.ndx <<END \n"
                mindist+="17\n"
                mindist+="18\n"
                os.system(mindist)

        column_1, column_2=[],[]   
        with open ('mindist.xvg','r') as mindist:
            line_25_to_end=mindist.readlines()[24:]
            for line in line_25_to_end:
                line=list(map(float,line.split()))
                column_1.append(line[0])
                column_2.append(line[1])
            m=statistics.mean(column_2)
            s=statistics.stdev(column_2)
            plt.plot(column_2)
            plt.title("MINDIST CHO_O1 & BSS_C22")
            plt.xlabel("Time(ns)")
            plt.ylabel("Distance (nm)")
            plt.savefig("Minimum_Distance.png")
            plt.close()
            mindist_file.write("KCATBH00")
            mindist_file.write(str(x))
            mindist_file.write("         ")
            mindist_file.write(str(m))
            mindist_file.write("         ")
            mindist_file.write(str(s))
            mindist_file.write("\n")
            mindist.close()
    
    run_gromacs_simulation()

rmsd_file.close()
mindist_file.close()

print("Programming is like Sex. One mistake you have to support it for the rest of your life. - Michael Sinz")
   
#...........................................................................................................................................................................
