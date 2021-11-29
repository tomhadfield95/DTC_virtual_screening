#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 14:27:18 2021

@author: hadfield
"""
import argparse
from ccdc.docking import Docker
from ccdc.io import MoleculeReader, EntryReader
import tempfile
from rdkit import Chem
import glob
from pathlib import Path
from os import mkdir, chdir
from time import time
from dataclasses import dataclass, field
from multiprocessing import Pool
import sys
import os


def dock_ligands_GOLD(ligandsToDockSDF, cavityLigandFile, protein, outputSDF, n_docks = 10, tmpOutputFile = 'docked_ligands.sdf', constraint = None, all_poses = False):
        
    #Method for the (rigid) docking of ligands, without Multiprocessing
    
    #Arguments:
    #ligandsToDockSDF   - an SDF file containing the ligands you wish to dock - should have been protonated in RDKit (or elsewhere)
    #cavityLigandFile   - an SDF file which allows GOLD to determine where the binding pocket is 
                     #  - if you are docking molecules into a structure which already had a ligand bound, use the SDF of that crystal ligand
    #protein            - a PDB file containing the protein you want to dock into
    #n_docks            - the number of poses you want GOLD to sample. More docks enables better conformational coverage at the expense of time
    #outputFile         - the name of the output file for the docked poses - the output of the function is the best pose of each ligand but they will be saved in tmp named as {outputFile} as well
    
    #constraint         - mol2 file for doing constrained docking in GOLD
    
    
    
    
    #Set up docking class
    docker = Docker()
    settings = docker.settings
    settings.add_protein_file(protein)
    protein = settings.proteins[0]

    crystal_ligand = MoleculeReader(cavityLigandFile)[0]
    crystal_ligand.identifier = 'crystalLigand' #Change the identifier of the cavity ligand
    settings.binding_site = settings.BindingSiteFromLigand(protein, crystal_ligand, 10.0) #Define the binding site

    #ligandsToDock = MoleculeReader(ligandsToDockSDF)
    settings.add_ligand_file(ligandsToDockSDF, n_docks) #Generate n_docks poses per ligand

    settings.fitness_function = 'plp'
    settings.autoscale = 10.
    settings.early_termination = False
    batch_tempd = tempfile.mkdtemp()
    settings.output_directory = batch_tempd
    settings.output_file = tmpOutputFile

    if constraint is not None:
        scaffold = MoleculeReader(constraint)[0]
        settings.add_constraint(settings.ScaffoldMatchConstraint(scaffold))

    #Dock the ligands
    results = docker.dock()

    #Get docked poses:
    docksLocation = results.ligands.file_name[0]
    topRanked = glob.glob(f'{docksLocation[0:docksLocation.index(tmpOutputFile)]}ranked*_m*_1.sdf')
    mols = [Chem.SDMolSupplier(fn)[0] for fn in topRanked]
    
    if not all_poses:
    #Now return the top ranked mols
        topRanked = sorted(glob.glob(f'{docksLocation[0:docksLocation.index(tmpOutputFile)]}ranked*_m*_1.sdf')) #sort the list for consistency                                                                                                        
        mols = [Chem.SDMolSupplier(fname)[0] for fname in topRanked]
    else:
        allPoses = sorted(glob.glob(f'{docksLocation[0:docksLocation.index(tmpOutputFile)]}ranked*_m*_*.sdf'))
        mols = [Chem.SDMolSupplier(fname)[0] for fname in allPoses]
    
    
    
    #Either return just the docks or the docks and associated fitness scores

    
        
    w = Chem.SDWriter(outputSDF)
    for m in mols:
        w.write(m)
    w.close()
    
    print(f'Written {len(mols)} docked ligands to file: {outputSDF}')
    
    
    return 0
    


def dock_ligands_GOLD_MP(ligandsToDockSDF, cavityLigandFile, protein, outputSDF, n_docks = 10, tmpOutputFile = 'docked_ligands.sdf', constraint = None,  all_poses = False, n_processes = 7):

    """
    Dock the molecules from the supplied input file in parallel.
    Adapted from the gold_multi_map.py file supplied by ccdc
    """
    
    protein = expand_path(protein)
    cavityLigandFile = expand_path(cavityLigandFile)
    ligandsToDockSDF = expand_path(ligandsToDockSDF)
    outputSDF = expand_path(outputSDF)
    

    t0 = time()  # Script start time


    if not n_processes > 0:
        print(f"Error! Number of processes must be an integer greater than zero.")
        sys.exit(1)

    #Get number of molecules to be docked
    with EntryReader(ligandsToDockSDF) as reader:
        n_molecules = len(reader)

    print(f"There are {n_molecules} molecules to dock on {n_processes} processes...")

    output_dir = tempfile.mkdtemp() #Create temporary output file

    #########################################

    # Determine the sets of parameters defining the chunks...
    chunks = []  # List of records that define the chunks

    # Work out the size of each chunk, and hence the start and finish indices in the input file...
    basic_size = n_molecules // n_processes  # Basic size of a chunk, which must obviously be integral
    remainder  = n_molecules %  n_processes  # Number of molecules that would not be included in basic-sized chunks

    finish = 0  # Finish index

    for n in range(1, n_processes + 1):  # Recall that the number of chunks is the same as the number of processes
        start = finish + 1  # Start index
        chunk_size = basic_size + 1 if n <= remainder else basic_size  # Add one to the basic chunk sizes until the remainder are taken care of
        finish = start + chunk_size - 1
        chunks.append(chunkForMP(n=n, start=start, finish=finish, output_dir=output_dir, crystal_ligand = cavityLigandFile, constraint = constraint, protein = protein, ligand_file = ligandsToDockSDF, n_docks=n_docks))
        print(f"chunk {n}: size: {chunk_size}; start: {start}, finish: {finish}.")

    ##########################################

    # Dock the chunks in parallel...
    with Pool(n_processes) as pool:
        _ = pool.map(do_chunk, chunks)  # No output; docks are saved in the output_dir directory

    #print(f"Finished docking in {time() - t0:.1f} seconds.")

    if not all_poses:
    #Now return the top ranked mols
        topRanked = sorted(glob.glob(f'{output_dir}/*/ranked*_m*_1.sdf')) #sort the list for consistency                                                                                                        
        mols = [Chem.SDMolSupplier(fname)[0] for fname in topRanked]
    else:
        allPoses = sorted(glob.glob(f'{output_dir}/*/ranked*_m*_*.sdf'))
        mols = [Chem.SDMolSupplier(fname)[0] for fname in allPoses]

    # All done.
    print(f"Finished in {time() - t0:.1f} seconds.")

    w = Chem.SDWriter(outputSDF)
    for m in mols:
        w.write(m)
    w.close()
    
    print(f'Written {len(mols)} docked ligands to file: {outputSDF}')

    
    return 0

####Additional Code for parallel docking####

@dataclass
class chunkForMP:

    #Contains the information to provide to GOLD for docking
    #Store a bunch of chunks in a list and use them to do the docking in parallel

    n: int       # Chunk number
    start: int   # Index of first molecule in chunk
    finish: int  # Index of last molecule in chunk
    n_docks: int  # Number of docking poses to generate
    #conf_file: Path   # GOLD configuration file
    output_dir: Path  # Output dir, in which the chunk sub-directory will be created
    crystal_ligand: Path #Location of the crystal ligand (used for binding site definition)
    constraint: Path #Location of mol2 file used to constrain the docking
    protein: Path #Location of protein pdb File
    ligand_file: Path #Location of ligands to be docked

    dir: Path = field(init=False) # Sub-directory for chunk, see __post_init__ below

    def __post_init__(self):
        self.dir = f'{self.output_dir}/chunk_{self.n:02d}'

def do_chunk(chunk):

    """
    Dock a chunk of the input file.

    :param chunk: a record holding the parameters defining the chunk

    As we can't return a GOLD results object from a pool process (it can't be pickled as it wraps C++ objects),
    we simply return a boolean recording whether GOLD exited with a 'success' status code.
    """


    #Set up docking class
    docker = Docker()
    settings = docker.settings
    settings.add_protein_file(chunk.protein)

    settings.fitness_function = 'plp'
    settings.autoscale = 10.
    settings.early_termination = False
#     batch_tempd = tempfile.mkdtemp()
#     settings.output_directory = batch_tempd
#     settings.output_file = outputFile


    # Create and enter the sub-directory for this chunk...
    mkdir(chunk.dir)
    chdir(chunk.dir)
    settings.output_directory = '.'  # Ensure GOLD writes output to the chunk sub-directory


    #Define protein, binding site, scaffold and ligands to dock
    protein = settings.proteins[0]
    crystal_ligand = MoleculeReader(chunk.crystal_ligand)[0]
    crystal_ligand.identifier = 'crystalLigand' #Change the identifier of the cavity ligand
    settings.binding_site = settings.BindingSiteFromLigand(protein, crystal_ligand, 10.0) #Define the binding site

    # Specify the chunk of molecules to dock...
#     ligand_file = settings.ligand_files[0]  # The ligand file info will be overwritten, so store for reference below
    settings.clear_ligand_files() #Clear just to be sure
    settings.add_ligand_file(chunk.ligand_file, ndocks=chunk.n_docks, start=chunk.start, finish=chunk.finish)

    #Set scaffold constraint:
    if chunk.constraint is not None:
        scaffold =  MoleculeReader(chunk.constraint)[0]
        settings.add_constraint(settings.ScaffoldMatchConstraint(scaffold))


    # Run docking...
    print(f"Starting chunk {chunk.n} (ligand indices {chunk.start} - {chunk.finish})...")
    docker = Docker(settings=settings)
    results = docker.dock()



    print(f"Finished chunk {chunk.n}.")
    # As we can't return the results (as they are not picklable) and the poses have already been written to disk, we just return the status code

    return results.ligands.file_name[0] #return docks location
                                                                     
def expand_path(p):
        return os.path.abspath(os.path.expanduser(p))



    
if __name__=='__main__':
    

    parser = argparse.ArgumentParser()
    parser.add_argument('ligands', type=str, help='SDF file containing ligands you wish to dock')
    parser.add_argument('protein', type = str, help = 'Protein PDB file')
    parser.add_argument('--output_file', '-o', type=str,
                        help='Location to save output (SDF file)')
    parser.add_argument('--cavity_lig_file', '-c',
                        help = 'Cavity Ligand File')
    parser.add_argument('--constraint_file', '-s', type = str, default = None,
                        help = 'Location of a mol2 file to constrain the docking')
    parser.add_argument('--n_docks', '-n', type = float, default = 10,
                        help = 'Number of poses per ligand')
    parser.add_argument('--multiprocessing', '-mp', type = int, default = -1,
                        help = 'Use multiprocessing: Specify the number of cores you wish to use. If -1, GOLD will just use a single core')
    parser.add_argument('--all', '-a', action = 'store_true',
                        help = 'Save all the poses generated by GOLD, rather than just the top-ranked one')
    


    arguments = parser.parse_args()
    
    num_cores = int(arguments.multiprocessing)
    
    if num_cores == -1:
        dock_ligands_GOLD(arguments.ligands, arguments.cavity_lig_file, arguments.protein, arguments.output_file, int(arguments.n_docks), constraint=arguments.constraint_file, all_poses=arguments.all)
    elif num_cores >= 1:
        dock_ligands_GOLD_MP(arguments.ligands, arguments.cavity_lig_file, arguments.protein, arguments.output_file, int(arguments.n_docks), constraint=arguments.constraint_file, all_poses=arguments.all, n_processes = num_cores)
    else:
        ValueError('Please specify a valid number of cores to use for the docking, or omit the multiprocessing argument to run the job on a single core')
    