import numpy as np
import matplotlib.pyplot as plt

def read_Etot(file):
    '''Extract the info about converged total energy from QE scf runs

    Parameters
    ----------
    file: str
        relative path to output file of scf run
    
    Returns
    ----------
    float
        converged total energy in eV
    '''
    E_tot=0
    read_obj = open(file,'r')
    for line in read_obj:
        if '!    total energy' in line:
            E_tot = float(line.split()[-2])
    if E_tot == 0:
        print('No converged total energy in this file :(')
    return E_tot

def read_fermi(file):
    ''' Extract the Fermi energy from QE scf run
    Note: In case of fixed occupations, no Fermi energy is estimated;
          highest occupied level is assumed as Fermi energy estimate.

    Parameters
    ----------
    file: str
        relative path to output file of scf run

    Returns
    ----------
    float
        Fermi energy in eV
    '''
    E_F=0
    read_obj = open(file,'r')
    for line in read_obj:
        if 'Fermi energy' in line:
            E_F = float(line.split()[-2])
        elif 'highest occupied level (ev):' in line:
            print('Fermi energy set to highest occupied level')
            E_F = float(line.split()[-1])
    return E_F

def read_bands(file_wo_ext):
    ''' Reads the .out.gnu file and formats the column of bands into a 2D-array
    
    Parameters
    ----------
    file: str
        realtive path to output of bands.x without extensions

    Returns
    ----------
    k_points: 1D-array, float
        reciprocal distance measured along the chosen k-path
    bands: 2D-array, float
        electronic bands in units of eV along the k-path
    '''
    head = np.loadtxt(file_wo_ext+'.out',max_rows=1,dtype='str')
    k_num = int(head[int(np.where(head=='nks=')[0][0])+1])
    k_points = np.loadtxt(file_wo_ext+'.out.gnu',usecols=0,max_rows=k_num)
    bands = np.loadtxt(file_wo_ext+'.out.gnu',usecols=1)
    bands=bands.reshape(-1,k_num).T
    return k_points,bands

def read_struct_in(in_file):
    ''' Extracts all given structural parameters from scf input file

    Parameters
    ----------
    file: str
        relative path to scf input file

    Returns
    ----------
    [A, B, C, cosAB, cosAC, cosBC]: list of float
        lattice parameters and spate angles
    atoms: list of str
        names of atom types
    positions: ndarray, float
        position of atoms in crystal coordinates
    '''
    content = open(in_file,'r')
    lines = content.readlines()
    content = open(in_file,'r')
    n_at = 0
    struct_i,l = 0,0
    A, B, C, cosAB, cosAC, cosBC = 0,0,0,0,0,0
    for line in content:
        l+=1
        if 'nat' in line:
            n_at = int(line.replace('\n','').replace(' ','').split('=')[-1])
        elif ('A' in line)&('=' in line):
            A = float(line.replace('A = ',''))
        elif ('B' in line)&('=' in line):
            B = float(line.replace('B = ',''))
        elif ('C' in line)&('=' in line):
            C = float(line.replace('C = ',''))
        elif 'cosAB = ' in line:
            cosAB = float(line.replace('cosAB = ',''))
        elif 'cosAC = ' in line:
            cosAC = float(line.replace('cosAC = ',''))
        elif 'cosBC = ' in line:
            cosBC = float(line.replace('cosBC = ',''))
        elif 'ATOMIC_POSITIONS' in line:
            struct_i = l
    atoms = []
    positions = []
    for l in lines[struct_i:struct_i+n_at]:
        p  = np.array(l.replace('\n','').split(' '))
        at = p[0]
        p  = p[[p!='' for p in p]][1:].astype(float)
        positions.append(np.array(p).astype(float))
        atoms.append(at)
    atoms = np.array(atoms)
    positions = np.array(positions)
    content.close()
    return [A,B,C,cosAB,cosAC,cosBC],atoms,positions

def legend_without_duplicate_labels(ax,**kwargs):
    ''' Reduces the labels put into legends to unique entries

    Parameters
    ----------
    ax: matplotlib.pyplot.ax
        axis object from matplotlib
    **kwargs: 
        further statements passed to legend function of matplotlib.pyplot
    '''
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique),**kwargs)
    return
