#!/usr/bin/env python3
import requests
import io


def atom_charge(symbol):
    periodic = ('X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
    periodic += ('Cs', 'Ba', 'La')
    periodic += ('Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu')
    periodic += ('Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')
    atoms = {v.lower():i for i,v in enumerate(periodic)}
    return atoms[symbol.lower()]

def core_charge(atom_charge):
    """
    :returns: valence charge of atom.
    """
    if atom_charge <= 2:     # H-He
        return 0
    elif atom_charge <= 10:  # Li-Ne
        return 2
    elif atom_charge <= 18:  # Na-Ar
        return 10
    elif atom_charge <= 26:  # K-Zn
        return 18
    elif atom_charge <= 36:  # Ga-Kr
        return 28
    elif atom_charge <= 48:  # Rb-Cd
        return 36
    elif atom_charge <= 54:  # In-Xe
        return 46
    elif atom_charge <= 71:  # Cs-Lu
        return 54
    elif atom_charge <= 80:  # Hf-Hg
        return 68
    else:
        raise NotImplementedError("AREP Trail & Needs PP didn't support elements after Hg")


base_path = 'http://www.tcm.phy.cam.ac.uk/~mdt26/pseudo_lib/'

def DiracFock_AREP():
    periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
    periodic += ('Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')
#    for atom in periodic:
#        source_path = base_path + atom.lower() + '/1/pp.data'
#        target_path = 'DiracFock_AREP/' + atom.lower() + '_pp.data'
#        response = requests.get(source_path, stream=True)
#        handle = open(target_path, "wb")
#        for chunk in response.iter_content(chunk_size=512):
#            if chunk:  # filter out keep-alive new chunks
#                handle.write(chunk)
    target_path = 'DiracFock_AREP_gaussian.bas'
    with open(target_path, 'w') as output_file:
        for atom in periodic:
            source_path = base_path + atom.lower() + '/2/pp_gaussian'
            response = io.StringIO(requests.get(source_path).text)
            for line in response:
                if line.startswith('Element'):
                    atomic_symbol = line.split()[2]
                    print('newecp {}'.format(atomic_symbol), file=output_file)
                    ncore = core_charge(atom_charge(atomic_symbol))
                    print(' N_core {}'.format(ncore), file=output_file)
                    line = response.readline()
                    line = response.readline()
                    for i in range(3):
                        lmax = response.readline()[1]
                        if i == 0:
                            print(' lmax {}'.format(lmax), file=output_file)
                        num = int(response.readline()[1])
                        print(' {} {}'.format(lmax, num), file=output_file)
                        for j in range(num):
                            print('  {i}  {data[1]:>20} {data[2]:>20} {data[0]:>4}'.format(i=j+1, data=response.readline().split()), file=output_file)
                    print('end', file=output_file)


def HartreeFock():
    periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    for atom in periodic:
        source_path = base_path + atom.lower() + '/4/pp.data'
        target_path = 'HartreeFock/' + atom.lower() + '_pp.data'
        response = requests.get(source_path, stream=True)
        handle = open(target_path, "wb")
        for chunk in response.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)

def SofterDF_AREP():
    periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',       'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
    periodic += ('Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')
    for atom in periodic:
        source_path = base_path + atom.lower() + '/1_2/pp.data'
        target_path = 'SofterDF_AREP/' + atom.lower() + '_pp.data'
        response = requests.get(source_path, stream=True)
        handle = open(target_path, "wb")
        for chunk in response.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)

def SmallCoreDF_AREP():
    periodic += ('Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn')
    periodic += ('Mo', )
    periodic += ('Hf', )
    for atom in periodic:
        source_path = base_path + atom.lower() + '/7/pp.data'
        target_path = 'SmallCoreDF_AREP/' + atom.lower() + '_pp.data'
        response = requests.get(source_path, stream=True)
        handle = open(target_path, "wb")
        for chunk in response.iter_content(chunk_size=512):
            if chunk:  # filter out keep-alive new chunks
                handle.write(chunk)

DiracFock_AREP()
