#!/usr/bin/env python3
import requests
import io


class BaseLoader:

    base_path = 'http://www.tcm.phy.cam.ac.uk/~mdt26/pseudo_lib/'

    @staticmethod
    def atom_charge(symbol):
        periodic = ('X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
        periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
        periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
        periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
        periodic += ('Cs', 'Ba', 'La')
        periodic += ('Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu')
        periodic += ('Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')
        atoms = {v.lower(): i for i, v in enumerate(periodic)}
        return atoms[symbol.lower()]

    @staticmethod
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
        elif atom_charge <= 30:  # K-Zn
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

    def load_casino(self):

        for atom in self.periodic:
            source_path = self.base_path + atom.lower() + self.casino
            target_path = self.name + '/' + atom.lower() + '_pp.data'
            response = requests.get(source_path, stream=True)
            handle = open(target_path, "wb")
            for chunk in response.iter_content(chunk_size=512):
                if chunk:  # filter out keep-alive new chunks
                    handle.write(chunk)

    def load_orca(self):

        target_path = self.name + '_gaussian.bas'
        with open(target_path, 'w') as output_file:
            for atom in self.periodic:
                source_path = self.base_path + atom.lower() + self.pp_gaussian
                response = io.StringIO(requests.get(source_path).text)
                for line in response:
                    if line.startswith('Element'):
                        atomic_symbol = line.split()[2]
                        print('newecp {}'.format(atomic_symbol),
                              file=output_file)
                        ncore = self.core_charge(self.atom_charge(atomic_symbol))
                        print(' N_core {}'.format(ncore), file=output_file)
                        response.readline()
                        response.readline()
                        for i in range(3):
                            lmax = response.readline()[1]
                            if i == 0:
                                print(' lmax {}'.format(lmax),
                                      file=output_file)
                            num = int(response.readline()[1])
                            print(' {} {}'.format(lmax, num), file=output_file)
                            for j in range(num):
                                print(
                                    '  {i}  {data[1]:>20} {data[2]:>20} {data[0]:>4}'.format(
                                        i=j + 1,
                                        data=response.readline().split()),
                                    file=output_file)
                        print('end', file=output_file)

    def load_qchem(self):

        target_path = self.name + '_gaussian.bas'
        with open(target_path, 'w') as output_file:
            for atom in self.periodic:
                source_path = self.base_path + atom.lower() + self.pp_gaussian
                response = io.StringIO(requests.get(source_path).text)
                for line in response:
                    if line.startswith('Element'):
                        atomic_symbol = line.split()[2]
                        print(' {}'.format(atomic_symbol), file=output_file)
                        ncore = self.core_charge(self.atom_charge(atomic_symbol))
                        print(' {} {} {}'.format(self.name, 2, ncore), file=output_file)
                        response.readline()
                        response.readline()
                        for i in range(3):
                            lmax = response.readline()[1]
                            num = int(response.readline()[1])
                            if lmax == 'd':
                                print(' d potential', file=output_file)
                                print(' {}'.format(num-1), file=output_file)
                            elif lmax == 's':
                                print(' s-d potential'.format(lmax), file=output_file)
                                print(' {}'.format(num), file=output_file)
                            elif lmax == 'p':
                                print(' p-d potential'.format(lmax), file=output_file)
                                print(' {}'.format(num), file=output_file)
                            for j in range(num):
                                r, e, c = response.readline().split()
                                r = int(r)
                                if r == 1:
                                    # Coulombic term Z_val/r, where Z_val is the number of valence electrons
                                    continue
                                print(
                                    '  {data[0]:>2} {data[1]:>20} {data[2]:>20}'.format(
                                        data=(r-2, e, c)),
                                    file=output_file)
                        print('****', file=output_file)


class DiracFock_AREP(BaseLoader):

    def __init__(self):
        self.name = 'DiracFock_AREP'
        self.casino = '/1/pp.data'
        self.pp_gaussian = '/2/pp_gaussian'
        self.periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
        self.periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
        self.periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
        self.periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
        self.periodic += ('Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')


class HartreeFock(BaseLoader):

    def __init__(self):
        self.name = 'HartreeFock'
        self.casino = '/4/pp.data'
        self.pp_gaussian = '/5/pp_gaussian'
        self.periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
        self.periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
        self.periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')


class SofterDF_AREP(BaseLoader):

    def __init__(self):
        self.name = 'SofterDF_AREP'
        self.casino = '/1_2/pp.data'
        self.pp_gaussian = ''
        self.periodic = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
        self.periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
        self.periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
        self.periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',       'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
        self.periodic += ('Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')


class SmallCoreDF_AREP(BaseLoader):

    @staticmethod
    def core_charge(atom_charge):
        """
        :returns: valence charge of atom.
        """
        if 21 <= atom_charge <= 30:  # Sc-Zn
            return 18
        elif atom_charge == 42:  # Mo
            return 30
        elif atom_charge == 72:  # Hf
            return 62
        else:
            raise NotImplementedError("Small Core AREP Trail & Needs PP didn't support this element")

    def __init__(self):
        self.name = 'SmallCoreDF_AREP'
        self.casino = '/7/pp.data'
        self.pp_gaussian = '/8/pp_gaussian'
        self.periodic = ('Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn')
        self.periodic += ('Mo', )
        self.periodic += ('Hf', )


loader = DiracFock_AREP()
loader.load_orca()
loader.load_qchem()

loader = HartreeFock()
loader.load_orca()
loader.load_qchem()

loader = SmallCoreDF_AREP()
loader.load_orca()
loader.load_qchem()
