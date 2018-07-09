import os
import csv
from math import sqrt
from operator import itemgetter
from datetime import timedelta


def atom_charge(symbol):
    periodic = ('X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
    # periodic += ('Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg')
    atoms = {v.lower():i for i,v in enumerate(periodic)}
    return atoms[symbol.lower()]

def get_atomic_symbols(molecule):
    """Get atomic symbols set"""
    with open(os.path.join('..', 'chem_database', molecule + '.xyz'), 'r') as input_geometry:
        natoms = int(input_geometry.readline())
        charge, mult = map(int, input_geometry.readline().split())
        atomic_symbols = set()
        for atom in range(natoms):
            symbol, x, y, z = input_geometry.readline().split()
            atomic_symbols.add(symbol)
    return atomic_symbols

def get_XYZ(molecule):
    """Load XYZ-geometry from file."""
    with open(os.path.join('..', 'chem_database', molecule + '.xyz'), 'r') as input_geometry:
        natoms = int(input_geometry.readline())
        charge, mult = map(int, input_geometry.readline().split())
        geometry = []
        for atom in range(natoms):
            symbol, x, y, z = input_geometry.readline().split()
            geometry.append((atom_charge(symbol), map(float, (x, y, z))))
    return geometry

def get_max_Z(molecule):
    """Get maximal Z for atoms in molecule."""
    return max(Z for Z, _ in get_XYZ(molecule))

def get_ae_cutoffs(molecule):
    """Create AE_cutoff initial values.
    Used for Backflow format.
    """
    result = []
    for i, _ in enumerate(get_XYZ(molecule)):
        result.append('{i}         {i}         0.5                          0'.format(i=i+1))
    return '\n  '.join(result)

def get_atom_labels(molecule):
    """Returns number of atoms in a set and list of labels for this set.
    Used for generic JASTROW format.
    """
    result = []
    for i, _ in enumerate(get_XYZ(molecule)):
        result.append('{i}'.format(i=i+1))
    return i+1, ' '.join(result)

def  casino_time(*path):
    """Get CASINO time.
     Total CASINO CPU time  : : :      378.0500
    """
    regexp = re.compile(' Total CASINO CPU time  : : :\s+(?P<energy_error>\d+\.\d+)')
    try:
        with open(os.path.join(*path, 'out'), 'r') as casino_out:
            # we are only interested in the last occurrence
            return float(re.findall(regexp, casino_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None

def vmc_energy(*path):
    """Get VMC energy without JASTROW optimisation.
     -152.988424660763 +/- 0.003047553900      Correlation time method
    """

    regexp = re.compile(' (?P<energy>[-+]?\d+\.\d+) \+/- (?P<energy_error>[-+]?\d+\.\d+)      Correlation time method')
    try:
        with open(os.path.join(*path, 'out'), 'r') as vmc_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, vmc_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def vmc_variance(*path):
    """Get VMC variance with JASTROW optimisation.
      Sample variance of E_L (au^2/sim.cell) : 3.169677109628 +- 0.034986257092
    """

    regexp = re.compile('Sample variance of E_L \(au\^2/sim.cell\) : (?P<variance>[-+]?\d+\.\d+) \+- (?P<variance_error>[-+]?\d+\.\d+)')
    try:
        with open(os.path.join(*path, 'out'), 'r') as vmc_opt_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, vmc_opt_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def dmc_energy(*path):
    """Get DMC energy.
          mean:   -153.795024411601 +/-       0.001346260888
    """

    dir = os.path.join(*path)
    try:
        open(os.path.join(dir, '.casino_finished'), 'r').close()
        regexp = re.compile('mean:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
        with open(os.path.join(dir, 'out'), 'r') as dmc_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, dmc_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def dmc_stderr(*path):
    """Get DMC standard error.
          stderr:      0.000906128433 +/-       0.000046917552
    """

    dir = os.path.join(*path)
    try:
        open(os.path.join(dir, '.casino_finished'), 'r').close()
        regexp = re.compile('stderr:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
        with open(os.path.join(dir, 'out'), 'r') as dmc_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, dmc_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def dmc_ncorr(*path):
    """Get DMC correlation N.
          N_corr:      0.000906128433 +/-       0.000046917552
    """

    dir = os.path.join(*path)
    try:
        open(os.path.join(dir, '.casino_finished'), 'r').close()
        regexp = re.compile('N_corr:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
        with open(os.path.join(dir, 'out'), 'r') as dmc_out:
            # we are only interested in the last occurrence
            return map(float, re.findall(regexp, dmc_out.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None, None

def dmc_stats_nstep(*path):
    """Get DMC statistic accumulation steps.
          dmc_stats_nstep   : 96000
    """

    dir = os.path.join(*path)
    try:
        regexp = re.compile('dmc_stats_nstep   :\s+(?P<nstep>\d+)')
        with open(os.path.join(dir, 'input'), 'r') as dmc_input:
            # we are only interested in the last occurrence
            return int(re.findall(regexp, dmc_input.read())[-1])
    except (FileNotFoundError, IndexError) as e:
        print(e)
        return None

def get_all_inputs():
    "get file names of all *.xyz input files"
    return sorted((os.path.splitext(filename)[0] for filename in os.listdir(config['INPUTS_DIR']) if os.path.splitext(filename)[1] == '.xyz'))


wildcard_constraints:
    molecule='[-\w+=.]+',
    method='[-\w()]+',
    basis='[-\w]+',
    jastrow_type='[_\w]+',
    jastrow_rank='[_\w]+',
    backflow_rank='[_\w]+',
    jastrow_opt_method='\w+',

####################################################################################################################

rule HF_RESULTS:
    output: 'hf_results.csv'
    run:
        with open(output[0], 'w', newline='') as result_file:
            energy_data = csv.writer(result_file, dialect=csv.unix_dialect, quoting=csv.QUOTE_NONE)
            fieldnames = [
                'method', 'basis', 'molecule', 'hf_energy', 'hf_time',
            ]
            energy_data.writerow(fieldnames)
            for method in METHODS:
                for basis in BASES:
                    for molecule in MOLECULES:
                        path = (method, basis, molecule)
                        try:
                            energy_data.writerow((
                                *path,
                                hf_energy(*path),
                                hf_time(*path),
                            ))
                        except (FileNotFoundError, IndexError) as e:
                            print(e)

rule RESULTS:
    output: 'results.csv'
    run:
        with open(output[0], 'w', newline='') as result_file:
            writer = csv.writer(result_file, dialect=csv.unix_dialect, quoting=csv.QUOTE_NONE)
            fieldnames = [
                'method', 'basis', 'molecule', 'hf_energy', 'hf_time',
                'vmc_energy', 'vmc_energy_error', 'vmc_variance', 'vmc_variance_error', 'vmc_time',
                'jastrow_rank',
                'vmc_opt_energy', 'vmc_opt_energy_error', 'vmc_opt_variance', 'vmc_opt_variance_error', 'vmc_opt_time', 'vmc_opt_energy_time',
                'dmc_energy', 'dmc_energy_error', 'dmc_stderr', 'dmc_stderr_error', 'dmc_ncorr', 'dmc_ncorr_error', 'dmc_time'
            ]
            writer.writerow(fieldnames)
            vmc_path = ('VMC', '10000000')
            for method in METHODS:
                for basis in BASES:
                    for molecule in MOLECULES:
                        path = (method, basis, molecule)
                        for jastrow_rank in JASTROW_RANKS:
                            vmc_opt_path = ('VMC_OPT', 'emin', jastrow_rank, '10000')
                            vmc_opt_energy_path = ('VMC_OPT', 'emin', jastrow_rank, '1000000_9')
                            dmc_path = ('VMC_DMC', 'emin', jastrow_rank, 'tmax_2_1024_1')
                            try:
                                writer.writerow((
                                    *path,
                                    hf_energy(*path),
                                    hf_time(*path),
                                    *vmc_energy(*path, *vmc_path),
                                    *vmc_variance(*path, *vmc_path),
                                    casino_time(*path, *vmc_path),
                                    jastrow_rank,
                                    *vmc_energy(*path, *vmc_opt_energy_path),
                                    *vmc_variance(*path, *vmc_opt_energy_path),
                                    casino_time(*path, *vmc_opt_path),
                                    casino_time(*path, *vmc_opt_energy_path),
                                    *dmc_energy(*path, *dmc_path),
                                    *dmc_stderr(*path, *dmc_path),
                                    *dmc_ncorr(*path, *dmc_path),
                                    casino_time(*path, *dmc_path),
                                ))
                            except (FileNotFoundError, IndexError) as e:
                                print(e)

rule VMC_DMC_RUN:
    input:      '{path}/VMC_DMC/{jastrow_opt_method}/{jastrow_rank}/tmax_2_{nconfig}_1/input',
                '{path}/VMC_DMC/{jastrow_opt_method}/{jastrow_rank}/tmax_2_{nconfig}_1/gwfn.data',
                '{path}/VMC_DMC/{jastrow_opt_method}/{jastrow_rank}/tmax_2_{nconfig}_1/parameters.casl',
    output:     '{path}/VMC_DMC/{jastrow_opt_method}/{jastrow_rank}/tmax_2_{nconfig}_1/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'


rule VMC_DMC_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/out',
                '{method}/{basis}/{molecule}/VMC/10000000/out'
    output:     '{method}/{basis}/{molecule}/VMC_DMC/{jastrow_opt_method}/{jastrow_rank}/tmax_2_{nconfig}_1/input'
    params:
        dt_relative_step = 2.0,
        stderr = 0.001,
        # roughly
        # nstep ~ 1.5 * (E(VMC) - E(DMC))/(nconfig*dtdmc*stderr*stderr)
        # (E(VMC) - E(DMC)) ~ 6 * (E(HF) - E(VMC))
        magic_const = 1.5/6.0
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        hf, _ = vmc_energy(wildcards.method, wildcards.basis, wildcards.molecule, *('VMC', '10000000'))
        vmc, _ = vmc_energy(wildcards.method, wildcards.basis, wildcards.molecule, *('VMC_OPT', wildcards.jastrow_opt_method, wildcards.jastrow_rank, '1000000_9'))
        dtdmc = 1.0/(get_max_Z(wildcards.molecule)**2 * 3.0 * params.dt_relative_step)
        nstep = params.magic_const*(hf - vmc)/(int(wildcards.nconfig)*dtdmc*params.stderr*params.stderr)
        nstep=max(50000, int(round(nstep, -3)))
        if wildcards.basis.endswith('_PP'):
            tmove = 'T'
        else:
            tmove = 'F'
        with open(output[0], 'w') as f:
            f.write(open('../vmc_dmc.tmpl').read().format(
                neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep // 1000,
                tmove=tmove, backflow='F'
            ))

rule VMC_DMC_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_DMC/{jastrow_opt_method}/{jastrow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_DMC/{jastrow_opt_method}/{jastrow_rank}/tmax_2_{nconfig}_1/parameters.casl'
    run:
        shell('ln -s ../../../../VMC_OPT/{wildcards.jastrow_opt_method}/{wildcards.jastrow_rank}/10000/parameters.9.casl "{output}"')
        # workaround in multireference case
        source_path = os.path.join(wildcards.method, wildcards.basis, wildcards.molecule, 'VMC_OPT', wildcards.jastrow_opt_method, wildcards.jastrow_rank, '10000', 'correlation.out.9')
        target_path = os.path.join(os.path.dirname(output[0]), 'correlation.data')
        shell('[[ -e "{source_path}" ]] && ln -s ../../../../VMC_OPT/{wildcards.jastrow_opt_method}/{wildcards.jastrow_rank}/10000/correlation.out.9 "{target_path}"; exit 0')
        # workaround in pseudopotential
        if wildcards.basis.endswith('_PP'):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')


rule VMC_DMC_GWFN:
    input:      '{path}/gwfn.data',
    output:     '{path}/VMC_DMC/{jastrow_opt_method}/{jastrow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_OPT_ENERGY_RUN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/input',
                '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/gwfn.data',
                '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/parameters.casl'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_OPT_ENERGY_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../vmc_opt_energy.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule, backflow='F'))

rule VMC_OPT_ENERGY_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/out'
    output:     '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/parameters.casl'
    run:
        shell('ln -s ../10000/parameters.9.casl "{output}"')
        # workaround in multireference case
        source_path = os.path.join(wildcards.method, wildcards.basis, wildcards.molecule, 'VMC_OPT', wildcards.jastrow_opt_method, wildcards.jastrow_rank, '10000', 'correlation.out.9')
        target_path = os.path.join(os.path.dirname(output[0]), 'correlation.data')
        shell('[[ -e "{source_path}" ]] && ln -s ../10000/correlation.out.9 "{target_path}"; exit 0')
        # workaround in pseudopotential
        if wildcards.basis.endswith('_PP'):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')


rule VMC_OPT_ENERGY_GWFN:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/1000000_9/gwfn.data'
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_OPT_RUN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/input',
                '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/gwfn.data',
                '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/parameters.casl',
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_OPT_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../opt_plan/{}.tmpl'.format(wildcards.jastrow_opt_method)).read().format(
                neu=neu, ned=ned, nconfig=10000, molecule=wildcards.molecule
            ))

rule VMC_OPT_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/parameters.casl'
    run:
        with open(output[0], 'w') as f:
            f.write(open('../casl/{}.tmpl'.format(wildcards.jastrow_rank)).read())
        # workaround in multireference case
        source_path = os.path.join(wildcards.method, wildcards.basis, wildcards.molecule, 'correlation.data')
        target_path = os.path.join(os.path.dirname(output[0]), 'correlation.data')
        shell('[[ -e "{source_path}" ]] && ln -s ../../../../correlation.data "{target_path}"; exit 0')
        # workaround in pseudopotential
        if wildcards.basis.endswith('_PP'):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')


rule VMC_OPT_GWFN:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/{jastrow_rank}/10000/gwfn.data'
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_RUN:
    input:      '{path}/VMC/10000000/input'
    output:     '{path}/VMC/10000000/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC/10000000/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC/10000000/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../vmc.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule))
        # workaround in multireference case
        source_path = os.path.join(wildcards.method, wildcards.basis, wildcards.molecule, 'correlation.data')
        target_path = os.path.join(os.path.dirname(output[0]), 'correlation.data')
        shell('[[ -e "{source_path}" ]] && ln -s ../../correlation.data "{target_path}"; exit 0')
        # workaround in pseudopotential
        if wildcards.basis.endswith('_PP'):
            for symbol in get_atomic_symbols(wildcards.molecule):
                symbol = symbol.lower()
                shell('cd "$(dirname "{output}")" && ln -s ../../../../../../ppotential/DiracFock_AREP/{symbol}_pp.data')

rule VMC_GWFN:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC/10000000/gwfn.data'
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_DMC_BF_RUN:
    input:      '{path}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/input',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/correlation.data',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/parameters.casl',
    output:     '{path}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_DMC_BF_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data',
                '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/out',
                '{method}/{basis}/{molecule}/VMC/10000000/out'
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/input'
    params:
        dt_relative_step = 2.0,
        stderr = 0.001,
        # roughly
        # nstep ~ 1.5 * (E(VMC) - E(DMC))/(nconfig*dtdmc*stderr*stderr)
        # (E(VMC) - E(DMC)) ~ 6 * (E(HF) - E(VMC))
        magic_const = 1.5/6.0
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        hf, _ = vmc_energy(wildcards.method, wildcards.basis, wildcards.molecule, *('VMC', '10000000'))
        vmc, _ = vmc_energy(wildcards.method, wildcards.basis, wildcards.molecule, *('VMC_OPT_BF', wildcards.jastrow_opt_method, wildcards.jastrow_rank + '__' + wildcards.backflow_rank, '1000000_9'))
        dtdmc = 1.0/(get_max_Z(wildcards.molecule)**2 * 3.0 * params.dt_relative_step)
        nstep = params.magic_const*(hf - vmc)/(int(wildcards.nconfig)*dtdmc*params.stderr*params.stderr)
        nstep=max(50000, int(round(nstep, -3)))
        if wildcards.basis.endswith('_PP'):
            tmove = 'T'
        else:
            tmove = 'F'
        with open(output[0], 'w') as f:
            f.write(open('../vmc_dmc.tmpl').read().format(
                neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep // 1000,
                tmove=tmove, backflow='T'
            ))

rule VMC_DMC_BF_DATA_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/out',
                '{method}/{basis}/{molecule}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/correlation.data'
    shell:      'ln -s "../../../../VMC_OPT_BF/{wildcards.jastrow_opt_method}/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/10000/correlation.out.9" "{output}"'

rule VMC_DMC_BF_CASL_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/out',
                '{method}/{basis}/{molecule}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/parameters.casl'
    shell:      'ln -s "../../../../VMC_OPT_BF/{wildcards.jastrow_opt_method}/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/10000/parameters.9.casl" "{output}"'

rule VMC_DMC_BF_GWFN:
    input:      '{method}/{basis}/{molecule}/gwfn.data',
    output:     '{method}/{basis}/{molecule}/VMC_DMC_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_OPT_BF_ENERGY_RUN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/input',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/gwfn.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/correlation.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/parameters.casl'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_OPT_BF_ENERGY_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../vmc_opt_energy.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule, backflow='T'))

rule VMC_OPT_BF_DATA_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/out'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/correlation.data'
    shell:      'ln -s ../10000/correlation.out.9 "{output}"'

rule VMC_OPT_BF_CASL_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/out'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/parameters.casl'
    shell:      'ln -s ../10000/parameters.9.casl "{output}"'

rule VMC_OPT_BF_ENERGY_GWFN:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/1000000_9/gwfn.data'
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'

####################################################################################################################

rule VMC_OPT_BF_CASL_RUN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/input',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/gwfn.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/correlation.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/parameters.casl',
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/out'
    shell:      'cd "$(dirname "{output}")" && runqmc'

rule VMC_OPT_BF_INPUT:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/input'
    run:
        neu, ned = get_up_down(wildcards.method, wildcards.basis, wildcards.molecule)
        with open(output[0], 'w') as f:
            f.write(open('../opt_plan/{}.tmpl'.format(wildcards.jastrow_opt_method)).read().format(
                neu=neu, ned=ned, nconfig=10000, molecule=wildcards.molecule
            ))

rule VMC_OPT_BF_DATA_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/correlation.data'
    run:
        jastrow = wildcards.jastrow_rank.split('_')
        with open(output[0], 'w') as f:
            backflow = wildcards.backflow_rank.split('_')
            ae_cutoffs = get_ae_cutoffs(wildcards.molecule)
            number, labels = get_atom_labels(wildcards.molecule)
            template = '../backflow_eta_mu.tmpl' if backflow[2] == '00' else '../backflow_eta_mu_phi.tmpl'
            f.write(open(template).read().format(
                eta_term=backflow[0],
                mu_number_of_atoms=number, mu_atom_labels=labels, mu_term=backflow[1],
                phi_number_of_atoms=number, phi_atom_labels=labels, phi_term_eN=backflow[2][0], phi_term_ee=backflow[2][1],
                ae_cutoffs=ae_cutoffs))
        source_path = os.path.join(wildcards.method, wildcards.basis, wildcards.molecule, 'correlation.data')
        shell('[[ -e "{source_path}" ]] && cat "{source_path}" >> "{output}"; exit 0')


rule VMC_OPT_BF_CASL_JASTROW:
    input:      '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/gwfn.data'
    output:     '{method}/{basis}/{molecule}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/parameters.casl'
    run:
        jastrow = wildcards.jastrow_rank.split('_')
        with open(output[0], 'w') as f:
            f.write(open('../casl/{}.tmpl'.format(wildcards.jastrow_rank)).read())


rule VMC_OPT_BF_GWFN:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/{jastrow_rank}__{backflow_rank}/10000/gwfn.data'
    shell:      'mkdir -p "$(dirname "{output}")" && ln -rs "{input}" "{output}"'
