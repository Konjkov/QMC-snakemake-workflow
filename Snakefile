import os
from math import sqrt
from operator import itemgetter
from datetime import timedelta


INPUTS_DIR = '../chem_database'

def atom_charge(symbol):
    periodic = ('X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne')
    periodic += ('Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar')
    periodic += ('K', 'Ca', 'Sc', 'Ti', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr')
    periodic += ('Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe')
    atoms = {v:i for i,v in enumerate(periodic)}
    return atoms[symbol]

def get_max_Z(gwfn_file):
    """Get maximal Z for atoms in molecule."""
    with open(gwfn_file, "r") as gwfn:
        result = ''
        line = gwfn.readline()
        while line and not line.startswith('Atomic numbers for each atom:'):
            line = gwfn.readline()
        if not line:
            raise Exception
        line = gwfn.readline()
        while line and not line.startswith('Valence charges for each atom:'):
            result += line
            line = gwfn.readline()
    return max(map(int, result.split()))

def get_atom_list(molecule):
    with open(os.path.join('..', 'chem_database', molecule + '.xyz'), 'r') as input_geometry:
        result = dict.fromkeys(ATOMS, 0)
        for line in input_geometry:
            if line.startswith(' '):
                atom_symbol = line.split()[0].lower()
                result[atom_symbol] += 1
        return result

def get_ae_cutoffs(molecule):
    with open(os.path.join('..', 'chem_database', molecule + '.xyz'), 'r') as input_geometry:
        i = 0
        result = []
        for line in input_geometry:
            if line.startswith(' '):
                i += 1
                result.append('{i}         {i}         0.2                          1'.format(i=i))
        return '\n  '.join(result)

def get_atom_labels(molecule):
    """Returns number of atoms in a set and list of labels for this set.
    Used for generic JASTROW format.
    """
    with open(os.path.join('..', 'chem_database', molecule + '.xyz'), 'r') as input_geometry:
        i = 0
        result = []
        for line in input_geometry:
            if line.startswith(' '):
                i += 1
                result.append('{i}'.format(i=i))
        return i, ' '.join(result)

def vmc_energy(molecule, method, basis):
    """Get VMC energy without JASTROW optimisation.
     -152.988424660763 +/- 0.003047553900      Correlation time method
    """

    regexp = re.compile(' (?P<energy>[-+]?\d+\.\d+) \+/- (?P<energy_error>[-+]?\d+\.\d+)      Correlation time method')
    with open(os.path.join(molecule, method, basis, 'VMC', '10000000', 'out'), 'r') as vmc_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, vmc_out.read())[-1])
    return value, error

def vmc_opt_energy(molecule, method, basis, opt_method, rank):
    """Get VMC energy with JASTROW optimisation.
     -153.693436512511 +/- 0.003006326588      Correlation time method
    """

    regexp = re.compile(' (?P<energy>[-+]?\d+\.\d+) \+/- (?P<energy_error>[-+]?\d+\.\d+)      Correlation time method')
    with open(os.path.join(molecule, method, basis, opt_method, 'emin', 'casl', rank, '1000000_9', 'out'), 'r') as vmc_opt_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, vmc_opt_out.read())[-1])
    return value, error

def vmc_opt_variance(molecule, method, basis):
    """Get VMC variance with JASTROW optimisation.
      Sample variance of E_L (au^2/sim.cell) : 3.169677109628 +- 0.034986257092
    """

    regexp = re.compile('Sample variance of E_L \(au\^2/sim.cell\) : (?P<variance>[-+]?\d+\.\d+) \+- (?P<variance_error>[-+]?\d+\.\d+)')
    with open(os.path.join(molecule, method, basis, 'VMC_OPT', 'emin', 'casl', '8_8_44', '1000000_9', 'out'), 'r') as vmc_opt_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, vmc_opt_out.read())[-1])
    return value, error


def dmc_energy(molecule, method, basis):
    """Get DMC energy.
          mean:   -153.795024411601 +/-       0.001346260888
    """

    dir = os.path.join(molecule, method, basis, 'VMC_DMC', 'emin', 'casl', '8_8_44', 'tmax_2_1024_1')
    open(os.path.join(dir, '.casino_finished'), 'r').close()
    regexp = re.compile('mean:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
    with open(os.path.join(dir, 'out'), 'r') as dmc_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, dmc_out.read())[-1])
    return value, error

def dmc_stderr(molecule, method, basis):
    """Get DMC standard error.
          stderr:      0.000906128433 +/-       0.000046917552
    """

    dir = os.path.join(molecule, method, basis, 'VMC_DMC', 'emin', 'casl', '8_8_44', 'tmax_2_1024_1')
    open(os.path.join(dir, '.casino_finished'), 'r').close()
    regexp = re.compile('stderr:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
    with open(os.path.join(dir, 'out'), 'r') as dmc_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, dmc_out.read())[-1])
    return value, error

def dmc_ncorr(molecule, method, basis):
    """Get DMC correlation N.
          N_corr:      0.000906128433 +/-       0.000046917552
    """

    dir = os.path.join(molecule, method, basis, 'VMC_DMC', 'emin', 'casl', '8_8_44', 'tmax_2_1024_1')
    open(os.path.join(dir, '.casino_finished'), 'r').close()
    regexp = re.compile('N_corr:\s+(?P<energy>[-+]?\d+\.\d+) \+/- \s+(?P<energy_error>[-+]?\d+\.\d+)')
    with open(os.path.join(dir, 'out'), 'r') as dmc_out:
        # we are only interested in the last occurrence
        value, error = map(float, re.findall(regexp, dmc_out.read())[-1])
    return value, error

def dmc_stats_nstep(molecule, method, basis):
    """Get DMC statistic accumulation steps.
          dmc_stats_nstep   : 96000
    """

    dir = os.path.join(molecule, method, basis, 'VMC_DMC', 'emin', 'casl', '8_8_44', 'tmax_2_1024_1')
    regexp = re.compile('dmc_stats_nstep   :\s+(?P<nstep>\d+)')
    with open(os.path.join(dir, 'input'), 'r') as dmc_input:
        # we are only interested in the last occurrence
        value = int(re.findall(regexp, dmc_input.read())[-1])
    return value

def get_all_inputs():
    "get file names of all *.in input files"
    return [os.path.splitext(filename)[0] for filename in os.listdir(INPUTS_DIR) if os.path.splitext(filename)[1] == '.xyz']


wildcard_constraints:
    molecule='[-\w]+',
    method='[-\w()]+',
    basis='[-\w]+',
    jastrow_type='[_\w]+',
    jastrow_rank='[_\w]+',
    backflow_rank='[_\w]+',
    jastrow_opt_method='\w+',
    order='_\d+'

####################################################################################################################

rule VMC_DMC_RUN:
    input:      '{path}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/input',
                '{path}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/gwfn.data',
                '{path}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/parameters.casl',
    output:     '{path}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/out'
    shell:      'cd {wildcards.path}/VMC_DMC/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}/tmax_2_1024_1 && runqmc'


rule VMC_DMC_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/.keep',
                '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/out',
                '{molecule}/{method}/{basis}/VMC/10000000/out'
    output:     '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/input'
    params:
        dt_relative_step = 2.0,
        stderr = 0.001,
        # roughly
        # nstep ~ 1.5 * (E(VMC) - E(DMC))/(nconfig*dtdmc*stderr*stderr)
        # (E(VMC) - E(DMC)) ~ 6 * (E(HF) - E(VMC))
        magic_const = 1.5/6.0
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            gwfn_file = os.path.join(wildcards.molecule, wildcards.method, wildcards.basis, 'gwfn.data')
            hf, _ = vmc_energy(wildcards.molecule, wildcards.method, wildcards.basis)
            vmc, _ = vmc_opt_energy(wildcards.molecule, wildcards.method, wildcards.basis, 'VMC_OPT', wildcards.jastrow_rank)
            dtdmc = 1.0/(get_max_Z(gwfn_file)**2 * 3.0 * params.dt_relative_step)
            nstep = params.magic_const*(hf - vmc)/(int(wildcards.nconfig)*dtdmc*params.stderr*params.stderr)
            nstep=max(50000, int(round(nstep, -3)))
            with open(file_name, 'w') as f:
                f.write(open('../vmc_dmc.tmpl').read().format(
                    neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep // 1000
                ))

rule VMC_DMC_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/out',
                '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/parameters.casl'
    shell:      'ln -s ../../../../../VMC_OPT/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}/10000/parameters.9.casl {output}'

rule VMC_DMC_GWFN:
    input:      '{molecule}/{method}/{basis}/gwfn.data',
                '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_DMC_DIRS:
    input:      '{molecule}/{method}/{basis}/gwfn.data'
    output:     '{molecule}/{method}/{basis}/VMC_DMC/{jastrow_opt_method}/casl/{jastrow_rank}/tmax_2_{nconfig}_1/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_OPT_ENERGY_RUN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/input',
                '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/gwfn.data',
                '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/parameters.casl'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/out'
    shell:      'cd {wildcards.path}/VMC_OPT/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}/1000000_9 && runqmc'

rule VMC_OPT_ENERGY_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc_opt_energy.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule))

rule VMC_OPT_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/out'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/parameters.casl'
    shell:      'ln -s ../10000/parameters.9.casl {output}'

rule VMC_OPT_ENERGY_GWFN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/.keep'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_OPT_ENERGY_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/1000000_9/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_OPT_RUN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/input',
                '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/gwfn.data',
                '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/parameters.casl',
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/out'
    shell:      'cd {wildcards.path}/VMC_OPT/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}/10000 && runqmc'

rule VMC_OPT_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc_opt.tmpl').read().format(neu=neu, ned=ned, nconfig=10000, method=wildcards.jastrow_opt_method, cycles=9, molecule=wildcards.molecule))

rule VMC_OPT_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/parameters.casl'
    run:
        for file_name in output:
            jastrow = wildcards.jastrow_rank.split('_')
            with open(file_name, 'w') as f:
                f.write(open('../casl.tmpl').read().format(term_2_0=jastrow[0], term_1_1=jastrow[1], term_2_1_1=jastrow[2][0], term_2_1_2=jastrow[2][1]))
            # workaround in multireference case
            source_path = os.path.join(wildcards.molecule, wildcards.method, wildcards.basis, 'correlation.data')
            target_path = os.path.join(os.path.dirname(file_name), 'correlation.data')
            shell('[[ -e {source_path} ]] && ln -s ../../../../../correlation.data {target_path}')


rule VMC_OPT_GWFN:
    input:      '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/.keep'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_OPT_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT/{jastrow_opt_method}/casl/{jastrow_rank}/10000/.keep'
    shell:      'touch {output}'

####################################################################################################################


rule VMC_RUN:
    input:      '{path}/VMC/10000000/input'
    output:     '{path}/VMC/10000000/out'
    shell:      'cd {wildcards.path}/VMC/10000000 && runqmc'

rule VMC_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC/10000000/gwfn.data'
    output:     '{molecule}/{method}/{basis}/VMC/10000000/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule))
            # workaround in multireference case
            source_path = os.path.join(wildcards.molecule, wildcards.method, wildcards.basis, 'correlation.data')
            target_path = os.path.join(os.path.dirname(file_name), 'correlation.data')
            shell('[[ -e {source_path} ]] && ln -s ../../correlation.data {target_path}')

rule VMC_GWFN:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC/10000000/gwfn.data'
    shell:      'ln -s ../../gwfn.data {output}'

rule VMC_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC/10000000/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_DMC_BF_RUN:
    input:      '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/input',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/correlation.data',
                '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/parameters.casl',
    output:     '{path}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/out'
    shell:      'cd {wildcards.path}/VMC_DMC_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/tmax_2_1024_1 && runqmc'

rule VMC_DMC_BF_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep',
                '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/out',
                '{molecule}/{method}/{basis}/VMC/10000000/out'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/input'
    params:
        dt_relative_step = 2.0,
        stderr = 0.001,
        # roughly
        # nstep ~ 1.5 * (E(VMC) - E(DMC))/(nconfig*dtdmc*stderr*stderr)
        # (E(VMC) - E(DMC)) ~ 6 * (E(HF) - E(VMC))
        magic_const = 1.5/6.0
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            gwfn_file = os.path.join(wildcards.molecule, wildcards.method, wildcards.basis, 'gwfn.data')
            hf, _ = vmc_energy(wildcards.molecule, wildcards.method, wildcards.basis)
            vmc, _ = vmc_opt_energy(wildcards.molecule, wildcards.method, wildcards.basis, 'VMC_OPT_BF', wildcards.jastrow_rank + '__' + wildcards.backflow_rank)
            dtdmc = 1.0/(get_max_Z(gwfn_file)**2 * 3.0 * params.dt_relative_step)
            nstep = params.magic_const*(hf - vmc)/(int(wildcards.nconfig)*dtdmc*params.stderr*params.stderr)
            nstep=max(50000, int(round(nstep, -3)))
            with open(file_name, 'w') as f:
                f.write(open('../vmc_dmc_bf.tmpl').read().format(
                    neu=neu, ned=ned, nconfig=wildcards.nconfig, dtdmc=dtdmc, molecule=wildcards.molecule, nstep=nstep, nblock=nstep // 1000
                ))

rule VMC_DMC_BF_DATA_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out',
                '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/correlation.data'
    shell:      'ln -s ../../../../../VMC_OPT_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/10000/correlation.out.9 {output}'

rule VMC_DMC_BF_CASL_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out',
                '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/parameters.casl'
    shell:      'ln -s ../../../../../VMC_OPT_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/10000/parameters.9.casl {output}'

rule VMC_DMC_BF_GWFN:
    input:      '{molecule}/{method}/{basis}/gwfn.data',
                '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_DMC_BF_DIRS:
    input:      '{molecule}/{method}/{basis}/gwfn.data'
    output:     '{molecule}/{method}/{basis}/VMC_DMC_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/tmax_2_{nconfig}_1/.keep'
    shell:      'touch {output}'

####################################################################################################################


rule VMC_OPT_BF_ENERGY_RUN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/input',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/gwfn.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/correlation.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/parameters.casl'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/out'
    shell:      'cd {wildcards.path}/VMC_OPT_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/1000000_9 && runqmc'

rule VMC_OPT_BF_ENERGY_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc_opt_energy_bf.tmpl').read().format(neu=neu, ned=ned, molecule=wildcards.molecule))

rule VMC_OPT_BF_DATA_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/correlation.data'
    shell:      'ln -s ../10000/correlation.out.9 {output}'

rule VMC_OPT_BF_CASL_ENERGY_JASTROW:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/parameters.casl'
    shell:      'ln -s ../10000/parameters.9.casl {output}'

rule VMC_OPT_BF_ENERGY_GWFN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/.keep'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_OPT_BF_ENERGY_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/1000000_9/.keep'
    shell:      'touch {output}'

####################################################################################################################

rule VMC_OPT_BF_CASL_RUN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/input',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/gwfn.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/correlation.data',
                '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/parameters.casl',
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/out'
    shell:      'cd {wildcards.path}/VMC_OPT_BF/{wildcards.jastrow_opt_method}/casl/{wildcards.jastrow_rank}__{wildcards.backflow_rank}/10000 && runqmc'

rule VMC_OPT_BF_INPUT:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/input'
    run:
        for file_name in output:
            neu, ned = get_up_down(wildcards.molecule, wildcards.method, wildcards.basis)
            with open(file_name, 'w') as f:
                f.write(open('../vmc_opt_bf.tmpl').read().format(
                    neu=neu, ned=ned, nconfig=10000,
                    method=wildcards.jastrow_opt_method, cycles=9, molecule=wildcards.molecule))

rule VMC_OPT_BF_DATA_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/correlation.data'
    run:
        for file_name in output:
            jastrow = wildcards.jastrow_rank.split('_')
            with open(file_name, 'w') as f:
                backflow = wildcards.backflow_rank.split('_')
                ae_cutoffs = get_ae_cutoffs(wildcards.molecule)
                number, labels = get_atom_labels(wildcards.molecule)
                template = '../backflow.tmpl' if backflow[2] == '00' else '../backflow.all'
                f.write(open(template).read().format(
                    eta_term=backflow[0],
                    mu_number_of_atoms=number, mu_atom_labels=labels, mu_term=backflow[1],
                    phi_number_of_atoms=number, phi_atom_labels=labels, phi_term_eN=backflow[2][0], phi_term_ee=backflow[2][1],
                    ae_cutoffs=ae_cutoffs))
            source_path = os.path.join(wildcards.molecule, wildcards.method, wildcards.basis, 'correlation.data')
            shell('[[ -e {source_path} ]] && cat {source_path} >> {file_name}')


rule VMC_OPT_BF_CASL_JASTROW:
    input:      '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    output:     '{molecule}/{method}/{basis}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/parameters.casl'
    run:
        for file_name in output:
            jastrow = wildcards.jastrow_rank.split('_')
            with open(file_name, 'w') as f:
                f.write(open('../casl.tmpl').read().format(term_2_0=jastrow[0], term_1_1=jastrow[1], term_2_1_1=jastrow[2][0], term_2_1_2=jastrow[2][1]))


rule VMC_OPT_BF_GWFN:
    input:      '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/gwfn.data'
    shell:      'ln -s ../../../../../gwfn.data {output}'

rule VMC_OPT_BF_DIRS:
    input:      '{path}/gwfn.data'
    output:     '{path}/VMC_OPT_BF/{jastrow_opt_method}/casl/{jastrow_rank}__{backflow_rank}/10000/.keep'
    shell:      'touch {output}'
