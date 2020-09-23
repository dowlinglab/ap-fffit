import flow
from flow import FlowProject, directives
import templates.ndcrc
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)


class Project(FlowProject):
    pass


@Project.operation
@Project.post.isfile("in.data")
def create_data(job):
    """Copy the lammps data file to the job directory"""
    import shutil
    shutil.copy("../data/APScaled-6-9-7.data", job.fn("in.data"))


@Project.operation
@Project.post.isfile("in.lammps")
def create_input(job):
    """Create the lammps input file"""

    content = _generate_input(job)

    with open(job.fn("in.lammps"), "w") as f:
        f.write(content)


@Project.label
def simulation_complete(job):
    if job.isfile("final.rst"):
        return True
    else:
        return False


@Project.operation
@Project.pre.after(create_data)
@Project.pre.after(create_input)
@Project.post(simulation_complete)
@flow.with_job
@flow.cmd
def simulate(job):
    """Run the lammps simulation!"""
    return "mpirun -np 2 lmp_intel_gmt -pk intel 0 omp 1 lrt yes -sf intel -in in.lammps"


@Project.operation
@Project.pre.after(simulate)
@Project.post(lambda job: "a" in job.doc)
@Project.post(lambda job: "b" in job.doc)
@Project.post(lambda job: "c" in job.doc)
@Project.post(lambda job: "a_unc" in job.doc)
@Project.post(lambda job: "b_unc" in job.doc)
@Project.post(lambda job: "c_unc" in job.doc)
def calculate_lattice_parameters(job):
    """Calculate the average lattice parameters"""

    import numpy as np
    import lammps_thermo
    from block_average import block_average

    n_unitcell_x = 6
    n_unitcell_y = 9
    n_unitcell_z = 7

    # Load the thermo data
    thermo = lammps_thermo.load(job.fn("log.prod"))
    lx = thermo.prop("lx") / n_unitcell_x
    ly = thermo.prop("ly") / n_unitcell_y
    lz = thermo.prop("lz") / n_unitcell_z

    # Save averages
    job.doc.a = np.mean(lx)
    job.doc.b = np.mean(ly)
    job.doc.c = np.mean(lz)

    # Compute an approx uncertainty with block averaging and save
    (means_est, vars_est, vars_err) = block_average(lx)
    job.doc.a_unc = np.max(np.sqrt(vars_est))
    (means_est, vars_est, vars_err) = block_average(ly)
    job.doc.b_unc = np.max(np.sqrt(vars_est))
    (means_est, vars_est, vars_err) = block_average(lz)
    job.doc.c_unc = np.max(np.sqrt(vars_est))


#####################################################################
################# HELPER FUNCTIONS BEYOND THIS POINT ################
#####################################################################


def _generate_input(job):

    contents = f"""
# Ammonium Perchlorate, 6-9-7, 1512 ClO4-, 1512 NH4+
# February 12th 2018
# Garrett Tow

# Initialization
units real
dimension 3
boundary p p p
atom_style full

read_data in.data

timestep 1.0

thermo 50
thermo_style custom step vol temp press pe ke etotal evdwl ecoul elong enthalpy density lx ly lz ebond eangle
thermo_modify line one

# Structural Definition
#
# Pair Coeffs
#
# 1  Cl
# 2  H
# 3  N
# 4  O

# Bond Coeffs
#
# 1  Cl-O
# 2  H-N

# Angle Coeffs
#
# 1  H-N-H
# 2  O-Cl-O

pair_style lj/cut/coul/long 15.0 15.0
kspace_style pppm 1.0e-5

pair_coeff 1 1 {job.sp.epsilon_Cl:10.4f} {job.sp.sigma_Cl:10.3f}
pair_coeff 2 2 {job.sp.epsilon_H:10.4f} {job.sp.sigma_H:10.3f}
pair_coeff 3 3 {job.sp.epsilon_N:10.4f} {job.sp.sigma_N:10.3f}
pair_coeff 4 4 {job.sp.epsilon_O:10.4f} {job.sp.sigma_O:10.3f}


# Bond Coeffs
#
# 1  Cl-O
# 2  H-N

bond_style harmonic

bond_coeff 1 426.42 1.4523
bond_coeff 2 413.55 1.0300

# Angle Coeffs
#
# 1  H-N-H
# 2  O-Cl-O

angle_style harmonic

angle_coeff 1 33.45 109.5
angle_coeff 2 107.60 109.5

neigh_modify one 8000

group perchlorate type 1 4
group ammonium type 2 3

velocity all create {job.sp.T} 293288 dist gaussian mom yes rot yes
fix 1 all npt temp {job.sp.T} {job.sp.T} 100.0 aniso {job.sp.P} {job.sp.P} 1000.0

run {job.sp.nsteps.eq}

dump 1 all custom 10000 out.lammpstrj id mol type xu yu zu
dump_modify 1 sort id
dump_modify 1 format line "%d %d %d %20.15g %20.15g %20.15g"

log log.prod

run {job.sp.nsteps.prod}
write_restart final.rst
"""

    return contents


if __name__ == "__main__":
    Project().main()
