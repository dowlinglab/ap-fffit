import numpy as np
import signac
import pandas as pd

from .ap import APConstants


def save_signac_results(project, param_names, csv_name):
    """Save the signac results to a CSV file.

    Parameters
    ----------
    project : signac.Project
        signac project to load
    param_names : tuple
        parameters that are varied
    csv_name : string
        name of csv file to save results
    """
    AP = APConstants()
    if type(param_names) not in (list, tuple):
        raise TypeError("param_names must be a list or tuple")
    job_groupby = tuple(param_names)

    # Store data here before converting to dataframe
    data = []
    umask = np.logical_and(AP.unit_cell_elements != "H", AP.unit_cell_elements != "O")
    umask_Cl = AP.unit_cell_elements == "Cl"
    umask_O = AP.unit_cell_elements == "O"
    umask_N = AP.unit_cell_elements == "N"
    umask_H = AP.unit_cell_elements == "H"


    # Loop over all jobs in project and group by parameter sets
    for params, job_group in project.groupby(job_groupby):

        for job in job_group:
            # Extract the parameters into a dict
            new_row = {name: param for (name, param) in zip(job_groupby, params)}

            # Extract the temperature for each job.
            # Assumes temperature increments >= 1 K
            temperature = round(job.sp.T)
            new_row["temperature"] = temperature

            # Some jobs in the early rounds fail. Save everything as np.nan
            try:
                fdata = []
                with open(job.fn("Sim_UC.xyz")) as f:
                    for line in f:
                        fdata.append(line.strip().split())
                ucell = np.array(fdata[2:])
                ucell = np.array(ucell[:,1:], dtype=np.float64)
            except FileNotFoundError:
                print(f"Job failed: {job}")
                new_row["uc_mean_distance"] = np.nan
                new_row["uc_mean_distance_Cl"] = np.nan
                new_row["uc_mean_distance_O"] = np.nan
                new_row["uc_mean_distance_N"] = np.nan
                new_row["uc_mean_distance_H"] = np.nan
                new_row["lattice_ape"] = np.nan
                new_row["abs(HB3-HB4)"] = np.nan
                new_row["log(abs(HB3-HB4))"] = np.nan
                new_row["abs(HA3-HA4)"] = np.nan
                new_row["log(abs(HA3-HA4))"] = np.nan
                continue

            if temperature > 10:

                scale = np.array([
                        AP.expt_lattice_a[temperature] / AP.expt_lattice_a[10],
                        AP.expt_lattice_b[temperature] / AP.expt_lattice_b[10],
                        AP.expt_lattice_c[temperature] / AP.expt_lattice_c[10],
                    ]
                )
                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz * scale, axis=1) * umask
                uc_mean_distance = np.sum(masked_result) / np.count_nonzero(masked_result)

                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz * scale, axis=1) * umask_Cl
                uc_mean_distance_Cl = np.sum(masked_result) / np.count_nonzero(masked_result)

                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz * scale, axis=1) * umask_O
                uc_mean_distance_O = np.sum(masked_result) / np.count_nonzero(masked_result)

                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz * scale, axis=1) * umask_N
                uc_mean_distance_N = np.sum(masked_result) / np.count_nonzero(masked_result)

                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz * scale, axis=1) * umask_H
                uc_mean_distance_H = np.sum(masked_result) / np.count_nonzero(masked_result)


            else:

                uc_mean_distance = np.mean(np.linalg.norm(ucell - AP.unit_cell_xyz, axis=1))

                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz, axis=1) * umask_Cl
                uc_mean_distance_Cl = np.sum(masked_result) / np.count_nonzero(masked_result)

                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz, axis=1) * umask_O
                uc_mean_distance_O = np.sum(masked_result) / np.count_nonzero(masked_result)

                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz, axis=1) * umask_N
                uc_mean_distance_N = np.sum(masked_result) / np.count_nonzero(masked_result)

                masked_result = np.linalg.norm(ucell - AP.unit_cell_xyz, axis=1) * umask_H
                uc_mean_distance_H = np.sum(masked_result) / np.count_nonzero(masked_result)

                # H-bonds/angles only well-defined at 10 K
                HBond1 = job.doc.HBond1
                HBond2 = job.doc.HBond2
                HBond3 = job.doc.HBond3
                HBond4 = job.doc.HBond4

                HAngle1 = job.doc.HAngle1
                HAngle2 = job.doc.HAngle2
                HAngle3 = job.doc.HAngle3
                HAngle4 = job.doc.HAngle4

                new_row["abs(HB3-HB4)"] = np.abs(HBond3-HBond4)
                new_row["log(abs(HB3-HB4))"] = np.log(np.abs(HBond3-HBond4))
                new_row["abs(HA3-HA4)"] = np.abs(HAngle3-HAngle4)
                new_row["log(abs(HA3-HA4))"] = np.log(np.abs(HAngle3-HAngle4))


            new_row["uc_mean_distance"] = uc_mean_distance
            new_row["uc_mean_distance_Cl"] = uc_mean_distance_Cl
            new_row["uc_mean_distance_O"] = uc_mean_distance_O
            new_row["uc_mean_distance_N"] = uc_mean_distance_N
            new_row["uc_mean_distance_H"] = uc_mean_distance_H

            new_row["lattice_ape"] = 100 * np.mean(
                [ np.abs(job.doc.a - AP.expt_lattice_a[temperature]) / AP.expt_lattice_a[temperature],
                  np.abs(job.doc.b - AP.expt_lattice_b[temperature]) / AP.expt_lattice_b[temperature],
                  np.abs(job.doc.c - AP.expt_lattice_c[temperature]) / AP.expt_lattice_c[temperature]
                ]
            )

            new_row["job_id"] = job.id
            data.append(new_row)


    # Save to csv file for record-keeping
    df = pd.DataFrame(data)
    df.to_csv(csv_name)

