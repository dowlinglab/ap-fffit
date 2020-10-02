import signac
import sys

from fffit.signac import save_signac_results
from utils.ap import AP


def main():
    if len(sys.argv) != 2:
        print("Usage: python extract_lattice.py [iteration number]")
        exit(1)
    else:
        iternum = sys.argv[1]

    run_path = "/scratch365/rdefever/ap-fffit/ap-fffit/runs/"
    itername = "lattice-iter" + str(iternum)
    project_path = run_path + itername
    csv_name = "csv/" + itername + "-results.csv"

    property_names = ["a", "b", "c"]
    project = signac.get_project(project_path)
    save_signac_results(project, AP.param_names, property_names, csv_name)


if __name__ == "__main__":
    main()
