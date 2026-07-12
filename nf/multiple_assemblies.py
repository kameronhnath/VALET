import argparse
import subprocess
import os

### Python script to run valet with multiple assemblies of the same reads

### ARGS:
##      -r1 and -r2         -> reads files
##      -p (--profile)      -> custom profile to use in the config file
##      -a (--assemblies)   -> list of assembly files separated by commas (i.e. -a a1.fasta,a2.fasta)
##      -n (--names)        -> list of names for each assembly (i.e. -n name1,name2), needs to be the same length of the list of assemblies
##      -o (--out_dir)      -> output directory for the results

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-r1",
        dest="reads1"
    )

    parser.add_argument(
        "-r2",
        dest="reads2"
    )

    parser.add_argument(
        "-p", "--profile",
        dest="profile", default=""
    )

    parser.add_argument(
        "-a", "--assemblies",
        dest="assemblies"
    )

    parser.add_argument(
        "-n", "--names",
        dest="names"
    )

    parser.add_argument(
        "-o", "--out_dir",
        dest="out_dir"
    )

    args = parser.parse_args()

    assemblies = args.assemblies.split(",")
    names = args.names.split(",")

    outdir = os.path.abspath(args.out_dir)

    summary_files = []

    for (assembly,name) in list(zip(assemblies,names)):
        nf_process = ["nextflow", "valet.nf", "--assembly", assembly, 
                        "--reads1", args.reads1, "--reads2", args.reads2,
                        "-output-dir", outdir + "/" + name]
        if args.profile != "":
            nf_process = nf_process + ["-profile", args.profile]
        subprocess.run(nf_process)
        summary_files.append(outdir + "/" + name + "/summary_table/summary.tsv")

    print(summary_files)
    R_process = ["Rscript", "../src/R/compare_assemblies.R", ','.join(summary_files), ','.join(names), outdir + '/comparison_plots']
    subprocess.run(R_process)


if __name__=="__main__":
    main()
