import argparse
import os

### Script to detect errors with a Valet summary file and an errors.bed file generated form test_data_generator_v2.py

### ARGS:
##      -g -> Ground truth error file from test_data_generator_v2.py
##      -v -> Valet summary.bed file
##      -l -> Label to give the current run (i.e. seed number used to gen the data)
##      -f -> Output filename; Can be provided to output the results in a tabular format

def read_bed(bedfile):
    errors = []

    with open(bedfile) as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.rstrip().split("\t")

            errors.append((cols[0],int(cols[1]),int(cols[2]),cols[3]))

    return errors


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-g", "--ground_truth",
        dest="ground_truth"
    )

    parser.add_argument(
        "-v", "--valet_bed",
        dest="valet_bed"
    )

    parser.add_argument(
        "-l", "--label",
        dest="label"
    )

    parser.add_argument(
        "-f", "--file",
        dest="filename", default=""
    )

    args = parser.parse_args()
    real_errors = read_bed(args.ground_truth)
    valet_errors = read_bed(args.valet_bed)

    valid_preds = set()
    detected_errors = set()
    d_i = set()
    t_i = 0
    d_lc = set()
    t_lc = -0
    d_du = set()
    t_du = 0
    d_de = set()
    t_de = 0


    # Count total number of each real error type
    for error in real_errors:
        match error[3]:

            case "inversion":
                t_i += 1

            case "low_coverage":
                t_lc += 1

            case "duplication":
                t_du += 1

            case "deletion":
                t_de += 1

    
    # Search for each predicted error in the real errors
    for i,pred_error in enumerate(valet_errors):
        for j,real_error in enumerate(real_errors):
            center = (real_error[1] + real_error[2]) / 2
            
            s = pred_error[1]
            e = pred_error[2]
            pred_center = (e+s) / 2

            # Check if the 'center' of the real error is within 1kb of the start/end of the predicted error (from valet)
            if pred_center - 1000 <= center and pred_center + 1000 >= center:

                # If the error is properly detected, add it to the valid preds and detected error sets
                valid_preds.add(i)
                detected_errors.add(j)

                # Also add it to each error-type set
                match real_error[3]:

                    case "inversion":
                        d_i.add(j)

                    case "low_coverage":
                        d_lc.add(j)

                    case "duplication":
                        d_du.add(j)

                    case "deletion":
                        d_de.add(j)

    # Print results
    print(str(len(valid_preds)) + " valid preds out of " + str(len(valet_errors)) + " total errors.")
    print(str(len(detected_errors)) + " detected errors out of " + str(len(real_errors)) + " real simulated errors.")
    print(str(len(d_i)) + " / " + str(t_i) + " inversions")
    print(str(len(d_lc)) + " / " + str(t_lc) + " low-coverage")
    print(str(len(d_du)) + " / " + str(t_du) + " duplications")
    print(str(len(d_de)) + " / " + str(t_de) + " deletions")

    safe = lambda a, b: (a / b * 100) if b else 0.0

    # Log results as a tsv row if a log file is provided
    if args.filename != "":

        # If the file does not already exist - add it and a header
        if not os.path.exists(args.filename):
            with open(args.filename, "w") as f:
                f.write(
                    "label\tvalid_preds\ttotal_preds\tvalid_pct\t"
                    "detected_errors\ttotal_errors\tdetected_pct\t"
                    "inversions\ttotal_inversions\tinv_pct\t"
                    "lowcov\ttotal_lowcov\tlowcov_pct\t"
                    "duplications\ttotal_duplications\tdup_pct\t"
                    "deletions\ttotal_deletions\tdel_pct\n"
                )

        with open(args.filename, "a") as f:
            f.write(
                f"{args.label}\t"
                f"{len(valid_preds)}\t{len(valet_errors)}\t{safe(len(valid_preds), len(valet_errors)):.3f}\t"
                f"{len(detected_errors)}\t{len(real_errors)}\t{safe(len(detected_errors), len(real_errors)):.3f}\t"
                f"{len(d_i)}\t{t_i}\t{safe(len(d_i), t_i):.3f}\t"
                f"{len(d_lc)}\t{t_lc}\t{safe(len(d_lc), t_lc):.3f}\t"
                f"{len(d_du)}\t{t_du}\t{safe(len(d_du), t_du):.3f}\t"
                f"{len(d_de)}\t{t_de}\t{safe(len(d_de), t_de):.3f}\n"
            )




if __name__=="__main__":
    main()