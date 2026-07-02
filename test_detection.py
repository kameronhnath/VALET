import argparse

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

    args = parser.parse_args()
    real_errors = read_bed(args.ground_truth)
    valet_errors = read_bed(args.valet_bed)

    valid_preds = set()
    detected_errors = set()
    for i,pred_error in enumerate(valet_errors):
        for j,real_error in enumerate(real_errors):
            center = (real_error[1] + real_error[2]) / 2
            
            s = pred_error[1]
            e = pred_error[2]

            if s - 1000 <= center and e + 1000 >= center:
                valid_preds.add(i)
                detected_errors.add(j)

    print(str(len(valid_preds)) + " valid preds out of " + str(len(valet_errors)) + " total errors.")
    print(str(len(detected_errors)) + " detected errors out of " + str(len(real_errors)) + " real simulated errors.")



if __name__=="__main__":
    main()