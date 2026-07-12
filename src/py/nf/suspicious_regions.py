""" Find regions that overlap and contain different mis-assembly signatures.

Args:
    options: Command line arguments.
    bed_filename: BED file to find suspicious regions.
    output_dir: Output directory of current assembly.

Returns:
    True if the bed_filename parameter contains any entries.
"""


def find_suspicious_regions(options, bed_filename, output_dir):
    

    # bedtools sort -i - | bedtools merge -d 1000 -o distinct -c 4 -d 1000
    with open(bed_filename, 'r') as bed_file:
        # Truncate all BED entries down to 4 fields.
        regions = []
        for line in bed_file:
            regions.append(line.strip().split()[0:4])

        if len(regions) <= 0:
            return False
        regions.sort(key = lambda region: (region[0], int(region[1]), int(region[2])))

        intermediate_bed_file = open(output_dir + "/tmp.bed", 'w')
        intermediate_bed_file.write('\n'.join(['\t'.join(region) for region in regions]))
        intermediate_bed_file.close()

        # Merge the BEDfile.
        merged_bed_file = open(output_dir + '/tmp.merged.bed','w')
        call_arr = ['bedtools', 'merge',\
                '-i', intermediate_bed_file.name,\
                '-d', str(options.suspicious_flank_size), '-o', 'distinct', '-c', '4']
        run(call_arr, stdout=merged_bed_file)
        merged_bed_file.close()

        # Only save flagged regions that meet the minimum number of signatures.
        with open(merged_bed_file.name, 'r') as merged_file,\
                open(output_dir + "/suspicious.bed", 'w') as suspicious_file:
            for entry in merged_file:
                if len(entry.split()[3].split(',')) >= options.min_suspicious_regions:
                    suspicious_file.write(entry)

    return True