import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Transfer leafcutter introns to bed file")
    # Set up arguments
    parser.add_argument("--leafcutter_intron_file", type=str, required=True,
                                                 help="Leafcutter intron file")
    parser.add_argument("--bed_file", type=str, required=True, help="Output bed file")

    # Load parameters
    args = parser.parse_args()
    leafcutter_intron_file = args.leafcutter_intron_file
    bed_file = args.bed_file

    # Set BED-file header
    bed_file_text = 'track name="Leafcutter Introns" description="Leafcutter Introns" visibility=2 itemRgb="On"\n'

    with open(leafcutter_intron_file) as f:
        leafcutter_intron_lines = f.readlines()     # returns list of strings

    line_count = 0
    for line in leafcutter_intron_lines:
        line_count += 1

        if line_count > 1:      # Skip first line -> headers
            intron_position = line.split(" ")[0]

            intron_array = intron_position.split(":")
            if len(intron_array) >= 3:  # Skips last empty line...
                chromosome = intron_array[0]
                start = intron_array[1]
                end = intron_array[2]

                bed_file_text += "\t".join([chromosome, start, end, "intron" + str(line_count),  "0", "+", start, end, "255,0,0\n"])

    with open(bed_file, "w") as f:
        f.write(bed_file_text)