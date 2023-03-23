import os
import argparse


def transform_gtf_chromosome_names(input_gtf_file, output_gtf_file, log_file):
    """
    Transforms the chromosome names in the gtf file to the ones used in aligned BAM-files.
    I.e.: 1 -> chr1, 2 -> chr2, ..., X -> chrX, Y -> chrY, MT -> chrM

    Explanation:
    # Use double braces as escape for single braces in awk command
    # Also <"> are escaped...
    # Explanation:
    # 1. BEGIN: Sets the input and output field separator to tab
    # 2. $1=($1~/^\#/?$1:sprintf(\"chr%s\",$1)): Sets first column to "chr" + first column if first column does not start with "#"
    # 2.1 /.../: Regular expression
    # 2.2 $1: First column
    # 2.3 ~: Matches
    # 2.4 /^\#/: Starts with "#"
    # 2.5 ?$1: If true, then $1
    # 2.6 :$1:sprintf(\"chr%s\",$1): If false, then "chr" + $1
    # 3. 1 is default action in awk to print full record
    """
    command = "awk 'BEGIN {{FS=OFS=\"\t\"}} {{$1=($1~/^\#/?$1:sprintf(\"chr%s\",$1))}} 1' {input_gtf_file} " \
              "> {output_gtf_file} 2> {log_file}".format(
                input_gtf_file=input_gtf_file,
                output_gtf_file=output_gtf_file,
                log_file=log_file
    )
    os.system(command)


if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Transforms the chromosome names in the gtf file to the ones used in "
                                                 "aligned BAM-files (i.e. transforms 1 to chr1 for example...)")
    parser.add_argument("--input_gtf_file", type=str, required=True, help="The input gtf file.")
    parser.add_argument("--output_gtf_file", type=str, required=True, help="The output gtf file.")
    parser.add_argument("--log_file", type=str, required=True, help="The log file.")

    # Parse arguments
    args = parser.parse_args()
    input_gtf_file = args.input_gtf_file
    output_gtf_file = args.output_gtf_file
    log_file = args.log_file

    # Transform chromosome names
    transform_gtf_chromosome_names(input_gtf_file, output_gtf_file, log_file)