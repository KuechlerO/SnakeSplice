# This script can be used to collect the Ensembl gene ID -> gene name mappings

import pybiomart


def get_ensembl_mappings(output_file):
    """
    Gets the Ensembl mappings: gene ID -> gene name
    :return:
    """
    # Set up connection to Ensembl Biomart server
    dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

    # List the types of data we want
    my_attributes = ['ensembl_gene_id', 'external_gene_name']
    my_filters = {}     # No filters -> get all genes

    # Query BioMart
    output_df = dataset.query(attributes=my_attributes, filters=my_filters)
    # Rename the columns
    output_df.rename(columns={'Gene stable ID': 'ensembl_gene_id', 'Gene name': 'gene_name'}, inplace=True)

    output_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    # Output file
    output_file = snakemake.output[0]
    # Get the mappings
    get_ensembl_mappings(output_file)
