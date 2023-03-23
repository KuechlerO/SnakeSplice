import os


def create_leafcutter_group_file(output_dir, condition, control_sample_ids, patient_sample_ids):
    """
    Creates a group file for leafcutter analysis
    :param output_dir:              Output directory
    :param condition:               Current condition's name
    :param control_sample_ids:      List of control sample ids
    :param patient_sample_ids:
    :return:
    """
    # Create group file
    group_file = os.path.join(output_dir, f"{condition}_group_file.txt")

    group_file_text = ""
    for sample_id in control_sample_ids:
        group_file_text += f"{sample_id}\tcontrol\n"
    for sample_id in patient_sample_ids:
        group_file_text += f"{sample_id}\tpatient\n"

    with open(group_file, "w") as f:
        f.write(group_file_text)


if __name__ == '__main__':
    # Get snakemake variables
    output_dir = snakemake.params.output_dir
    control_samples_ids = snakemake.params.control_samples["sample_name"].tolist()
    condition_samples_array = snakemake.params.condition_samples_array

    for condition_samples in condition_samples_array:
        condition_samples_ids = condition_samples["sample_name"].tolist()
        current_condition = condition_samples["condition"][0]
        print("condition: ", current_condition)
        print("control samples: ", control_samples_ids)
        print("condition samples: ", condition_samples)
        print("condition samples IDs: ", condition_samples_ids)
        create_leafcutter_group_file(output_dir, current_condition, control_samples_ids, condition_samples_ids)
