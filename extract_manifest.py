import os

# Define the directory containing the data files
data_dir = '/Users/calebchin/Documents/met_wasserstein/cs4775_project_data'

# Define the output TSV file
output_tsv = '/Users/calebchin/Documents/met_wasserstein/hmp2_16s_manifest.tsv'

# Open the output TSV file for writing
with open(output_tsv, 'w') as tsv_file:
    # Write the header
    tsv_file.write("sample-id\tabsolute-filepath\n")
    
    # Iterate through the files in the data directory
    for filename in os.listdir(data_dir):
        # Construct the absolute file path
        absolute_filepath = os.path.join(data_dir, filename)
        
        # Extract the sample ID from the filename (assuming the sample ID is the part before the first underscore)
        sample_id = filename.split('.')[0]
        
        # Write the sample ID and absolute file path to the TSV file
        tsv_file.write(f"{sample_id}\t{absolute_filepath}\n")