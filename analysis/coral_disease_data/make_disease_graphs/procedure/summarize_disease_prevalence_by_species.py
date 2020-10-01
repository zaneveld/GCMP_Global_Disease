import pandas as pd
import argparse

# Parser for arguments

parser = argparse.ArgumentParser(description="Creates a .csv table with species, clade, total % with disease")

parser.add_argument('-i', type=str,
                          dest="input_data",
                          help="Enter path to a comma seperated values file (.csv)"
                    )

parser.add_argument('-x', '--index_column_name', type=str,
                                     help='The name of a column in the csv files. Counts will be summed for each unique entry in this column. Example: If you pick species,'
                    )

parser.add_argument('-c', '--count_columns', type=str,
                                       nargs='*',
                                       help="Enter columns you would like to add to the sum table"
                    )
parser.add_argument('-o', '--output_file', type=str, help="Enter a path to write an output .csv; default is %default",default='../outputs/species_counts.csv')


opts = parser.parse_args()

# read csv
data_in = pd.read_csv('%s' % opts.input_data)


columns_to_keep = [opts.index_column_name] + opts.count_columns

data_in = data_in.loc[:,columns_to_keep]

#pull names from species list given by user and format into a list
species_list = data_in[opts.index_column_name].unique().tolist()

# assign sum of all values to a list for the final dataframe
sum_row = data_in.sum(axis=1).tolist()

# create dictionary
d = {'Species':species_list,'Total':sum_row}

#creat a new dataframe using the dictionary created previously
summary_df = pd.DataFrame(d)

print(summary_df)

#create csv file 
summary_df.to_csv(opts.output_file, index=False)

