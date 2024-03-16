import os
import pandas as pd


def reformat_finemap_to_caviarbf(infile_name, outfile_name):
    df = pd.read_csv(infile_name, sep=' ')
    df[['rsid', 'beta', 'se']].to_csv(outfile_name, index=False, header=False, sep=' ')


data_dir = '/media/zhaoyang-new/workspace/westonelison_finemap_simulator/CSE_284_Finemapping/simulations/simulation_summary_stats/'
output_dir = '/example/CV2/'

directories = [d for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d))]

for dir_itr in directories:
    folder_path = os.path.join(data_dir, dir_itr)
    files_in_folder = os.listdir(folder_path)
    z_file = [f for f in files_in_folder if f.endswith('.z')][0]

    output_name = z_file.split('.finemap')[0] + '.caviarbf.z'
    reformat_finemap_to_caviarbf(os.path.join(folder_path, z_file), os.path.join(output_dir, output_name))
