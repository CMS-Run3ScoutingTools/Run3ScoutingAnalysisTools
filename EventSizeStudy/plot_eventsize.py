####################################################################################
# A script to plot pie charts of the average event size as result of the command:  #
# > edmEventSize -v rootfilename.root > Ref.txt                                    #
# Example for running (Ref.txt is provided as example file):                       #
# > python3 plot_eventsize.py Ref.txt --output Ref.pdf                             #
# Contact person: Sanu Varghese, Elisa Fontanesi                                   #
####################################################################################

import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep 
import argparse
import sys
import warnings
warnings.filterwarnings("ignore", message="The value of the smallest subnormal for <class 'numpy.float64'> type is zero.")

# Set the style to CMS
# --------------------
hep.style.use('CMS')

def read_data(filename, column_index):
    branches = []
    sizes = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Skip the first 2 lines which contains a space and general file info and the second line which is the header
        for line in lines[3:]:
            parts = line.strip().split()
            branch_name = ' '.join(parts[:-2]).replace('__MYHLT', '')
            try:
                size = float(parts[-column_index])  # Make sure this accesses the correct index
            except IndexError:
                # print(f"IndexError for line: {line.strip()}")
                continue  # Skip lines that do not have enough parts
            branches.append(branch_name)
            sizes.append(size)
    return branches, sizes

def main():
    parser = argparse.ArgumentParser(description='Plot the breakdown of the edm event size')
    parser.add_argument('filename', type=str, help='Filename to process')
    parser.add_argument('--uncompressed', action='store_true', help='Use uncompressed sizes instead of compressed')
    parser.add_argument('--output', type=str, help='Output file name for the plot', default='output.pdf')

    args = parser.parse_args()

    column_index = 2 if args.uncompressed else 1  # 1 for compressed, 2 for uncompressed

    branches, sizes = read_data(args.filename, column_index)

    total_size_bytes = sum(sizes)
    percentages = [(size / total_size_bytes) * 100 for size in sizes]

    # NOTE: define threshold below which everything gets clubbed to others
    # --------------------------------------------------------------------
    threshold_percentage = 3
    new_branches = []
    new_sizes = []
    other_size = 0

    for branch, size, percentage in zip(branches, sizes, percentages):
        if percentage < threshold_percentage:
            other_size += size
        else:
            new_branches.append(branch)
            new_sizes.append(size)

    if other_size > 0:
        new_branches.append('Others')
        new_sizes.append(other_size)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', 
              '#7f7f7f', '#bcbd22', '#17becf'] * 2
    colors = colors[:len(new_branches)]

    def make_autopct(sizes):
        def my_autopct(pct):
            total = sum(sizes)
            val = int(round(pct * total / 100.0))
            return '{:.1f}%\n({} B)'.format(pct, val)
        return my_autopct

    fig, ax = plt.subplots(figsize=(10, 6), dpi=650)
    wedges, texts, autotexts = ax.pie(new_sizes, labels=new_branches, colors=colors, startangle=90,
                                      counterclock=False, radius=1.0, autopct=make_autopct(new_sizes), pctdistance=0.85,
                                      textprops={'fontsize': 8})

    centre_circle = plt.Circle((0, 0), 0.50, fc='white')
    fig.gca().add_artist(centre_circle)

    # Unit conversion: convention that 1 kB = 1024 bytes
    # --------------------------------------------------
    total_size_kB = total_size_bytes / 1024

    plt.text(0, 0, 'Total Size:\n{:.1f} kB'.format(total_size_kB), ha='center', va='center', fontsize=12)

    hep.cms.label(ax=ax, data=True, loc=0, fontsize=10 ,com =13.6)
    #plt.title(args.filename.split('.')[0])

    # Save the figure as specified by the output argument
    # ---------------------------------------------------
    fig.savefig(args.output, format=args.output.split('.')[-1], dpi=650)
    #plt.show()

if __name__ == '__main__':
    main()
