import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def vanilla_plot(gene_index=None, data_dir=None, save_dir=None, plot_args=None):
    print("Loading data......")

    ## LOAD DATA
    data = pd.read_csv(data_dir)
    print("Loading finished!")

    ## TAKE NEEDED DATA
    t = data.iloc[:, 1]
    y1 = np.floor(data.iloc[:, gene_index])
    gene_name = data.columns[gene_index]
    log_data = np.log(y1 + 1)

    plt.figure(figsize=(10, 8))
    plt.scatter(t, log_data, s=10, c=log_data, cmap=plt.get_cmap(plot_args['cmap']))
    plt.xlabel("Pseudotime", fontsize=18)
    plt.ylabel("Expression log(FPKM+1)", fontsize=18)
    #plt.show()
    plt.savefig(save_dir + str(gene_index) + ".png")