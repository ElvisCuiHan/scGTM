from scGTM_Valley_Only import *

import argparse

parser = argparse.ArgumentParser(description='Sing-cell Kinetics Generalized Additive Model')

parser.add_argument('--gene.start', type=int, default=2, metavar='GENESTART',
                    help="index of gene to start (default: {:d})".format(2))
parser.add_argument('--gene.end', type=int, default=3, metavar='GENEEND',
                    help="index of gene to end (default: {:d})".format(3))
parser.add_argument('--model.iter', type=int, default=100, metavar='ITERNUM',
                    help="Number of iteration in PSO (default: 100)")
parser.add_argument('--data.dir', type=str, default=None, metavar='DATADIR',
                    help="directory to data (default: None)")
parser.add_argument('--model.marginal', type=str, default="ZIP", metavar='MARGINAL',
                    help="Marginal distribution (default: ZIP)")
parser.add_argument('--model.save_dir', type=str, default="Results/ZIP/", metavar='SAVEDIR',
                    help="directory to save results (default: Results/ZIP)")
args = vars(parser.parse_args())

if __name__ == '__main__':
    parallel(args)
