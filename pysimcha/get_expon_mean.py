import pandas as pd
import argparse
import os
from scipy.stats import expon


def get_expon_mean(a, b, c, n_cutoff):
	path_to_file = os.path.abspath(__file__)
	this_path = os.path.dirname(path_to_file)
	data_path = os.path.join(this_path, "../out/pcawg_filtered_500/")
	data = pd.read_csv(os.path.join(data_path, "clones.tsv"), sep="\t")
	params = {"stress": float(a), "tsg": float(b), "og": float(b), "ess": float(c)}
	myfitfunc = lambda x:  (params[x] * data[x]).to_numpy()
	mySum = sum(map(myfitfunc, params.keys()))

	data['delta_fitness'] = mySum
	
	dFit_dist = data["delta_fitness"].to_numpy()
	# remove the negative values
	dFit_dist = dFit_dist[dFit_dist > 0]
    # Aim for no fitness advantage when there aren't enough samples to create an exponential fit
	if len(dFit_dist) < n_cutoff:
		return 0, 0
		
	# fit the exponential distribution
	loc, scale = expon.fit(dFit_dist)
	return loc, scale


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="Get the scale of the exponential distribution of the filtered PCAWG dataset, given a set of input alpha, beta, gamma values")
	parser.add_argument("-a", default=1, help="Cellular stress weight")
	parser.add_argument("-b", default=1, help="TSG/OG weight")
	parser.add_argument("-c", default=1, help="Essentiality weight")
	parser.add_argument("-n", default=100, help="Minimum number of samples needed to fit the exponential distribution")
	args = parser.parse_args()

	loc, scale =  get_expon_mean(args.a, args.b, args.c, args.n)
	print(scale, end='')
