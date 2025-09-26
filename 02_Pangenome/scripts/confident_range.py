import sys
from scipy.stats import nbinom

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

model_file = sys.argv[1]

for line in openfile(model_file):
    if len(line.split()) > 0:
        if line.split()[0] == "kmercov":
            kmercov = float(line.split()[1])
        if line.split()[0] == "bias":
            bias = max(0.001, float(line.split()[1]))

mu_1 = kmercov
size_1 = mu_1 / bias
p_1 = size_1 / (size_1 + mu_1)
lower_bound = nbinom.ppf(0.0005, n=size_1, p=p_1)

mu_2 = 2 * kmercov
size_2 = mu_2 / bias
p_2 = size_2 / (size_2 + mu_2)
upper_bound = nbinom.ppf(0.9995, n=size_2, p=p_2)

print(str(lower_bound) + "\t" + str(upper_bound))
