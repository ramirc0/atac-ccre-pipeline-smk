#Jill E. Moore
#Moore Lab - UMass Chan
#ENCODE4 cCRE Pipeline
#December 2024

import numpy, sys, math

# init lists to collect values calculated in for loop

# values after for loop calculations:
# [[raw signal val], [log10(signal)]] - unless val missing/zero, then assigned "Zero" instead in for loop below
sig=[[],[]]
# peak name (anchor ID) for each row to preserve original order
masterPeak=[]
# log10 vals for non-zero entries only
calculate=[]

for line in open(sys.argv[1]):
    line=line.rstrip().split("\t")
    # index col 4 = mean bigWig signal, col 2 = # bases covered (double check: in case value is missing instead of 0)
    if float(line[4]) == 0 or float(line[2]) == 0:
        # if there is no signal, set "Zero" to assign -10 later
        sig[1].append("Zero")
        # set sig[0] (raw signal) to 0 so that the index stays aligned across lists
        # if sig[0] is shorter (no 0), you'd get an index error when combined with other entries
        sig[0].append((float(line[4])))
        masterPeak.append(line[0])
    else:
        # for non-zero bigWig signal, take log10 of mean signal and save for z-score calculation
        sig[1].append(math.log(float(line[4]),10)) # 10 to specificy log10
        sig[0].append((float(line[4])))
        # only non-zeros contribute
        calculate.append(math.log(float(line[4]),10))
        masterPeak.append(line[0])

# compute mean and std of log10 signal
# non-zero peaks only: using zero peaks would pull the distribution down
lmean=numpy.mean(calculate)
lstd=numpy.std(calculate)

# output z-scores in original peak order
i=0
for entry in sig[1]:
    if entry != "Zero":
        # calculate standard z-score: (log10_signal - mean) / std
        # output peak_name, zscore, raw_signal, log10_signal
        print(masterPeak[i], "\t", (entry-lmean)/lstd, "\t", sig[0][i], "\t", sig[1][i])
    else:
        # if peak signal = 0, hardcode to -10
        print(masterPeak[i], "\t", -10, "\t", 0, "\t", 0)
    i+=1
