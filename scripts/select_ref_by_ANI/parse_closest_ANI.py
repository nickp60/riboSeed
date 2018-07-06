import sys
counter=0
matrix = []
with open(sys.argv[1], 'r') as inf:
    for line in inf:
        thisline = line.strip().split("\t")
        matrix.append(thisline)
        # for i,x in enumerate(thisline):
        #     if x == "contigs" or x == "scaffolds":
                # print(i + 2) # unix indexes at 1, and line starts with a blak tab
#                 sys.exit()

# raise ValueError("contigs not found in header")
# just start by getting the first forward comparison (ie, the vertical part of the matrix
# we assume rthere is just one entry

# get the index of the row with id "contigs"
contigs_idx = [i for i, x in enumerate(matrix) if x[0].startswith("contigs")][0]

# select our headers
headers = matrix[0]
# note here that we have to trim the first column from the row of interest, we have one less header column than the rest of the rows, cause the first one is empty
line_of_interest = sorted(zip(headers, matrix[contigs_idx][1:]), key=lambda x: x[1], reverse=True)

# bam
# we need to do this stupid loop thing because if there is a score tie, the contigs
# entry won't neccesdarily be first.  So, we go through the line, printing the first entry that doesnt start with "contigs" and isn't identical (probably an artifact)
for pair in line_of_interest:
    # the new version appends _## to the names
    if not pair[0].startswith("contigs") and pair[1] != "1":
        print("_".join(pair[0].split("_")[0:-1]))
        sys.exit(0)
    
# print(line_of_interest[1][0])
