import sys
counter=0
with open(sys.argv[1], 'r') as inf:
    for line in inf:
        for i,x in enumerate(line.strip().split("\t")):
            if x == "contigs":
                print(i + 2) # unix indexes at 1, and line starts with a blak tab
                sys.exit()
raise ValueError("contigs not found in header")
