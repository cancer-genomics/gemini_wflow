import sys
import pysam

bamFile = sys.argv[1]
region = sys.argv[2]
minPhred = sys.argv[3]
minMAPQ = sys.argv[4]
outFile = sys.argv[5]

def getReadPileup(read, region, minPhred, minMAPQ):

  if read == 1:
    flag_filter = 3468
  elif read == 2:
    flag_filter = 3404
  else:
    print("Error: read must be 1 or 2.")
    exit()

  iter = alnFile.pileup(region = region,
                        max_depth = 1e9,
                        ignore_overlaps = False,
                        min_base_quality = int(minPhred),
                        min_mapping_quality = int(minMAPQ),
                        flag_require = 3,
                        flag_filter = flag_filter,
                        truncate = False)

  for i in iter:

     qnames = i.get_query_names()
     chrom = i.reference_name
     pos = i.reference_pos + 1 # Convert to 1-based coordinate system
     obs = i.get_query_sequences(add_indels = False) # Do not include insertions in the pileup, deletions are an empty string

     for j in range(i.get_num_aligned()):
       obs_j = obs[j]
       if obs_j != "": # Exclude deletions
         if obs_j.islower():
           strand = "-"
           obs_j = obs_j.upper()
         elif obs_j.isupper():
           strand = "+"
         else:
           print("Error: the observed base is neither uppercase nor lowercase..")
           exit()

         if obs_j != "N":
           f.write(qnames[j] + "\t" + str(read) + "\t" + strand + "\t" + chrom + "\t" + str(pos) + "\t" + obs_j + "\n")

alnFile = pysam.AlignmentFile(bamFile, "rb")

f = open(outFile, "w")
f.write("fragment" + "\t" "read" "\t" + "strand" + "\t" + "chrom" + "\t" + "pos" + "\t" + "obs" + "\n")
getReadPileup(read = 1, region = region, minPhred = minPhred, minMAPQ = minMAPQ)
getReadPileup(read = 2, region = region, minPhred = minPhred, minMAPQ = minMAPQ)
f.close()

alnFile.close()

