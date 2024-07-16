#############################################################################
# Simple pileup parser to count alleles at each reference site
#
# Gregg Thomas, May 2024
#############################################################################

import sys
import re

#############################################################################

def splitBases(bases):
    bases = re.sub(r'\^.', '', bases);
    # Remove ^ and the following character from the bases string
    # ^ indicates the start of a read and the following character is the mapping quality

    indels = re.findall(r'(\+|-)(\d+)', bases);
    # For the current column of base information, split out the indel strings

    parsed_bases = [];
    # Initialize an empty list to store the parsed bases

    for indel in indels:
        sign, length = indel;
        length = int(length);
        # For each indel, extract the sign and length of the indel

        start = bases.index(sign + str(length));
        end = start + length + len(str(length)) + 1;
        # Find the start and end positions of the indel in the base string
        # The start will always be the first occurrence of the sign and length
        # since characters are removed from the bases string as it is parsed

        parsed_bases.extend(list(bases[:start]));
        # Add the bases before the indel to the result list

        parsed_bases.append(bases[start:end]);
        # Add the indel to the result list

        bases = bases[end:];
        # Remove all the characters up to the end of the indel from the bases string

    parsed_bases.extend(list(bases));
    # Add the remaining bases to the result list

    return parsed_bases;

#############################################################################

if len(sys.argv) != 3:
    print("Usage: python pileup_parse.py <input_pileup_file> <output_file>");
    sys.exit(1);
# Check if the correct number of arguments were provided    

pileup_file = sys.argv[1];
output_file = sys.argv[2];
# Extract the input and output file names from the command line arguments

with open(pileup_file, 'r') as pileup, open(output_file, 'w') as parsed:
## Open the input and output files

    headers = ['ref.id', 'ref.pos', 'ref.base', 'total.reads', 'A', 'T', 'C', 'G', 'N', 'ins', 'del'];
    parsed.write('\t'.join(headers) + '\n');
    # Write the headers to the output file

    for line in pileup:
    ## Iterate over the lines in the pileup file

        cols = line.strip().split('\t');
        # Split the line into columns

        ref_base = cols[2].upper();
        # Extract the reference base

        total_reads = int(cols[3]);
        # Extract the total number of reads at the site

        counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'N': 0, 'ins': 0, 'del': 0, '*' : 0};
        # Initialize a dictionary to store the counts of each base

        parsed_bases = splitBases(cols[4]);
        # Parse the bases in the pileup column

        for base in parsed_bases:
        ## Iterate over the parsed bases

            if base.upper() in counts:
                counts[base.upper()] += 1;
            # If the base is a standard base, increment the count for that base
            
            elif base in ['.', ',']:
                counts[ref_base] += 1;
            # If the base is a period or comma, increment the count for the reference base

            elif base.startswith('+'):
                counts['ins'] += 1
            # If the base is an insertion, increment the count for insertions
                            
            elif base.startswith('-'):
                counts['del'] += 1
            # If the base is a deletion, increment the count for deletions

            elif base == '*':
                counts['*'] += 1

            #else:
            #    print(base, "skipped")
            # No else statement, so all other characters are skipped
        ## End base loop
       
        sum_counts = sum(v for k, v in counts.items() if k not in ['del', 'ins']);
        #sum_counts = sum(v for k, v in counts.items() if k != 'del');
         # Calculate the sum of the counts, excluding deletions

        if sum_counts != total_reads:
            print(f"{line}\n\n{parsed_bases}\n\n{counts}\n\nERROR: sum of counts ({sum_counts}) does not match total number of reads ({total_reads})");
            sys.exit(1);
        # Check if the sum of the counts is equal to the total number of reads


        outline = cols[:4] + [str(counts[base]) for base in ['A', 'T', 'C', 'G', 'N', 'ins', 'del']];
        parsed.write('\t'.join(outline) + '\n');
        # Write the counts to the output file
        #print("----")
    ## End line loop
## Close files

#############################################################################