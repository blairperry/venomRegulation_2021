import sys
import os

infile = sys.argv[1]
outfile =  os.path.splitext(infile)[0] + '.singleLine.fa'

with open(infile) as f_input, open(outfile, 'w') as f_output:
    block = []

    for line in f_input:
        if line.startswith('>'):
            if block:
                f_output.write(''.join(block) + '\n')
                block = []
            f_output.write(line)
        else:
            block.append(line.strip())

    if block:
        f_output.write(''.join(block) + '\n')
