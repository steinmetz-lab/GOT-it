## Off-Spotter Workflow


_The following document gives an overview of the required steps for
the prediction of possible (off-) targets of given gRNAs within a genome
by the program Off-Spotter (download link: https://cm.jefferson.edu/downloads/off-spotter-code/)_

\- Jonas Weidenhausen, 05th October 2015

### I. Make the executables

After downloading the zipped files, simply type _make_ in the terminal
in the same location.

### II. Downloading/Creation of Genome Data

_This only has to be done once per genome!_

For the lab yeast strains BY(...) the yeast reference sequence S288C_R64_2_1  has been
downloaded from http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/.
However, the .fasta or .txt genome file(s) has to use a specific format to be accepted by the
program.

The requirements are:

* Chromosomes need to start with a '>'
* The number of the chromosome has to be directly behind '>'
* After the number, any character (excepts a new number of course) is excepted
* Bases must be in **one** line (without new line \n)

_example:_
\>1\_chr1\_S288C\_R64\_
ATAGCTGTCGTAGCTGA(...)GCTGTGACTAGCTAGCATCG
\>2\_chr2\_S288C\_R64\_
ATACGTCGACTGGCTAG(...)CCCGTATGCGAGCTACGCGC

The new line characters can be deleted by the following
terminal command:
tr -d '\n' < genome_sequence.txt/fasta > genome_sequence_without_new_lines.txt
However, new lines after chromosome headers need to be updated either manually or
by an appropriate command. Also, the numbers after '>' have to be written manually
or by another command. In this case, it was done manually because there was only one
genome sequence. Additionally, the information about the mitochondrial chromosome
was deleted (manually) leading to a file consisting of 16 chromosomes (in our case).


### III. Table Creation of all PAM Findings within the Genome Data

_This only has to be done once per genome!_

The program searches for all occurring PAMs within the given genome file. It creates two
tables of all possible hits separated by the according PAM (NGG, NAG, NNNNACA or NNGRRT) with
an assigned value called data.bin and index.bin. Those files can be used in future if using the
same genome.
In order to create those files, the following command is used:
./Table_Creation -i /input_path/full_genome.txt -o /same_folder/Off-Spotter_bin/
_This progress will need something around 24GB of memory!_

### IV. Loading the Memory and Start a Query

 
