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
tr -d '\n' < genome\_sequence.txt/fasta \> genome\_sequence\_without\_new\_lines.txt  
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
./Table\_Creation -i /input\_path/full\_genome.txt -o /same\_folder/Off-Spotter\_bin/  
_This progress will need something around 24GB of memory!_

### IV. Loading the Memory

./Load\_Memory -t /folder\_name\_of\_bin\_files/ -g w303  
w303 is the name of the used genome. So far, the program doesn't allow any other name than
one of the four implemented genomes. It doesn't matter which name is used since we are using 
our own genome with the bin files. w303 has been chosen because it is also a yeast genome.  

### V. Searching for Targets  

Once the tables are loaded, we can start searching for (mis)matches of given gRNA within the genome.
There are different options that are explained as follows:  

* -i: input gRNA(s) or file of gRNA(s) (mandatory)
* -f: flag indicating that we use a file (optional)
* -p: used PAM (mandatory): G for NGG, A for NAG, C for NNNNACA and R for NNGRRT
* -n: number of mismatches: \[0-5\] (mandatory)
* -g: genome name: w303 (mandatory)
* -t: path to bin table files (optional)
* -o: path of output file (optional - shows in terminal if not given)
* -m: memory option that has been used (optional)
* -a: path of annotation files if used (optional)  

If giving gRNA(s) with the command, they have to be separated by two '-' characters:  
_ACTGTACGGTACTGCTGCTA--AGCTGTAACAACGTAACGATC_  
Otherwise, it must be a list of gRNA(s) separated by new line characters within the file.  
Example command:  
./Results -i gRNA\_SEQUENCE(S) -pG -n5 -gw303 -o /output\_file.txt  
or:  
./Results -i /gRNA(s)\_file.txt -f -pG -n5 -gw303 -o /output\_file.txt  

### VI. Output  

The output file has different information separated by tabs. Each line represents one hit.
It will return the following format:  
* chromosome name (number)
* strand
* coordinate start
* coordinate end
* gRNA sequence (without PAM)
* genomic hit with PAM and in lowercases: mismatches
* number of mismatches
* annotation info (optional)  

### VII. Detaching Memory

**_Important!_**  
One has to detach the permanently loaded memory by using this command:  
./Detach_Table -gw303




 
