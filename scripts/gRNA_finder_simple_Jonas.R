##################################################################################
#
#   simple R coding for finding possible gRNAs within a desired region
#   
#   input: -genomeName: e.g. "BSgenome.Scerevisiae.UCSC.sacCer3"
#          -chr: chromosome, e.g. "chr02"
#          -start: start of region, e.g. 4000
#          -end: end of region, e.g. 4500
#
#   result: data frame with list and aditional information about the possible gRNAs
#
#  -Jonas, October 2015
#
###################################################################################

design_simple_gRNA = function(genomeName, chr, start, end) {
  
  if(missing(genomeName)){
    stop("Name of the genome is missing.")
  }
  if(missing(chr)){
    stop("Chromosome was not chosen (correctly).")
  }
  if(missing(start)){
    stop("Start of your region is missing.")
  }
  if(missing(end)){
    stop("End of your region is missing.")
  }
  
  require(Biostrings)
  ## genome is the genome sequence as DNAStringSet object
  ## load genome data, for yeast: BSgenome.Scerevisiae.UCSC.sacCer3
  ## chr is the chromosome where my target sequence is, use format 'chr01', 'chr14' etc
  ## start and end define the target region
  library(genomeName, character.only=TRUE)
  genome <- get(genomeName)
  
  ##the following loop checks and changes the chromosome names 
  nameSpace = vector()
  for(y in seqnames(genome)){
    ##get roman number as String
    romanNumber = gsub("chr", "", y)
    ##check if it is NOT a number! else use the provided chromosome name
    if((grepl("[0-9]", romanNumber)) == FALSE){
      number = 0
      name = ""
      ##filter mitochondrial chromosome
      if (romanNumber!="M") {
        number = as.numeric(as.roman(romanNumber))
        if(nchar(number) == 1){
          number = paste("0", number, sep="")
        }
        name = paste("chr", number, sep="")
      } else {
        name = "chrM"
      }
      nameSpace = append(nameSpace, name, (length(nameSpace)+1))
    } else {
      nameSpace = append(nameSpace, y, (length(nameSpace)+1))
    }
  }
  ##rename seqnames with the created vector
  seqnames(genome) <- nameSpace
  
  ##obtain the desired sequence as an X-String object
  ##forward strand: take into account that start is 17 bases earlier and the end 6 bases later due to the cutting site
  forStrand_SubSeq <- DNAStringSet(getSeq(genome, chr, start=(start-17), end=(end+6), strand="+"))
  ##reverse strand: take into account that the end is 17 bases later and the start 6 bases earlier due to the cutting site
  revStrand_SubSeq <- DNAStringSet(getSeq(genome, chr, start=(start-6), end=(end+17), strand="-"))
  
  sequence <- c(forStrand_SubSeq, revStrand_SubSeq)
  
  ##create vectors holding the correct positions
  sequenceFor = c((start-17):(end+6))
  sequenceRev = c((end+17):(start-6))
  
  ## given a sequence, find the PAM site (with PAM.MOTIF) and retrieve the sequence as gRNA
  ## if the given sequence is not long enough (at the boundary), it ignores these PAM sites
  GRNA.LENGTH=20
  ## define PAM WITHOUT N, later: make it an option
  PAM.MOTIF ="GG"

  mysubseq = function(x,ir) {
    ## given ir, IRanges -> overgiven from gRNA.seq = (...)
    ## x, DNAStringSet -> overgiven from gRNA.seq = (...)
    ## retrieve sequences
    ## output will looke like:
    ################################################
    ##     pStart pEnd                  seq
    ##  1      24   43 TCTCGATTTGTTACTTGATT
    ##  2      54   73 AATCCCACCACATCCATCCA
    ##  3     190  209 CCGTCAAGCGTCTTTAGTCG
    ##  4     230  249 TGTTGTAATCAGTGTACAGT
    ##  5     271  290 TCGACGTATACGTTGCGCTG
    ################################################
    if(length(ir)>0 & (sum((start(ir) -GRNA.LENGTH)>0)>0)) {
      sel = which((start(ir) - 1 - GRNA.LENGTH)>0)
      
      rv = data.frame(
        pStart =  start(ir)[sel] - 1 - GRNA.LENGTH,
        pEnd   =  end(ir)[sel] - nchar(PAM.MOTIF) - 1, ## minus the N site
        stringsAsFactors=FALSE)
      rv$seq = as.character(subseq(rep(x, length(sel)),
                                   start=rv$pStart, end=rv$pEnd))
      rv
      
    } else {
      rv = NULL
    }
    rv
  }
  #    if(is.null(names(seqs))) {
  #        names(seqs) = paste0("target",1:length(seqs))
  #    }
  #seqs = fsa
  
  ##in order to count the GC content of a given String
  gcContent <-function(x) {
    if(is.character(x))
      x = DNAStringSet(x)
    alf <- alphabetFrequency(x, as.prob=TRUE)
    rowSums(alf[,c("G", "C"),drop=FALSE])
  }
  
  rv = do.call(rbind, lapply(1:2, function(count) {
    
    strand <- sequence[count]
    pam_pos = vmatchPattern(PAM.MOTIF,strand,fixed=FALSE)
    ##result of pam_pos looks like the matrix of the next step with fewer rows/hits
    
    ## merge the polymers (e.g. if there are hits like GGG)
    ## pam_pos returns a matrix with start and end of the found PAM sequence within the searched strand
    ## result looks like this:
    ###########################
    ##    IRanges of length 15
    ##    start end width
    ##    [1]     45  46     2
    ##    [2]     75  76     2
    ##    [3]    211 212     2
    ##    [4]    251 252     2
    ##    [5]    292 293     2
    ##    ...    ... ...   ...
    ##    [11]   353 354     2
    ##    [12]   370 371     2
    ##    [13]   420 421     2
    ##    [14]   451 452     2
    ##    [15]   478 479     2
    #############################
    pam_pos = lapply(pam_pos, function(x){
      rv = reduce(x)
      sel = width(rv)>nchar(PAM.MOTIF)
      start(rv[sel]) = start(rv[sel]) + 1
      end(rv[sel]) = start(rv[sel]) + 1
      rv
      
    })
    ##here, length is always 1 since pam_pos and strand are lists with one element
    ##don't know why we do it that way
    gRNA.seq = lapply(1:length(pam_pos),  function(i){
      mysubseq(strand[i], pam_pos[[i]])
    }
    )
    names(gRNA.seq) = names(strand)
    
    # desired output
    # seqnames  start	end	width	strand	seq	gc
    # chr04	411862	411881	20	-	TCTACCAGCATTCAAGTGGC	0.5
    
    ##get the right strand sign
    ##and sort for correct sorting vector
    strandSign = ""
    sequenceVector = c()
    if(count == 1) {
      strandSign = "+"
      sequenceVector = sequenceFor
    }
    else if (count == 2) {
      strandSign = "-"
      sequenceVector = sequenceRev
    }
    
    ## handle if there are no results on the strand
    if(!is.null(gRNA.seq[[1]])) {
     
      finalDataFrame = data.frame(
        chro=rep(chr, sapply(gRNA.seq,nrow)),
        start=
          if(count==1){
            sequenceVector[unlist(sapply(gRNA.seq, "[","pStart"))]
          } 
        else if (count == 2){
          sequenceVector[unlist(sapply(gRNA.seq, "[","pEnd"))]
        },
        end=
          if(count==1){
            sequenceVector[unlist(sapply(gRNA.seq, "[","pEnd"))]
          } 
        else if (count == 2){
          sequenceVector[unlist(sapply(gRNA.seq, "[","pStart"))]
        },
        width=abs(sequenceVector[unlist(sapply(gRNA.seq, "[","pStart"))] - sequenceVector[unlist(sapply(gRNA.seq, "[","pEnd"))]) + 1,
        strand=rep(strandSign, sapply(gRNA.seq,nrow)),
        gRNA_Seq=unlist(sapply(gRNA.seq, "[","seq")),
        gc=gcContent(unlist(sapply(gRNA.seq, "[","seq"))),
        stringsAsFactors=FALSE)
      
      return(finalDataFrame)
      
    }
    
  }))
  
  ##process structure of rv with indices
  indicesVector = c(1:nrow(rv))
  for(i in 1:nrow(rv)){
    indicesVector[i] = paste("gRNA_", i, sep="")
  }
  rownames(rv) <- indicesVector
  
  print(rv)
}
