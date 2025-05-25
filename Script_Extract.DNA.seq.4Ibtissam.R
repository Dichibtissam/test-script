require(rentrez)


# Load the tsv table exported from the Organnelle NCBI web page https://www.ncbi.nlm.nih.gov/datasets/organelle
MetaData.MitoGenome = read.delim("Data/ncbi_dataset_Mitochondrion.tsv", sep = "\t", h=T)
head(MetaData.MitoGenome)


### Données de Ibtissam
Classification.Metadata = read.delim("Data/out_reformat_taxids_with_headers.csv", sep = ";", h=T)

# merge the meta.data table with the accession with the table including the taxonomy of the taxon provided by Ibtissam using tghe TaxonKit tool.

Meta.mata.mitogenome.Metazoa.test = merge(MetaData.MitoGenome, Classification.Metadata, by.x=6 , by.y = 1, all.x = T)

### Select the accession number for the Colepotera only.
Test.Accession.coleo = Meta.mata.mitogenome.Metazoa.test[which(Meta.mata.mitogenome.Metazoa.test$order == "Coleoptera"), "GenBank.accession"]





#######################
### Extract data coleoptera

# Load the function to extract the mitogenome in a fasta file, as well as the meta data and the sequence annotations

source("Scripts/Get.Seq.Metadata.NCBI.R")


cpu0 = Sys.time()
Data.coleo = Get.Seq.Metadata.NCBI(Accession.NB = Test.Accession.coleo, output.Path = "/home/idich/internship/ibtissam/alignement/coleoptera", output.Name = "Test.Extraction.Coleo.bis")
cpu1 = Sys.time()
cpu1 - cpu0 # 13.95946 mins

names(Data.coleo)


## Homogeneization of the gene name

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = tolower(Data.coleo$Annotation.Sequence.DF$gene_name)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub("ox", "o", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub("iii", "3", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub("ii", "2", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub("i", "1", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub("nadh", "nd", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub("nad", "nd", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub(" ", "", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub("atpase", "atp", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.coleo$Annotation.Sequence.DF$gene_name.corrected = gsub("cob", "cytb", Data.coleo$Annotation.Sequence.DF$gene_name.corrected)

Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected = gsub("co21", "co3", Data.Diplostraca$Annotation.Sequence.DF$gene_name.corrected)

table(Data.coleo$Annotation.Sequence.DF$gene_name.corrected)




##########################################################################################
#### Export the 13 genes of interests in 13 fsata file : emprical example with the coleoptera.

require(seqinr)
# Load the Mitogenome exported by the function "Get.Seq.Metadata.NCBI.R" into R in as an
# "alignment" object
Mitogenomes.fasta = seqinr::read.alignment("/home/idich/internship/ibtissam/alignement/coleoptera/Test.Extraction.Coleo.And.outgroup.Align.fas",format = "fasta")

# source the function tosplit the mitogeoes into 13 fasta alignment
source("Scripts/Split.MitoGenome.R")

mkdir Separate.Align

# Run the function, export the 13 gene alignments in fasta format.
# see directly n the file of the function to get the information about the different arguments.ø
Annot.DF.Coleo.test = Split.MitoGenome(input.fasta = Mitogenomes.fasta, Gene.to.get = c("co1", "co2", "co3", "cytb", "atp6", "atp8", "nd1", "nd2", "nd3","nd4", "nd4l", "nd5", "nd6"), Annotations = Data.coleo$Annotation.Sequence.DF, Annotation.col = 6, Access.NB.col = 1, Starting.col = 3, Ending.col = 4, output.folder = "Results/Separate.Align")




