##install dependency packages
#process Go terms
install.packages("ontologyIndex")
#process KEGG terms
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("KEGGREST")
#encrichment analysis
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")

##make Go-Name related file
#download go.obo which includes all Go information at https://geneontology.org/docs/download-ontology/go.obo
library(ontologyIndex)
ont<-get_ontology("go.obo")
GO_db<-data.frame(GO=names(ont$name),Name=as.character(ont$name),stringsAsFactors = FALSE)
write.csv(GO_db,"GO_db.csv",row.names = F)

OR
# Read the GO OBO file as a character vector
obo_content <- readLines("go.obo")
# Initialize vectors to store GO ID, Name, and Namespace
go_id <- character()
go_name <- character()
go_namespace <- character()
# Initialize a flag to indicate whether a new term is being processed
new_term <- FALSE
# Process each line of the OBO file
for (line in obo_content) {
  if (startsWith(line, "[Term]")) {
    new_term <- TRUE
  } else if (new_term) {
    if (startsWith(line, "id:")) {
      go_id <- c(go_id, substring(line, 5))
    } else if (startsWith(line, "name:")) {
      go_name <- c(go_name, substring(line, 7))
    } else if (startsWith(line, "namespace:")) {
      go_namespace <- c(go_namespace, substring(line, 12))
    } else if (line == "") {
      new_term <- FALSE
    }
  }
}

# Create a data frame with GO ID, Name, and Namespace
go_df <- data.frame(
  id = go_id,
  name = go_name,
  namespace = go_namespace
)

# Display the first few rows of the result
head(go_df)



##make KEGG-Name related file, the same as Go-Name
library(KEGGREST)
kegg_path <- keggList("pathway")
kegg_path_db <- data.frame(Pathway = names(kegg_path), Name = as.character(kegg_path), stringsAsFactors = FALSE)
write.csv(kegg_path_db,"KEGG_db.csv",row.names = F)

##make Go-Gene and Kegg-Gene related file
library(dplyr)
emapper<-read.delim("MM_1uvpq7b3.emapper.annotations.tsv",sep="\t",skip=4)
emapper <- emapper %>% dplyr::filter(!grepl("#", X.query))

GO_gene <- emapper %>% dplyr::select("X.query", "GOs")
GO_gene <- GO_gene[GO_gene[,2] != "-", ]
GO_gene <- GO_gene %>% tidyr::separate_rows(GOs, sep = ",") %>% dplyr::rename(GID = X.query, GO = GOs)

Kegg_gene<- emapper %>% dplyr::select("X.query", "KEGG_Pathway")
Kegg_gene <- Kegg_gene[Kegg_gene[,2] != "-", ]
Kegg_gene <- Kegg_gene %>% tidyr::separate_rows(KEGG_Pathway, sep = ",") %>% dplyr::rename(GID = X.query, Kegg_path = KEGG_Pathway)

##enrichment analysis
enrich.res <- enricher(gene = geneList, # gene list
         pvalueCutoff = 0.05, 
         pAdjustMethod = "fdr",
         qvalueCutoff = 0.05, 
         TERM2GENE = TERM2GENE, 
         TERM2NAME = TERM2NAME)


gene<-read.table("S2.significant.constrict.genes",header = F)
gene<-as.character(gene[,1])
GOdb<-read.csv("GO_db.csv",header = T)
GOgene<-read.csv("S2_Gene_Go.csv",header = T)
term2gene<-GOgene[,c(2,1)]
term2name<-GOdb[,c(1,2)]
cons_go_all<-enricher(gene=gene,pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "fdr",TERM2GENE = term2gene,TERM2NAME = term2name)
cons_go<-enricher(gene=gene,pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "fdr",TERM2GENE = term2gene,TERM2NAME = term2name)
write.csv(cons_go_all,"S2_cons_go_all.csv",row.names = F)
barplot(cons_go)
dotplot(cons_go)


gene<-read.table("S2.significant.constrict.genes",header = F)
gene<-as.character(gene[,1])
KEGGdb<-read.csv("KEGG_db.csv",header = T)
KEGGgene<-read.csv("S2_Gene_kegg_new.csv",header = T)
term2gene<-KEGGgene[,c(2,1)]
term2name<-KEGGdb[,c(1,2)]
kegg_all<-enricher(gene=gene,pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "fdr",TERM2GENE = term2gene,TERM2NAME = term2name)
kegg<-enricher(gene=gene,pvalueCutoff = 0.05,qvalueCutoff = 0.05,pAdjustMethod = "fdr",TERM2GENE = term2gene,TERM2NAME = term2name)
write.csv(kegg_all,"S2_constrict_kegg_all.csv",row.names = F)
dotplot(kegg)
barplot(kegg)
