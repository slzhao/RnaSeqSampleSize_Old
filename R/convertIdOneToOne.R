##' convertId
##' 
##' A function to convert ID based on the biomaRt package.
##' 
##' A function to convert ID based on the biomaRt package..
##' 
##' @param x the Ids need to be converted.
##' @param verbose Logical. Indicate report extra information on progress or not.
##' @inheritParams biomaRt::getBM
##' @inheritParams biomaRt::useMart
##' @return A converted ID character with the same order of parameter x.
##' @importFrom biomaRt getBM useMart
##' @export
##' @examples x<-c("Q04837","P0C0L4","P0C0L5","O75379","Q13068","A2MYD1","P60709","P30462","P30475","P30479")
##' convertIdOneToOne(x,filters="uniprot_swissprot",verbose=TRUE)
convertIdOneToOne<-function(x,dataset="hsapiens_gene_ensembl",filters="uniprot_swissprot",attributes =c(filters,"entrezgene"),verbose=FALSE) {
	if (verbose) {
		cat("Now conectting with ensembl. Internet acess is needed and it may use 30 seconds.\n")
		flush.console()
	}
	
	ensembl = useMart("ensembl",dataset=dataset)
	newIdTable<-getBM(attributes =attributes,filters=filters,values=x,mart = ensembl)
	newIdTable<-newIdTable[which(newIdTable[,1]!="" & newIdTable[,2]!=""),]
	result<-newIdTable[match(x, newIdTable[,1]),2]
	if (verbose & any(is.na(result))) {
		cat(paste0(length(which(is.na(result)))," Ids cant't be converted to ",attributes[2],", will be set as NA\n"))
	}
	names(result)<-x
	return(result)
}



