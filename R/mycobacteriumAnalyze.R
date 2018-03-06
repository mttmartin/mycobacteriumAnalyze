#' import clusterProfiler
#' import UniProt.ws

#' Reads in data from the disk. This data should be in CSV format with the
#' first column as the protein names and the second column as the gene names.
#'
#' @param filename The filename of the CSV file you want to read in.
#' @return A dataframe with two columns. The first labeled "protein" the second "gene"
#' @export
get_data <- function(filename)
{
	data <- read.csv(filename, header=TRUE)
	names(data) <- c("protein", "gene")
	return(data)
}

#' Retrieves the UniProt IDs for a given species that are compatible for 
#' KEGG enrichment.  
#'
#' @param species The species you are analyzing (either 'avium' or 'abscessus')
#' @param data The orignal protein data composed of a 'gene' column and a 'protein' column
#' @export
get_uniprot_IDs <- function(data, species)
{
	if (species == "avium")
	{
		uniprot_IDs <- data$protein
	}
	else if (species == "abscessus")
	{
		up <- UniProt.ws::UniProt.ws(taxId=36809)
		entrez_IDs <- clusterProfiler::bitr(data$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mabscessus.eg.db")$ENTREZID
		uniprot_IDs <- select(up, entrez_IDs, "UNIPROTKB", "ENTREZ_GENE")$UNIPROTKB
	}
	else
	{
		stop(paste("get_uniprot_IDs: unrecognized species", species))
	}

	return(uniprot_IDs)
}

#' Retrieves entrez IDs for a given species that are compatible for
#' GO enrichment analysis.
#'
#' @param data The orignal protein data composed of a 'gene' column and a 'protein' column
#' @param species The species you are analyzing (either 'avium' or 'abscessus')
#' @return A vector of entrez IDs
#' @export
get_entrez_IDs <- function(data, species)
{
	if (species == "avium")
	{
		stop("get_entrez_IDs: Avium not supported")
	}
	else if (species == "abscessus")
	{
		up <- UniProt.ws(taxId=36809)
		entrez_IDs <- bitr(data$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mabscessus.eg.db")$ENTREZID
	}
	else
	{
		stop(paste("get_entrez_IDs: unrecognized species", species))
	}
	
	return (entrez_IDs)
}


#' Calculates and returns results for KEGG enrichment given a species
#' and corresponding UniProt IDs.
#'
#' @param species the species you are analyzing (either 'avium' or 'abscessus') 
#' @param uniprot_IDs The corresponding UniProt IDs you want to test
#' @param pvalueCutoff The p-value threshold that is considered significant (optional, default=0.05)
#' @return A results object from clusterProfiler
#' @export
get_KEGG_enrichment <- function(uniprot_IDs, species, pvalueCutoff=0.05)
{
	if (species == "avium")
	{
		kegg_results <- clusterProfiler::enrichKEGG(gene=as.character(uniprot_IDs), organism='mav', keyType='uniprot', pvalueCutoff=pvalueCutoff)
	}
	else if (species == "abscessus")
	{
		kegg_results <- clusterProfiler::enrichKEGG(gene=as.character(uniprot_IDs), organism='mab', keyType='uniprot', pvalueCutoff=pvalueCutoff)
	}
	else
	{
		stop(paste("get_KEGG_enrichment: unrecognized species", species))
	}

	return (kegg_results)
}


write_table <- function(table, filename)
{
	write.csv(table, file=filename, quote=TRUE, row.names=FALSE)
}

#' Writes the results from the KEGG enrichment object to disk as a comma
#' separated value format file.
#' 
#' @param KEGG_enrichment_object The clusterProfiler KEGG enrichment results
#' @param filename The location where the file should be written
#' @return No return value
#' @export
write_KEGG_enrichment <- function(KEGG_enrichment_object, filename)
{
	write_table(KEGG_enrichment_object@result, filename)
}


get_GO_enrichment <- function(species, entrez_IDs, ontology, pvalueCutoff=0.05)
{
	if (species == "avium")
	{
		stop("get_GO_enrichment: M. avium not implemented")
	}
	else if (species == "abscessus")
	{
		library(org.Mabscessus.eg.db)
		GO_results <- enrichGO(gene=entrez_IDs, OrgDb=org.Mabscessus.eg.db, ont=ontology, pAdjustMethod="BH", pvalueCutoff=pvalueCutoff, readable=TRUE, keyType="ENTREZID")

	}
	else
	{
		stop(paste("get_GO_enrichment: unrecognized species", species))
	}

	return(GO_results)
}

write_GO_enrichment <- function(GO_enrichment_object, filename)
{
	write_table(GO_enrichment_object@result)
}

#' Analyzes the data(CSV file with protien and gene information) from a
#' user specified file and carries out GO and KEGG enrichment analysis and
#' finally outputs the results as tables at the user-specified location. 
#' This function is essentially a wrapper around the rest of the functions
#' provided, so that standard analysis can be done easily.
#'
#' @param input_filename The location of file you want to analyze. It should be in CSV format with a protein column and gene column
#' @param output_prefix The output prefix you want added to all output files.
#' @param species The species you are analyzing (either 'avium' or 'abscessum')
#' @param do_KEGG Whether you want KEGG enrichment analysis (optional; default TRUE)
#' @param do_GO Whether you want to do GO enrichment analysis (optional; default FALSE)
#' @export
analyze_mycobacterium_data <- function(input_filename, output_prefix, species, do_KEGG=TRUE, do_GO=FALSE)
{
	data <- get_data(input_filename)
	
	if (do_KEGG)
	{
		uniprot_IDs <- get_uniprot_IDs(data, species)
		KEGG_enrichment_object <- get_KEGG_enrichment(uniprot_IDs, species)
		KEGG_filename <- paste(output_prefix, "_KEGG.csv", sep="")
		write_KEGG_enrichment(KEGG_enrichment_object, KEGG_filename)
	}

	if (do_GO)
	{
		if (species == "avium")
		{
			print("Warning: Avium GO enrichment not supported, skipping.")
		}
		else if (species == "abscessus")
		{
			entrez_IDs <- get_entrez_IDs(data, species)
			GO_enrichment_object <- get_GO_enrichment(entrez_IDs, species)
			GO_filename <- paste(output_prefix, "_GO.csv", sep="")
			write_GO_enrichment(GO_enrichment_object, GO_filename)
		}

	}
}

