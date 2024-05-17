
library(tidyverse)
library(DECIPHER)
library(bioseq)
library(Biostrings)
library(ape)



read_fasta_df <- function (file = "") {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 
                                                       1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  tib <- tibble(label = gsub(">", "", fasta[ind]), sequence = seqs)
  return(tib)
}

read_fasta_df <- function (file = "") {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 
                                                       1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  tib <- tibble(label = gsub(">", "", fasta[ind]), sequence = seqs)
  return(tib)
}

write_fasta <- function(sequences, names, filename) {
  con <- file(filename, "w")
  for (i in seq_along(sequences)) {
    writeLines(paste0(">", names[i]), con)
    if (!is.character(sequences[[i]])) {
      sequences[[i]] <- as.character(sequences[[i]])
    }
    writeLines(sequences[[i]], con)
  }
  close(con)
}


write_fasta_df <- function (data, filename) 
{
  fastaLines = c()
  for (rowNum in 1:nrow(data)) {
    fastaLines = c(fastaLines, as.character(paste(">", 
                                                  data[rowNum, "label"], sep = "")))
    fastaLines = c(fastaLines, as.character(data[rowNum, 
                                                 "sequence"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

dna_to_DNAbin <- function (dna){
  DNAbin <- as_DNAbin(dna)
  names(DNAbin) <- names(dna)
  return(DNAbin)
}
dna_to_DNAStringset <- function(x) 
{
  bioseq:::check_dna(x)
  DNAstr <- DNAStringSet(paste(x))
  names(DNAstr) <- names(x)
  return(DNAstr)
}

DNAStringSet_to_dna <- function(x){
  x_dna <- as_dna(paste(x))
  names(x_dna) <- names(x)
  res <- tibble(label = names(x), sequence = x_dna)
  return(res)
}

# Convert DNAstringset to DNAbin
DNAStringSet_to_DNAbin <- function(DNAStringSet){
  DNAbin <- as.DNAbin(DNAStringSet)
  return(DNAbin)
}

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2315-y
palette <- c("A" = "#46ff2d", 
             "G" = "#ffae01", 
             "C" = "#f24641", 
             "T" = "#4294fa", 
             "K" = "#8b4816",
             "M" = "#83831f",
             "R" = "#ffff81",
             "S" = "#ff9d80",
             "Y" = "#e381f2",
             "W" = "#80fff2",
             "V" = "#fde4b8",
             "B" = "#f9c1bf",
             "H" = "#c0d9f9",
             "D" = "#c7ffba",
             "U" = "#8989fb",
             "N" = "black", 
             "-" = "white",
             "+" = "deeppink")


pal_df <- data.frame(names = names(palette), col = palette)


setwd("/Users/bross/Desktop/AIMS/Analysis/markers_alignment")

########### Align and pull together all sequences from different species based on marker

fasta_file1 <- "consensus_sequences.fasta"
sequences1 <- readDNAStringSet(fasta_file1)
fasta_file2 <- "all_markers.fasta"
sequences2 <- readDNAStringSet(fasta_file2)

fasta_files <- c("consensus_sequences.fasta", "all_markers.fasta")
combined_sequences <- DNAStringSet()
for (file in fasta_files) {
  combined_sequences <- c(combined_sequences, readDNAStringSet(file))
}

marker_list <- c("Cob", "Cox1", "cp23s")

aligned_sequences_list <- list()
# Iterate over each marker
for (substring in marker_list) {
  # Initialize a DNAStringSet object to store sequences for the current marker
  marker_sequences <- DNAStringSet()
  
  # Iterate over combined sequences
  for (seq_index in 1:length(combined_sequences)) {
    # Check if the sequence name contains the current marker
    if (grepl(substring, names(combined_sequences)[seq_index])) {
      # If it contains the marker, add the sequence to marker_sequences
      marker_sequences <- c(marker_sequences, combined_sequences[seq_index])
    }
  }
  
  # Perform multiple sequence alignment
  aligned_sequences <- AlignSeqs(marker_sequences)
  
  # Store aligned sequences in the list
  aligned_sequences_list[[substring]] <- aligned_sequences
  
  # Write the aligned sequences to a file
  writeXStringSet(aligned_sequences, paste0("aligned_sequences_", substring, ".fasta"))
}

#######################
# View the alignments
aligned_df <- data.frame()
for(i in 1:length(aligned_sequences_list)){
  a_pair <- aligned_sequences_list[[i]] %>% writeXStringSet("temp_file.fasta")
  a_df <- read_fasta_df("temp_file.fasta")
  aligned_df <- rbind(aligned_df, a_df)
}

aligned_plotting <- aligned_df %>%
  mutate(sample_id = str_replace(label, "(Cob|Cox1|cp23s)", ""),
         label = str_extract(label, "(Cob|Cox1|cp23s)"))


# Create profile key
key <- aligned_plotting %>%
  tibble::rownames_to_column(var = "id")


# Create long dataframe for ggplot
long_sequences <- str_split(aligned_plotting$sequence, "") %>%
  reshape2::melt() %>%
  group_by(L1) %>%
  mutate(x = row_number(),
         L1 = as.character(L1)) %>%
  left_join(., key, by = c("L1" = "id")) %>%
  ungroup()

# Plot alignment
p1 <- ggplot(long_sequences, aes(y = sample_id, x = x)) +
  geom_tile(aes(fill = value), size = 1, name = "base") +
  facet_wrap(~ label, nrow = 3, scales = "free") +
  geom_vline(xintercept = seq(50, max(long_sequences$x), by = 50), color = "black", linetype = "dashed", size = 0.05) +
  scale_fill_manual(values = palette) +
  theme(aspect.ratio = 0.3,
        axis.title.y = element_blank(),
        strip.text = element_text(margin = margin(0, 0, 10, 0)), # Add space below strip text
        strip.background = element_blank(), # Remove strip background
        strip.placement = "outside",
        axis.text.x = element_text(size = 0.5),
        axis.ticks.width = unit(0.1, "cm"))+
  scale_x_continuous(expand = c(0, 0), breaks = seq(min(long_sequences$x), max(long_sequences$x), by = 100)) +
  xlab("Position")
p1

ggsave("marker_alignment_all_Cladocopium.pdf", plot = p1, device="pdf", width = 10, height = 4)


####### Concatenate based on species

fasta_files_prex <- c("aligned_sequences_Cob.fasta", "aligned_sequences_Cox1.fasta", "aligned_sequences_cp23s.fasta")
fasta_files_2 <- readDNAStringSet(fasta_files_prex)

# Define substrings to search for
species <- c("SCF049-1", "SCF049-2", "SCF049-3", "madreporum", "sodalum", "patulum", "proliferum", "vulgare")

concatenated_sequences <- list()

# Iterate over each substring
for (substring in species) {
  # Initialize a vector to store concatenated sequences
  concatenated_sequences[[substring]] <- DNAStringSet()
  
  # Iterate over each aligned fasta sequence
  for (seq_index in 1:length(fasta_files_2)) {
    # Find sequences with the current substring in their names
    matching_indices <- grep(substring, names(fasta_files_2[seq_index]))
    
    # Extract matching sequences
    matching_sequences <- fasta_files_2[seq_index][matching_indices]
    
    # Concatenate matching sequences to the list
    concatenated_sequences[[substring]] <- c(concatenated_sequences[[substring]], matching_sequences)
  }
  
  # Convert concatenated sequences to character strings
  concatenated_sequences_as_strings <- lapply(concatenated_sequences[[substring]], as.character)
  
  # Combine concatenated sequences into a single string
  concatenated_string <- paste(concatenated_sequences_as_strings, collapse = "")
  
  # Write the concatenated sequences to a new fasta file
  writeLines(c(paste0(">", substring), concatenated_string), paste0("concatenated_sequences_", substring, ".fasta"))
}


all_sequences <- list.files(pattern="concatenated_sequences_.*\\.fasta", full.names = TRUE)
fasta_contents <- lapply(all_sequences, readLines)
all_sequences_content <- unlist(fasta_contents)
writeLines(all_sequences_content, "all_sequences.fasta")





tree <- read.tree(text = "((patulum:0.000001,(madreporum:0.000001,(vulgare:0.000407,sodalum:0.000001):0.000407):0.000509):0.001528,(SCF049-1:0.000618,(SCF049-3:0.000416,SCF049-2:0.000001):0.004063):0.001494,proliferum:0.000001);")

pdf("first_tree.pdf")
plot(tree)
dev.off()
