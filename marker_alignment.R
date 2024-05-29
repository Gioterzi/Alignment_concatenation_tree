
library(tidyverse)
library(DECIPHER)
library(bioseq)
library(Biostrings)
library(ape)
library(sangerseqR)



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

replace_plus_with_minus <- function(seq) {
  return(chartr("+", "-", seq))
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

read_ab1_df <- function(file) {
  ab1_data <- read.abif(file)
  if (inherits(ab1_data, "abif")) {
    seq <- as.character(sangerseqR::primarySeq(sangerseqR::sangerseq(ab1_data)))
    return(DNAStringSet(seq, names = basename(file)))
  } else {
    stop("The file does not seem to be a valid AB1 file.")
  }
}


setwd("/Users/bross/Desktop/AIMS/Analysis/markers_alignment")

########### Align and pull together all sequences from different species based on marker
psba_all <- sort(list.files(path = "/Users/bross/Desktop/AIMS/Analysis/markers_alignment/psbA", pattern = "psbA", full.names = TRUE))

fasta_files <- c("final_sequences.fasta", "all_markers.fasta")
combined_sequences <- DNAStringSet()
for (file in fasta_files) {
  combined_sequences <- c(combined_sequences, readDNAStringSet(file))
}

for (file in psba_all) {
  combined_sequences <- c(combined_sequences, readDNAStringSet(file))
}

keep_indices <- !grepl("SCF049-2|SCF049-3|SCF055-4|SCF055-5", names(combined_sequences))
filtered_sequences <- combined_sequences[keep_indices]



marker_list <- c("Cob", "Cox1", "cp23s", "LSU", "psbA")

aligned_sequences_list <- list()
# Iterate over each marker
for (substring in marker_list) {
  # Initialize a DNAStringSet object to store sequences for the current marker
  marker_sequences <- DNAStringSet()
  
  # Iterate over combined sequences
  for (seq_index in 1:length(filtered_sequences)) {
    # Check if the sequence name contains the current marker
    if (grepl(substring, names(filtered_sequences)[seq_index])) {
      # If it contains the marker, add the sequence to marker_sequences
      marker_sequences <- c(marker_sequences, filtered_sequences[seq_index])
    }
  }
  
  # Perform multiple sequence alignment
  aligned_sequences <- marker_sequences %>%
    OrientNucleotides() %>%
    AlignSeqs()
  
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
  mutate(sample_id = str_replace(label, "(Cob|Cox1|cp23s|LSU|psbA)", ""),
         label = str_extract(label, "(Cob|Cox1|cp23s|LSU|psbA)"))


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


fasta_files_prex <- c("aligned_sequences_Cob.fasta", "aligned_sequences_Cox1.fasta", "aligned_sequences_cp23s.fasta", "aligned_sequences_LSU.fasta", "aligned_sequences_psbA.fasta")

filtered_df <- fasta_files_prex %>%
  map(read_fasta_df) %>% # read in all the files
  purrr::reduce(rbind) # reduce with rbind into one dataframe

final_sequences <- filtered_df %>%
  mutate(sequence = case_when(str_detect(label, "Cox1") ~ str_sub(sequence, start = 20, end = 900),
                              str_detect(label, "LSU") ~ str_sub(sequence, start = 70, end = 500),
                              str_detect(label, "cp23s") ~ str_sub(sequence, start = 50, end = 500),
                              str_detect(label, "Cob") ~ str_sub(sequence, start = 20, end = 850),
                              str_detect(label, "psbA") ~ str_sub(sequence, start = 1001, end = 1050),
                              TRUE ~ str_sub(sequence, start = 10, end = 900)))

aligned_plotting <- final_sequences %>%
  mutate(sample_id = str_replace(label, "(Cob|Cox1|cp23s|LSU|psbA)", ""),
         label = str_extract(label, "(Cob|Cox1|cp23s|LSU|psbA)"))


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

ggsave("marker_alignment_all_Cladocopium_trimmed.pdf", plot = p1, device="pdf", width = 10, height = 4)

trimmed_fasta_filenames <- c("trimmed_Cox1.fasta", "trimmed_LSU.fasta", "trimmed_cp23s.fasta", "trimmed_Cob.fasta")

# Write trimmed sequences to respective trimmed FASTA files
for (i in seq_along(trimmed_fasta_filenames)) {
  # Filter dataframe for sequences with corresponding labels
  trimmed_sequences <- final_sequences[grepl(gsub("trimmed_", "", gsub(".fasta", "", trimmed_fasta_filenames[i])), final_sequences$label), ]
  # Write sequences to trimmed FASTA file using your custom function
  write_fasta(trimmed_sequences$sequence, trimmed_sequences$label, trimmed_fasta_filenames[i])
}

fasta_files_prex2 <- c("trimmed_Cox1.fasta", "trimmed_LSU.fasta", "trimmed_cp23s.fasta", "trimmed_Cob.fasta") #, "aligned_sequences_psbA.fasta")
fasta_files_2 <- lapply(fasta_files_prex2, readDNAStringSet)
# Define substrings to search for
species <- c("SCF049","SCF055", "madreporum", "sodalum", "patulum", "proliferum", "vulgare")


concatenate_sequences <- function(species_name, fasta_list) {
  concatenated_seq <- DNAString("")
  for (fasta in fasta_list) {
    # Find the sequence for the given species in the current fasta file
    seq_index <- which(sapply(names(fasta), function(x) grepl(species_name, x)))
    if (length(seq_index) == 1) {
      concatenated_seq <- c(concatenated_seq, fasta[[seq_index]])
    } else {
      stop(paste("Sequence for species", species_name, "not found or duplicated in one of the FASTA files"))
    }
  }
  concatenated_seq
}

# Loop over each species and generate concatenated sequences
for (sp in species) {
  cat("Processing species:", sp, "\n")
  concat_seq <- concatenate_sequences(sp, fasta_files_2)
  
  # Create a DNAStringSet and set the name to the species name
  concat_seq_set <- DNAStringSet(concat_seq)
  names(concat_seq_set) <- sp
  
  # Write the concatenated sequence to a new FASTA file
  writeXStringSet(concat_seq_set, paste0(sp, "_concatenated.fasta"))
}

all_sequences <- list.files(pattern="concatenated.fasta", full.names = TRUE)
fasta_contents <- lapply(all_sequences, readLines)
all_sequences_content <- unlist(fasta_contents)
writeLines(all_sequences_content, "all_sequences.fasta")


clean <- readDNAStringSet("all_sequences.fasta")

clean <- lapply(clean, replace_plus_with_minus)
clean <- DNAStringSet(clean)
writeXStringSet(clean, "tree.fasta")


#raxml-ng -msa tree.fasta -prefix output_tree -model GTR -seed 12345 -all -threads 4 

tree <- read.tree(text = "(sodalum:0.001542,vulgare:0.000385,(madreporum:0.001156,(patulum:0.000771,((SCF055:0.000001,SCF049:0.000001):0.000772,proliferum:0.000001):0.001029):0.000514):0.000001);
")

pdf("first_tree.pdf")
plot(tree)
dev.off()











##########################


       #psbA - now included in the all.inclusive code up here


##########################
psba_all <- sort(list.files(path = "/Users/bross/Desktop/AIMS/Analysis/markers_alignment/psbA", pattern = "psbA", full.names = TRUE))

psba_all_fasta <- psba_all %>%
  map(read_fasta_df) %>% # read in all the files
  purrr::reduce(rbind) # reduce with rbind into one dataframe

psba_strings <- psba_all_fasta %>%
  deframe() %>%
  as_dna() %>%
  dna_to_DNAStringset()

alligned_psba <- psba_strings %>%
  OrientNucleotides() %>%
  AlignSeqs()

writeXStringSet(alligned_psba, paste0("aligned_sequences_psbA.fasta"))

alligned_df <- read_fasta_df("aligned_sequences_psbA.fasta")

key <- alligned_df %>%
  tibble::rownames_to_column(var = "id")

long_sequences <- str_split(key$sequence, "") %>%
  reshape2::melt() %>%
  group_by(L1) %>%
  mutate(x = row_number(),
         L1 = as.character(L1)) %>%
  left_join(., key, by = c("L1" = "id")) %>%
  ungroup()

# Plot alignment
p2 <- ggplot(long_sequences, aes(y = label, x = x)) +
  geom_tile(aes(fill = value), size = 1, name = "base") +
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
p2

trim <- alligned_df %>%
  mutate(sequence = case_when(str_detect(label, "psbA") ~ str_sub(sequence, start = 900, end = 1150),
                              TRUE ~ str_sub(sequence, start = 30, end = 600)))

key <- trim %>%
  tibble::rownames_to_column(var = "id")

long_sequences <- str_split(key$sequence, "") %>%
  reshape2::melt() %>%
  group_by(L1) %>%
  mutate(x = row_number(),
         L1 = as.character(L1)) %>%
  left_join(., key, by = c("L1" = "id")) %>%
  ungroup()

p2 <- ggplot(long_sequences, aes(y = label, x = x)) +
  geom_tile(aes(fill = value), size = 1, name = "base") +
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
p2


ggsave("psbA_alignment_all_Cladocopium.pdf", plot = p2, device="pdf", width = 10, height = 4)

tree <- read.tree(text = "((madreporum-psbA:0.033327,patulum-psbA:0.019071):0.041766,((proliferum-psbA:0.000001,SCF049-psbA:0.022572):0.268198,vulgare-psbA:0.028671):0.032945,sodalum-psbA:0.030775);")
plot(tree)




