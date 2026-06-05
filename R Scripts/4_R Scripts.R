# Title: Multidimensional Beta Diversity Dominated by Replacement Reveals Hidden Complementarity of Moss-Associated Oribatid Mites in Alpine Petrifying Springs
# Authors: Yuanyuan Zhou, Qiang Wei, Yan Shen, Chengyi Li, Hu Chen, Zhihui Wang, Zhaohui Zhang*
# R code implementation and inspection: Yuanyuan Zhou (1181559954@qq.com); and Qiang Wei (wq_tendermoments@163.com)
# Date: June 2026
# R 4.4.1
# Correspondence: Zhaohui Zhang (zhaozhang9@hotmail.com)




library(RColorBrewer)
library(tidyverse)
library(rotl)
library(ape)
library(ggtree)
library(ggplot2)
library(ggtreeExtra)
library(ggnewscale)
library(BAT)
library(ggtern)
library(cowplot)
library(grid)
library(vegan)
library(lme4)
library(car)
library(performance)
library(pheatmap)
library(rdacca.hp)
citation("rotl")



#### 0. Global colour palettes for all visualizations ####

library(RColorBrewer)

# display.brewer.all()
# brewer.pal.info
# row.names(brewer.pal.info)

palette_set3 <- brewer.pal(12, "Set3")   # Discrete colours for landscape types  
palette_prgn <- brewer.pal(11, "PRGn")   # Base diverging palette for correlation plots  

landscape_cols <- c("Spr_TL" = palette_set3[5], "Cas_TL" = palette_set3[7],
                    "Cus_TL" = palette_set3[6], "Rap_TL" = palette_set3[4],
                    "Rim_TL" = palette_set3[10])
palette_corr <- colorRampPalette(palette_prgn)(100)






#### 1. Community and trait data preparation ####

#### 1.1 Community abundance matrix construction: comm_abund ####

comm <- read.csv("oribatida.csv", check.names = FALSE); names(comm)

library(tidyverse)

sample_info <- comm %>%
  dplyr::select(Type, Site) %>%
  dplyr::mutate(SampleID = paste(Type, Site, sep = "_")); head(sample_info)

comm_abund <- comm %>%
  dplyr::select(-Type, -Site) %>%
  dplyr::mutate(dplyr::across(everything(), ~ as.numeric(as.character(.)))) %>%
  as.data.frame()

rownames(comm_abund) <- sample_info$SampleID

## Optional checks: inspect matrix structure, abundance totals, and missing values. 
# head(comm_abund); str(comm_abund); summary(comm_abund); rowSums(comm_abund); colSums(comm_abund); sum(is.na(comm_abund))




#### 1.2 Trait matrix preparation, type conversion, and functional tree construction ####

#### 1.2.1 Trait matrix preparation and type conversion: traits_final ####

traits <- read.csv("traits.csv", check.names = FALSE); summary(traits)

traits_final <- traits %>%
  tibble::column_to_rownames("Genus") %>%
  dplyr::mutate(
    body_length = as.numeric(body_length), body_width = as.numeric(body_width),
    sclerotization = factor(sclerotization, levels = c(1, 2, 3), labels = c("low", "medium", "high"), ordered = TRUE),
    reproduction = factor(reproduction, levels = c(0, 1), labels = c("sexual", "parthenogenetic")),
    pteromorph = factor(pteromorph, levels = c("none", "humeral_projection", "fixed_pteromorph", "movable_pteromorph")),
    claw_number = factor(claw_number, levels = c(1, 2, 3), labels = c("monodactylous", "bidactylous", "tridactylous"), ordered = TRUE),
    life_history = factor(life_history, levels = c(1, 2, 3), labels = c("r_selected", "intermediate", "K_selected"), ordered = TRUE),
    trophic_guild = factor(trophic_guild, levels = c("herbivore", "primary_decomposer", "fungal_feeder", "predator"))
  )

sp_common <- intersect(colnames(comm_abund), rownames(traits_final)); sp_common

## Optional checks: verify taxon matching between community and trait matrices. 
# length(sp_common); setdiff(colnames(comm_abund), rownames(traits_final)); setdiff(rownames(traits_final), colnames(comm_abund)); identical(colnames(comm_abund), rownames(traits_final))
# str(traits_final); summary(traits_final); colSums(is.na(traits_final)); sum(is.na(traits_final))




#### 1.2.2 Functional distance calculation and functional tree construction: fd_tree ####

comm_pa <- comm_abund
comm_pa[comm_pa > 0] <- 1
comm_pa <- as.data.frame(comm_pa)

fd_taxa <- intersect(colnames(comm_pa), rownames(traits_final))
traits_fd <- traits_final[fd_taxa, , drop = FALSE]
comm_pa_fd <- comm_pa[, fd_taxa, drop = FALSE]


## Functional distance among genera based on mixed functional traits.
fd_dist <- BAT::gower(traits_fd, st = "range")


## Build and evaluate candidate functional trees. 
fd_tree  <- BAT::tree.build(fd_dist, func = "nj")
fd_tree1 <- BAT::tree.build(fd_dist, func = "upgma")
fd_tree2 <- BAT::tree.build(fd_dist, func = "best")

BAT::tree.quality(fd_dist, fd_tree)
BAT::tree.quality(fd_dist, fd_tree1)
BAT::tree.quality(fd_dist, fd_tree2)

plot(fd_tree, cex = 0.6)
# ape::write.tree(fd_tree, file = "02 Oribatida_functional_tree_NJ.tre")

comm_pa_fd <- comm_pa_fd[, fd_tree$tip.label, drop = FALSE]
identical(colnames(comm_pa_fd), fd_tree$tip.label)




#### 1.2.3 Functional tree visualization: fd_tree ####

library(ape)
library(ggtree)
library(ggplot2)

fd_tree_plot <- ape::ladderize(fd_tree, right = FALSE)

p_fd_tree_rect <- ggtree(fd_tree_plot) +
  geom_tree(linewidth = 0.27, color = "grey10") +
  geom_tiplab(size = 3.5, color = "grey10", align = TRUE, linetype = "longdash", linesize = 0.20, offset = 0.02) +
  geom_tippoint(size = 1.2, color = "#803F8F") +
  xlab("Functional distance") +
  theme_tree2() +
  theme(legend.position = "none",
        plot.margin = margin(5, 15, 5, 5),
        axis.text.x = element_text(size = 10, color = "grey10"),
        axis.title.x = element_text(size = 10, color = "grey10"),
        axis.line.x = element_line(linewidth = 0.27, color = "grey10"),
        axis.ticks.x = element_line(linewidth = 0.27, color = "grey10"),
        axis.ticks.length = unit(2, "pt")); p_fd_tree_rect

ggsave(filename = "02 Functional tree_rectangular.pdf", plot = p_fd_tree_rect, width = 12, height = 18, units = "cm", device = cairo_pdf)








#### 2. Phylogenetic tree construction and visualization ####

#### 2.1 Phylogenetic tree construction and taxonomic correction: phylo_tree ####

library(rotl)

clean_genera <- stringr::str_replace(colnames(comm_abund), "^X?\\d+_", ""); clean_genera


## Match genus names against the Open Tree of Life taxonomy. 
resolved_names <- rotl::tnrs_match_names(names = clean_genera, context_name = "Animals"); print(resolved_names)
# write.csv(resolved_names, file = "01 resolved_names_OTT_matching_results.csv", row.names = FALSE)


## Identify taxonomic matches requiring manual inspection. 
problem_matches <- resolved_names %>%
  dplyr::filter(is.na(ott_id) | is_synonym == TRUE | approximate_match == TRUE | flags != "" |
                  number_matches > 1 | tolower(unique_name) != tolower(search_string)); print(problem_matches)

## Eight genera require manual inspection, especially Trimalaconothrus, Hermannia, and Protokalumma. 
inspect(resolved_names, taxon_name = "hypochthonius")      # Correct match. 
inspect(resolved_names, taxon_name = "acrotritia")         # Correct match. 
inspect(resolved_names, taxon_name = "nothrus")            # Correct match.
inspect(resolved_names, taxon_name = "epidamaeus")         # Correct match. 
inspect(resolved_names, taxon_name = "Hermannia")          # Correct match, but extra text needs removal. 
inspect(resolved_names, taxon_name = "trimalaconothrus")   # Treated as Tyrphonothrus based on taxonomic references. 
inspect(resolved_names, taxon_name = "achipteria")         # Correct match. 
inspect(resolved_names, taxon_name = "protokalumma")       # Removed from phylogenetic tree construction. 

ids_to_check <- c(891290, 6165554, 709418, 357575, 856695, 7086843, 882337, 99284)

info_list2 <- lapply(ids_to_check, function(x) rotl::taxonomy_taxon_info(x, include_lineage = TRUE))

for (i in seq_along(info_list2)) {
  cat("\n====================\n"); cat("OTT ID:", ids_to_check[i], "\n")
  print(rotl::tax_name(info_list2[[i]])); print(rotl::unique_name(info_list2[[i]])); print(rotl::tax_rank(info_list2[[i]]))
  print(rotl::synonyms(info_list2[[i]])); print(rotl::tax_lineage(info_list2[[i]]))
}


## Build the taxonomic correction table for phylogenetic tree construction. 
name_map <- tibble::tibble(
  Original_Col = colnames(comm_abund), Clean_Input = clean_genera, Matched_Name = resolved_names$unique_name,
  Final_Name = dplyr::case_when(tolower(clean_genera) == "trimalaconothrus" ~ "Tyrphonothrus",
                                tolower(clean_genera) == "protokalumma" ~ NA_character_,
                                tolower(clean_genera) == "hermannia" ~ "Hermannia",
                                TRUE ~ stringr::str_replace(resolved_names$unique_name, " \\(.*\\)", "")),
  ott_id = dplyr::case_when(tolower(clean_genera) == "protokalumma" ~ NA_real_, TRUE ~ resolved_names$ott_id)
) %>%
  dplyr::filter(!is.na(Final_Name), !is.na(ott_id)); name_map

# write.csv(name_map, "name_map_phylogeny.csv", row.names = FALSE)


## Construct the genus-level phylogenetic tree based on Open Tree of Life. 
rotl::is_in_tree(unique(name_map$ott_id))

phylo_tree_raw <- rotl::tol_induced_subtree(ott_ids = unique(name_map$ott_id), label_format = "name"); phylo_tree_raw

phylo_tree_raw$tip.label <- as.character(phylo_tree_raw$tip.label)
phylo_tree_raw$tip.label <- phylo_tree_raw$tip.label %>%
  stringr::str_replace("_ott.*", "") %>%
  stringr::str_replace(" \\(.*\\)", "") %>%
  stringr::str_replace("_in_Opisthokonta.*", "")

phylo_tree_raw$tip.label[phylo_tree_raw$tip.label == "Hermannia_(genus_in_Opisthokonta)"] <- "Hermannia"
phylo_tree_raw$tip.label[phylo_tree_raw$tip.label == "Hermannia_(genus"] <- "Hermannia"

length(phylo_tree_raw$tip.label)
setequal(sort(phylo_tree_raw$tip.label), sort(name_map$Final_Name))
setdiff(sort(phylo_tree_raw$tip.label), sort(name_map$Final_Name))
setdiff(sort(name_map$Final_Name), sort(phylo_tree_raw$tip.label))


## Resolve polytomies and assign branch lengths using Grafen's method. 
set.seed(123)
phylo_tree <- phylo_tree_raw
phylo_tree <- ape::multi2di(phylo_tree)
phylo_tree <- ape::compute.brlen(phylo_tree, method = "Grafen", power = 1)

ape::is.binary(phylo_tree)
ape::is.rooted(phylo_tree)
all(phylo_tree$edge.length > 0)

# ape::write.tree(phylo_tree, file = "01 Oribatida_phylo_tree_final_Grafen.tre")




#### 2.2 Phylogenetic tree visualization with community annotations ####

#### 2.2.1 Alignment between phylogenetic tree and community matrix ####

standardize_taxon_names <- function(x) {
  x %>%
    stringr::str_replace("^X?\\d+_", "") %>%
    stringr::str_replace("_ott.*$", "") %>%
    stringr::str_replace(" \\(.*\\)$", "") %>%
    dplyr::recode("Trimalaconothrus" = "Tyrphonothrus")
}

comm_abund_phylo <- comm_abund %>%
  as.data.frame()

colnames(comm_abund_phylo) <- colnames(comm_abund_phylo) %>%
  standardize_taxon_names(); colnames(comm_abund_phylo)

comm_pa_phylo <- comm_abund_phylo
comm_pa_phylo[comm_pa_phylo > 0] <- 1


## Check taxon matching between the phylogenetic tree and community matrix. 
tree_tips <- phylo_tree$tip.label; tree_tips
comm_taxa <- colnames(comm_abund_phylo); comm_taxa

setdiff(tree_tips, comm_taxa)
setdiff(comm_taxa, tree_tips)   # Protokalumma is excluded. 


## Retain taxa shared by the phylogenetic tree and community matrix. 
common_taxa <- intersect(phylo_tree$tip.label, colnames(comm_abund_phylo)); common_taxa
length(tree_tips); length(comm_taxa); length(common_taxa)

phylo_tree_plot <- ape::keep.tip(phylo_tree, common_taxa)

comm_abund_phylo <- comm_abund_phylo %>%
  dplyr::select(dplyr::all_of(phylo_tree_plot$tip.label))

comm_pa_phylo <- comm_pa_phylo %>%
  dplyr::select(dplyr::all_of(phylo_tree_plot$tip.label))

identical(colnames(comm_abund_phylo), phylo_tree_plot$tip.label)
identical(colnames(comm_pa_phylo), phylo_tree_plot$tip.label)
setdiff(phylo_tree_plot$tip.label, colnames(comm_abund_phylo))
setdiff(colnames(comm_abund_phylo), phylo_tree_plot$tip.label)
dim(comm_abund_phylo); dim(comm_pa_phylo); length(phylo_tree_plot$tip.label)




#### 2.2.2 Outer-ring annotation data for the phylogenetic tree ####

## Relative distribution of each genus across landscape types.
comm_relabund <- vegan::decostand(comm_abund_phylo, method = "total", MARGIN = 2) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleID")

sample_info2 <- sample_info %>%
  dplyr::mutate(SampleID = as.character(SampleID),
                Type = factor(Type, levels = c("Spr_TL", "Cas_TL", "Cus_TL", "Rap_TL", "Rim_TL")))

abund_type <- comm_relabund %>%
  dplyr::left_join(sample_info2 %>% dplyr::select(SampleID, Type), by = "SampleID") %>%
  dplyr::group_by(Type) %>%
  dplyr::summarise(dplyr::across(.cols = dplyr::all_of(phylo_tree_plot$tip.label), .fns = ~ sum(.x, na.rm = TRUE)),
                   .groups = "drop") %>%
  tidyr::pivot_longer(cols = dplyr::all_of(phylo_tree_plot$tip.label),
                      names_to = "Genus", values_to = "Landscape_relative_abundance") %>%
  dplyr::mutate(Genus = factor(Genus, levels = phylo_tree_plot$tip.label),
                Type = factor(Type, levels = c("Spr_TL", "Cas_TL", "Cus_TL", "Rap_TL", "Rim_TL")))


## Occurrence frequency of each genus across samples. 
occ_freq <- comm_pa_phylo %>%
  as.data.frame() %>%
  dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~ mean(.x > 0, na.rm = TRUE) * 100)) %>%
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Genus", values_to = "Occurrence_frequency") %>%
  dplyr::mutate(Genus = factor(Genus, levels = phylo_tree_plot$tip.label), occ_group = "Occurrence")

range(occ_freq$Occurrence_frequency)


## Classify genera by total abundance. 
total_abund <- comm_abund_phylo %>%
  as.data.frame() %>%
  dplyr::summarise(dplyr::across(.cols = dplyr::everything(), .fns = ~ sum(.x, na.rm = TRUE))) %>%
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Genus", values_to = "Total_abundance") %>%
  dplyr::mutate(Genus = factor(Genus, levels = phylo_tree_plot$tip.label), abund_group = "Total abundance")

range(total_abund$Total_abundance)

abundance_status <- total_abund %>%
  dplyr::mutate(Total_relative_abundance = Total_abundance / sum(Total_abundance, na.rm = TRUE) * 100,
                Abundance_status = dplyr::case_when(Total_relative_abundance >= 10 ~ "Dominant",
                                                    Total_relative_abundance >= 1 ~ "Common",
                                                    Total_relative_abundance < 1 ~ "Rare"),
                Abundance_status = factor(Abundance_status, levels = c("Dominant", "Common", "Rare"))) %>%
  dplyr::select(Genus, Total_abundance, Total_relative_abundance, Abundance_status)

abundance_status %>%
  dplyr::count(Abundance_status)

abundance_ring <- abundance_status %>%
  dplyr::mutate(Genus = as.character(Genus),
                Abundance_status = factor(Abundance_status, levels = c("Dominant", "Common", "Rare")),
                value = 1) %>%
  dplyr::filter(Genus %in% phylo_tree_plot$tip.label)




#### 2.2.3 Infraorder annotation for matched genera ####

infraorder_cols <- c("Brachypylina" = "#B3DE69", "Desmonomata" = "#FDB462",
                     "Enarthronota" = "#80B1D3", "Mixonomata" = "#FB8072")


## Extract infraorder information from the OTT lineage. 
get_infraorder <- function(ott_id) {
  tax_info <- rotl::taxonomy_taxon_info(ott_ids = ott_id, include_lineage = TRUE)
  lineage_list <- rotl::tax_lineage(tax_info)
  if (is.null(lineage_list) || length(lineage_list) == 0) return(NA_character_)
  lineage_df <- lineage_list[[1]] %>%
    tibble::as_tibble() %>%
    dplyr::rename_with(tolower)
  if (!all(c("rank", "name") %in% names(lineage_df))) return(NA_character_)
  infraorder <- lineage_df %>%
    dplyr::filter(rank == "infraorder") %>%
    dplyr::pull(name) %>%
    unique()
  if (length(infraorder) == 0 || is.na(infraorder[1])) return(NA_character_)
  return(infraorder[1])
}


## Match each genus in the plotted phylogeny to infraorder information. 
lineage_df <- tibble::tibble(Genus = phylo_tree_plot$tip.label) %>%
  dplyr::left_join(
    name_map %>%
      dplyr::filter(!is.na(ott_id), Final_Name %in% phylo_tree_plot$tip.label) %>%
      dplyr::distinct(Final_Name, ott_id) %>%
      dplyr::mutate(Infraorder = purrr::map_chr(ott_id, get_infraorder)) %>%
      dplyr::select(Genus = Final_Name, ott_id, Infraorder),
    by = "Genus"
  ) %>%
  dplyr::mutate(Infraorder = dplyr::if_else(is.na(Infraorder), "Unresolved", Infraorder),
                Genus = factor(Genus, levels = phylo_tree_plot$tip.label))

lineage_df %>%
  dplyr::count(Infraorder, sort = TRUE)

lineage_df %>%
  dplyr::arrange(Infraorder, Genus) %>%
  print(n = Inf)


## Prepare infraorder ring annotation data. 
lineage_ring <- lineage_df %>%
  dplyr::mutate(Genus = as.character(Genus),
                Infraorder = factor(Infraorder, levels = names(infraorder_cols)),
                value = 1) %>%
  dplyr::filter(Genus %in% phylo_tree_plot$tip.label)

setdiff(lineage_ring$Genus, phylo_tree_plot$tip.label)
setdiff(phylo_tree_plot$tip.label, lineage_ring$Genus)
setdiff(levels(lineage_ring$Infraorder), names(infraorder_cols))

# Structure of the phylogenetic tree and outer-ring annotations. 
# Inner tree: genus-level phylogenetic relationships, using phylo_tree_plot. 
# Ring 1: infraorder colour band, using lineage_ring. 
# Ring 2: genus labels.
# Ring 3: abundance-status markers, using abundance_ring.
# Ring 4: landscape-level relative abundance heatmap, using abund_type. 
# Ring 5: occurrence-frequency bar plot, using occ_freq. 
# Ring 6: total-abundance bar plot, using total_abund.




#### 2.2.4 Circular phylogenetic tree plotting ####

if (!exists("is.waive")) is.waive <- function(x) inherits(x, "waiver")

library(ggplot2)
library(ggtreeExtra)
library(ggnewscale)

## Inner phylogenetic tree. 
P_tree_inner <- ggtree(phylo_tree_plot, layout = "fan", open.angle = 90) +
  geom_tree(linewidth = 0.25, color = "grey20") +
  geom_aline(linetype = "longdash", color = "grey85", linewidth = 0.1, show.legend = FALSE) +
  theme(plot.margin = margin(5, 5, 5, 5)); p_tree_inner


## Ring 1: infraorder colour band. 
p_tree_infraorder <- p_tree_inner +
  geom_fruit(data = lineage_ring %>% dplyr::mutate(Genus = as.character(Genus)),
             geom = geom_bar, stat = "identity", width = 1,
             mapping = aes(y = Genus, x = value, fill = Infraorder),
             offset = 0.03, pwidth = 1.45) +
  scale_fill_manual(name = "Infraorder", values = infraorder_cols) +
  new_scale_fill(); p_tree_infraorder

# Ring 2: genus labels. 
p_tree_infraorder_label <- p_tree_infraorder +
  geom_tiplab2(aes(label = label), offset = 0.3, size = 2.5, color = "grey15",
               align = FALSE, show.legend = FALSE); p_tree_infraorder_label


## Ring 3: abundance status. 
p_tree_abundance_status <- p_tree_infraorder_label +
  geom_fruit(data = abundance_ring, geom = geom_bar, stat = "identity", width = 0.6,
             mapping = aes(y = Genus, x = value, fill = Abundance_status, color = Abundance_status),
             pwidth = 0.12, offset = 0.08) +
  scale_fill_manual(name = "Abundance status",
                    values = c("Dominant" = "#1B7837", "Common" = "#A6DBA0", "Rare" = "white")) +
  scale_color_manual(guide = "none",
                     values = c("Dominant" = "#1B7837", "Common" = "#A6DBA0", "Rare" = "#A6DBA0")) +
  new_scale_fill(); p_tree_abundance_status


## Ring 4: landscape-level relative abundance heatmap. 
p_tree_abund_type <- p_tree_abundance_status +
  geom_fruit(data = abund_type %>%
               dplyr::mutate(Genus = as.character(Genus),
                             Type = factor(Type, levels = c("Spr_TL", "Cas_TL", "Cus_TL", "Rap_TL", "Rim_TL"))),
             geom = geom_tile, mapping = aes(y = Genus, x = Type, fill = Landscape_relative_abundance),
             offset = 0.10, pwidth = 0.85, color = "grey100", linewidth = 0.5) +
  scale_fill_gradient(name = "Mean relative abundance", low = "#F6F5F6", high = "#803F8F") +
  new_scale_fill(); p_tree_abund_type


# Ring 5: occurrence-frequency bar plot. 
p_tree_occ_freq <- p_tree_abund_type +
  geom_fruit(data = occ_freq %>% dplyr::mutate(Genus = as.character(Genus)),
             geom = geom_bar, stat = "identity", width = 0.6,
             mapping = aes(y = Genus, x = Occurrence_frequency, fill = occ_group),
             pwidth = 0.6, offset = 0.15) +
  scale_fill_manual(name = "Occurrence frequency", values = c("Occurrence" = "#8DD3C7")) +
  new_scale_fill(); p_tree_occ_freq


## Ring 6: total-abundance bar plot. 
p_tree_total_abund <- p_tree_occ_freq +
  geom_fruit(data = total_abund %>% dplyr::mutate(Genus = as.character(Genus)),
             geom = geom_bar, stat = "identity", width = 0.6,
             mapping = aes(y = Genus, x = Total_abundance, fill = abund_group),
             pwidth = 0.6, offset = 0.15) +
  scale_fill_manual(name = "Total abundance", values = c("Total abundance" = "#FDB462")) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "right", legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)); p_tree_total_abund

# ggsave(filename = "01 Phylogenetic tree.pdf", plot = p_tree_total_abund, width = 13, height = 13, units = "in", device = cairo_pdf)






#### 3. Multidimensional β-diversity decomposition based on presence-absence data ####

library(BAT)

#### 3.1 Taxonomic β-diversity decomposition: TD  ####

## Build genus-level presence-absence community matrix. 
comm_pa <- comm_abund
comm_pa[comm_pa > 0] <- 1
comm_pa <- as.data.frame(comm_pa)

#table(as.matrix(comm_pa)); rowSums(comm_pa); colSums(comm_pa); dim(comm_pa); head(comm_pa[, 1:6])


## Pairwise taxonomic β-diversity and its replacement and richness-difference components. 
td_pair_pa <- BAT::beta(comm = as.matrix(comm_pa), func = "Jaccard", abund = FALSE, comp = TRUE)

td_total_pa <- as.matrix(td_pair_pa$Btotal)
td_repl_pa  <- as.matrix(td_pair_pa$Brepl)
td_rich_pa  <- as.matrix(td_pair_pa$Brich)


## Multiple-site taxonomic β-diversity decomposition. 
td_multi_pa <- BAT::beta.multi(comm = as.matrix(comm_pa), func = "Jaccard", abund = FALSE); td_multi_pa

td_repl_ratio_pa <- td_multi_pa["Brepl", "Average"] / td_multi_pa["Btotal", "Average"]; td_repl_ratio_pa
td_rich_ratio_pa <- td_multi_pa["Brich", "Average"] / td_multi_pa["Btotal", "Average"]; td_rich_ratio_pa




#### 3.2 Functional β-diversity decomposition: FD ####

## Pairwise functional β-diversity and its replacement and richness-difference components. 
fd_pair_pa <- BAT::beta(comm = as.matrix(comm_pa_fd), tree = fd_tree, func = "Jaccard", abund = FALSE, comp = TRUE)

fd_total_pa <- as.matrix(fd_pair_pa$Btotal)
fd_repl_pa  <- as.matrix(fd_pair_pa$Brepl)
fd_rich_pa  <- as.matrix(fd_pair_pa$Brich)


## Multiple-site functional β-diversity decomposition. 
fd_multi_pa <- BAT::beta.multi(comm = as.matrix(comm_pa_fd), tree = fd_tree, func = "Jaccard", abund = FALSE); fd_multi_pa

fd_repl_ratio_pa <- fd_multi_pa["Brepl", "Average"] / fd_multi_pa["Btotal", "Average"]; fd_repl_ratio_pa
fd_rich_ratio_pa <- fd_multi_pa["Brich", "Average"] / fd_multi_pa["Btotal", "Average"]; fd_rich_ratio_pa





#### 3.3 Phylogenetic β-diversity decomposition: PD ####

## Use the phylogenetic tree and presence-absence matrix aligned in Chapter 2. 
phylo_tree_pd <- phylo_tree_plot
comm_pd_pa <- comm_pa_phylo

identical(colnames(comm_pd_pa), phylo_tree_pd$tip.label)
dim(comm_pd_pa); table(as.matrix(comm_pd_pa)); rowSums(comm_pd_pa); colSums(comm_pd_pa)


## Pairwise phylogenetic β-diversity and its replacement and richness-difference components. 
pd_pair_pa <- BAT::beta(comm = as.matrix(comm_pd_pa), tree = phylo_tree_pd, func = "Jaccard", abund = FALSE, comp = TRUE)

pd_total_pa <- as.matrix(pd_pair_pa$Btotal)
pd_repl_pa  <- as.matrix(pd_pair_pa$Brepl)
pd_rich_pa  <- as.matrix(pd_pair_pa$Brich)


## Multiple-site phylogenetic β-diversity decomposition. 
pd_multi_pa <- BAT::beta.multi(comm = as.matrix(comm_pd_pa), tree = phylo_tree_pd, func = "Jaccard", abund = FALSE); pd_multi_pa

pd_repl_ratio_pa <- pd_multi_pa["Brepl", "Average"] / pd_multi_pa["Btotal", "Average"]; pd_repl_ratio_pa
pd_rich_ratio_pa <- pd_multi_pa["Brich", "Average"] / pd_multi_pa["Btotal", "Average"]; pd_rich_ratio_pa




#### 3.4 Organizing β-diversity decomposition results ####

## Pairwise β-diversity matrices.

beta_list <- list(
  TD = list(total = td_total_pa, repl = td_repl_pa, rich = td_rich_pa),
  FD = list(total = fd_total_pa, repl = fd_repl_pa, rich = fd_rich_pa),
  PD = list(total = pd_total_pa, repl = pd_repl_pa, rich = pd_rich_pa)
)


## Pairwise β-diversity data frame.

extract_upper_triangle <- function(dist_obj, sample_info, value_name = "distance") {
  dist_mat <- as.matrix(dist_obj)
  site_group <- setNames(sample_info$Type, sample_info$SampleID)
  if (!all(rownames(dist_mat) %in% names(site_group))) stop("Some SampleIDs in dist_obj are not found in sample_info.")
  idx <- which(upper.tri(dist_mat), arr.ind = TRUE)
  df <- data.frame(sample1 = rownames(dist_mat)[idx[, 1]], sample2 = colnames(dist_mat)[idx[, 2]], stringsAsFactors = FALSE)
  df$group1 <- site_group[df$sample1]
  df$group2 <- site_group[df$sample2]
  df$pair_type <- ifelse(df$group1 == df$group2, "within", "between")
  df$group_pair <- apply(df[, c("group1", "group2")], 1, function(x) paste(sort(x), collapse = " vs "))
  df[[value_name]] <- dist_mat[idx]
  rownames(df) <- NULL
  return(df)
}

site_group <- setNames(sample_info$Type, sample_info$SampleID); site_group

pa_dfs <- list(
  extract_upper_triangle(beta_list$TD$total, sample_info, "TD_total"),
  extract_upper_triangle(beta_list$TD$repl,  sample_info, "TD_repl"),
  extract_upper_triangle(beta_list$TD$rich,  sample_info, "TD_rich"),
  extract_upper_triangle(beta_list$FD$total, sample_info, "FD_total"),
  extract_upper_triangle(beta_list$FD$repl,  sample_info, "FD_repl"),
  extract_upper_triangle(beta_list$FD$rich,  sample_info, "FD_rich"),
  extract_upper_triangle(beta_list$PD$total, sample_info, "PD_total"),
  extract_upper_triangle(beta_list$PD$repl,  sample_info, "PD_repl"),
  extract_upper_triangle(beta_list$PD$rich,  sample_info, "PD_rich")
)

beta_df_pa <- Reduce(function(x, y) {
  dplyr::full_join(x, y, by = c("sample1", "sample2", "group1", "group2", "pair_type", "group_pair"))
}, pa_dfs)

dim(beta_df_pa); head(beta_df_pa); summary(beta_df_pa)
colSums(is.na(beta_df_pa)); table(beta_df_pa$group_pair); table(beta_df_pa$pair_type)

# write.csv(beta_df_pa, file = "02 beta_pairwise_presence_absence.csv", row.names = FALSE)


## Multiple-site β-diversity summary

beta_multi_df <- dplyr::bind_rows(
  data.frame(Data_type = "Presence-absence", Dimension = "TD", Btotal = td_multi_pa["Btotal", "Average"],
             Replacement = td_multi_pa["Brepl", "Average"], Richness_difference = td_multi_pa["Brich", "Average"]),
  data.frame(Data_type = "Presence-absence", Dimension = "FD", Btotal = fd_multi_pa["Btotal", "Average"],
             Replacement = fd_multi_pa["Brepl", "Average"], Richness_difference = fd_multi_pa["Brich", "Average"]),
  data.frame(Data_type = "Presence-absence", Dimension = "PD", Btotal = pd_multi_pa["Btotal", "Average"],
             Replacement = pd_multi_pa["Brepl", "Average"], Richness_difference = pd_multi_pa["Brich", "Average"])
) %>%
  dplyr::mutate(Data_type = factor(Data_type, levels = "Presence-absence"),
                Dimension = factor(Dimension, levels = c("TD", "FD", "PD")),
                Replacement_ratio = Replacement / Btotal,
                Richness_ratio = Richness_difference / Btotal,
                Replacement_percent = Replacement_ratio * 100,
                Richness_percent = Richness_ratio * 100)

head(beta_multi_df); dim(beta_multi_df)
colSums(is.na(beta_multi_df)); table(beta_multi_df$Data_type); table(beta_multi_df$Dimension)

# write.csv(beta_multi_df, "03 beta_multi_decomposition_summary_presence_absence.csv", row.names = FALSE)




#### 3.5 Visualization of multidimensional β-diversity decomposition ####

#### 3.5.1 Ternary plots of pairwise β-diversity components ####

library(ggtern)
library(cowplot)
library(grid)

plot_beta_ternary <- function(data, total_col, repl_col, rich_col, point_col = "grey80",
                              density_low = "#E7D4E8", density_high = "#762A83",
                              x_lab = "Similarity", y_lab = "Replacement", z_lab = "Richness difference") {
  tern_df <- data %>%
    dplyr::transmute(Similarity = 1 - .data[[total_col]],
                     Replacement = .data[[repl_col]],
                     Richness = .data[[rich_col]])
  
  tern_df_density <- tern_df %>%
    dplyr::mutate(Similarity = pmax(Similarity, 0.001),
                  Replacement = pmax(Replacement, 0.001),
                  Richness = pmax(Richness, 0.001))
  
  centroid_df <- tern_df %>%
    dplyr::summarise(Similarity = mean(Similarity, na.rm = TRUE),
                     Replacement = mean(Replacement, na.rm = TRUE),
                     Richness = mean(Richness, na.rm = TRUE))
  
  arrow_df <- tibble::tibble(axis = c("Similarity", "Replacement", "Richness"),
                             x = centroid_df$Similarity,
                             y = centroid_df$Replacement,
                             z = centroid_df$Richness,
                             xend = c(centroid_df$Similarity, 0, 1 - centroid_df$Richness),
                             yend = c(1 - centroid_df$Similarity, centroid_df$Replacement, 0),
                             zend = c(0, 1 - centroid_df$Replacement, centroid_df$Richness))
  
  axis_cols <- c("Similarity" = "grey40", "Replacement" = "#FDB462", "Richness" = "#8DD3C7")
  
  p <- ggtern(data = tern_df, aes(x = Similarity, y = Replacement, z = Richness)) +
    geom_point(colour = point_col, alpha = 0.6, size = 2, shape = 16) +
    stat_density_tern(data = tern_df_density,
                      aes(x = Similarity, y = Replacement, z = Richness,
                          fill = after_stat(level), alpha = after_stat(level)),
                      geom = "polygon", bdl = 0.001, na.rm = TRUE) +
    scale_fill_gradient(low = density_low, high = density_high, guide = "none") +
    scale_alpha_continuous(range = c(0.2, 0.8), guide = "none") +
    geom_segment(data = arrow_df,
                 aes(x = x, y = y, z = z, xend = xend, yend = yend, zend = zend, colour = axis),
                 linewidth = 0.27, arrow = grid::arrow(type = "closed", length = grid::unit(0.18, "cm"), angle = 15),
                 show.legend = FALSE) +
    geom_point(data = centroid_df, aes(x = Similarity, y = Replacement, z = Richness),
               inherit.aes = FALSE, colour = "#E64B35", size = 3.2, shape = 16) +
    scale_colour_manual(values = axis_cols, guide = "none") +
    labs(x = x_lab, y = y_lab, z = z_lab) +
    theme_bw(base_size = 10) +
    theme(tern.plot.background = element_rect(fill = "white", colour = NA),
          tern.panel.background = element_rect(fill = "white", colour = NA),
          tern.panel.grid.major.T = element_line(colour = scales::alpha(axis_cols["Replacement"], 0.6), linewidth = 0.27, linetype = 2),
          tern.panel.grid.major.L = element_line(colour = scales::alpha(axis_cols["Similarity"], 0.6), linewidth = 0.27, linetype = 2),
          tern.panel.grid.major.R = element_line(colour = scales::alpha(axis_cols["Richness"], 0.6), linewidth = 0.27, linetype = 2),
          tern.panel.grid.minor.T = element_blank(),
          tern.panel.grid.minor.L = element_blank(),
          tern.panel.grid.minor.R = element_blank(),
          tern.axis.line.T = element_line(colour = axis_cols["Replacement"], linewidth = 0.27),
          tern.axis.line.L = element_line(colour = axis_cols["Similarity"], linewidth = 0.27),
          tern.axis.line.R = element_line(colour = axis_cols["Richness"], linewidth = 0.27),
          tern.axis.arrow.show = TRUE,
          tern.axis.arrow.T = element_line(colour = axis_cols["Replacement"], linewidth = 0.27),
          tern.axis.arrow.L = element_line(colour = axis_cols["Similarity"], linewidth = 0.27),
          tern.axis.arrow.R = element_line(colour = axis_cols["Richness"], linewidth = 0.27),
          tern.axis.ticks.major.T = element_line(colour = axis_cols["Replacement"], linewidth = 0.27),
          tern.axis.ticks.major.L = element_line(colour = axis_cols["Similarity"], linewidth = 0.27),
          tern.axis.ticks.major.R = element_line(colour = axis_cols["Richness"], linewidth = 0.27),
          tern.axis.ticks.minor.T = element_blank(),
          tern.axis.ticks.minor.L = element_blank(),
          tern.axis.ticks.minor.R = element_blank(),
          tern.axis.ticks.length.major = ggplot2::rel(0.1),
          tern.axis.ticks.length.minor = ggplot2::rel(0.1),
          tern.axis.title.T = element_blank(),
          tern.axis.title.L = element_blank(),
          tern.axis.title.R = element_blank(),
          tern.axis.text.T = element_text(colour = axis_cols["Replacement"], size = 10),
          tern.axis.text.L = element_text(colour = axis_cols["Similarity"], size = 10),
          tern.axis.text.R = element_text(colour = axis_cols["Richness"], size = 10),
          legend.position = "none")
  
  return(list(plot = p, data = tern_df, centroid = centroid_df))
}


p_td_pa <- plot_beta_ternary(data = beta_df_pa,
                             total_col = "TD_total", repl_col = "TD_repl", rich_col = "TD_rich",
                             point_col = "grey80", density_low = "#E7D4E8", density_high = "#762A83",
                             x_lab = "Similarity (1 - TD total)",
                             y_lab = "TD replacement",
                             z_lab = "TD richness difference"); p_td_pa$plot

p_fd_pa <- plot_beta_ternary(data = beta_df_pa,
                             total_col = "FD_total", repl_col = "FD_repl", rich_col = "FD_rich",
                             point_col = "grey80", density_low = "#D9F0D3", density_high = "#1B7837",
                             x_lab = "Similarity (1 - FD total)",
                             y_lab = "FD replacement",
                             z_lab = "FD richness difference"); p_fd_pa$plot

p_pd_pa <- plot_beta_ternary(data = beta_df_pa,
                             total_col = "PD_total", repl_col = "PD_repl", rich_col = "PD_rich",
                             point_col = "grey80", density_low = "#E0F3F8", density_high = "#3B55A4",
                             x_lab = "Similarity (1 - PD total)",
                             y_lab = "PD replacement",
                             z_lab = "PD richness difference"); p_pd_pa$plot

g_td_pa <- ggplotGrob(p_td_pa$plot)
g_fd_pa <- ggplotGrob(p_fd_pa$plot)
g_pd_pa <- ggplotGrob(p_pd_pa$plot)

p_beta_ternary_pa <- cowplot::plot_grid(g_td_pa, g_fd_pa, g_pd_pa,
                                        ncol = 3, nrow = 1, align = "none"); p_beta_ternary_pa

# ggsave(filename = "03 Beta ternary_presence_absence.pdf", plot = p_beta_ternary_pa, width = 18.5, height = 6.5, units = "cm", device = cairo_pdf)




#### 3.5.2 Stacked bar plot of multiple-site β-diversity components ####

beta_multi_plot_long <- beta_multi_df %>%
  dplyr::filter(Data_type == "Presence-absence") %>%
  dplyr::select(Data_type, Dimension, Replacement, Richness_difference, Replacement_ratio, Richness_ratio) %>%
  tidyr::pivot_longer(cols = c(Replacement, Richness_difference),
                      names_to = "Component", values_to = "Value") %>%
  dplyr::mutate(Ratio = dplyr::case_when(Component == "Replacement" ~ Replacement_ratio,
                                         Component == "Richness_difference" ~ Richness_ratio),
                Component = dplyr::recode(Component, "Richness_difference" = "Richness difference"),
                Component = factor(Component, levels = c("Replacement", "Richness difference")),
                Dimension = factor(Dimension, levels = c("PD", "FD", "TD")),
                label = paste0(round(Value, 3), "\n(", sprintf("%.1f%%", Ratio * 100), ")")); beta_multi_plot_long


p_beta_stack_pa <- ggplot(beta_multi_plot_long, aes(x = Dimension, y = Value, fill = Component)) +
  geom_col(width = 0.8, color = "grey100", linewidth = 0.27) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            hjust = 0.5, size = 3.2, color = "grey100", lineheight = 1.1) +
  scale_fill_manual(values = c("Replacement" = "#80B1D3", "Richness difference" = "#BC80BD")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.03))) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  coord_flip() +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.27),
        axis.text.x = element_text(colour = "black", size = 10),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.line = element_line(linewidth = 0.27, colour = "black"),
        axis.ticks = element_line(linewidth = 0.27, colour = "black"),
        legend.position = "top", legend.text = element_text(size = 10),
        legend.key.size = grid::unit(0.45, "cm")); p_beta_stack_pa

# ggsave(filename = "03 Beta stacked_components_presence_absence.pdf", plot = p_beta_stack_pa, width = 9, height = 6.5, units = "cm", device = cairo_pdf)






#### 4. Relationships among taxonomic, functional and phylogenetic β diversity ####

#### 4.1 Mantel and partial Mantel tests among β-diversity dimensions ####

library(vegan)

run_mantel_pa <- function(beta_list, permutations = 9999) {
  set.seed(123)
  
  beta_info <- data.frame(beta_key = c("total", "repl", "rich"),
                          Beta_component = c("Total beta diversity", "Replacement component", "Richness-difference component"))
  dimension_pairs <- list(c("TD", "FD"), c("TD", "PD"), c("FD", "PD"))
  mantel_list <- list(); partial_list <- list(); k1 <- 1; k2 <- 1
  
  for (j in seq_len(nrow(beta_info))) {
    for (pair in dimension_pairs) {
      beta_key <- beta_info$beta_key[j]; x_name <- pair[1]; y_name <- pair[2]
      z_name <- setdiff(c("TD", "FD", "PD"), pair)
      
      m <- vegan::mantel(as.dist(beta_list[[x_name]][[beta_key]]),
                         as.dist(beta_list[[y_name]][[beta_key]]),
                         method = "spearman", permutations = permutations)
      pm <- vegan::mantel.partial(as.dist(beta_list[[x_name]][[beta_key]]),
                                  as.dist(beta_list[[y_name]][[beta_key]]),
                                  as.dist(beta_list[[z_name]][[beta_key]]),
                                  method = "spearman", permutations = permutations)
      
      mantel_list[[k1]] <- data.frame(Test = "Mantel", Data_type = "Presence-absence",
                                      Beta_component = beta_info$Beta_component[j],
                                      Comparison = paste0(x_name, " vs ", y_name),
                                      Mantel_r = as.numeric(m$statistic), P_value = m$signif)
      partial_list[[k2]] <- data.frame(Test = "Partial Mantel", Data_type = "Presence-absence",
                                       Beta_component = beta_info$Beta_component[j],
                                       Comparison = paste0(x_name, " vs ", y_name),
                                       Mantel_r = as.numeric(pm$statistic), P_value = pm$signif)
      k1 <- k1 + 1; k2 <- k2 + 1
    }
  }
  
  mantel_df <- dplyr::bind_rows(mantel_list)
  partial_df <- dplyr::bind_rows(partial_list)
  
  all_df <- dplyr::bind_rows(mantel_df, partial_df) %>%
    dplyr::mutate(P_label = dplyr::case_when(P_value <= 0.001 ~ "< 0.001", P_value < 0.01 ~ "< 0.01",
                                             P_value < 0.05 ~ "< 0.05", TRUE ~ "ns"),
                  Significance = dplyr::case_when(P_value <= 0.001 ~ "***", P_value < 0.01 ~ "**",
                                                  P_value < 0.05 ~ "*", TRUE ~ "ns"),
                  Data_type = factor(Data_type, levels = "Presence-absence"),
                  Beta_component = factor(Beta_component, levels = c("Total beta diversity", "Replacement component", "Richness-difference component")),
                  Comparison = factor(Comparison, levels = c("TD vs FD", "TD vs PD", "FD vs PD")),
                  Test = factor(Test, levels = c("Mantel", "Partial Mantel"))) %>%
    dplyr::arrange(Beta_component, Comparison, Test)
  
  return(list(Mantel = mantel_df, Partial_Mantel = partial_df, All_results = all_df))
}

mantel_results <- run_mantel_pa(beta_list = beta_list, permutations = 9999); mantel_results

mantel_results_df <- mantel_results$Mantel; mantel_results_df
partial_mantel_results_df <- mantel_results$Partial_Mantel; partial_mantel_results_df
mantel_all_results_df <- mantel_results$All_results; mantel_all_results_df

# write.csv(mantel_all_results_df, "05 Mantel and Partial Mantel test results_presence_absence.csv", row.names = FALSE)




#### 4.2 Linear mixed models for pairwise β-diversity relationships ####

library(lme4)
library(car)
library(performance)

plot_beta_relation <- function(data, xvar, yvar, xlab = NULL, ylab = NULL, add_interaction = TRUE) {
  df <- data %>%
    dplyr::select(sample1, sample2, pair_type, group_pair, dplyr::all_of(c(xvar, yvar))) %>%
    dplyr::rename(x = dplyr::all_of(xvar), y = dplyr::all_of(yvar)) %>%
    dplyr::mutate(pair_type = factor(pair_type, levels = c("within", "between")))
  
  form <- if (add_interaction) {
    y ~ x * pair_type + (1 | group_pair)
  } else {
    y ~ x + pair_type + (1 | group_pair)
  }
  
  fit <- lme4::lmer(form, data = df)
  
  print(summary(fit)); print(car::Anova(fit, test.statistic = "F")); print(performance::r2(fit))
  print(performance::check_normality(fit)); print(performance::check_heteroscedasticity(fit))
  
  pred_df <- df %>%
    dplyr::group_by(pair_type) %>%
    dplyr::summarise(x_min = min(x, na.rm = TRUE), x_max = max(x, na.rm = TRUE), .groups = "drop") %>%
    dplyr::rowwise() %>%
    dplyr::do(data.frame(pair_type = .$pair_type, x = seq(.$x_min, .$x_max, length.out = 100))) %>%
    dplyr::ungroup()
  
  mm <- model.matrix(if (add_interaction) ~ x * pair_type else ~ x + pair_type, data = pred_df)
  beta <- lme4::fixef(fit); vcov_mat <- as.matrix(vcov(fit))
  pred_df$predicted <- as.numeric(mm %*% beta)
  pred_se <- sqrt(diag(mm %*% vcov_mat %*% t(mm)))
  pred_df$conf.low <- pred_df$predicted - 1.96 * pred_se
  pred_df$conf.high <- pred_df$predicted + 1.96 * pred_se
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df, ggplot2::aes(x = x, y = y, colour = pair_type),
                        shape = 21, fill = NA, size = 2, stroke = 0.5, alpha = 0.85) +
    ggplot2::geom_ribbon(data = pred_df, ggplot2::aes(x = x, ymin = conf.low, ymax = conf.high, group = pair_type, fill = pair_type),
                         alpha = 0.20, color = NA) +
    ggplot2::geom_line(data = pred_df, ggplot2::aes(x = x, y = predicted, group = pair_type, colour = pair_type), linewidth = 0.75) +
    ggplot2::scale_color_manual(values = c("within" = "#1E88E5", "between" = "#E53935")) +
    ggplot2::scale_fill_manual(values = c("within" = "#1E88E5", "between" = "#E53935")) +
    ggplot2::labs(x = if (is.null(xlab)) xvar else xlab, y = if (is.null(ylab)) yvar else ylab) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.27),
                   axis.title = ggplot2::element_text(size = 10),
                   axis.text = ggplot2::element_text(colour = "black", size = 10),
                   axis.line = ggplot2::element_line(linewidth = 0.27, colour = "black"),
                   axis.ticks = ggplot2::element_line(linewidth = 0.27, colour = "black"),
                   legend.position = "none")
  
  return(list(model = fit, plot = p, pred = pred_df, data = df))
}


beta_relation_specs <- data.frame(
  Object = c("res_pa_td_fd_total", "res_pa_td_pd_total", "res_pa_fd_pd_total",
             "res_pa_td_fd_repl",  "res_pa_td_pd_repl",  "res_pa_fd_pd_repl",
             "res_pa_td_fd_rich",  "res_pa_td_pd_rich",  "res_pa_fd_pd_rich"),
  Beta_component = c(rep("Total beta diversity", 3), rep("Replacement component", 3), rep("Richness-difference component", 3)),
  Comparison = rep(c("TD vs FD", "TD vs PD", "FD vs PD"), 3),
  xvar = c("TD_total", "TD_total", "FD_total", "TD_repl", "TD_repl", "FD_repl", "TD_rich", "TD_rich", "FD_rich"),
  yvar = c("FD_total", "PD_total", "PD_total", "FD_repl", "PD_repl", "PD_repl", "FD_rich", "PD_rich", "PD_rich"),
  xlab = c("Taxonomic beta diversity", "Taxonomic beta diversity", "Functional beta diversity",
           "Taxonomic replacement", "Taxonomic replacement", "Functional replacement",
           "Taxonomic richness difference", "Taxonomic richness difference", "Functional richness difference"),
  ylab = c("Functional beta diversity", "Phylogenetic beta diversity", "Phylogenetic beta diversity",
           "Functional replacement", "Phylogenetic replacement", "Phylogenetic replacement",
           "Functional richness difference", "Phylogenetic richness difference", "Phylogenetic richness difference"),
  stringsAsFactors = FALSE
); beta_relation_specs

beta_relation_pa <- setNames(vector("list", nrow(beta_relation_specs)), beta_relation_specs$Object)

for (i in seq_len(nrow(beta_relation_specs))) {
  beta_relation_pa[[i]] <- plot_beta_relation(data = beta_df_pa,
                                              xvar = beta_relation_specs$xvar[i],
                                              yvar = beta_relation_specs$yvar[i],
                                              xlab = beta_relation_specs$xlab[i],
                                              ylab = beta_relation_specs$ylab[i],
                                              add_interaction = TRUE)
}

list2env(beta_relation_pa, envir = .GlobalEnv)

beta_relation_pa$res_pa_td_fd_total; beta_relation_pa$res_pa_td_pd_total

p_beta_relation_pa_all <- cowplot::plot_grid(
  res_pa_td_fd_total$plot, res_pa_td_pd_total$plot, res_pa_fd_pd_total$plot,
  res_pa_td_fd_repl$plot,  res_pa_td_pd_repl$plot,  res_pa_fd_pd_repl$plot,
  res_pa_td_fd_rich$plot,  res_pa_td_pd_rich$plot,  res_pa_fd_pd_rich$plot,
  ncol = 3, nrow = 3, align = "hv"); p_beta_relation_pa_all

# ggsave(filename = "04 Beta relation_presence_absence.pdf", plot = p_beta_relation_pa_all, width = 18, height = 17.5, units = "cm", device = cairo_pdf)




#### 4.3 Extracting linear mixed model results ####

extract_beta_lmm_pa <- function(beta_relation_pa, beta_relation_specs) {
  extract_one_lmm <- function(res_object, beta_component, comparison) {
    fit <- res_object$model
    anova_df <- as.data.frame(car::Anova(fit, test.statistic = "F"))
    anova_df$Term <- rownames(anova_df)
    f_col <- grep("^F", names(anova_df), value = TRUE)[1]
    p_col <- grep("Pr", names(anova_df), value = TRUE)[1]
    get_val <- function(term, col_name) if (term %in% anova_df$Term) as.numeric(anova_df[anova_df$Term == term, col_name][1]) else NA_real_
    
    beta <- lme4::fixef(fit)
    slope_within <- as.numeric(beta["x"])
    slope_between <- if ("x:pair_typebetween" %in% names(beta)) as.numeric(beta["x"] + beta["x:pair_typebetween"]) else NA_real_
    r2_res <- suppressWarnings(performance::r2_nakagawa(fit))
    
    data.frame(Data_type = "Presence-absence", Beta_component = beta_component, Comparison = comparison,
               Slope_within = slope_within, Slope_between = slope_between,
               F_x = get_val("x", f_col), P_x = get_val("x", p_col),
               F_pair_type = get_val("pair_type", f_col), P_pair_type = get_val("pair_type", p_col),
               F_interaction = get_val("x:pair_type", f_col), P_interaction = get_val("x:pair_type", p_col),
               R2_marginal = as.numeric(r2_res$R2_marginal),
               R2_conditional = as.numeric(r2_res$R2_conditional))
  }
  
  beta_lmm_results_df <- dplyr::bind_rows(lapply(seq_len(nrow(beta_relation_specs)), function(i) {
    extract_one_lmm(beta_relation_pa[[beta_relation_specs$Object[i]]],
                    beta_relation_specs$Beta_component[i],
                    beta_relation_specs$Comparison[i])
  })) %>%
    dplyr::mutate(P_x_label = dplyr::case_when(P_x <= 0.001 ~ "< 0.001", P_x < 0.01 ~ "< 0.01", P_x < 0.05 ~ "< 0.05", TRUE ~ "ns"),
                  Significance_x = dplyr::case_when(P_x <= 0.001 ~ "***", P_x < 0.01 ~ "**", P_x < 0.05 ~ "*", TRUE ~ "ns"),
                  P_pair_type_label = dplyr::case_when(P_pair_type <= 0.001 ~ "< 0.001", P_pair_type < 0.01 ~ "< 0.01", P_pair_type < 0.05 ~ "< 0.05", TRUE ~ "ns"),
                  Significance_pair_type = dplyr::case_when(P_pair_type <= 0.001 ~ "***", P_pair_type < 0.01 ~ "**", P_pair_type < 0.05 ~ "*", TRUE ~ "ns"),
                  P_interaction_label = dplyr::case_when(P_interaction <= 0.001 ~ "< 0.001", P_interaction < 0.01 ~ "< 0.01", P_interaction < 0.05 ~ "< 0.05", TRUE ~ "ns"),
                  Significance_interaction = dplyr::case_when(P_interaction <= 0.001 ~ "***", P_interaction < 0.01 ~ "**", P_interaction < 0.05 ~ "*", TRUE ~ "ns"),
                  Data_type = factor(Data_type, levels = "Presence-absence"),
                  Beta_component = factor(Beta_component, levels = c("Total beta diversity", "Replacement component", "Richness-difference component")),
                  Comparison = factor(Comparison, levels = c("TD vs FD", "TD vs PD", "FD vs PD"))) %>%
    dplyr::select(Data_type, Beta_component, Comparison, Slope_within, Slope_between,
                  F_x, P_x, P_x_label, Significance_x,
                  F_pair_type, P_pair_type, P_pair_type_label, Significance_pair_type,
                  F_interaction, P_interaction, P_interaction_label, Significance_interaction,
                  R2_marginal, R2_conditional) %>%
    dplyr::arrange(Beta_component, Comparison)
  
  return(beta_lmm_results_df)
}


beta_lmm_results_df <- extract_beta_lmm_pa(beta_relation_pa = beta_relation_pa, beta_relation_specs = beta_relation_specs); beta_lmm_results_df

# write.csv(beta_lmm_results_df, "06 Beta LMM results_presence_absence.csv", row.names = FALSE)






#### 5. Effects of landscape type on multidimensional β diversity ####

#### 5.1 Custom functions for landscape-effect tests ####

pairwise_adonis <- function(dist_mat, group, nperm = 9999) {
  group <- as.factor(group); lv <- levels(group); dist_mat <- as.matrix(dist_mat); res <- list()
  for (i in 1:(length(lv) - 1)) for (j in (i + 1):length(lv)) {
    idx <- which(group %in% c(lv[i], lv[j]))
    fit <- vegan::adonis2(as.dist(dist_mat[idx, idx]) ~ droplevels(group[idx]), permutations = nperm)
    res[[paste(lv[i], lv[j], sep = "_vs_")]] <- data.frame(group1 = lv[i], group2 = lv[j],
                                                           F = fit$F[1], R2 = fit$R2[1],
                                                           p = fit$`Pr(>F)`[1])
  }
  out <- dplyr::bind_rows(res); out$p_adj <- p.adjust(out$p, method = "BH"); out
}


pairwise_anosim <- function(dist_mat, group, nperm = 9999) {
  group <- as.factor(group); lv <- levels(group); dist_mat <- as.matrix(dist_mat); res <- list()
  for (i in 1:(length(lv) - 1)) for (j in (i + 1):length(lv)) {
    idx <- which(group %in% c(lv[i], lv[j]))
    fit <- vegan::anosim(as.dist(dist_mat[idx, idx]), grouping = droplevels(group[idx]), permutations = nperm)
    res[[paste(lv[i], lv[j], sep = "_vs_")]] <- data.frame(group1 = lv[i], group2 = lv[j],
                                                           R = unname(fit$statistic), p = fit$signif)
  }
  out <- dplyr::bind_rows(res); out$p_adj <- p.adjust(out$p, method = "BH"); out
}

pairwise_mrpp <- function(dist_mat, group, nperm = 9999) {
  group <- as.factor(group); lv <- levels(group); dist_mat <- as.matrix(dist_mat); res <- list()
  for (i in 1:(length(lv) - 1)) for (j in (i + 1):length(lv)) {
    idx <- which(group %in% c(lv[i], lv[j]))
    fit <- vegan::mrpp(as.dist(dist_mat[idx, idx]), grouping = droplevels(group[idx]), permutations = nperm)
    res[[paste(lv[i], lv[j], sep = "_vs_")]] <- data.frame(group1 = lv[i], group2 = lv[j],
                                                           delta = fit$delta, A = fit$A, p = fit$Pvalue)
  }
  out <- dplyr::bind_rows(res); out$p_adj <- p.adjust(out$p, method = "BH"); out
}


analyze_beta_landscape <- function(beta_mat, site_ids = NULL, sample_data = sample_info, type_col = "Type",
                                   type_levels = c("Spr_TL", "Cas_TL", "Cus_TL", "Rap_TL", "Rim_TL"),
                                   my_cols = landscape_cols, nperm = 9999) {
  beta_mat <- as.matrix(beta_mat)
  if (is.null(site_ids)) site_ids <- rownames(beta_mat)
  if (is.null(site_ids)) site_ids <- seq_len(nrow(beta_mat))
  group <- factor(sample_data[[type_col]], levels = type_levels)
  dist_obj <- as.dist(beta_mat)
  
  adonis_res <- vegan::adonis2(dist_obj ~ group, permutations = nperm)
  anosim_res <- vegan::anosim(dist_obj, grouping = group, permutations = nperm)
  mrpp_res <- vegan::mrpp(dist_obj, grouping = group, permutations = nperm)
  overall_table <- data.frame(Method = c("PERMANOVA", "ANOSIM", "MRPP"), Statistic = c("F", "R", "A"),
                              Value = c(adonis_res$F[1], unname(anosim_res$statistic), mrpp_res$A),
                              R2 = c(adonis_res$R2[1], NA, NA),
                              P_value = c(adonis_res$`Pr(>F)`[1], anosim_res$signif, mrpp_res$Pvalue))
  
  bd <- vegan::betadisper(dist_obj, group = group)
  bd_anova <- anova(bd); bd_perm <- vegan::permutest(bd, permutations = nperm)
  bd_anova_df <- as.data.frame(bd_anova)
  dispersion_table <- data.frame(Method = c("betadisper_anova", "betadisper_permutest"),
                                 Statistic = "F",
                                 Value = c(bd_anova_df[1, grep("^F", names(bd_anova_df))],
                                           bd_perm$tab[1, "F"]),
                                 P_value = c(bd_anova_df[1, grep("Pr", names(bd_anova_df))],
                                             bd_perm$tab[1, "Pr(>F)"]))
  
  pcoa_res <- stats::cmdscale(dist_obj, eig = TRUE, k = 2)
  scores <- as.data.frame(pcoa_res$points); colnames(scores) <- c("PCoA1", "PCoA2")
  scores$Type <- factor(group, levels = type_levels); scores$Site <- site_ids
  eig <- pcoa_res$eig; var_explained <- eig / sum(eig[eig > 0])
  centroids <- scores %>% dplyr::group_by(Type) %>% dplyr::summarise(PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2), .groups = "drop")
  
  p_pcoa <- ggplot(scores, aes(PCoA1, PCoA2)) +
    stat_ellipse(aes(fill = Type, color = Type), geom = "polygon", level = 0.65, alpha = 0.10, linewidth = 0.27, show.legend = FALSE) +
    geom_point(aes(color = Type), shape = 16, size = 2, alpha = 0.50) +
    geom_point(data = centroids, aes(PCoA1, PCoA2, color = Type), shape = 4, size = 2.5, stroke = 0.75, show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey85", linewidth = 0.27) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey85", linewidth = 0.27) +
    scale_color_manual(values = my_cols) + scale_fill_manual(values = my_cols) +
    labs(x = paste0("PCoA1 (", round(var_explained[1] * 100, 2), "%)"),
         y = paste0("PCoA2 (", round(var_explained[2] * 100, 2), "%)")) +
    coord_equal() + theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.27),
          axis.title = element_text(size = 10), axis.text = element_text(colour = "black", size = 10),
          axis.line = element_line(linewidth = 0.27, colour = "black"),
          axis.ticks = element_line(linewidth = 0.27, colour = "black"))
  
  list(overall = list(adonis = adonis_res, anosim = anosim_res, mrpp = mrpp_res, summary_table = overall_table),
       pairwise = list(adonis = pairwise_adonis(beta_mat, group, nperm), anosim = pairwise_anosim(beta_mat, group, nperm), mrpp = pairwise_mrpp(beta_mat, group, nperm)),
       dispersion = list(betadisper = bd, anova = bd_anova, permutest = bd_perm, summary_table = dispersion_table),
       pcoa = list(raw = pcoa_res, scores = scores, centroids = centroids, var_explained = var_explained, plot = p_pcoa))
}




#### 5.2 Landscape effects on taxonomic β diversity ####

res_td_total_pa <- analyze_beta_landscape(beta_mat = td_total_pa, site_ids = rownames(comm_pa)); res_td_total_pa
res_td_repl_pa  <- analyze_beta_landscape(beta_mat = td_repl_pa,  site_ids = rownames(comm_pa)); res_td_repl_pa
res_td_rich_pa  <- analyze_beta_landscape(beta_mat = td_rich_pa,  site_ids = rownames(comm_pa)); res_td_rich_pa


#### 5.3 Landscape effects on functional β diversity ####

res_fd_total_pa <- analyze_beta_landscape(beta_mat = fd_total_pa, site_ids = rownames(comm_pa_fd)); res_fd_total_pa
res_fd_repl_pa  <- analyze_beta_landscape(beta_mat = fd_repl_pa,  site_ids = rownames(comm_pa_fd)); res_fd_repl_pa
res_fd_rich_pa  <- analyze_beta_landscape(beta_mat = fd_rich_pa,  site_ids = rownames(comm_pa_fd)); res_fd_rich_pa


#### 5.4 Landscape effects on phylogenetic β diversity ####

res_pd_total_pa <- analyze_beta_landscape(beta_mat = pd_total_pa, site_ids = rownames(comm_pd_pa)); res_pd_total_pa
res_pd_repl_pa  <- analyze_beta_landscape(beta_mat = pd_repl_pa,  site_ids = rownames(comm_pd_pa)); res_pd_repl_pa
res_pd_rich_pa  <- analyze_beta_landscape(beta_mat = pd_rich_pa,  site_ids = rownames(comm_pd_pa)); res_pd_rich_pa


#### 5.5 PCoA visualization of landscape effects ####

p_pcoa_pa_all <- cowplot::plot_grid(
  res_td_total_pa$pcoa$plot, res_td_repl_pa$pcoa$plot, res_td_rich_pa$pcoa$plot,
  res_fd_total_pa$pcoa$plot, res_fd_repl_pa$pcoa$plot, res_fd_rich_pa$pcoa$plot,
  res_pd_total_pa$pcoa$plot, res_pd_repl_pa$pcoa$plot, res_pd_rich_pa$pcoa$plot,
  ncol = 3, nrow = 3, align = "hv"); p_pcoa_pa_all

# ggsave(filename = "05 PCoA_landscape_effects_presence_absence.pdf", plot = p_pcoa_pa_all, width = 18, height = 18, units = "cm", device = cairo_pdf)




#### 5.6 Extracting landscape-effect test results ####

extract_beta_landscape_results <- function(res_obj, dimension, component) {
  overall_df <- data.frame(Data_type = "Presence-absence", Dimension = dimension, Component = component,
                           Method = c("PERMANOVA", "ANOSIM", "MRPP"), Statistic = c("F", "R", "A"),
                           Value = c(res_obj$overall$adonis$F[1], unname(res_obj$overall$anosim$statistic), res_obj$overall$mrpp$A),
                           R2 = c(res_obj$overall$adonis$R2[1], NA, NA),
                           P_value = c(res_obj$overall$adonis$`Pr(>F)`[1], res_obj$overall$anosim$signif, res_obj$overall$mrpp$Pvalue))
  
  pairwise_df <- dplyr::bind_rows(
    res_obj$pairwise$adonis %>% dplyr::transmute(group1, group2, Method = "PERMANOVA", Statistic = "F", Value = F, R2 = R2, P_value = p, P_adj = p_adj),
    res_obj$pairwise$anosim %>% dplyr::transmute(group1, group2, Method = "ANOSIM", Statistic = "R", Value = R, R2 = NA_real_, P_value = p, P_adj = p_adj),
    res_obj$pairwise$mrpp %>% dplyr::transmute(group1, group2, Method = "MRPP", Statistic = "A", Value = A, R2 = NA_real_, P_value = p, P_adj = p_adj)
  ) %>%
    dplyr::mutate(Data_type = "Presence-absence", Dimension = dimension, Component = component) %>%
    dplyr::select(Data_type, Dimension, Component, group1, group2, Method, Statistic, Value, R2, P_value, P_adj)
  
  list(overall = overall_df, pairwise = pairwise_df)
}

beta_landscape_result_list <- list(
  extract_beta_landscape_results(res_td_total_pa, "TD", "Total beta diversity"),
  extract_beta_landscape_results(res_td_repl_pa,  "TD", "Replacement component"),
  extract_beta_landscape_results(res_td_rich_pa,  "TD", "Richness-difference component"),
  extract_beta_landscape_results(res_fd_total_pa, "FD", "Total beta diversity"),
  extract_beta_landscape_results(res_fd_repl_pa,  "FD", "Replacement component"),
  extract_beta_landscape_results(res_fd_rich_pa,  "FD", "Richness-difference component"),
  extract_beta_landscape_results(res_pd_total_pa, "PD", "Total beta diversity"),
  extract_beta_landscape_results(res_pd_repl_pa,  "PD", "Replacement component"),
  extract_beta_landscape_results(res_pd_rich_pa,  "PD", "Richness-difference component")
)

beta_landscape_overall_df <- dplyr::bind_rows(lapply(beta_landscape_result_list, function(x) x$overall)); beta_landscape_overall_df
beta_landscape_pairwise_df <- dplyr::bind_rows(lapply(beta_landscape_result_list, function(x) x$pairwise)); beta_landscape_pairwise_df

# write.csv(beta_landscape_overall_df, "09 Beta landscape overall results_presence_absence.csv", row.names = FALSE)
# write.csv(beta_landscape_pairwise_df, "10 Beta landscape pairwise results_presence_absence.csv", row.names = FALSE)




#### 5.7 Extracting dispersion and PCoA1 results ####

extract_dispersion_results <- function(res_obj, dimension, component) {
  res_obj$dispersion$summary_table %>%
    dplyr::transmute(Data_type = "Presence-absence", Dimension = dimension, Component = component,
                     Method, Statistic, Value, P_value)
}

dispersion_results_list <- list(
  extract_dispersion_results(res_td_total_pa, "TD", "Total beta diversity"),
  extract_dispersion_results(res_td_repl_pa,  "TD", "Replacement component"),
  extract_dispersion_results(res_td_rich_pa,  "TD", "Richness-difference component"),
  extract_dispersion_results(res_fd_total_pa, "FD", "Total beta diversity"),
  extract_dispersion_results(res_fd_repl_pa,  "FD", "Replacement component"),
  extract_dispersion_results(res_fd_rich_pa,  "FD", "Richness-difference component"),
  extract_dispersion_results(res_pd_total_pa, "PD", "Total beta diversity"),
  extract_dispersion_results(res_pd_repl_pa,  "PD", "Replacement component"),
  extract_dispersion_results(res_pd_rich_pa,  "PD", "Richness-difference component")
)

dispersion_overall_df <- dplyr::bind_rows(dispersion_results_list); dispersion_overall_df
# write.csv(dispersion_overall_df, "11 Beta landscape dispersion results_presence_absence.csv", row.names = FALSE)


extract_pcoa1 <- function(res_obj, col_name) {
  out <- res_obj$pcoa$scores %>% dplyr::select(Site, PCoA1)
  colnames(out)[2] <- col_name
  out
}

pcoa1_list <- list(
  extract_pcoa1(res_td_total_pa, "TD_total_PA"),
  extract_pcoa1(res_td_repl_pa,  "TD_repl_PA"),
  extract_pcoa1(res_td_rich_pa,  "TD_rich_PA"),
  extract_pcoa1(res_fd_total_pa, "FD_total_PA"),
  extract_pcoa1(res_fd_repl_pa,  "FD_repl_PA"),
  extract_pcoa1(res_fd_rich_pa,  "FD_rich_PA"),
  extract_pcoa1(res_pd_total_pa, "PD_total_PA"),
  extract_pcoa1(res_pd_repl_pa,  "PD_repl_PA"),
  extract_pcoa1(res_pd_rich_pa,  "PD_rich_PA")
)

pcoa1_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Site"), pcoa1_list); pcoa1_df
# write.csv(pcoa1_df, "12 Beta landscape PCoA1 scores_presence_absence.csv", row.names = FALSE)






#### 6. Drivers of multidimensional beta diversity ####

#### 6.1 Correlation screening between environmental factors and beta-diversity axes ####

#### 6.1.1 Correlation heatmap function ####

correlation_heatmap <- function(community, environment, method = c("pearson", "spearman"),
                                adjust = c("BH", "none"), cluster_rows = FALSE,
                                cluster_cols = FALSE, cutree_rows = NA, cutree_cols = NA) {
  method <- match.arg(method)
  adjust <- match.arg(adjust)
  
  if (!is.matrix(community) && !is.data.frame(community)) stop("community must be a matrix or data.frame.")
  if (!is.matrix(environment) && !is.data.frame(environment)) stop("environment must be a matrix or data.frame.")
  
  community <- as.data.frame(community)
  environment <- as.data.frame(environment)
  
  cor_result <- psych::corr.test(x = community, y = environment, method = method,
                                 adjust = adjust, use = "pairwise", ci = FALSE)
  r_mat <- cor_result$r
  p_mat <- cor_result$p
  
  sig_mat <- matrix("", nrow = nrow(p_mat), ncol = ncol(p_mat), dimnames = dimnames(p_mat))
  sig_mat[p_mat <= 0.001] <- "***"
  sig_mat[p_mat > 0.001 & p_mat <= 0.01] <- "**"
  sig_mat[p_mat > 0.01 & p_mat <= 0.05] <- "*"
  
  my_cols <- colorRampPalette(c("#08519C", "#F7F7F7", "#A50F15"))(100)
  
  heatmap <- pheatmap::pheatmap(mat = r_mat, color = my_cols, breaks = seq(-1, 1, length.out = 101),
                                cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                                clustering_method = "average", cutree_rows = cutree_rows,
                                cutree_cols = cutree_cols, treeheight_row = 25,
                                treeheight_col = 15, border_color = "grey95",
                                cellwidth = 13, cellheight = 13, fontsize_row = 10,
                                fontsize_col = 10, display_numbers = sig_mat,
                                fontsize_number = 10, number_color = "black",
                                legend_breaks = c(-1, -0.5, 0, 0.5, 1),
                                legend_labels = c("-1", "-0.5", "0", "0.5", "1"))
  
  return(list(correlation_matrix = r_mat, p_value_matrix = p_mat,
              significance_matrix = sig_mat, heatmap = heatmap))
}




#### 6.1.2 Correlations between environmental factors and beta-diversity PCoA1 axes ####

res_cor_pa <- correlation_heatmap(community = env_raw[, 3:26], environment = pcoa1_df[, 2:10],
                                  method = "spearman", adjust = "none",
                                  cluster_rows = TRUE, cluster_cols = TRUE,
                                  cutree_rows = 4, cutree_cols = 2); res_cor_pa

res_cor_environment <- correlation_heatmap(community = env_raw[, 3:26], environment = env_raw[, 3:26],
                                           method = "spearman", adjust = "none",
                                           cluster_rows = TRUE, cluster_cols = TRUE,
                                           cutree_rows = 2, cutree_cols = 2); res_cor_environment

res_cor_pcoa1 <- correlation_heatmap(community = pcoa1_df[, 2:10], environment = pcoa1_df[, 2:10],
                                     method = "spearman", adjust = "none",
                                     cluster_rows = TRUE, cluster_cols = TRUE,
                                     cutree_rows = 2, cutree_cols = 2); res_cor_pcoa1

res_cor_predictors <- correlation_heatmap(community = vp_predictors_df[, 2:8], environment = vp_predictors_df[, 2:8],
                                          method = "spearman", adjust = "none",
                                          cluster_rows = TRUE, cluster_cols = TRUE,
                                          cutree_rows = 2, cutree_cols = 2); res_cor_predictors


ggplot2::ggsave("06 Correlation heatmap_environment_vs_beta_PCoA1.pdf", plot = res_cor_pa$heatmap$gtable,
                width = 18.5, height = 22, units = "cm", device = cairo_pdf)

ggplot2::ggsave("06 Correlation heatmap_environment_variables.pdf", plot = res_cor_environment$heatmap$gtable,
                width = 18.5, height = 18.5, units = "cm", device = cairo_pdf)

ggplot2::ggsave("06 Correlation heatmap_beta_PCoA1_axes.pdf", plot = res_cor_pcoa1$heatmap$gtable,
                width = 16, height = 16, units = "cm", device = cairo_pdf)

ggplot2::ggsave("06 Correlation heatmap_predictor_axes.pdf", plot = res_cor_predictors$heatmap$gtable,
                width = 14, height = 14, units = "cm", device = cairo_pdf)




#### 6.2 Ecological latent variables for driver analysis ####

#### 6.2.1 Data import and core object preparation ####

moss_raw <- read.csv("moss.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
env_raw <- read.csv("env.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

mite_abund <- comm_abund
mite_pa <- comm_pa

mite_sample_info <- sample_info %>%
  dplyr::mutate(SampleID = as.character(SampleID),
                Type = factor(Type, levels = c("Spr_TL", "Cas_TL", "Cus_TL", "Rap_TL", "Rim_TL")))

mite_beta_pa <- beta_df_pa

moss_comm <- moss_raw %>%
  dplyr::mutate(SampleID = paste(Type, Site, sep = "_")) %>%
  tibble::column_to_rownames("SampleID") %>%
  dplyr::select(-Type, -Site)

moss_pa <- moss_comm
moss_pa[moss_pa > 0] <- 1

env_data <- env_raw %>%
  dplyr::mutate(SampleID = paste(Type, Site, sep = "_"),
                Type = factor(Type, levels = c("Spr_TL", "Cas_TL", "Cus_TL", "Rap_TL", "Rim_TL"))) %>%
  tibble::column_to_rownames("SampleID")

identical(rownames(mite_abund), rownames(mite_pa))
identical(rownames(mite_abund), rownames(moss_comm))
identical(rownames(mite_abund), rownames(env_data))
identical(rownames(mite_abund), mite_sample_info$SampleID)
setdiff(rownames(mite_abund), rownames(moss_comm))
setdiff(rownames(mite_abund), rownames(env_data))




#### 6.2.2 Ecological latent-variable matrices ####

lv_moss_comm <- moss_comm

lv_hydro_substrate <- env_data %>%
  dplyr::select(Temp, RH, pH, Cond, OC, TN, TP, `C:N`)

lv_moss_function <- env_data %>%
  dplyr::select(MTC, MTN, MTP, `TC:TN`, Chl_a, Chl_b, Biomass, WC, MDA, TSS, POD, SOD, CAT)

lv_space <- env_data %>%
  dplyr::select(Lng, Lat, Alt); lv_space




#### 6.2.3 Distance matrices of ecological latent variables ####

lv_hydro_substrate_scaled <- lv_hydro_substrate %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(scale(.x))))

lv_moss_function_scaled <- lv_moss_function %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(scale(.x))))

lv_space_scaled <- lv_space %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(scale(.x))))

dist_moss_comm <- vegan::vegdist(lv_moss_comm, method = "bray")
dist_hydro_substrate <- dist(lv_hydro_substrate_scaled, method = "euclidean")
dist_moss_function <- dist(lv_moss_function_scaled, method = "euclidean")
dist_space <- dist(lv_space_scaled, method = "euclidean")

latent_dist_list <- list(Moss_community = dist_moss_comm,
                         Hydro_substrate = dist_hydro_substrate,
                         Moss_function = dist_moss_function,
                         Spatial_factors = dist_space)

sapply(latent_dist_list, length)




#### 6.3 Mantel tests between beta-diversity matrices and ecological drivers ####

#### 6.3.1 Beta-diversity distance matrices based on presence-absence data ####

mite_beta_dist <- list(
  TD_total = as.dist(td_total_pa),
  TD_repl = as.dist(td_repl_pa),
  TD_rich = as.dist(td_rich_pa),
  FD_total = as.dist(fd_total_pa),
  FD_repl = as.dist(fd_repl_pa),
  FD_rich = as.dist(fd_rich_pa),
  PD_total = as.dist(pd_total_pa),
  PD_repl = as.dist(pd_repl_pa),
  PD_rich = as.dist(pd_rich_pa)
)

sapply(mite_beta_dist, length)
names(mite_beta_dist)




#### 6.3.2 Mantel test function ####

run_mantel <- function(beta_dist_list, latent_dist_list, data_type_name = "Presence-absence") {
  set.seed(123)
  
  purrr::map_dfr(names(beta_dist_list), function(beta_name) {
    purrr::map_dfr(names(latent_dist_list), function(driver_name) {
      res <- vegan::mantel(xdis = beta_dist_list[[beta_name]],
                           ydis = latent_dist_list[[driver_name]],
                           method = "spearman", permutations = 9999)
      data.frame(Data_type = data_type_name, Beta_metric = beta_name, Driver = driver_name,
                 Mantel_r = unname(res$statistic), p_value = res$signif)
    })
  }) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"),
                  Significance = dplyr::case_when(p_adj <= 0.001 ~ "***",
                                                  p_adj <= 0.01 ~ "**",
                                                  p_adj <= 0.05 ~ "*",
                                                  TRUE ~ ""))
}

mantel_pa <- run_mantel(beta_dist_list = mite_beta_dist,
                        latent_dist_list = latent_dist_list,
                        data_type_name = "Presence-absence"); mantel_pa

# write.csv(mantel_pa, "11 Mantel_PA_results.csv", row.names = FALSE)




#### 6.3.3 Mantel heatmap visualization ####

mantel_plot_pa_horizontal <- mantel_pa %>%
  dplyr::mutate(Driver = factor(Driver, levels = c("Spatial_factors", "Hydro_substrate", "Moss_function", "Moss_community")),
                Beta_metric = factor(Beta_metric, levels = c("TD_total", "TD_repl", "TD_rich",
                                                             "FD_total", "FD_repl", "FD_rich",
                                                             "PD_total", "PD_repl", "PD_rich")),
                x_pos = dplyr::case_when(Beta_metric == "TD_total" ~ 1, Beta_metric == "TD_repl" ~ 2, Beta_metric == "TD_rich" ~ 3,
                                         Beta_metric == "FD_total" ~ 4.4, Beta_metric == "FD_repl" ~ 5.4, Beta_metric == "FD_rich" ~ 6.4,
                                         Beta_metric == "PD_total" ~ 7.8, Beta_metric == "PD_repl" ~ 8.8, Beta_metric == "PD_rich" ~ 9.8),
                y_pos = dplyr::case_when(Driver == "Moss_community" ~ 1, Driver == "Moss_function" ~ 2,
                                         Driver == "Hydro_substrate" ~ 3, Driver == "Spatial_factors" ~ 4),
                Significance_raw = dplyr::case_when(p_value <= 0.001 ~ "***",
                                                    p_value <= 0.01 ~ "**",
                                                    p_value <= 0.05 ~ "*",
                                                    TRUE ~ ""),
                label_text = ifelse(Significance_raw == "", sprintf("%.2f", Mantel_r),
                                    paste0(sprintf("%.2f", Mantel_r), "\n", Significance_raw))); mantel_plot_pa_horizontal


p_mantel_pa_horizontal <- ggplot2::ggplot(mantel_plot_pa_horizontal,
                                          ggplot2::aes(x = x_pos, y = y_pos, fill = Mantel_r)) +
  ggplot2::geom_tile(width = 0.75, height = 0.75, color = "grey95") +
  ggplot2::geom_text(ggplot2::aes(label = label_text), size = 3, color = "black") +
  ggplot2::scale_fill_gradient2(low = "#3182BD", mid = "white", high = "#8C6BB1", midpoint = 0) +
  ggplot2::scale_x_continuous(breaks = c(1, 2, 3, 4.4, 5.4, 6.4, 7.8, 8.8, 9.8),
                              labels = c("Taxonomic total", "Taxonomic repl.", "Taxonomic rich.",
                                         "Functional total", "Functional repl.", "Functional rich.",
                                         "Phylogenetic total", "Phylogenetic repl.", "Phylogenetic rich."),
                              expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
  ggplot2::scale_y_continuous(breaks = 1:4,
                              labels = c("Moss community", "Moss function", "Hydro-substrate", "Spatial factors"),
                              expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_classic() +
  ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.75),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black", size = 10),
                 axis.text.y = ggplot2::element_text(colour = "black", size = 10),
                 axis.line = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_line(linewidth = 0.75, colour = "black"),
                 legend.title = ggplot2::element_text(size = 10, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10),
                 legend.key.height = grid::unit(0.5, "cm"),
                 legend.key.width = grid::unit(0.35, "cm")); p_mantel_pa_horizontal

# ggplot2::ggsave("06 Mantel_PA_horizontal.pdf", p_mantel_pa_horizontal, width = 18, height = 6.5, units = "cm")




#### 6.4 Mantel correlations among ecological latent variables ####

set.seed(123)

latent_cor_result <- purrr::map_dfr(names(latent_dist_list), function(driver1) {
  purrr::map_dfr(names(latent_dist_list), function(driver2) {
    res <- vegan::mantel(xdis = latent_dist_list[[driver1]],
                         ydis = latent_dist_list[[driver2]],
                         method = "spearman", permutations = 9999)
    data.frame(Driver1 = driver1, Driver2 = driver2,
               Mantel_r = unname(res$statistic), p_value = res$signif)
  })
}) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"),
                Significance = dplyr::case_when(p_adj <= 0.001 ~ "***",
                                                p_adj <= 0.01 ~ "**",
                                                p_adj <= 0.05 ~ "*",
                                                TRUE ~ ""),
                x_pos = as.numeric(factor(Driver2, levels = names(latent_dist_list))),
                y_pos = as.numeric(factor(Driver1, levels = rev(names(latent_dist_list)))),
                label_text = ifelse(Significance == "", sprintf("%.2f", Mantel_r),
                                    paste0(sprintf("%.2f", Mantel_r), "\n", Significance))); latent_cor_result

# write.csv(latent_cor_result, "13 Latent Mantel results.csv", row.names = FALSE)


p_latent_mantel <- ggplot2::ggplot(latent_cor_result, ggplot2::aes(x = x_pos, y = y_pos, fill = Mantel_r)) +
  ggplot2::geom_tile(width = 0.95, height = 0.95, color = "grey95") +
  ggplot2::geom_text(ggplot2::aes(label = label_text), size = 3.5, color = "black", lineheight = 1) +
  ggplot2::scale_fill_gradient2(low = "#3182BD", mid = "white", high = "#8C6BB1", midpoint = 0, name = "Mantel r") +
  ggplot2::scale_x_continuous(breaks = 1:4, labels = names(latent_dist_list),
                              expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
  ggplot2::scale_y_continuous(breaks = 1:4, labels = rev(names(latent_dist_list)),
                              expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_classic() +
  ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black", size = 10),
                 axis.text.y = ggplot2::element_text(colour = "black", size = 10),
                 axis.line = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_line(linewidth = 0.5, colour = "black"),
                 legend.title = ggplot2::element_text(size = 10, face = "bold"),
                 legend.text = ggplot2::element_text(size = 9),
                 legend.key.height = grid::unit(0.5, "cm"),
                 legend.key.width = grid::unit(0.35, "cm")); p_latent_mantel

ggplot2::ggsave("06 Latent Mantel.pdf", plot = p_latent_mantel, width = 4, height = 3, units = "in")




#### 6.5 Multiple regression on distance matrices ####

#### 6.5.1 MRM models based on presence-absence beta diversity ####

set.seed(123)

mrm_result_pa_full <- purrr::map_dfr(names(mite_beta_dist), function(beta_name) {
  mrm_res <- ecodist::MRM(mite_beta_dist[[beta_name]] ~ latent_dist_list$Moss_community +
                            latent_dist_list$Hydro_substrate + latent_dist_list$Moss_function +
                            latent_dist_list$Spatial_factors, nperm = 9999)
  
  as.data.frame(mrm_res$coef) %>%
    tibble::rownames_to_column("Term") %>%
    stats::setNames(c("Term", "Estimate", "p_value")) %>%
    dplyr::mutate(Data_type = "Presence-absence", Beta_metric = beta_name,
                  R2 = as.numeric(mrm_res$r.squared["R2"]),
                  R2_p_value = as.numeric(mrm_res$r.squared["pval"]),
                  F_value = as.numeric(mrm_res$F.test["F"]),
                  F_p_value = as.numeric(mrm_res$F.test["F.pval"]))
}) %>%
  dplyr::mutate(Term = dplyr::recode(Term,
                                     "Int" = "Intercept",
                                     "latent_dist_list$Moss_community" = "Moss_community",
                                     "latent_dist_list$Hydro_substrate" = "Hydro_substrate",
                                     "latent_dist_list$Moss_function" = "Moss_function",
                                     "latent_dist_list$Spatial_factors" = "Spatial_factors")); mrm_result_pa_full

mrm_result_pa <- mrm_result_pa_full %>%
  dplyr::mutate(Significance = dplyr::case_when(p_value <= 0.001 ~ "***",
                                                p_value <= 0.01 ~ "**",
                                                p_value <= 0.05 ~ "*",
                                                TRUE ~ "")) %>%
  dplyr::select(Beta_metric, Data_type, Term, Estimate, p_value, Significance); mrm_result_pa


mrm_model_fit_pa <- mrm_result_pa_full %>%
  dplyr::distinct(Data_type, Beta_metric, R2, R2_p_value, F_value, F_p_value) %>%
  dplyr::mutate(Beta_metric = factor(Beta_metric, levels = c("TD_total", "TD_repl", "TD_rich",
                                                             "FD_total", "FD_repl", "FD_rich",
                                                             "PD_total", "PD_repl", "PD_rich"))) %>%
  dplyr::arrange(Beta_metric); mrm_model_fit_pa

# write.csv(mrm_result_pa, "14 MRM coefficients PA.csv", row.names = FALSE)
# write.csv(mrm_model_fit_pa, "15 MRM model fit PA.csv", row.names = FALSE)




#### 6.6 Variation partitioning and hierarchical partitioning ####

#### 6.6.1 Predictor axes for variation partitioning ####

get_pcoa_axes <- function(dist_obj, prefix, cum_var = 0.80, max_axes = 2) {
  dist_obj <- as.dist(dist_obj)
  attr(dist_obj, "Labels") <- rownames(mite_abund)
  pcoa_res <- stats::cmdscale(dist_obj, eig = TRUE, k = nrow(mite_abund) - 1)
  eig_values <- pcoa_res$eig
  pos_eig <- eig_values[eig_values > 0]
  eig_prop <- pos_eig / sum(pos_eig)
  eig_cum <- cumsum(eig_prop)
  n_axes <- max(min(which(eig_cum >= cum_var)[1], max_axes), 1)
  axes_df <- as.data.frame(pcoa_res$points[, 1:n_axes, drop = FALSE])
  colnames(axes_df) <- paste0(prefix, "_PCoA", seq_len(n_axes))
  rownames(axes_df) <- rownames(mite_abund)
  return(list(axes = axes_df, eig_prop = eig_prop, eig_cum = eig_cum, n_axes = n_axes))
}


get_pca_axes <- function(data_scaled, prefix, cum_var = 0.80, max_axes = 3) {
  pca_res <- stats::prcomp(data_scaled, center = FALSE, scale. = FALSE)
  eig_prop <- (pca_res$sdev^2) / sum(pca_res$sdev^2)
  eig_cum <- cumsum(eig_prop)
  n_axes <- max(min(which(eig_cum >= cum_var)[1], max_axes), 1)
  axes_df <- as.data.frame(pca_res$x[, 1:n_axes, drop = FALSE])
  colnames(axes_df) <- paste0(prefix, "_PC", seq_len(n_axes))
  rownames(axes_df) <- rownames(data_scaled)
  loading_df <- as.data.frame(pca_res$rotation[, 1:n_axes, drop = FALSE]) %>%
    tibble::rownames_to_column("Variable")
  colnames(loading_df)[-1] <- paste0(prefix, "_PC", seq_len(n_axes))
  return(list(axes = axes_df, loadings = loading_df, eig_prop = eig_prop,
              eig_cum = eig_cum, n_axes = n_axes, pca = pca_res))
}

pcoa_moss_comm <- get_pcoa_axes(latent_dist_list$Moss_community, "MossComm", cum_var = 0.80, max_axes = 2); pcoa_moss_comm
pca_hydro_substrate <- get_pca_axes(lv_hydro_substrate_scaled, "HydroSub", cum_var = 0.80, max_axes = 2); pca_hydro_substrate
pca_moss_function <- get_pca_axes(lv_moss_function_scaled, "MossFunc", cum_var = 0.80, max_axes = 2); pca_moss_function
pca_space <- get_pca_axes(lv_space_scaled, "Space", cum_var = 0.80, max_axes = 1); pca_space

axis_summary_formal <- data.frame(
  Predictor_set = c("Moss_community", "Hydro_substrate", "Moss_function", "Spatial_factors"),
  Ordination_method = c("PCoA", "PCA", "PCA", "PCA"),
  Retained_axes = c(pcoa_moss_comm$n_axes, pca_hydro_substrate$n_axes,
                    pca_moss_function$n_axes, pca_space$n_axes)
); axis_summary_formal

vp_predictors_formal_2axes <- list(Moss_community = pcoa_moss_comm$axes,
                                   Hydro_substrate = pca_hydro_substrate$axes,
                                   Moss_function = pca_moss_function$axes,
                                   Spatial_factors = pca_space$axes)

sapply(vp_predictors_formal_2axes, ncol)
all(rownames(vp_predictors_formal_2axes$Moss_community) == rownames(mite_abund),
    rownames(vp_predictors_formal_2axes$Hydro_substrate) == rownames(mite_abund),
    rownames(vp_predictors_formal_2axes$Moss_function) == rownames(mite_abund),
    rownames(vp_predictors_formal_2axes$Spatial_factors) == rownames(mite_abund))

vp_predictors_df <- dplyr::bind_cols(SampleID = rownames(mite_abund),
                                     vp_predictors_formal_2axes$Moss_community,
                                     vp_predictors_formal_2axes$Hydro_substrate,
                                     vp_predictors_formal_2axes$Moss_function,
                                     vp_predictors_formal_2axes$Spatial_factors); vp_predictors_df

# write.csv(vp_predictors_df, file = "18 Variation partitioning predictor axes.csv", row.names = FALSE)




#### 6.6.2 dbRDA-based hierarchical partitioning ####

rdacca_pa_list <- purrr::map(names(mite_beta_dist), function(beta_name) {
  rdacca.hp::rdacca.hp(dv = mite_beta_dist[[beta_name]], iv = vp_predictors_formal_2axes,
                       method = "dbRDA", type = "adjR2", var.part = TRUE)
})

names(rdacca_pa_list) <- names(mite_beta_dist)

rdacca_hp_pa_df <- purrr::map_dfr(names(rdacca_pa_list), function(beta_name) {
  rdacca_pa_list[[beta_name]]$Hier.part %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Predictor_set") %>%
    dplyr::mutate(Data_type = "Presence-absence", Beta_metric = beta_name,
                  Total_explained_variation = rdacca_pa_list[[beta_name]]$Total_explained_variation)
}) %>%
  dplyr::select(Data_type, Beta_metric, Predictor_set, Unique, Average.share,
                Individual, `I.perc(%)`, Total_explained_variation); rdacca_hp_pa_df

# write.csv(rdacca_hp_pa_df, "19 rdacca hp PA.csv", row.names = FALSE)




#### 6.6.3 Significance tests for unique fractions using partial dbRDA  ####

test_unique_fraction <- function(response_dist, predictor_list, nperm = 9999) {
  X_all <- dplyr::bind_cols(predictor_list)
  
  purrr::map_dfr(names(predictor_list), function(test_group) {
    test_vars <- colnames(predictor_list[[test_group]])
    cond_vars <- setdiff(colnames(X_all), test_vars)
    
    mod <- vegan::dbrda(stats::as.formula(paste0("response_dist ~ ",
                                                 paste(test_vars, collapse = " + "),
                                                 " + Condition(",
                                                 paste(cond_vars, collapse = " + "),
                                                 ")")), data = X_all)
    perm_res <- vegan::anova.cca(mod, permutations = nperm)
    
    data.frame(Predictor_set = test_group,
               Adj_R2_unique = vegan::RsquareAdj(mod)$adj.r.squared,
               F_value = perm_res$F[1],
               p_value = perm_res$`Pr(>F)`[1])
  }) %>%
    dplyr::mutate(Significance = dplyr::case_when(p_value <= 0.001 ~ "***",
                                                  p_value <= 0.01 ~ "**",
                                                  p_value <= 0.05 ~ "*",
                                                  TRUE ~ ""))
}

set.seed(123)

target_beta_metrics_pa <- c("TD_total", "TD_repl", "TD_rich",
                            "FD_total", "FD_repl", "FD_rich",
                            "PD_total", "PD_repl", "PD_rich")

varpart_unique_pa_formal <- purrr::map_dfr(target_beta_metrics_pa, function(beta_name) {
  test_unique_fraction(response_dist = mite_beta_dist[[beta_name]],
                       predictor_list = vp_predictors_formal_2axes,
                       nperm = 9999) %>%
    dplyr::mutate(Data_type = "Presence-absence", Beta_metric = beta_name)
}) %>%
  dplyr::select(Data_type, Beta_metric, Predictor_set, Adj_R2_unique,
                F_value, p_value, Significance); varpart_unique_pa_formal

# write.csv(varpart_unique_pa_formal, "20 partial dbRDA unique PA.csv", row.names = FALSE)




#### 6.6.4 Overall dbRDA model fit and explained variation ####

set.seed(123)

varpart_model_fit_pa_formal <- purrr::map_dfr(target_beta_metrics_pa, function(beta_name) {
  X_all <- dplyr::bind_cols(vp_predictors_formal_2axes)
  mod <- vegan::dbrda(mite_beta_dist[[beta_name]] ~ ., data = X_all)
  r2_res <- vegan::RsquareAdj(mod)
  perm_res <- vegan::anova.cca(mod, permutations = 9999)
  
  data.frame(Data_type = "Presence-absence", Beta_metric = beta_name,
             Total_axes = ncol(X_all),
             R2 = r2_res$r.squared,
             Adj_R2 = r2_res$adj.r.squared,
             Residual = 1 - r2_res$adj.r.squared,
             F_value = perm_res$F[1],
             p_value = perm_res$`Pr(>F)`[1])
}) %>%
  dplyr::mutate(Significance = dplyr::case_when(p_value <= 0.001 ~ "***",
                                                p_value <= 0.01 ~ "**",
                                                p_value <= 0.05 ~ "*",
                                                TRUE ~ "")); varpart_model_fit_pa_formal

# write.csv(varpart_model_fit_pa_formal, "21 overall dbRDA PA.csv", row.names = FALSE)


p_dbRDA_fit_pa <- varpart_model_fit_pa_formal %>%
  dplyr::mutate(
    Beta_metric = factor(Beta_metric, levels = c("TD_total", "TD_repl", "TD_rich",
                                                 "FD_total", "FD_repl", "FD_rich",
                                                 "PD_total", "PD_repl", "PD_rich")),
    Dimension = factor(stringr::str_extract(Beta_metric, "TD|FD|PD"), levels = c("TD", "FD", "PD")),
    Beta_label = dplyr::recode(as.character(Beta_metric),
                               "TD_total" = "TD total", "TD_repl" = "TD repl", "TD_rich" = "TD rich",
                               "FD_total" = "FD total", "FD_repl" = "FD repl", "FD_rich" = "FD rich",
                               "PD_total" = "PD total", "PD_repl" = "PD repl", "PD_rich" = "PD rich"),
    Beta_label = factor(Beta_label, levels = c("TD total", "TD repl", "TD rich",
                                               "FD total", "FD repl", "FD rich",
                                               "PD total", "PD repl", "PD rich")),
    label = paste0(round(Adj_R2 * 100, 1), "%", Significance)
  ) %>%
  ggplot2::ggplot(ggplot2::aes(x = Beta_label, y = Adj_R2, fill = Dimension)) +
  ggplot2::geom_col(width = 0.72, color = "grey100", linewidth = 0.5) +
  ggplot2::geom_text(ggplot2::aes(label = label), vjust = -0.35, size = 3.5, color = "black") +
  ggplot2::facet_wrap(~ Dimension, nrow = 1, scales = "free_x") +
  ggplot2::scale_fill_manual(values = c("TD" = "#80B1D3", "FD" = "#BC80BD", "PD" = "#FDB462")) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                              expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
  ggplot2::labs(x = NULL, y = "Adjusted R²", fill = NULL) +
  ggplot2::theme_classic() +
  ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
                 axis.title = ggplot2::element_text(face = "bold", size = 11),
                 axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, colour = "black", size = 10),
                 axis.text.y = ggplot2::element_text(colour = "black", size = 10),
                 axis.line = ggplot2::element_blank(),
                 legend.position = "none",
                 strip.background = ggplot2::element_rect(fill = "grey95", colour = "black", linewidth = 0.5),
                 strip.text = ggplot2::element_text(face = "bold", size = 11)); p_dbRDA_fit_pa

ggplot2::ggsave(filename = "09 rdacca_fit_PA.pdf", plot = p_dbRDA_fit_pa, width = 18, height = 8, units = "cm")




#### 6.7 Visualization of hierarchical partitioning results ####

plot_rdacca_hp <- function(rdacca_obj, beta_title) {
  hp_df <- rdacca_obj$Hier.part %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Predictor") %>%
    dplyr::mutate(Predictor = factor(Predictor, levels = c("Moss_community", "Hydro_substrate", "Moss_function", "Spatial_factors")),
                  Predictor_label = dplyr::recode(Predictor,
                                                  "Moss_community" = "Moss community",
                                                  "Hydro_substrate" = "Hydro-substrate",
                                                  "Moss_function" = "Moss function",
                                                  "Spatial_factors" = "Spatial factors"))
  
  p_bar <- ggplot2::ggplot(hp_df, ggplot2::aes(x = reorder(Predictor_label, -`I.perc(%)`),
                                               y = `I.perc(%)`, fill = Predictor)) +
    ggplot2::geom_col(width = 0.75, color = "grey100", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(`I.perc(%)`, 1), "%")),
                       vjust = -0.4, size = 3.8, color = "black") +
    ggplot2::scale_fill_manual(values = c("Moss_community" = "#23793A",
                                          "Hydro_substrate" = "#436FA9",
                                          "Moss_function" = "#EE7C6E",
                                          "Spatial_factors" = "#6A3906")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
    ggplot2::labs(x = NULL, y = "Individual contribution (%)", title = beta_title) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.5),
                   axis.title = ggplot2::element_text(face = "bold", size = 11),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, colour = "black", size = 10),
                   axis.text.y = ggplot2::element_text(colour = "black", size = 10),
                   axis.line = ggplot2::element_blank(),
                   legend.position = "none",
                   plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0.5))
  
  p_rdacca <- plot(rdacca_obj, plot.perc = TRUE)
  return(list(bar = p_bar, rdacca = p_rdacca, data = hp_df))
}

beta_title_names <- c("TD_total" = "Taxonomic total", "TD_repl" = "Taxonomic replacement", "TD_rich" = "Taxonomic richness difference",
                      "FD_total" = "Functional total", "FD_repl" = "Functional replacement", "FD_rich" = "Functional richness difference",
                      "PD_total" = "Phylogenetic total", "PD_repl" = "Phylogenetic replacement", "PD_rich" = "Phylogenetic richness difference")

rdacca_pa_plot_list <- purrr::map(names(rdacca_pa_list), function(beta_name) {
  plot_rdacca_hp(rdacca_obj = rdacca_pa_list[[beta_name]],
                 beta_title = beta_title_names[beta_name])
})

names(rdacca_pa_plot_list) <- names(rdacca_pa_list)

p_rdacca_pa_bar_all <- cowplot::plot_grid(
  rdacca_pa_plot_list$TD_total$bar, rdacca_pa_plot_list$TD_repl$bar, rdacca_pa_plot_list$TD_rich$bar,
  rdacca_pa_plot_list$FD_total$bar, rdacca_pa_plot_list$FD_repl$bar, rdacca_pa_plot_list$FD_rich$bar,
  rdacca_pa_plot_list$PD_total$bar, rdacca_pa_plot_list$PD_repl$bar, rdacca_pa_plot_list$PD_rich$bar,
  ncol = 3, nrow = 3, align = "hv"); p_rdacca_pa_bar_all


# ggplot2::ggsave(filename = "08 rdacca_contribution_bars_PA.pdf", plot = p_rdacca_pa_bar_all, width = 18, height = 20, units = "cm")

pdf("08 rdacca hp PA built-in plots.pdf", width = 8, height = 8)
par(mfrow = c(3, 3), mar = c(3, 3, 3, 1))
plot(rdacca_pa_list$TD_total, plot.perc = TRUE); title("Taxonomic total")
plot(rdacca_pa_list$TD_repl, plot.perc = TRUE); title("Taxonomic replacement")
plot(rdacca_pa_list$TD_rich, plot.perc = TRUE); title("Taxonomic richness difference")
plot(rdacca_pa_list$FD_total, plot.perc = TRUE); title("Functional total")
plot(rdacca_pa_list$FD_repl, plot.perc = TRUE); title("Functional replacement")
plot(rdacca_pa_list$FD_rich, plot.perc = TRUE); title("Functional richness difference")
plot(rdacca_pa_list$PD_total, plot.perc = TRUE); title("Phylogenetic total")
plot(rdacca_pa_list$PD_repl, plot.perc = TRUE); title("Phylogenetic replacement")
plot(rdacca_pa_list$PD_rich, plot.perc = TRUE); title("Phylogenetic richness difference")
dev.off()





