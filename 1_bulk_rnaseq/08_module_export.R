#!/usr/bin/env Rscript

# Purpose: collect the GO terms of the up- and down-regulated modules into one
#          workbook, one sheet per module.
# Input:   results/deg_go/Aged_vs_Middle_GO_enrichment_pvalue.csv
# Output:  results/modules/Up_Down_Group_GO_Results.xlsx
# Usage:   Rscript 08_module_export.R
#
# The module definitions are fixed: they group the enriched GO terms reported
# in the manuscript into G1-G6.

suppressPackageStartupMessages({
    library(openxlsx)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
    stop("cannot locate this script, run it as: Rscript 08_module_export.R")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1])))
setwd(script_dir)
source("config.R")

go_file <- file.path(DIR_RESULTS, "deg_go",
                     paste0(DEG_TEST, "_vs_", DEG_CONTROL,
                            "_GO_enrichment_pvalue.csv"))
output_dir <- file.path(DIR_RESULTS, "modules")
output_file <- file.path(output_dir, "Up_Down_Group_GO_Results.xlsx")

if (!file.exists(go_file)) {
    stop("GO result not found: ", go_file,
         "\n  run 07_deg_go.R first")
}

up_pathways <- list(
    "UP_G1_Ion Transport & Excitability" = c(
        "inorganic cation import across plasma membrane",
        "inorganic ion import across plasma membrane",
        "import across plasma membrane",
        "import into cell",
        "potassium ion transmembrane transport",
        "potassium ion transport",
        "potassium ion import across plasma membrane",
        "sodium ion transmembrane transport",
        "sodium ion transport",
        "multicellular organismal signaling",
        "cardiac conduction",
        "heart process",
        "action potential",
        "regulation of action potential",
        "cardiac muscle cell action potential"
    ),
    "UP_G2_Calcium-Triggered Secretion" = c(
        "regulation of catecholamine secretion",
        "catecholamine secretion",
        "regulation of amine transport",
        "amine transport",
        "regulation of calcium ion-dependent exocytosis",
        "response to calcium ion",
        "cellular response to calcium ion",
        "cellular response to metal ion",
        "cellular response to inorganic substance",
        "response to metal ion",
        "regulation of regulated secretory pathway",
        "calcium-ion regulated exocytosis",
        "regulation of exocytosis",
        "regulated exocytosis",
        "exocytosis",
        "signal release",
        "neurotransmitter secretion",
        "signal release from synapse",
        "neurotransmitter transport",
        "vesicle-mediated transport in synapse",
        "synaptic vesicle cycle"
    ),
    "UP_G3_Immune Activation & Cytotox" = c(
        "leukocyte activation involved in immune response",
        "cell activation involved in immune response",
        "leukocyte degranulation",
        "natural killer cell mediated immunity",
        "natural killer cell degranulation",
        "leukocyte mediated cytotoxicity",
        "cell killing",
        "lymphocyte differentiation",
        "T cell differentiation",
        "mononuclear cell differentiation",
        "alpha-beta T cell activation",
        "T cell differentiation in thymus",
        "immune response-regulating cell surface receptor signaling pathway",
        "T cell receptor signaling pathway",
        "T cell selection",
        "positive T cell selection",
        "positive regulation of leukocyte activation",
        "positive regulation of T cell differentiation in thymus"
    ),
    "UP_G5_Chemokine Recruitment & IFN-g" = c(
        "neutrophil chemotaxis",
        "granulocyte chemotaxis",
        "lymphocyte chemotaxis",
        "response to interferon-gamma"
    )
)

down_pathways <- list(
    "DN_G4_Signal Reception HW Decline" = c(
        "sensory organ morphogenesis",
        "camera-type eye morphogenesis",
        "eye morphogenesis",
        "sensory system development",
        "visual system development",
        "camera-type eye development",
        "eye development",
        "embryonic organ development"
    ),
    "DN_G5_Guidance & Quiescence Loss" = c(
        "axon guidance",
        "neuron projection guidance",
        "axon development",
        "negative regulation of cell development",
        "negative regulation of developmental growth",
        "negative regulation of growth",
        "negative regulation of cell growth",
        "regulation of axon extension involved in axon guidance",
        "semaphorin-plexin signaling pathway"
    ),
    "DN_G6_ECM Assembly Collapse" = c(
        "external encapsulating structure organization",
        "extracellular matrix organization",
        "extracellular structure organization",
        "collagen fibril organization"
    )
)

selected_pathways <- c(up_pathways, down_pathways)

go_raw <- read.csv(go_file, header = TRUE, stringsAsFactors = FALSE)
cols_keep <- c("Ontology", "ID", "Description", "GeneRatio", "BgRatio",
               "pvalue", "p.adjust", "qvalue", "Count", "geneID")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

wb <- createWorkbook()
header_style <- createStyle(textDecoration = "bold", halign = "center",
                            fgFill = "#F2F2F2", border = "bottom")

for (gname in names(selected_pathways)) {
    g_df <- go_raw[go_raw$Description %in% selected_pathways[[gname]],
                   cols_keep, drop = FALSE]
    g_df <- g_df[order(g_df$p.adjust), ]
    g_df$Direction <- if (grepl("^UP_", gname)) "Up" else "Down"

    sheet_name <- substr(gname, 1, 31)
    addWorksheet(wb, sheetName = sheet_name, gridLines = TRUE)
    writeData(wb, sheet = sheet_name, x = g_df)
    addStyle(wb, sheet_name, header_style, rows = 1,
             cols = seq_len(ncol(g_df)), gridExpand = TRUE)
    setColWidths(wb, sheet_name, cols = seq_len(ncol(g_df)), widths = "auto")
}

saveWorkbook(wb, output_file, overwrite = TRUE)
message("wrote ", output_file)
