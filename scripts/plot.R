print("starting")

library(data.table)
library(ggplot2)
library(purrr)
library(tidyr)
library(dplyr)
library(tibble)
library(magrittr)


if (exists("snakemake")){
    pileup_path <- snakemake@input$pileup
    png_path <- snakemake@output$png
    start_idx <- snakemake@params$start
    end_idx <- snakemake@params$end
}else{
    setwd("/home/vy/Documents/Sarah")
    pileup_path <- "mpileup/FAW82928_pass_barcode01_2687e1a7_f06207b4_1.mpileup"
    png_path <- "plots/barplot.png"
}
start_idx <- 1
end_idx <- 3000

print("loading pileup")
pileup <- fread(pileup_path) %>% 
    set_colnames(c("name",
                   "pos",
                   "ref",
                   "n_reads",
                   "matches",
                   "qualities",
                   "mappingQualities"))

print("processing pileup")
melted_pileup <- pileup %>% 
    select(pos, qualities, ref) %>% 
    melt(id=c("pos", "ref")) %>% 
    select(pos, ref, value) %>% 
    set_colnames(c("pos", "ref", "quality"))

string_to_column <- function(df_row){
    convert_qual <- function(x){
        x %>% charToRaw %>%
            as.numeric %>%
            return
    }
    quality <- df_row$quality %>% strsplit("") %>% unlist
    numericQuality <- (quality %>% sapply(convert_qual))
    data.table(position=rep(df_row$pos), base=rep(df_row$ref), quality, numericQuality) %>%
        return
}

final_df <- list()
for (row in 1:nrow(melted_pileup)){
    final_df[[row]] <- string_to_column(melted_pileup[row])
}
plottable_pileup <- final_df %>% rbindlist


print("plotting pileup")
every_nth = function(n) {
  return(function(x) {x[c(TRUE, rep(FALSE, n - 1))]})
}
filtered_pileup <- plottable_pileup %>% 
              filter(position>=start_idx) %>% 
              filter(position<end_idx)


# boxplot
#plt <- ggplot(filtered_pileup, aes(x=factor(position))) + 
#    geom_boxplot(size=0.1, outlier.size=0.1, aes(y=numericQuality)) + 
#    theme_bw() + 
#    theme(axis.text=element_text(size=3),
#          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#    xlab("position") +
#    ylab("score") + 
#    scale_x_discrete(breaks=every_nth(n=10))
#ggsave(png_path, plt, height=3, width=250, units='cm', limitsize=FALSE)



mean_pileup <- plottable_pileup[, list(mean=mean(numericQuality)),by=position] %>% 
    merge(melted_pileup %>% select(pos,ref), 
          by.x='position', by.y="pos")

plt <- ggplot(mean_pileup, 
              aes(x=factor(position), y=mean)) + 
    geom_point(size=0.1) + 
    theme_bw() + 
    xlab("position") +
    ylab("score") +
    scale_x_discrete(breaks=every_nth(n=10)) +
    theme(axis.text=element_text(size=4),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
ggsave(png_path, plt, height=5, width=250, units='cm', limitsize=FALSE)

#y <- apply(x, 2, string_to_column)
#
#x[5] %>% string_to_column
#
#
#x$V7 <- lapply(x$V6, convert_qual) %>% rbindlist
#
#x %>% mutate(V7=as.numeric(charToRaw(V6)))
