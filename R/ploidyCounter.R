library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
db <- args[1]
cnf <- "--defaults-file=/usr/local/share/stacks/sql/mysql.cnf"

mypipe <- pipe(paste("mysql", cnf, "-N -e", "'select * from alleles where length(allele)<2'", db))
data  <- read.table(mypipe)
colnames(data) <- c("id", "sample_id", "tag_id", "allele", "read_pct", "read_cnt")

mypipe  <- pipe(paste("mysql", cnf, "-N -e", "'select sample_id,file from samples'", db))
samples <- read.table(mypipe)
colnames(samples) <- c("id", "file")

data$depth <- round(data$read_cnt * 100 / data$read_pct)
data2 <- data[data$depth > 30 & data$depth < 50,]

for (s in 1:nrow(samples)) {
	sample <- samples[s,]
	svg(paste0("density/", sample$file, ".svg"))
	print(qplot(read_pct, colour = factor(sample_id), data = data2[data2$sample_id == sample$id,], geom="density", binwidth = 1) + theme(legend.position="none"))
	invisible(dev.off())
	svg(paste0("histo/", sample$file, ".svg"))
	print(qplot(read_pct, colour = factor(sample_id), data = data2[data2$sample_id == sample$id,], geom="histogram", binwidth = 1) + theme(legend.position="none"))
	invisible(dev.off())
}
