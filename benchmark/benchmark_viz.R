library(tidyverse)
library(ggsignif)

# data ----
filenames <- c("reference_true", "reference_dash", "reference_total",
               "learn_true", "learn_dash", "learn_gray", "learn_red", "learn_total")
for (filename in filenames){
  tmp <- read_csv(paste0("~/ciber/benchmark/", filename, ".csv"))
  tmp <- tmp %>% column_to_rownames(., var = "Algorithm")
  assign(filename, tmp)
}; rm(tmp, filename, filenames)

# reference_true+reference_dash-reference_total
# learn_true+learn_dash+learn_gray+learn_red-learn_total

recall <- round(learn_true/reference_true, 2)
precision <- round(learn_true/learn_total, 2)
relaxed_precision <- round((learn_true+learn_gray)/learn_total, 2)
recall_wd <- round((learn_true+learn_dash)/(reference_true+reference_dash), 2)
precision_wd <- round((learn_true+learn_dash)/learn_total, 2)
relaxed_precision_wd <- round((learn_true+learn_gray+learn_dash)/learn_total, 2)

df <- recall %>% rownames_to_column(., var = "Algorithm") %>% 
  pivot_longer(., cols = colnames(.)[-1])
colnames(df) <- c("Algorithm", "dataset", "recall")
# df <- read.csv("./benchmark_metrics.csv")
df$Algorithm <- factor(df$Algorithm, levels = c("CIBER","MMHC","RSMAX2","H2PC"))
df$dataset <- factor(df$dataset, levels = c("HMMA","HMR","HHR","HHATAC-peak","HHATAC-TF","IntestD","TcellD"))

df$precision <- precision %>% rownames_to_column(., var = "Algorithm") %>% 
  pivot_longer(., cols = colnames(.)[-1]) %>% .[["value"]]
df$relaxed_precision <- relaxed_precision %>% rownames_to_column(., var = "Algorithm") %>% 
  pivot_longer(., cols = colnames(.)[-1]) %>% .[["value"]]

df$recall_wd <- recall_wd %>% rownames_to_column(., var = "Algorithm") %>% 
  pivot_longer(., cols = colnames(.)[-1]) %>% .[["value"]]
df$precision_wd <- precision_wd %>% rownames_to_column(., var = "Algorithm") %>% 
  pivot_longer(., cols = colnames(.)[-1]) %>% .[["value"]]
df$relaxed_precision_wd <- relaxed_precision_wd %>% rownames_to_column(., var = "Algorithm") %>% 
  pivot_longer(., cols = colnames(.)[-1]) %>% .[["value"]]

rm(list = ls(pattern = "reference")); rm(list = ls(pattern = "learn"))
rm(recall, precision, relaxed_precision,
   recall_wd, precision_wd, relaxed_precision_wd)

df_l <- df %>% pivot_longer(., cols = colnames(.)[-c(1:2)])
colnames(df_l) <- gsub("name", "metric", colnames(df_l))
df_l$metric <- factor(df_l$metric, levels = c("recall", "precision", "relaxed_precision",
                                             "recall_wd", "precision_wd", "relaxed_precision_wd"))
write.csv(df, file = "./benchmark_metrics2.csv", row.names = FALSE, quote = FALSE)

# data summary ----
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

data_summary(df_l, "value", groupnames = c("Algorithm", "metric"))
# df_l %>% group_by(Algorithm, metric) %>% summarise_at(vars(value), list(mean = mean, sd=sd))


# plot ----
df %>%  ggplot(., aes(y=recall, x=precision, col=Algorithm, shape=dataset)) +
  # scale_shape_manual(values = c(0,1,2,3,4,5,15,16,17)) +
  # scale_shape_manual(values = c(21:25,3,4)) +
  scale_shape_manual(values = c(24,25,23,21,22,3,4)) +
  geom_point(position=position_jitter(h=0.005, w=0.005), alpha = 0.8, stroke=1, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(0,1)

df %>% ggplot(., aes(y=recall, x=relaxed_precision, col=Algorithm, shape=dataset, fill=Algorithm)) +
  scale_shape_manual(values = c(24,25,23,21,22,10,11)) +
  # scale_shape_manual(values = c(0,1,2,3,4,5,15,16,17)) +
  geom_point(position=position_jitter(h=0.01, w=0.01), alpha = 0.8, stroke=1, size=3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(0,1)

df_l %>% 
  subset(., metric %in% c("recall", "precision", "relaxed_precision")) %>%
  # subset(., metric %in% c("recall_wd", "precision_wd", "relaxed_precision_wd")) %>%
  ggplot(., aes(x=Algorithm, y=value, fill=Algorithm)) +
  stat_summary(fun=mean, geom="bar", col="black", width=0.8) +
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", width=.2) +
  geom_point(position = position_jitter(width = 0.3), color = "black", size = 0.6) + 
  geom_signif(comparisons=list(c("CIBER", "MMHC")),
              y_position = 1.04, tip_length = 0.02, vjust=0.5, map_signif_level = T) +
  geom_signif(comparisons=list(c("CIBER", "RSMAX2")),
              y_position = 1.09, tip_length = 0.02, vjust=0.5, map_signif_level = T) +
  geom_signif(comparisons=list(c("CIBER", "H2PC")),
              y_position = 1.14, tip_length = 0.02, vjust=0.5, map_signif_level = T) +
  ggtitle("Statistics without unconventional edges (dashed lines)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.3)) + 
  facet_wrap( . ~ metric, ncol = 3)
