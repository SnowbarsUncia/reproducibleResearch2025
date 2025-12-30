setwd ("C:\\Users\\User\\Desktop\\VosprIssled") # Установка рабочей директории

# Загрузка необходимых пакетов:
if (!("openxlsx" %in% installed.packages())) install.packages("openxlsx")
library (openxlsx)
if (!("ggplot2" %in% installed.packages())) install.packages("ggplot2")
library (ggplot2)
if (!("ggpubr" %in% installed.packages())) install.packages("ggpubr")
library(ggpubr)
if (!("scales" %in% installed.packages())) install.packages("scales")
library(scales)

# Построение графиков:

yeast_growth_inhibition <- read.xlsx("Table_yeast_growth_DLA0500.xlsx", sheet = 1)
ggplot(yeast_growth_inhibition, aes(x=Condition, y=RGR2, fill=Condition)) +
  geom_boxplot(show.legend = FALSE) +
  geom_pwc(method="wilcox_test", label = "p.adj") +
  scale_fill_manual(values = c("forestgreen", "#ff1")) +
  scale_y_continuous(labels=comma_format(decimal.mark=",")) +
  geom_jitter (width = 0.1, show.legend = FALSE) +
  geom_pwc(method = "wilcox_test", label="p.adj") +
  ylab("Relative growth rate") +
  xlab("") +
  theme_bw(base_size = 16)
ggsave("yeast_pH_compensation", device=png, width=20, height=12, units="cm") # Сохранение получившихся графиков в рабочей директории