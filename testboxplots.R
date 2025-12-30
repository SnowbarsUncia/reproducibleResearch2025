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

tbl <- read.xlsx("Test_table_2.xlsx", sheet = 2) # Прочтение таблицы с данными

str(tbl) # Проверка прочитанной таблицы
hist(tbl$PO.activity)
hist(tbl$GST.activity)
hist(tbl$CAT.activity)

# Построение графиков:

plot.PO <- ggplot(data=tbl, aes(x=Group, y=PO.activity, fill=Species)) +
  expand_limits(y=0) +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap(~Species) +
  ylab("PO activity, a.u.") +
  xlab("") +
  theme_bw(base_size = 16) +
  theme(strip.text = element_text(face="italic")) +
  scale_fill_manual(values=c("#D2AA6D", "forestgreen")) +
  geom_jitter (width = 0.1, show.legend = FALSE)
plot.PO # Отрисовка получившегося графика
ggsave("PO_with_stats.png", device=png, width=20, height=12, units="cm") # Сохранение полученного графика в рабочей директории
