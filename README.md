# Воспроизведение анализа статьи "Transcriptional response of Saccharomyces cerevisiae to lactic acid enantiomers"

## Резюме

В оригинальной статье рассматривается транскрипционная реакция дрожжей _Saccharomyces cerevisiae_ на энантиомеры молочной кислоты (лактата), то есть изменения транскриптома дрожжей в ответ на введение энантиомеров лактата. В данном воспроизведении вышеупомянутого исследования будет рассматриваться лишь реакция _Saccharomyces cerevisiae_ на D-лактата в концентрации 5 mM.

## Методы

### Использованные программы:

`R v4.3.2`;
`PuTTY v0.81`;
`NCBI SRA Toolkit v3.0.0`;
`Subread v2.0.4`

### Пакеты R, использованные в ходе работы:

`openxlsx`;
`ggplot2`;
`ggpubr`;
`scales`;
`BiocManager`;
`EnhancedVolcano`;
`DESeq2`

### В процессе работы на удалённых серверах bash, использовались следующие команды:

`pwd` (print working directory) - отобразить текущее расположение

`cd` (change directory) - перейти в другое расположение

`mkdir` (make directory) - Создать пустую папку в текущем расположении

`ls` (list) - просмотреть содержимое папки (по умолчанию - той, что является текущим расположением); `ls -lh` - для отображения содержимого определённой папки, например, `ls -ls /home/user/test`

`head`, `tail` - Просмотреть начало или конец файла; `head -n` или `tail -n` - для отображения конкретного количества строк сначала или с конца файла, например, `tail -n 20 file.txt`

`cat` - Просмотреть файл целиком

`grep` - Поиск текста в файле, например, `grep "seq"`

`|` - Передача вывода одной команды другой команде (pipeline)

`>` - Запись вывода команды в файл, например, `> file.txt`

`mv`, `cp` (move, copy) - Переместить или переименовать; копировать. `mv -r` и `cp -r` - для работы с конкретной папкой

`rm` (remove) - Удалить; `rm -r` - для удаления папки

`sudo` - Для запуска команд от имени администратора

`screen`- Создание экранов для запуска длительных процессов; screen `-S "ИМЯ"` - Для создания экрана; Сочетание клавиш Ctrl+A+D - для выхода из экрана; `screen -r "ИМЯ"` - для возврата к существующему экрану; `screen -ls` - для показа полного списка существующих экранов

`wget` - Для скачивания файлов из интернета на сервер с указанием гиперссылки, например `wget https://example.com/file.zip`

`ps` (process status) - Для отображения списка активных процессов; `ps aux` - для отображения полного списка процессов с подробной информацией о них

`kill` - Прервать активный процесс с указанием ID этого процесса, например, `kill 1234`


## Скачивание исходных данных:

Сначала нужно указать, где находится исполняемый файл программы sratoolkit:
`export PATH=$PATH:/media/secondary/apps/sratoolkit.3.0.0-ubuntu64/bin/`

Затем, можно начинать, собственно, скачивание данных:
`fasterq-dump --threads 2 -A --progress SRR24466389; fasterq-dump --threads 2 -A --progress SRR24466390; fasterq-dump --threads 2 -A --progress SRR24466391; fasterq-dump --threads 2 -A --progress SRR24466380; fasterq-dump --threads 2 -A --progress SRR24466381; fasterq-dump --threads 2 -A --progress SRR24466382`

## Выравнивание первичных прочтений на референс:

### Скачивание референсной последовательности:
`wget https://ftp.ensembl.org/pub/release-108/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz`;
`wget https://ftp.ensembl.org/pub/release-108/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.tople
vel.fa.gz`

### Распаковка архивов:
`gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz`;
`gunzip Saccharomyces_cerevisiae.R64-1-1.108.gtf.gz`;

### Построение индекса и подготовка файла с данными сплайсинга в hisat2:
`hisat2-build Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa yeast_index`;
`hisat2_extract_splice_sites.py Saccharomyces_cerevisiae.R64-1-1.108.gtf > yeast_splice_sites.txt`

### Выравнивание с помощью hisat2 и сортировка bam-файла с помощью samtools:
`for sample in `ls *_1.fastq`; do base=$(basename $sample "_1.fastq"); hisat2 -x yeast_index --known-splicesite-infile yeast_splice_sites.txt -p 8 -1 ${base}_1.fastq -2 ${base}_2.fastq | samtools view --threads 2 -bS | samtools sort --threads 2 -o $base.bam; done`

## Построение графиков в среде R

___

## Анализ дифференциальной экспрессии генов

### Указываем где находится исполняемый файл subread:

`export PATH=$PATH:/media/secondary/apps/subread-2.0.4-Linux-x86_64/bin`

### Подсчёт экспрессии:

`featureCounts -s 2 -T 2 -p -a Saccharomyces_cerevisiae.R64-1-1.108.gtf \-o allSamples.featureCounts.txt $(ls *.bam)`

### Скачивание данных с удалённого сервера для их дальнейшей визуализации в R (Windows OS):

`.\pscp -P 627 2025_RR_St1@bioinformatics.isu.ru:~/allSamples.featureCounts.txt`

### Анализ дифференциальной экспрессии в R:

___



