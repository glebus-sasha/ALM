# ALMetevsk
Анализ метагенома в нефти. Есть подозрение, что данные загрязнены человеком.
## Подключение
```
cd /home/alexandr/Downloads/063_annotator1
sudo openvpn externalwork3-client.conf
ssh oxkolpakova@pbx3
source activate alm
export PATH=$PATH:/home/oxkolpakova/programs/miniconda3/envs/alm/bi
scp -r oxkolpakova@pbx3:/home/oxkolpakova/scripts/BWAINDEX_FASTP_BWAMEM.nf /home/alexandr/Documents/ALM/scripts
screen -XS <session-id> quit
```
## Загрузка референса и создание индекса

Загрузить не получилось( Пришлось скачать вручную.
```
datasets-cli download genome accession GCF_000001405.40 --include gff3,genome --filename GCF_000001405.40.zip

bwa mem index human.fna
```
## Качество
Предарительные анализ показал большое количество polyG по всех образцах, а так же адаптеры в некоторых образцах. 
Это пример команда для очистки и создания отчета fastp.
```
fastp -q 20 -l 50  --trim_poly_g --thread 12 -h /home/oxkolpakova/data/result/fastp/report_html/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00.html \
    --in1 /home/oxkolpakova/data/raw/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00_R1.fq.gz \
    --in2 /home/oxkolpakova/data/raw/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00_R2.fq.gz \
    --out1 /home/oxkolpakova/data/result/fastp/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00_R1.fq.gz \
    --out2 /home/oxkolpakova/data/result/fastp/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00_R2.fq.gz

```

Для выполнения был написан скрипт fastp_do.sh с помощью GPT.
После тримминга результат значительно улучшился.

## Выравнивание
Разведочный анализ показал загрязнение человеком. Поэтому выполним выравнивание на геном человека и отфильтруем невыравненные с помощью bwa mem и samtools.

Принт для GPT:
Я хочу выравнить мои файлы с парными ридами на геном человека с помощью bwa, можешь написать команду?

Путь до референса fna:
/home/oxkolpakova/data/references/human.fna

Путь до файла разметки gtf:
/home/oxkolpakova/data/references/genomic.gtf

Путь до первого рида R1:
/home/oxkolpakova/data/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00_R1.fq.gz

Путь до второго рида R2:
/home/oxkolpakova/data/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00_R2.fq.gz

Путь до выходных файлов:
/home/oxkolpakova/data/result/human

Путь до невыравненных ридов:
/home/oxkolpakova/data/result/unmapped

## bwa mem

```
bwa mem -t 12 /home/oxkolpakova/data/references/human.fna \
  /home/oxkolpakova/data/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00_R1.fq.gz \
  /home/oxkolpakova/data/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00_R2.fq.gz \
  | samtools view -Sb - \
  | samtools sort -o /home/oxkolpakova/data/result/bwa/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00.sorted.bam

samtools view -b -o /home/oxkolpakova/data/result/bwa/unmapped/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00.bam -f 4 \
  /home/oxkolpakova/data/result/bwa/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00.sam
```

Для запуска используем скрипт do_bwa.sh


## Статистика после bwa mem

Считаем количество ридов, которые не были выравнены на геном человека с помощью скрипта
read_count_after_bwa.sh
Результат в файле num_of_reads_after_bwa.txt говорит о том, что все образцы загрязнены человеком: уровень чистоты в диапазоне от 0,42% до 27,96%

## bam2fasta для очищенных от человека

```
for file in ./*.bam; do
    output=$(basename "$file" .bam).fasta
    samtools fasta "$file" > "$output"
done
```

## Анализ с помощью kraken2
Загружаем датабазы с помощью скрипта
databases.sh

```
./kraken2-build --build --threads n --db kraken2_db
```
_Дальше мы делаем твист_
_И переходим ко второй попытке_
_Делаем всё заново_

## nextflow

Для этого написан скрипт nextflow так же его можно опционально докеризировать.
Доступны 3 инструмента
BWAINDEX
FASTP
BWAMEM

```
./BWAINDEX_FASTP_BWAMEM.nf with-report report.html -with-dag -resume
```
В этот раз используем только fastp с параметрами -q 20 -l 140 --trim_poly_g для очистки качества 20 и длинне ридов 140, так же обрежим поли G.

## kraken2
Для метагеномного анализа используем kraken2
базу данных скачаем в 
/srv/50f56420-22fa-4043-91a0-7d2a1709438f/oxkolpakova/kraken2_DB
с помощью скрипта kraken2-build

```
/home/oxkolpakova/programs/miniconda3/envs/alm/bin/kraken2-build --standard --threads 20 --db /srv/50f56420-22fa-4043-91a0-7d2a1709438f/oxkolpakova/kraken2_DB
```

