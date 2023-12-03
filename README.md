# ALMetevsk
Анализ метагенома образцов нефти. Есть подозрение, что данные загрязнены фрагментами генома человека.
## Подключение
```
cd /home/alexandr/Downloads/063_annotator1 && sudo openvpn externalwork3-client.conf
ssh oxkolpakova@pbx3
source activate alm
export PATH=$PATH:/home/oxkolpakova/programs/miniconda3/envs/alm 

scp -r oxkolpakova@pbx3:/home/oxkolpakova/data/results/kraken2* /home/alexandr/Documents/ALM/data/results/kraken2
 
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
В этот раз используем только fastp с параметрами -q 20 -l 140 --trim_poly_g для очистки качества 20 и длинне ридов 140, так же обрежем поли G.

## kraken2
Для метагеномного анализа используем kraken2
базу данных скачаем в 
/srv/50f56420-22fa-4043-91a0-7d2a1709438f/oxkolpakova/kraken2_DB
с помощью скрипта kraken2-build

```
/home/oxkolpakova/programs/miniconda3/envs/alm/bin/kraken2-build --standard --max-db-size 68719476736 --threads 20 --db /srv/50f56420-22fa-4043-91a0-7d2a1709438f/oxkolpakova/kraken2_DB
```

Проанализируем с помощью стандартной базы kraken2

Не получилось загрузить базы стандартно из за ограничений по RAM, скачиваем прекомпиленные
https://benlangmead.github.io/aws-indexes/k2
```
wget -c https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20231009.tar.gz
tar -xzf k2_standard_20231009.tar.gz
```
kraken2 пробный запуск
```
read1='/home/oxkolpakova/data/results/fastp/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_10_L00_R1_P.fastq.gz'
read2='/home/oxkolpakova/data/results/fastp/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_10_L00_R2_P.fastq.gz'
out='/home/oxkolpakova/data/results/kraken2/kraken2_result.txt'
report='/home/oxkolpakova/data/results/kraken2/kraken2_report.txt'
database='/srv/50f56420-22fa-4043-91a0-7d2a1709438f/oxkolpakova/kraken2_DB'
NCPUS=20
mkdir -p "$out"

kraken2 \
--db $database \
--threads $NCPUS \
--report $report \
--report-zero-counts \
--use-names \
--memory-mapping \
--minimum-base-quality 20 \
--gzip-compressed \
--paired \
$read1 $read2 
```

Запускаем kraken2 с помощью nextflow
```
./BWAINDEX_FASTP_BWAMEM_KRAKEN2.nf with-report report.html -with-dag -with-singularity -resume
```
А теперь bracken

```
input='/home/oxkolpakova/data/results/kraken2/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_22_L00_kraken2_report.txt'
out='/home/oxkolpakova/data/results/bracken/'
database='/srv/50f56420-22fa-4043-91a0-7d2a1709438f/oxkolpakova/kraken2_DB'
bracken -d $database -i $input -o $out/bracken_result.txt -r 100 -l S
```

тестируем test.fn
./test.nf with-report report.html -with-dag -with-singularity -resume


## Galaxy https://usegalaxy.org/workflows/list
данные секвенирования от 13.11 (4 образца)

1. Fastp (-q 20, -l 50)
```
ln -s '/jetstream2/scratch/main/jobs/54067920/inputs/dataset_5061d16a-b100-43c1-bb38-68a30e1afe22.dat' '202311131707_210701004_A_neft_13_11_2023_3_L00_R1_fq_gz.fastq.gz' &&  ln -s '/jetstream2/scratch/main/jobs/54067920/inputs/dataset_c37a61dd-9ddf-4195-a377-8443e35189e3.dat' '202311131707_210701004_A_neft_13_11_2023_3_L00_R2_fq_gz_R2.fastq.gz' &&    fastp  --thread ${GALAXY_SLOTS:-1} --report_title 'fastp report for 202311131707_210701004_A_neft_13_11_2023_3_L00_R1_fq_gz.fastq.gz'   -i '202311131707_210701004_A_neft_13_11_2023_3_L00_R1_fq_gz.fastq.gz' -o first.fastq.gz  -I '202311131707_210701004_A_neft_13_11_2023_3_L00_R2_fq_gz_R2.fastq.gz' -O second.fastq.gz       --detect_adapter_for_pe                 -q 20      -l 50                     &&  mv first.fastq.gz '/jetstream2/scratch/main/jobs/54067920/outputs/dataset_b01c9a77-43bc-4a52-be24-9158b1ec63c9.dat' && mv second.fastq.gz '/jetstream2/scratch/main/jobs/54067920/outputs/dataset_ee300622-ac1c-49e4-a431-2b96854e77df.dat'
```
2. Bowtie2 (почистить от контаминации чеоовеком, не получилось запустить, ошибка)
3. Assembly with MEGAHIT
```
if [[ -n "$GALAXY_MEMORY_MB" ]]; then MEMORY="-m $((GALAXY_MEMORY_MB * 1024))"; fi;  megahit --num-cpu-threads ${GALAXY_SLOTS:-4} -1 '/scratch4/nekrut/galaxy/main/staging/54068902/inputs/dataset_f46a33e4-418c-449f-945b-8db8872ee1ea.dat' -2 '/scratch4/nekrut/galaxy/main/staging/54068902/inputs/dataset_b4664554-7b5a-4626-9b3b-eb9d1316ca13.dat' --min-count '2' --k-list '21,29,39,59,79,99,119,141'  --bubble-level '2' --merge-level '20,0.95' --prune-level '2' --prune-depth '2' --disconnect-ratio '0.1' --low-local-ratio '0.2' --cleaning-rounds '5'   --min-contig-len '200' $MEMORY && cat megahit_out/log
```
  
4. Assembly	metaSPAdes
```
mkdir -p paired_reads1 && ln -s '/jetstream2/scratch/main/jobs/54069093/inputs/dataset_23ac917e-2891-47e0-8b6a-aeeafbd3c343.dat' 'paired_reads1/fastp_on_data_2236_and_data_2235:_Read_1_output.fastq.gz' &&  ln -s '/jetstream2/scratch/main/jobs/54069093/inputs/dataset_43fdc7f1-bf77-44b6-9c3f-2cb7dc144808.dat' 'paired_reads1/fastp_on_data_2236_and_data_2235:_Read_2_output.fastq.gz' &&       export OMP_THREAD_LIMIT=${GALAXY_SLOTS:-4} &&  metaspades.py -o 'output'  -t ${GALAXY_SLOTS:-4} -m $((${GALAXY_MEMORY_MB:-8192}/1024))   --pe-1 1 'paired_reads1/fastp_on_data_2236_and_data_2235:_Read_1_output.fastq.gz' --pe-2 1 'paired_reads1/fastp_on_data_2236_and_data_2235:_Read_2_output.fastq.gz' --pe-or 1 fr
```

Prokka on Contigs
