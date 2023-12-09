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
  
2.1 Kraken2(-confidence '0.1', db Standard)

```
kraken2 --threads ${GALAXY_SLOTS:-1} --db '/cvmfs/data.galaxyproject.org/managed/kraken2_databases/k2_standard_20210517'    --paired '/scratch4/nekrut/galaxy/main/staging/54212028/inputs/dataset_84304d8e-b002-4828-9009-3bb90db5e19a.dat' '/scratch4/nekrut/galaxy/main/staging/54212028/inputs/dataset_865cb6ff-d23a-4553-80e2-02247ab57ebb.dat'   --confidence '0.1' --minimum-base-quality '0' --minimum-hit-groups '2'    --report '/scratch4/nekrut/galaxy/main/staging/54212028/outputs/dataset_76b004a6-ab6b-49c3-9684-f03229aedf0c.dat'     > '/scratch4/nekrut/galaxy/main/staging/54212028/outputs/dataset_71208db0-67d4-4602-9fe2-9b16d6d35da5.dat'
```

2.1.1 Kraken taxonomic report

```
ln -s "/jetstream2/scratch/main/jobs/54212148/inputs/dataset_76b004a6-ab6b-49c3-9684-f03229aedf0c.dat" "Report: Kraken2 on data 2744 and data 2743" && ln -s "/jetstream2/scratch/main/jobs/54212148/inputs/dataset_38ca50e0-a39e-407f-8199-be23bede1781.dat" "Report: Kraken2 on data 2748 and data 2747" && ln -s "/jetstream2/scratch/main/jobs/54212148/inputs/dataset_d4d54f55-52a9-48fe-bf8d-f11ab8369dd5.dat" "Report: Kraken2 on data 2754 and data 2753" && ln -s "/jetstream2/scratch/main/jobs/54212148/inputs/dataset_b08cc66c-7d5d-409a-b9a8-854dd189dbf1.dat" "Report: Kraken2 on data 2758 and data 2757" &&  export KRAKEN_DB_PATH='/cvmfs/data.galaxyproject.org/managed/kraken_database/bacteria' && python '/jetstream2/scratch/main/jobs/54212148/tool_files/kraken_taxonomy_report.py'  --db 'Bacteria'  --header-line    --intermediate --sanitize-names   --output '/jetstream2/scratch/main/jobs/54212148/outputs/dataset_c7a9e954-a06d-493b-9220-45eeaa24a24b.dat'   'Report: Kraken2 on data 2744 and data 2743' 'Report: Kraken2 on data 2748 and data 2747' 'Report: Kraken2 on data 2754 and data 2753' 'Report: Kraken2 on data 2758 and data 2757'
```

2.2. BWA-MEM2

```
set -o | grep -q pipefail && set -o pipefail;  ln -s '/jetstream2/scratch/main/jobs/54212240/inputs/dataset_c1887965-2ff3-44ce-b942-b13563a87253.dat' 'localref.fa' && bwa-mem2 index 'localref.fa' &&    bwa-mem2 mem -t "${GALAXY_SLOTS:-1}" -v 1                 'localref.fa' '/jetstream2/scratch/main/jobs/54212240/inputs/dataset_7adb0438-470d-4694-abf4-8beca62555b5.dat' '/jetstream2/scratch/main/jobs/54212240/inputs/dataset_f78b0059-53d6-46ab-9af2-715b33b2c05a.dat'  | samtools sort -@${GALAXY_SLOTS:-2} -T "${TMPDIR:-.}" -O bam -o '/jetstream2/scratch/main/jobs/54212240/outputs/dataset_3fd22041-5f27-4e4b-a6ec-f68bd54a2abf.dat'
```

2.3 Assembly	metaSPAdes
```
mkdir -p paired_reads1 && ln -s '/jetstream2/scratch/main/jobs/54212306/inputs/dataset_b433568d-2829-4be3-b8f6-d4b0afec0d5e.dat' 'paired_reads1/fastP_S4_R1.fastq.gz.fastq.gz' &&  ln -s '/jetstream2/scratch/main/jobs/54212306/inputs/dataset_aa3beb4a-b89a-4d14-9d1f-86b5cf0b3d7d.dat' 'paired_reads1/fastP_S4_R2.fastq.gz.fastq.gz' &&       export OMP_THREAD_LIMIT=${GALAXY_SLOTS:-4} &&  metaspades.py -o 'output'  -t ${GALAXY_SLOTS:-4} -m $((${GALAXY_MEMORY_MB:-8192}/1024))   --pe-1 1 'paired_reads1/fastP_S4_R1.fastq.gz.fastq.gz' --pe-2 1 'paired_reads1/fastP_S4_R2.fastq.gz.fastq.gz' --pe-or 1 fr
```

Prokka on Contigs


## krona after kraken2 (Galaxy)

Перед запуском крона, нужно объединить файлы с помощью

```
combine_kreports -r Report_Kraken2_S1* -0 Galaxy.kreport
```

Затем строим krona

```
./kreport2krona.py -r Galaxy.kreport -o Galaxy.krona.txt --no-intermediate-ranks
ktImportText Galaxy.krona.txt -o Galaxy.krona.html
```


##GAPSEQ with nextflow

Не получилось отладить
```
./gapseq.nf with-report report.html -with-dag -resume -with-docker
./gapseq.nf -with-docker
```


Поэтому используем докер на компе
```
sudo docker run -it --rm --name gapseq_example1 -h gapseq_container -v /media/alexandr/KINGSTON/meta_SPAdes/:/opt/static cdiener/gapseq /bin/bash 
cd opt/static/
./for_gapseq_1.sh
```

```
sudo docker run -it --rm --name gapseq_example2 -h gapseq_container -v /media/alexandr/KINGSTON/meta_SPAdes/:/opt/static cdiener/gapseq /bin/bash 
cd opt/static/
./for_gapseq_2.sh
```

```
sudo docker run -it --rm --name gapseq_example3 -h gapseq_container -v /media/alexandr/KINGSTON/meta_SPAdes/:/opt/static cdiener/gapseq /bin/bash 
cd opt/static/
./for_gapseq_3.sh
```

```
sudo docker run -it --rm --name gapseq_example4 -h gapseq_container -v /media/alexandr/KINGSTON/meta_SPAdes/:/opt/static cdiener/gapseq /bin/bash 
cd opt/static/
./for_gapseq_4.sh
```
