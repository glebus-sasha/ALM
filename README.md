# ALMetevsk

## Подключение
```{}
cd /home/alexandr/Downloads/063_annotator1
sudo openvpn externalwork3-client.conf
ssh oxkolpakova@pbx3
source activate alm
export PATH=$PATH:/home/oxkolpakova/programs/miniconda3/envs/alm/bin
scp -r /home/alexandr/Documents/ALM/data/raw/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_[5-6]_* oxkolpakova@pbx3:/home/oxkolpakova/data/raw
scp -r oxkolpakova@pbx3:/home/oxkolpakova/data/result/fastqc*.html /home/alexandr/Documents/ALM/data/results/fastqc_after
```
## Загрузка референса и создание индекса

Загрузить не получилось( Пришлось скачать вручную.
```
datasets-cli download genome accession GCF_000001405.40 --include gff3,genome --filename GCF_000001405.40.zip

bwa mem index human.fna
```
## Качество

Это пример команда для очистки и создания отчета fastp.
```
fastp -q 20 -l 50  --trim_poly_g --thread 12 -h /home/oxkolpakova/data/result/fastp/report_html/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00.html \
    --in1 /home/oxkolpakova/data/raw/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00_R1.fq.gz \
    --in2 /home/oxkolpakova/data/raw/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00_R2.fq.gz \
    --out1 /home/oxkolpakova/data/result/fastp/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00_R1.fq.gz \
    --out2 /home/oxkolpakova/data/result/fastp/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_6_L00_R2.fq.gz

```

Для выполнения был написан скрипт fastp_do.sh с помощью GPT.

## Выравнивание

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
```
bwa mem -t 12 /home/oxkolpakova/data/references/human.fna /home/oxkolpakova/data/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00_R1.fq.gz /home/oxkolpakova/data/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00_R2.fq.gz > /home/oxkolpakova/data/result/human/alignment.sam
```
## Сохраняем невыравненные
```
samtools view -Sb alignment.sam | samtools sort -o alignment.sorted.bam
samtools view -b -o /home/oxkolpakova/data/result/unmapped/unmapped.bam -f 4 /home/oxkolpakova/data/result/human/alignment.sam
samtools index samtools index alignment.sorted.bam 
samtools idxstats alignment.sorted.bam 
```

## umgap
```
umgap-analyse.sh -1 202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00_R1.fq.gz \
-2 /home/oxkolpakova/data/202309251627_220601009_2P230329071US2S2721BX_B_neft250923_1_L00_R2.fq.gz -o home/oxkolpakova/data/result/umgap/output.fa

```

## bwa mem

Привет кода для bwa mem
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
Можно ли упростить этот код?




