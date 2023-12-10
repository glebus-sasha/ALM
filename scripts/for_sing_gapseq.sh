#!/bin/bash

# Полный путь к .sif файлу
sif_path='/home/oxkolpakova/programs/gapseq_latest.sif'

# Директория для вывода логов
log_dir='/home/oxkolpakova/ALM/data/results/meta_SPAdes/logs/'

# Функция для обработки одного файла
process_file() {
    local input_file=$1
    local output_dir="/home/oxkolpakova/ALM/data/results/meta_SPAdes/"
    
    # Имя файла без расширения
    xbase=$(basename $input_file)
    shname=${xbase%.*}
    
    # Журнал для вывода
    log_file="$log_dir/$shname.log"

    # Запустить gapseq команды внутри контейнера и логгировать
    singularity run -B $output_dir:/mnt $sif_path gapseq find -p all -b 100 -m Bacteria /mnt/$xbase > $log_file 2>&1
    singularity run -B $output_dir:/mnt $sif_path gapseq find-transport -b 100 /mnt/$xbase >> $log_file 2>&1
    singularity run -B $output_dir:/mnt $sif_path gapseq draft -r /mnt/$shname-all-Reactions.tbl -t /mnt/$shname-Transporter.tbl -p /mnt/$shname-all-Pathways.tbl -c /mnt/$xbase -u 100 -l 50 >> $log_file 2>&1
    singularity run -B $output_dir:/mnt $sif_path gapseq fill -m /mnt/$shname-draft.RDS -n /opt/gapseq/dat/media/meerwasser.csv -c /mnt/$shname-rxnWeights.RDS -b 100 -g /mnt/$shname-rxnXgenes.RDS >> $log_file 2>&1

    echo "Обработка файла $input_file завершена. Журнал: $log_file"
}

# Массив файлов для обработки
files=('/home/oxkolpakova/ALM/data/results/meta_SPAdes/metaSPAdes_S1_Contigs.fasta'
       '/home/oxkolpakova/ALM/data/results/meta_SPAdes/metaSPAdes_S2_Contigs.fasta'
       '/home/oxkolpakova/ALM/data/results/meta_SPAdes/metaSPAdes_S3_Contigs.fasta'
       '/home/oxkolpakova/ALM/data/results/meta_SPAdes/metaSPAdes_S4_Contigs.fasta')

# Создать директорию для логов, если она не существует
mkdir -p $log_dir

# Запустить обработку каждого файла в фоне
for file in "${files[@]}"; do
    process_file "$file" &
done

# Ожидать завершения всех фоновых процессов
wait

echo "Все файлы обработаны."
#!/bin/bash
gapseq find -p all -b 200 -m Bacteria /opt/static/metaSPAdes_S4_Contigs.fasta
filename='/opt/static/tren.fasta'
gapseq find -p all -b 100 -m Bacteria $filename
gapseq find-transport -b 100 $filename
gapseq draft -r $shname-all-Reactions.tbl -t $shname-Transporter.tbl -p $shname-all-Pathways.tbl -c $filename -u 100 -l 50
gapseq fill -m $shname-draft.RDS -n /opt/gapseq/dat/media/meerwasser.csv -c $shname-rxnWeights.RDS -b 100 -g $shname-rxnXgenes.RDS

