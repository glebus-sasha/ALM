#!/bin/bash

input_dir="/home/oxkolpakova/data/result/fastp/"
reference="/home/oxkolpakova/data/references/human.fna"
output_dir="/home/oxkolpakova/data/result/bwa/"
unmapped_dir="/home/oxkolpakova/data/result/bwa/unmapped/"

# Перейдите в папку результатов BWA, если она не существует
mkdir -p "$output_dir"
mkdir -p "$unmapped_dir"

for R1_file in "$input_dir"/*_R1.fq.gz; do
  R2_file="${R1_file%_R1.fq.gz}_R2.fq.gz"

  # Создайте имя для выходного файла на основе имени R1
  output_prefix=$(basename "$R1_file" _R1.fq.gz)

  # Выполните bwa mem и последующие операции для пары файлов
  bwa mem -t 12 "$reference" "$R1_file" "$R2_file" \
    | samtools view -Sb - \
    | samtools sort -o "$output_dir${output_prefix}_sorted.bam"

  samtools view -b -o "$unmapped_dir${output_prefix}.bam" -f 4 "$output_dir${output_prefix}_sorted.bam"
done

