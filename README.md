# Germinativo_GATK
Pipeline em Python e Bash para análise de dados germinativos utilizando a ferramenta GATK. O repositório irá abordar comandos básicos para a sequência do passo-a-passo

## Neste pipeline utilizaremos o google colab como nosso ambiente de programação.

# Parte 1 - Revisão dos comandos do Linux

Atenção: Para executar comandos Linux no Jupyter Notebook, é necessário utilizar uma exclamação (!) no início da linha ou colocar a linha `%%bash` no começo da célula.

O comando man permite acessar o manual referente a um comando do terminal. Esses manuais podem ser muito técnicos para o usuário iniciante, mas são bem convenientes para conferir os parâmetros que cada comando aceita.

Inicialmente iremos realizar o link do nosso ambiente teste para o drive para que possamos utilizar os arquivos que lá estão:

```pyton
from google.colab import drive
drive.mount('/content/drive')
```

Para mostrar qual é o diretório atual, usa-se o comando `pwd`:

```bash
%%bash
pwd
```
Para criar um diretório chamado teste, existe o comando `mkdir`:

```bash
%%bash
mkdir /content/drive/Shareddrives/T4-2022/LeonardoMedeiros/Agosto2022/logs
```
Confira que o diretório foi criado:

```bash
%%bash
ls -l /content/drive/Shareddrives/T4-2022/LeonardoMedeiros/Agosto2022
```

Nota: Para nos mantermos no diretório de interesse sem precisar entrar nele a todo momento, utilizaremos o seguinte comando:

```bash
%cd /content/drive/Shareddrives/T4-2022/LeonardoMedeiros/Agosto2022
```

# Parte 2 - Preparo do ambiente

## 2.1 - Instalação dos programas necessários

A maior parte dos programas usados nesse notebook podem ser instalados com o comando `apt install`:

```bash
%%bash
sudo apt install tree fastqc bwa samtools bedtools \
1>/content/drive/Shareddrives/T4-2022/LeonardoMedeiros/Agosto2022/logs/log_instalacao.txt \
2>/content/drive/Shareddrives/T4-2022/LeonardoMedeiros/Agosto2022/logs/erros_instalacao.txt

rm /usr/bin/cutadapt
sudo python3 -m pip install --upgrade cutadapt
```

O GATK e o Picard precisam ser instalados separadamente. Vamos baixar (`wget`) e descompactar (`unzip`) os arquivos do GATK.

```bash
%%bash
wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip \
1>logs/log_down_gatk.txt \
2>logs/erros_down_gatk.txt
```

O Picard já vem como executável, então podemos apenas baixá-lo (`wget`).

```bash
%%bash
wget https://github.com/broadinstitute/picard/releases/download/2.24.2/picard.jar \
1>logs/log_down_picard.txt 2>logs/erro_down_picard.txt
```

## 2.2 - Download dos dados em formato FASTQ

```bash
%%bash
wget https://github.com/Varstation/POS-BIOINFO/blob/master/dados/NextSeq_550_Nextera_Capture_NA12878_TruSight_One/FASTQ/NA12878/NA12878TSOENRCC-E_S5_L001_R1_001.fastq.gz?raw=true \
  -O AMOSTRA01_R1.fastq.gz 1>logs/log_down_amostra-r1.txt 2>logs/erro_down_amostra-r1.txt
wget https://github.com/Varstation/POS-BIOINFO/blob/master/dados/NextSeq_550_Nextera_Capture_NA12878_TruSight_One/FASTQ/NA12878/NA12878TSOENRCC-E_S5_L001_R2_001.fastq.gz?raw=true \
  -O AMOSTRA01_R2.fastq.gz 1>logs/log_down_amostra-r2.txt 2>logs/erro_down_amostra-r2.txt
```
Podemos verificar os dados que possuimos utilizando os seguintes comandos:

```bash
%%bash
zcat /content/drive/Shareddrives/T4-2022/LeonardoMedeiros/Agosto2022/AMOSTRA01_R1.fastq.gz | wc -l
zcat /content/drive/Shareddrives/T4-2022/LeonardoMedeiros/Agosto2022/AMOSTRA01_R2.fastq.gz | wc -l
```
Agora que já realizamos o preparo do nosso ambiente de trabalho podemos iniciar processamento dos dados.

# Parte 3 - Processamento dos dados

## 3.1 - Organizar estrutura de diretórios

Mesmo as análises bioinformáticas mais simples acabam por gerar uma quantidade grande de arquivos. É importante mantê-los organizados por data, amostra e tipo de arquivo, de modo que seja simples encontrar esses arquivos caso eles sejam necessários mais tarde. Como nesse caso temos apenas uma amostra, vamos criar uma pasta dados e, dentro dela, pastas para cada etapa de processamento:

```bash
%%bash
mkdir dados
mkdir dados/fastq
mkdir dados/bwa
mkdir dados/picard
mkdir dados/fastqc
mkdir dados/bedtools
mkdir dados/annovar
mkdir dados/freebayes
mkdir dados/gatk
```

O comando `tree` permite mostrar nossa estrutura de diretórios como uma árvore:
```bash
%bash
pwd
tree dados
# mkdir /content/drive/Shareddrives/T4-2022/reference
# mkdir /content/drive/Shareddrives/T4-2022/reference/hg38
```

Vamos mover nossos dados brutos (FASTQ) para a pasta designada na nossa estrutura:

```bash
%%bash
mv AMOSTRA01_R1.fastq.gz ./dados/fastq
mv AMOSTRA01_R2.fastq.gz ./dados/fastq
```
Podemos checar utilizando o comando `tree` como está nossa pasta com os novos arquivos:

```bash
%%bash
tree dados
```
O comando `ls -lh` nos permite ver os arquivos em uma pasta com alguns detalhes, como tamanho do arquivo e data de modificação.

```bash
%%bash
ls -lh dados/fastq
```
## 3.2 - Controle de qualidade dos reads

A seguir vamos usar o programa FASTQC para gerar relatórios de qualidade do sequenciamento:

```bash
%%bash
time fastqc -o dados/fastqc \
  dados/fastq/AMOSTRA01_R1.fastq.gz \
  dados/fastq/AMOSTRA01_R2.fastq.gz
```
# 3.3 - Baixando e indexando o genoma de referência

Para encontrar as variantes genéticas no nosso dado, precisamos compará-lo com a referência do genoma humano. Vamos criar uma estrutura de diretórios para guardar essa referência:

```bash
%%bash
mkdir referencia
mkdir referencia/hg38
```
A seguir, vamos baixá-la do UCSC. Será baixado apenas um cromossomo, para acelerar as análises. Em uma aplicação real, seria necessário baixar o arquivo com todos os cromossomos humanos.

```bash
%%bash
cd /content/drive/Shareddrives/T4-2022/reference/hg38/
ls -l
cat hg38.fa | grep "^>"
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz
```
Também é necessário indexar esse genoma. O objetivo é, de forma simples, criar um arquivo que permita aos algoritmos de mapeamento localizar rapidamente qualquer sequência nesse genoma. Sem essa etapa, não seria possível fazer o mapeamento e chamada de variantes.

```bash
%%bash
time bwa index -a bwtsw referencia/hg19/chr12.fa
```

Criar o indice da referencia (chr12)

```bash
%%bash
samtools faidx referencia/hg19/chr12.fa
```

Em resumo, o código a seguir utiliza o Picard Tools para criar um dicionário de sequência genômica a partir de um arquivo de referência genômica no formato FASTA 

```bash
%%bash
java -jar picard.jar CreateSequenceDictionary \
REFERENCE=referencia/hg19/chr12.fa \
OUTPUT=referencia/hg19/chr12.dict
```
# 3.4 - Alinhamento no genoma de referência

Agora usaremos o BWA para alinhar os reads ao genoma de referência. O arquivo gerado será do tipo SAM.

```bash
%%bash
NOME=NA12878
Biblioteca=TruSightOne
Plataforma=Illumina

time bwa mem -M -R "@RG\tID:CAP\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma" \
/content/drive/Shareddrives/T4-2022/reference/hg38/hg38.fa \
dados/fastq/AMOSTRA01_R1.fastq.gz \
dados/fastq/AMOSTRA01_R2.fastq.gz > dados/bwa/AMOSTRA01.sam
```
Para melhor entendimento do código anterior segue o detalhamento por linha de comando.

>Primeira linha

NOME=NA12878: Isso define uma variável chamada NOME com o valor "NA12878". Essa variável será usada posteriormente no script para fornecer um identificador de amostra.

>Segunda linha

Biblioteca=TruSightOne: Isso define uma variável chamada Biblioteca com o valor "TruSightOne". Essa variável representa o nome da biblioteca de sequenciamento usada na amostra.

>Terceira linha

Plataforma=Illumina: Isso define uma variável chamada Plataforma com o valor "Illumina". Essa variável representa a plataforma de sequenciamento usada para gerar os dados.

>Quinta linha

time: O comando time é usado para medir o tempo de execução do comando que segue. Ele será usado para medir quanto tempo leva para executar o comando bwa mem.

bwa mem -M -R ...: Este é o comando principal. Ele chama o programa bwa mem para realizar o alinhamento de sequência. Segue explicação das operações:

Operações utilizadas:
Operação|Descrição
---|---
-M|Esta opção é usada para marcar fragmentos de leitura como secundários para permitir um melhor processamento posterior, como o ordenamento com ferramentas como o Picard Tools.

-R|"@RG\tID:CAP\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma": Esta opção define um cabeçalho RG (Read Group) que será incluído no arquivo de saída SAM. Ele fornece informações sobre a amostra, biblioteca e plataforma usadas. As variáveis $NOME, $Biblioteca e $Plataforma são substituídas pelos valores definidos anteriormente.

/content/drive/Shareddrives/T4-2022/reference/hg38/hg38.fa: Este é o caminho para o arquivo de referência genômica usado para o alinhamento. O arquivo hg38.fa contém a sequência de referência do genoma humano.

dados/fastq/AMOSTRA01_R1.fastq.gz e dados/fastq/AMOSTRA01_R2.fastq.gz: Esses são os caminhos para os arquivos de leitura de sequência, que estão no formato FASTQ. O arquivo AMOSTRA01_R1.fastq.gz contém as leituras do par de extremidade 1 e o arquivo AMOSTRA01_R2.fastq.gz contém as leituras do par de extremidade 2.

>Oitava linha

dados/bwa/AMOSTRA01.sam: Isso redireciona a saída do comando para um arquivo chamado AMOSTRA01.sam no diretório dados/bwa. O formato SAM é frequentemente usado para representar o resultado do alinhamento.


### O comando bwa mem realiza o alinhamento das leituras de sequência às sequências de referência genômica e cria um arquivo SAM que contém as informações de alinhamento. O tempo de execução do comando é medido usando time.



Agora, o arquivo SAM vai passar por algumas etapas:

  - Correção de erros de pareamento dos reads
  - Compactação do arquivo (SAM → BAM)
  - Ordenação por coordenadas cromossômicas
  - Geração de índice

```bash
%%bash
time samtools fixmate dados/bwa/AMOSTRA01.sam dados/bwa/AMOSTRA01.bam
time samtools sort -O bam -o dados/bwa/AMOSTRA01_sorted.bam dados/bwa/AMOSTRA01.bam
time samtools index dados/bwa/AMOSTRA01_sorted.bam
```

Podemos ver o cabeçalho e algumas entradas do arquivo BAM, com o comando `samtools view`:

```bash
%%bash
time samtools view -H dados/bwa/AMOSTRA01_sorted.bam
```
A opção `-H` especifica que queremos visualizar apenas o cabeçalho do arquivo BAM, não as leituras alinhadas.

Vamos verificar como está o nosso arquivo alinhado

```bash
%%bash
time samtools view dados/bwa/AMOSTRA01_sorted.bam | head
```

## 3.5 - Conferir cobertura

Agora, vamos gerar um arquivo BED a partir do BAM. Isso é necessário porque o programa BedTools só é capaz de calcular cobertura utilizando um arquivo BED, para informar qual é a região na qual a cobertura deve ser contabilizada:

```bash
%%bash
bedtools bamtobed -i dados/bwa/AMOSTRA01_sorted.bam >dados/bedtools/AMOSTRA01_sorted.bed
bedtools merge -i dados/bedtools/AMOSTRA01_sorted.bed >dados/bedtools/AMOSTRA01_merged.bed
bedtools sort -faidx /content/drive/Shareddrives/T4-2022/reference/hg38/hg38.fa.fai -i dados/bedtools/AMOSTRA01_merged.bed >dados/bedtools/AMOSTRA01_merged_sorted.bed
```
Agora usa-se o comando `bedtools coverage` para encontrar a cobertura em cada região:

```bash
%%bash
bedtools coverage -a dados/bedtools/AMOSTRA01_merged_sorted.bed \
-b dados/bwa/AMOSTRA01_sorted.bam -mean \
>dados/bedtools/AMOSTRA01_coverageBed.bed
```

É possível usar o awk para realizar filtros nesse arquivo. Por exemplo, o código a seguir encontra as regiões com cobertura acima de 10:

```bash
%%bash
cat dados/bedtools/AMOSTRA01_coverageBed.bed | \
awk -F "\t" '{if($4>=10){print}}' \
> dados/bedtools/AMOSTRA01_coverageBed10x.bed
```

## 3.6 - Chamada de variantes com o GATK

Em geral recomendamos utilizar mais de um algoritmo para chamada de variantes. Isso traz maior sensibilidade para a detecção de variantes. Nesse caso, utilizaremos o GATK:

```bash
%%bash
sudo chmod 777 gatk-4.1.8.1/gatk HaplotypeCaller
time gatk-4.1.8.1/gatk HaplotypeCaller \
-R /content/drive/Shareddrives/T4-2022/reference/hg38/hg38.fa \
-I dados/bwa/AMOSTRA01_sorted.bam \
-O dados/gatk/AMOSTRA01_sorted.vcf \
-bamout dados/bwa/AMOSTRA01_sorted_bamout.bam
```
Uma explicação melhor detalhada sobre o código

>Primeira linha

sudo chmod 777 gatk-4.1.8.1/gatk HaplotypeCaller: Este comando altera as permissões do arquivo "gatk" e "HaplotypeCaller" para conceder permissões de leitura, escrita e execução para todos os usuários.

>Segunda linha
time gatk-4.1.8.1/gatk HaplotypeCaller ...: Este é o comando principal. Ele chama a ferramenta GATK HaplotypeCaller para realizar o chamado de variantes a partir dos dados de sequenciamento. Segue explicação das operações:

gatk-4.1.8.1/gatk HaplotypeCaller: Isso especifica o caminho para a ferramenta GATK HaplotypeCaller que será executada.

-R /content/drive/Shareddrives/T4-2022/reference/hg38/hg38.fa: Esta opção especifica o arquivo de referência genômica em que o chamado de variantes será baseado. O arquivo "hg38.fa" contém a sequência de referência do genoma humano.

-I dados/bwa/AMOSTRA01_sorted.bam: Esta opção especifica o arquivo BAM de entrada contendo as leituras de sequência genômica alinhadas.

-O dados/gatk/AMOSTRA01_sorted.vcf: Esta opção especifica o nome e localização do arquivo de saída no formato VCF (Variant Call Format), que conterá as variantes chamadas.

-bamout dados/bwa/AMOSTRA01_sorted_bamout.bam: Esta opção especifica o nome e localização do arquivo de saída BAM que conterá as leituras reformatadas para incluir informações sobre o chamado de variantes.


Usamos o comando a seguir para verificar nosso novo arquivo:

```bash
%%bash
cat dados/gatk/AMOSTRA01_sorted.vcf | grep "chr2" | head -n 15
```
Utilizamos o comando `grep` para identificar padrões de texto e nos retornar as linhas que contenham "chr2"

# Parte 4 - Anotação das Variantes

## 4.1 - Download do Annovar do Github

```bash
%%bash
wget https://github.com/Varstation/POS-BIOINFO/raw/master/annovar/annovar.zip
```

## 4.2 - Descompactar o Annovar

```bash
%%bash
#Descompartar 
unzip annovar.zip
#Remover o arquivo zipado; 
rm annovar.zip
```

## 4.3 - Download das bases de dados do Annovar para o hg38 (genoma de referência)

  - refGene
  - exac03
  - avsnp147
  - dbnsfp30a
  - clinvar_20210123

```bash
%%bash
#hg38 - refGene
perl annovar/annotate_variation.pl -buildver hg38 \
-downdb -webfrom annovar refGene annovar/humandb/
```
```bash
%%bash
#hg38 - exac03
perl annovar/annotate_variation.pl -buildver hg38 \
-downdb -webfrom annovar exac03 annovar/humandb/
```
```bash
%%bash
#hg38 - avsnp147
perl annovar/annotate_variation.pl -buildver hg38 \
-downdb -webfrom annovar avsnp147 annovar/humandb/
```
```bash
%%bash
#hg38 - dbnsfp30a
perl annovar/annotate_variation.pl -buildver hg38 \
-downdb -webfrom annovar dbnsfp30a annovar/humandb/
```
```bash
%%bash
#hg38 - clinvar_20210123
perl annovar/annotate_variation.pl -buildver hg38 \
-downdb -webfrom annovar clinvar_20210123 annovar/humandb/
```

## 4.4 - Annotar as variantes do VCF com o ANNOVAR
