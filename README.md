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

bwa mem -M -R ...: Este é o comando principal. Ele chama o programa bwa mem para realizar o alinhamento de sequência. Aqui estão os principais detalhes:

        -M: Esta opção é usada para marcar fragmentos de leitura como secundários para permitir um melhor processamento posterior, como o ordenamento com ferramentas como o Picard Tools.

        -R "@RG\tID:CAP\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma": Esta opção define um cabeçalho RG (Read Group) que será incluído no arquivo de saída SAM. Ele fornece informações sobre a amostra, biblioteca e plataforma usadas. As variáveis $NOME, $Biblioteca e $Plataforma são substituídas pelos valores definidos anteriormente.

        /content/drive/Shareddrives/T4-2022/reference/hg38/hg38.fa: Este é o caminho para o arquivo de referência genômica usado para o alinhamento. O arquivo hg38.fa contém a sequência de referência do genoma humano.

        dados/fastq/AMOSTRA01_R1.fastq.gz e dados/fastq/AMOSTRA01_R2.fastq.gz: Esses são os caminhos para os arquivos de leitura de sequência, que estão no formato FASTQ. O arquivo AMOSTRA01_R1.fastq.gz contém as leituras do par de extremidade 1 e o arquivo AMOSTRA01_R2.fastq.gz contém as leituras do par de extremidade 2.

>Oitava linha

dados/bwa/AMOSTRA01.sam: Isso redireciona a saída do comando para um arquivo chamado AMOSTRA01.sam no diretório dados/bwa. O formato SAM é frequentemente usado para representar o resultado do alinhamento.


### O comando bwa mem realiza o alinhamento das leituras de sequência às sequências de referência genômica e cria um arquivo SAM que contém as informações de alinhamento. O tempo de execução do comando é medido usando time.
 






