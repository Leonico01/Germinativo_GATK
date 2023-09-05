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



























