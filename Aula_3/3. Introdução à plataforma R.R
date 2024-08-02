# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações
# Aula: Introdução à plataforma R
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Local de Trabalho ---------------------------------------------------------
# Mostrar o diretório de trabalho atual
getwd()

# Listar os arquivos do diretório
dir() 

# Mudar o diretório de trabalho (substitua "novo_diretorio" pelo caminho do seu diretório)
setwd("novo_diretorio") 

# Instalação de Pacotes -----------------------------------------------------
# Instalar pacotes do CRAN
install.packages(c("tidyverse", "ggplot2", "scales", "msigdbr", "pheatmap", "rstatix", "ggpubr"))

# Verifica e instala o BiocManager, se necessário
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar pacotes do Bioconductor
BiocManager::install(c("limma", "fgsea", "CEMiTool", "GO.db", "HDO.db", 
                       "EnhancedVolcano"))

# Carregar Pacotes ----------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Operações Básicas e Tipos de Dados ----------------------------------------

# Operações Aritméticas
1 + 1  # Adição
2 - 1  # Subtração
2 * 3  # Multiplicação
6 / 2  # Divisão


# Tipos de Dados
## Numéricos
num <- 42
Num <- 42
Num == num

## Caracteres
char <- "Olá, mundo!"

## Lógicos
log <- FALSE

## Fatores
fact <- factor(c("Macho", "Fêmea", "Fêmea", "Macho"))
fact

# Vetores e Matrizes
## Vetores
vetor <- c(1, 2, 3, 4)

## Matrizes
matriz <- matrix(1:6, nrow = 2, ncol = 3)

# Manipulação de Dados ------------------------------------------------------
## Criando um Data Frame
df <- data.frame(
  id = 1:4,
  nome = c("Ana", "Pedro", "Bianca", "Carlos"),
  idade = c(28, 35, 23, 40)
)
df

## Escrita de Dados
write.csv(df, "dados.csv")
write.table(df, "dados.txt")

## Leitura de Dados
dados <- read.csv("dados.csv", row.names = 1)
dados_2 <- read.table("dados.txt", sep = " ")

# Funções Básicas de Manipulação de Dados -----------------------------------
## Subset
up_30 <- subset(df, idade > 30)

## Merge
df2 <- data.frame(
  id = 1:4,
  salario = c(4000, 5000, 4500, 5500)
)

df <- merge(df, df2, by = "id", all = FALSE)

## dplyr
df3 <- df |>
  dplyr::filter(idade > 30) |> 
  dplyr::select(nome, idade)

# Visualização de Dados -----------------------------------------------------
# Gráfico de dispersão
ggplot(data = df, aes(x = idade, y = salario)) +
  geom_point()

# Gráfico de barras
ggplot(data = df, aes(x = nome, y = salario)) +
  geom_bar(stat = "identity")

