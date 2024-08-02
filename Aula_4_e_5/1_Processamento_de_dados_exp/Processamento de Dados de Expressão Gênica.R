# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações
# Aula: Processamento de Dados de Expressão Gênica
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br


# Carregar o pacote readr para ler arquivos de texto

library(readr)
library(tidyverse)

# Definir o diretório de trabalho onde os arquivos estão localizados
# setwd()

# Ler os dados de expressão gênica do arquivo "gene_exp_data.txt" 

gene_exp_data <- read.delim("gene_exp_data.txt", row.names=1)

# Ler as informações dos pacientes do arquivo "sample_info.txt"

sample_annot <- read.delim("sample_info.txt")

# Selecionar as colunas relevantes dos dados de expressão gênica de acordo com os nomes dos pacientes
column_names <- sample_annot$SampleName[sample_annot$SampleName %in% names(gene_exp_data)]
gene_exp_data_ordered <- gene_exp_data[, column_names]

# Verificar se a ordem das colunas está correta comparando com os nomes dos pacientes

identical(names(gene_exp_data_ordered), sample_annot$SampleName)

# Ler o arquivo de anotação que contém informações sobre os genes

annotation <- read_delim("GPL13667.txt", delim = "\t", comment = "#")

# Selecionar apenas as colunas 'ID' e 'Gene Symbol', renomear 'Gene Symbol' para 'Gene_Symbol'
# e remover linhas onde 'Gene_Symbol' contém '/' ou '---'

annotation <- annotation %>% 
  select(ID, `Gene Symbol`) %>%
  rename(Gene_Symbol = `Gene Symbol`) %>%
  filter(!grepl("/", Gene_Symbol) & !grepl("---", Gene_Symbol))

# Converter gene_exp_data_ordered em dataframe e adicionar uma coluna 'ID' com os IDs das linhas

gene_exp_data_ordered <- as.data.frame(gene_exp_data_ordered)
gene_exp_data_ordered <- rownames_to_column(gene_exp_data_ordered, var = "ID")

# Juntar os dados de expressão gênica ordenados com as informações de anotação usando os IDs

gene_exp_data_ordered <- merge(annotation, gene_exp_data_ordered,
                               by = "ID", all = FALSE)

# Adicionar uma coluna 'AveExpr' com a média da expressão de cada sonda (probe)

gene_exp_data_ordered$AveExpr <- rowMeans(gene_exp_data_ordered[, c(3:201)])

# Ordenar os genes pela média de expressão em ordem decrescente

gene_exp_data_ordered <- gene_exp_data_ordered[order(gene_exp_data_ordered$AveExpr, decreasing = TRUE),]

# Remover sondas duplicadas para o mesmo símbolo, mantendo apenas aquela com o maior valor de expressão geral

gene_exp_data_ordered <- gene_exp_data_ordered[!duplicated(gene_exp_data_ordered$Gene_Symbol),]

# Remover as colunas 'ID' e 'AveExpr', que não são mais necessárias

gene_exp_data_ordered$ID <- NULL
gene_exp_data_ordered$AveExpr <- NULL

# Remover linhas com valores ausentes (NA)

gene_exp_data_ordered <- na.omit(gene_exp_data_ordered) 

# Escrever a tabela final em um arquivo chamado "gene_expression_final.txt"

write.table(gene_exp_data_ordered, file="gene_expression_final.txt", 
            sep="\t", row.names = FALSE, quote = FALSE)
