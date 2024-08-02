# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações
# Aula: Heatmap dos DEGs
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Carregar os pacotes necessários para leitura e manipulação de dados
library(readr)         # Para ler arquivos de texto
library(tidyverse)     # Para manipulação e visualização de dados
library(pheatmap)      # Para gerar heatmaps
library(RColorBrewer)  # Para paletas de cores

# Definir o diretório de trabalho onde os arquivos de dados estão localizados
# Altere para o diretório apropriado antes de rodar o script
#setwd("caminho/para/seu/diretorio")

# Leitura e manipulação dos dados ---------------------------------------------
# Ler os dados de expressão gênica
gene_exp_data <- read.delim("gene_expression_final.txt", row.names=1)

# Ler as informações dos pacientes do arquivo "sample_info.txt"
sample_annot <- read.delim("sample_info.txt")

# Selecionar as colunas relevantes dos dados de expressão gênica de acordo com os nomes dos pacientes
# Isto garante que os dados de expressão gênica estão na mesma ordem que os nomes dos pacientes
column_names <- sample_annot$SampleName[sample_annot$SampleName %in% names(gene_exp_data)]
gene_exp_data_ordered <- gene_exp_data[, column_names]

# Verificar se a ordem das colunas está correta comparando com os nomes dos pacientes
# Deve retornar TRUE se as colunas estiverem na ordem correta
identical(names(gene_exp_data_ordered), sample_annot$SampleName)

# Ler a lista de DEGs (genes diferencialmente expressos) do arquivo "alive_dead_DEGs.txt"
alive_dead_DEGs <- read.delim("alive_dead_DEGs.txt", row.names=1)

# Filtrar os genes que atendem aos critérios de significância estatística (adj.P.Val <= 0.05)
# e mudança de expressão (logFC >= 1 ou logFC <= -1)
# Adicionar os nomes dos genes como uma coluna
alive_dead_DEGs_filtered <- alive_dead_DEGs %>%
  rownames_to_column(var = "Gene_name") %>%
  filter(adj.P.Val <= 0.05 & (logFC >= 1 | logFC <= -1))

# Salvar os DEGs filtrados em um novo arquivo
write.table(alive_dead_DEGs_filtered, file = "alive_dead_DEGs_filtered.txt")

# Adicionar os nomes dos genes como uma coluna no conjunto de dados de expressão gênica
gene_exp_data_ordered <- rownames_to_column(gene_exp_data_ordered, var = "Gene_name")

# Filtrar os dados de expressão gênica para incluir apenas os genes diferencialmente expressos (DEGs)
filtered_gene_exp_data <- gene_exp_data_ordered %>%
  semi_join(alive_dead_DEGs_filtered, by = "Gene_name") %>%
  column_to_rownames(var = "Gene_name")

# Salvar os dados de expressão gênica filtrados em um arquivo de texto
write.table(filtered_gene_exp_data, file = "filtered_gene_exp_data.txt")

# Ajustar o conjunto de dados de anotações dos pacientes
sample_annot <- column_to_rownames(sample_annot, var = "SampleName")

# Criar um heatmap dos dados de expressão gênica filtrados ---------------------

pheatmap(filtered_gene_exp_data,
               color = colorRampPalette(c("#1111fd", "white", "red"))(50),
               cluster_rows = TRUE,   # Clustering de linhas
               cluster_cols = TRUE,  # Sem clustering de colunas
               scale = "row",         # Normaliza as linhas
               show_rownames = TRUE, # Mostrar nomes das linhas (genes)
               show_colnames = FALSE, # Não mostrar nomes das colunas (amostras)
               annotation_col = sample_annot) # Anotação das colunas (amostras)

# Salvar o heatmap como um arquivo PDF
dev.copy2pdf(file="heatmap_alive_HVs.pdf", width = 6, height = 6)


