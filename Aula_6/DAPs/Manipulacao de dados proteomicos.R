# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações 
# Aula: Manipulação de dados proteômicos
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Carregar as bibliotecas necessárias
library(tidyverse)  # Inclui ggplot2, dplyr, tidyr, readr, etc.
library(randomForest)  # Para imputação usando RandomForest

# Carregar o dataset inicial
prot_data <- read.delim("proteinGroups.txt")

# Selecione as colunas desejadas
prot_data_novo <- prot_data |>
  select(
    Protein.IDs,  # Identificadores de proteínas
    Peptides,  # Contagem de peptídeos
    Unique.peptides,  # Contagem de peptídeos únicos
    Potential.contaminant,  # Potenciais contaminantes
    Only.identified.by.site,  # Apenas identificados por site
    Reverse,  # Sequências reversas (controle)
    starts_with("LFQ.intensity.")  # Intensidades LFQ
  )

# Filtrar as linhas com Peptides >= 2 e Unique.peptides >= 1
prot_data_novo <- prot_data_novo |>
  filter(Peptides >= 2, Unique.peptides >= 1)

# Remover as linhas onde a coluna Protein.IDs contém "_BOVIN"
prot_data_novo <- prot_data_novo |>
  filter(!str_detect(Protein.IDs, "_BOVIN"))

# Filtrar as linhas que não contêm '+' nas colunas especificadas
prot_data_novo <- prot_data_novo |>
  filter(
    Potential.contaminant != "+",  # Remove proteínas potenciais contaminantes
    Only.identified.by.site != "+",  # Remove proteínas identificadas apenas por um único sítio de peptídeo
    Reverse != "+"  # Remove sequências de proteínas reversas usadas como controle negativo
  )

# Adicionar os rownames como uma coluna chamada "Protein.IDs"
prot_data_novo <- column_to_rownames(prot_data_novo, var = "Protein.IDs")

# Selecionar as colunas 6 até 82, substituir 0 por NA e aplicar a transformação logarítmica
prot_data_novo <- prot_data_novo |>
  select(6:82) |>
  mutate_all(~ ifelse(. == 0, NA, .)) |>
  mutate_at(vars(1:77), log2)

# Remover linhas com mais de 50% de NA
prot_data_novo <- prot_data_novo %>%
  filter(rowSums(is.na(.)) <= ncol(.) * 0.5)

# preparar os dados para imputação ---------------------------------------------
# Carregar a anotação
sample_annot <- read.delim("Phenodata.txt")  # Ler o arquivo de anotações

# Selecionar as colunas relevantes dos dados de expressão gênica de acordo com os nomes dos pacientes
column_names <- sample_annot$sample[sample_annot$sample %in% names(prot_data_novo)]
prot_data_novo_ordered <- prot_data_novo[, column_names]

# Verificar se a ordem das colunas está correta comparando com os nomes dos pacientes
identical(names(prot_data_novo_ordered), sample_annot$sample)

# Rotacionar o dataframe
prot_data_novo_ordered <- prot_data_novo_ordered |>
  rownames_to_column(var = "protein_ID") |>
  sjmisc::rotate_df(rn = "sample", cn = TRUE)  # Rotacionar o dataframe, mantendo "sample" como nomes das linhas

# Mesclar as anotações com o dataset rotacionado
prot_data_novo_ordered_annot <- merge(sample_annot, prot_data_novo_ordered, by = "sample")  # Mesclar pelo identificador "sample"

# Salvar o dataset anotado
write.table(prot_data_novo_ordered_annot, 
            file = "prot_data_novo_ordered_annot.txt", 
            sep = "\t", dec = ".", 
            row.names = FALSE, col.names = TRUE)  # Salvar o dataset mesclado

# Selecionar apenas as colunas de interesse
prot_data_novo_ordered_annot <- prot_data_novo_ordered_annot[2:357]  # Manter colunas 2 a 357 para a imputação


# Imputation with the "RandomForest" package -----------------------------------

# A função rfImpute é usada para preencher valores faltantes no dataset.
# A imputação é feita usando uma RandomForest, que é um conjunto de árvores de decisão.

# A imputação é feita com os seguintes parâmetros:
# as.factor(Group) ~ . : Especifica que a variável 'Group' é o fator de interesse e todas as outras variáveis (.) são preditores.
# ntree=1000 : Define que 1000 árvores de decisão serão usadas na floresta.
# iter = 10 : Define que 10 iterações de imputação serão realizadas.
# data = prot_data_novo_ordered_annot : O dataset de entrada para a imputação.

# Executar a imputação
prot_data_novo_ordered_annot <- rfImpute(as.factor(Group) ~ ., 
                                         ntree = 1000,  # Usar 1000 árvores na RandomForest
                                         iter = 10,  # Realizar 10 iterações para melhorar a precisão da imputação
                                         data = prot_data_novo_ordered_annot)  # Usar o dataset anotado como entrada

# Salvar o dataset imputado
write.table(prot_data_novo_ordered_annot, 
            file = "prot_data_novo_ordered_annot_imputation_12_03.txt", 
            sep = "\t", dec = ".", 
            row.names = FALSE, col.names = TRUE)  # Salvar o dataset imputado em um novo arquivo
