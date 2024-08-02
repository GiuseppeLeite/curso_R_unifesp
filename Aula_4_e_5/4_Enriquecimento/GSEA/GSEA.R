# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações
# Aula: Gene Set Enrichment Analysis (GSEA)
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Carregar as bibliotecas necessárias
library(tidyverse)  # Conjunto de pacotes para manipulação e visualização de dados
library(fgsea)      # Pacote para análise de enriquecimento de conjuntos de genes (GSEA)
library(msigdbr)    # Pacote para obter assinaturas de conjuntos de genes do MSigDB
set.seed(123)       # Definir semente para reprodutibilidade dos resultados



# Alive vs Control -------------------------------------------------------------

# Ler o arquivo de dados contendo genes diferencialmente expressos
alive_Control <- read.delim("alive_Control_DEGs.txt")
alive_Control <- rownames_to_column(alive_Control, var = "Gene_name")

# Extrair as classificações (ranks) para a análise de enriquecimento
ranks <- alive_Control |>
  select(Gene_name, t)  # Selecionar apenas as colunas de interesse: nomes dos genes e estatísticas t

# Ordenar as classificações em ordem decrescente
ranks <- ranks[order(-ranks$t, decreasing = TRUE), ]

# Converter as classificações para um vetor nomeado
ranks <- deframe(ranks)


# Definir a lista de assinaturas de genes Hallmark usando o msigdbr
pathways.hallmark <- msigdbr(species = "Homo sapiens", category = "H") |>
  distinct(gs_name, gene_symbol) |>
  nest(genes = c(gene_symbol)) |>    
  mutate(genes = map(genes, compose(as_vector, unname))) |> 
  deframe()  

# Realizar a análise de enriquecimento para o conjunto de genes Hallmark
fgseaRes <- fgsea(pathways = pathways.hallmark, stats = ranks, 
                  nPerm = 10000, # Realiza 10.000 permutações para calcular valores de p mais precisos
                  minSize = 15,  # Tamanho mínimo dos conjuntos de genes
                  maxSize = 500) # Tamanho máximo dos conjuntos de genes

# Remover a informação de borda principal dos resultados (leading edge)
fgseaRes$leadingEdge <- NULL

# Salvar os resultados em um arquivo
write.table(fgseaRes, file = "alive_Control_hallmark_MSigDB.txt")

# Definir a lista de assinaturas de genes Reactome usando o msigdbr
pathways.reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME") |>
  distinct(gs_name, gene_symbol) |>
  nest(genes = c(gene_symbol)) |>
  mutate(genes = map(genes, compose(as_vector, unname))) |>
  deframe()

# Realizar a análise de enriquecimento para o conjunto de genes Reactome
fgseaRes <- fgsea(pathways = pathways.reactome, stats = ranks,  
                  nperm=10000, 
                  minSize = 15,  
                  maxSize = 500) 

# Remover a informação de borda principal dos resultados (leading edge)
fgseaRes$leadingEdge <- NULL

# Salvar os resultados em um arquivo
write.table(fgseaRes, file = "alive_Control_reactome_MSigDB.txt")


# Dead vs Control --------------------------------------------------------------

# Ler o arquivo de dados contendo genes diferencialmente expressos
dead_Control <- read.delim("dead_Control_DEGs.txt")
dead_Control <- rownames_to_column(dead_Control, var = "Gene_name")  

# Extrair as classificações (ranks) para a análise de enriquecimento
ranks <- dead_Control |>
  select(Gene_name, t)

# Ordenar as classificações em ordem decrescente
ranks <- ranks[order(-ranks$t, decreasing = TRUE), ]

# Converter as classificações para um vetor nomeado
ranks <- deframe(ranks)

# Definir a lista de assinaturas de genes Hallmark usando o msigdbr
pathways.hallmark <- msigdbr(species = "Homo sapiens", category = "H") |>
  distinct(gs_name, gene_symbol) |>
  nest(genes = c(gene_symbol)) |> 
  mutate(genes = map(genes, compose(as_vector, unname))) |> 
  deframe() 

# Realizar a análise de enriquecimento para o conjunto de genes Hallmark
fgseaRes <- fgsea(pathways = pathways.hallmark, stats = ranks, 
                  nPerm = 10000,
                  minSize = 15, 
                  maxSize = 500)

# Remover a informação de borda principal dos resultados (leading edge)
fgseaRes$leadingEdge <- NULL

# Salvar os resultados em um arquivo
write.table(fgseaRes, file = "dead_Control_hallmark_MSigDB.txt")

# Definir a lista de assinaturas de genes Reactome usando o msigdbr
pathways.reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME") |>
  distinct(gs_name, gene_symbol) |>
  nest(genes = c(gene_symbol)) |>
  mutate(genes = map(genes, compose(as_vector, unname))) |>
  deframe()

# Realizar a análise de enriquecimento para o conjunto de genes Reactome
fgseaRes <- fgsea(pathways = pathways.reactome, stats = ranks,  
                  nperm=10000, 
                  minSize = 15,  
                  maxSize = 500) 

# Remover a informação de borda principal dos resultados (leading edge)
fgseaRes$leadingEdge <- NULL

# Salvar os resultados em um arquivo
write.table(fgseaRes, file = "dead_Control_reactome_MSigDB.txt")


# Alive vs Dead -----------------------------------------------------------------

# Ler o arquivo de dados contendo genes diferencialmente expressos
alive_dead <- read.delim("alive_dead_DEGs.txt")
alive_dead <- rownames_to_column(alive_dead, var = "Gene_name") 

# Extrair as classificações (ranks) para a análise de enriquecimento
ranks <- alive_dead |>
  select(Gene_name, t)

# Ordenar as classificações em ordem decrescente
ranks <- ranks[order(-ranks$t, decreasing = TRUE), ]

# Converter as classificações para um vetor nomeado
ranks <- deframe(ranks)

# Definir a lista de assinaturas de genes Hallmark usando o msigdbr
pathways.hallmark <- msigdbr(species = "Homo sapiens", category = "H") |>
  distinct(gs_name, gene_symbol) |> 
  nest(genes = c(gene_symbol)) |> 
  mutate(genes = map(genes, compose(as_vector, unname))) |> 
  deframe() 

# Realizar a análise de enriquecimento para o conjunto de genes Hallmark
fgseaRes <- fgsea(pathways = pathways.hallmark, stats = ranks, 
                  nPerm = 10000,
                  minSize = 15,
                  maxSize = 500)

# Remover a informação de borda principal dos resultados (leading edge)
fgseaRes$leadingEdge <- NULL

# Salvar os resultados em um arquivo
write.table(fgseaRes, file = "alive_dead_hallmark_MSigDB.txt")

# Definir a lista de assinaturas de genes Reactome usando o msigdbr
pathways.reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME") |>
  distinct(gs_name, gene_symbol) |>
  nest(genes = c(gene_symbol)) |>
  mutate(genes = map(genes, compose(as_vector, unname))) |>
  deframe()

# Realizar a análise de enriquecimento para o conjunto de genes Reactome
fgseaRes <- fgsea(pathways = pathways.reactome, stats = ranks,  
                  nperm=10000, 
                  minSize = 15,  
                  maxSize = 500) 

# Remover a informação de borda principal dos resultados (leading edge)
fgseaRes$leadingEdge <- NULL

# Salvar os resultados em um arquivo
write.table(fgseaRes, file = "alive_dead_reactome_MSigDB.txt")
