# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações 
# Aula: CEMiTool - Gene co-expression analysis
# Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Ref.: https://bioconductor.org/packages/release/bioc/vignettes/CEMiTool/inst/doc/CEMiTool.html


# Instalar o pacote CEMiTool (descomente a linha abaixo se ainda não tiver o pacote instalado)
# BiocManager::install("CEMiTool")

# Carregar os pacotes necessários
library(CEMiTool) # Pacote para análise de co-expressão gênica
library(tidyverse) # Pacote para manipulação e visualização de dados
library(fgsea)

# Preparação dos dados --------------------------------------------------------

# Carregar os dados de expressão gênica - mesmos dados usados para CIBERSORTx
# 'gene_expression_final.txt' contém os dados de expressão gênica
gene_exp_data <- read.delim("gene_expression_final.txt", row.names = 1)

# Carregar a anotação das amostras
# 'sample_info.txt' contém informações sobre as amostras dos pacientes
sample_annot <- read.delim("sample_info.txt")

# Selecionar as colunas relevantes dos dados de expressão gênica de acordo com os nomes dos pacientes
# Filtrar as colunas do gene_exp_data para incluir apenas os pacientes presentes em sample_annot
column_names <- sample_annot$SampleName[sample_annot$SampleName %in% names(gene_exp_data)]
gene_exp_data_ordered <- gene_exp_data[, column_names]

# Verificar se a ordem das colunas está correta comparando com os nomes dos pacientes
# Isso assegura que a ordem das colunas em gene_exp_data_ordered é a mesma que em sample_annot
identical(names(gene_exp_data_ordered), sample_annot$SampleName)

# Carregar o arquivo GMT (Gene Matrix Transposed) contendo informações de vias (pathways)
# O arquivo pathways.gmt é fornecido pelo pacote CEMiTool
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

#ou você pode baixar direito os arquivos do msigdb
# https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
gmt_in <- read_gmt("h.all.v2023.2.Hs.symbols.gmt")

# Carregar o arquivo de interações entre genes
# O arquivo interactions.tsv é fornecido pelo pacote CEMiTool
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)

# Executar CEMiTool -----------------------------------------------------------

# Realizar a análise de co-expressão gênica usando o CEMiTool
# Parâmetros:
# - filter: Filtrar genes com baixa variabilidade
# - filter_pval: Valor-p para o filtro de variabilidade
# - min_ngen: Número mínimo de genes por módulo
# - apply_vst: Aplicar Transformação de Variância Estabilizadora (FALSE)
# - eps: Tolerância de erro para a rede de co-expressão
# - plot: Gerar gráficos (TRUE)
# - verbose: Mostrar informações detalhadas no console (TRUE)


cem <- cemitool(gene_exp_data_ordered, sample_annot, gmt_in, interactions=int_df,
                filter=TRUE, filter_pval = 0.1, min_ngen = 50,
                apply_vst = FALSE, eps = 0.1,
                plot=TRUE, verbose=TRUE)

# Mostrar o gráfico de Análise de Enriquecimento de Vias (GSEA)
show_plot(cem, "gsea")

# Criar um relatório em formato HTML com os resultados da análise
generate_report(cem, directory="./Report", force=TRUE)

# Escrever os resultados da análise em arquivos
write_files(cem, directory="./Tables", force=TRUE)

# Salvar todos os gráficos gerados pela análise
save_plots(cem, "all", directory="./Plots", force=TRUE)
