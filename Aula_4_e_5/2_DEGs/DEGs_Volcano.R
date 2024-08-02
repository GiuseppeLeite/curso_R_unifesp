# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações
# Aula: DEGs & Volcano
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Carregar os pacotes necessários para leitura e manipulação de dados
library(readr)         # Para ler arquivos de texto
library(tidyverse)     # Para manipulação e visualização de dados
library(limma)         # Para análise de expressão gênica diferencial
library(EnhancedVolcano) # Para criar gráficos de vulcão

# Definir o diretório de trabalho onde os arquivos de dados estão localizados
#setwd()

# Ler os dados de expressão gênica do arquivo "gene_expression_final.txt"
gene_exp_data <- read.delim("gene_expression_final.txt", row.names=1)

# Ler as informações dos pacientes do arquivo "sample_info.txt"
sample_annot <- read.delim("sample_info.txt")

# Selecionar as colunas relevantes dos dados de expressão gênica de acordo com os nomes dos pacientes
column_names <- sample_annot$SampleName[sample_annot$SampleName %in% names(gene_exp_data)]
gene_exp_data_ordered <- gene_exp_data[, column_names]

# Verificar se a ordem das colunas está correta comparando com os nomes dos pacientes
identical(names(gene_exp_data_ordered), sample_annot$SampleName)

# Identificação de genes diferencialmente expressos (DEGs) ---------------------

# Criar a matriz de design 
f <- factor(sample_annot$groups, levels = unique(sample_annot$groups))
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)

# Aplicar os valores de intensidade ao modelo linear
fit <- lmFit(gene_exp_data_ordered, design)

# Criar uma matriz de contrastes
contrast.matrix <- makeContrasts("alive-Control",
                                 "dead-Control",
                                 "alive-dead",
                                 levels=design)

# Aplicar a matriz de contrastes aos dados modelados e calcular as estatísticas para os dados
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Produzir estatísticas para o conjunto de dados e salvar no computador

## Comparação alive_Control -------------------------------

# Extrair a tabela de genes diferencialmente expressos para a comparação "alive_Control"
alive_Control <- topTable(fit2, coef=1, number=Inf, adjust.method="BH")
write.table(alive_Control, file="alive_Control_DEGs.txt", sep="\t", quote=FALSE)

# Filtrar genes com valor p ajustado <= 0.05
alive_Control_sig <- alive_Control %>% filter(adj.P.Val <= 0.05)

# Ordenar genes pelo logFC (fold change)
alive_Control_sig <- alive_Control_sig[order(alive_Control_sig$logFC),]

# Selecionar os 5 genes mais regulados positivamente e negativamente
top_genes <- rbind(head(alive_Control_sig, 5), tail(alive_Control_sig, 5))

### Gráfico de Vulcão ###############################################################

# Definir cores para os pontos no gráfico de vulcão
keyvals <- ifelse(alive_Control$logFC < -0.37 & alive_Control$adj.P.Val < 0.05, "#20207d",
                  ifelse(alive_Control$logFC > 0.37 & alive_Control$adj.P.Val < 0.05, "#7d2020",
                         'gray'))

# Tratar valores NA (não disponíveis) atribuindo a cor cinza
keyvals[is.na(keyvals)] <- 'gray'

# Nomear as cores para a legenda do gráfico
names(keyvals)[keyvals == "#7d2020"] <- 'Up-regulated genes'  # Genes regulados positivamente
names(keyvals)[keyvals == 'gray'] <- 'NS'                 # Genes não significativos
names(keyvals)[keyvals == "#20207d"] <- 'Down-regulated genes' # Genes regulados negativamente

# Criar o gráfico de vulcão usando o pacote EnhancedVolcano
myvolcano <- EnhancedVolcano(alive_Control, 
                             lab = rownames(alive_Control),  # Rotulos dos pontos (nomes dos genes)
                             x = "logFC",  # Eixo X será o logFC (fold change)
                             y = "adj.P.Val",  # Eixo Y será o valor p ajustado
                             title = 'Alive vs. Control',  # Título do gráfico
                             selectLab = rownames(top_genes),  # Selecionar os genes de maior interesse
                             pCutoff = 0.05,  # Valor de corte para o valor p ajustado
                             FCcutoff = 0.37,  # Valor de corte para o logFC
                             pointSize = 3.0,  # Tamanho dos pontos
                             labSize = 4.0,  # Tamanho das labels
                             colCustom = keyvals,  # Cores personalizadas definidas anteriormente
                             colAlpha = 4/5,  # Transparência das cores
                             boxedLabels = TRUE,  # Labels em caixas
                             legendPosition = 'right',  # Posição da legenda
                             legendLabSize = 10,  # Tamanho da fonte da legenda
                             legendIconSize = 4.0,  # Tamanho dos ícones na legenda
                             drawConnectors = TRUE,  # Desenhar conectores entre pontos e labels
                             widthConnectors = 0.5,  # Largura dos conectores
                             colConnectors = 'black',  # Cor dos conectores
                             gridlines.major = TRUE,  # Exibir grades principais
                             gridlines.minor = TRUE,  # Exibir grades secundárias
                             border = 'full',  # Bordas completas ao redor do gráfico
                             borderWidth = 1,  # Largura das bordas
                             borderColour = 'black') +  # Cor das bordas
  ggplot2::coord_cartesian(xlim=c(-7, 7)) +  # Definir limites do eixo X
  ggplot2::scale_x_continuous(breaks=seq(-7, 7, 1)) +  # Definir os intervalos do eixo X
  ylab("-log10(adjusted p value)") +  # Rótulo do eixo Y
  theme(axis.text = element_text(color = 'black'),  # Definir cor do texto dos eixos
        panel.grid.minor = element_blank(),  # Remover grades secundárias
        panel.grid.major = element_blank(), # Remover grades principais
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))  

# Mostrar o gráfico de vulcão
myvolcano


# Salvar o gráfico de vulcão como PDF
dev.copy2pdf(file="Volcano_plot_alive_Control.pdf", width = 8, height = 8)
dev.off()

## Comparação dead_Control -------------------------------

# Extrair a tabela de genes diferencialmente expressos para a comparação "dead-Control"
dead_Control <- topTreat(fit2, coef=2, number=Inf, adjust.method="BH", confint=TRUE)
write.table(dead_Control, file="dead_Control_DEGs.txt", sep="\t", quote=FALSE)





## Comparação alive_dead -------------------------------

# Extrair a tabela de genes diferencialmente expressos para a comparação "alive-dead"
alive_dead <- topTreat(fit2, coef=3, number=Inf, adjust.method="BH", confint=TRUE)
write.table(alive_dead, file="alive_dead_DEGs.txt", sep="\t", quote=FALSE)

