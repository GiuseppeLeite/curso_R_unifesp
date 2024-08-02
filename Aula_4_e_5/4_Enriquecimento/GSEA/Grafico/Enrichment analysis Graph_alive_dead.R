# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações
# Aula: Gráfico do hallmark baseado em resultados filtrados
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Carregar a biblioteca necessária
library(ggplot2)  # Para criar gráficos

# Ler o arquivo de dados contendo os resultados filtrados dos hallmark gene sets
hallmark_filtered <- read.delim("alive_dead_hallmark_filtered.txt")

# Criar um vetor com os nomes dos pathways na ordem desejada
pathway_order <- unique(hallmark_filtered$pathway)

# Converter a coluna 'pathway' para um fator com a ordem de níveis desejada
hallmark_filtered$pathway <- factor(hallmark_filtered$pathway, levels = pathway_order)

# Criar o gráfico de barras usando ggplot2
p1 <- ggplot(hallmark_filtered, aes(x = pathway, y = NES)) +
  geom_bar(stat = "identity", aes(fill = NES >= 0), width = 0.8, color = "black", linewidth = 0.25) +  # Adicionar barras preenchidas com cores baseadas no valor de NES
  scale_fill_manual(values = c("#20207d", "#7d2020"), guide = FALSE) +  # Definir as cores do preenchimento manualmente
  geom_hline(yintercept = 0) +  # Adicionar uma linha horizontal no eixo y = 0
  coord_flip() +  # Inverter os eixos x e y para um gráfico horizontal
  labs(x = "Hallmark gene sets ", y = "NES") +  # Adicionar rótulos aos eixos
  theme_bw() +  # Usar o tema de fundo branco
  theme(
    axis.line = element_line(color = 'black'),  # Cor das linhas dos eixos
    plot.background = element_blank(),  # Remover o fundo do gráfico
    panel.grid.minor = element_blank(),  # Remover as grades menores
    panel.grid.major = element_blank(),  # Remover as grades maiores
    legend.position = "none",  # Remover a legenda
    axis.text = element_text(color = "black")  # Definir a cor do texto dos eixos
  )

# Exibir o gráfico
p1

# Salvar o gráfico em um arquivo PDF
dev.copy2pdf(file="bar_plot_enriched_graph.pdf", width = 7, height = 4)
