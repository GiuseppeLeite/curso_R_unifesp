# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações
# Aula: Gráfico do GO BP baseado em resultados filtrados (ORA)
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Carregar as bibliotecas necessárias
library(ggplot2)  # Para criar gráficos
library(tidyverse)  # Para manipulação de dados

# Ler o arquivo de dados contendo termos de GO BP (Biological Process) e valores associados
go_bp <- read.delim("go_bp_DEAD_Control.txt")

# Filtrar os termos de GO BP com FDR (False Discovery Rate) menor que 0.05
go_bp <- go_bp |>
  dplyr::filter(FDR < 0.05) 


# Criar o gráfico de barras usando ggplot2
p1 <- ggplot(go_bp, aes(x = reorder(Term, Fold.Enrichment), y = Fold.Enrichment)) +
  geom_bar(stat = "identity", aes(fill = -log10(FDR)), width = 0.8, color = "black", linewidth = 0.25) +  # Adicionar barras preenchidas com cores baseadas em -log10(FDR)
  scale_fill_gradient(low = "#d8bcbc", high = "#7d2020") +  # Definir a escala de cores do preenchimento
  geom_hline(yintercept = 0) +  # Adicionar uma linha horizontal no eixo y = 0
  coord_flip() +  # Inverter os eixos x e y para um gráfico horizontal
  labs(x = "Biological process (GO)", y = "Fold Enrichment") +  # Adicionar rótulos aos eixos
  theme_bw() +  # Usar o tema de fundo branco
  theme(
    axis.line = element_line(color = 'black'),  # Cor das linhas dos eixos
    plot.background = element_blank(),  # Remover o fundo do gráfico
    panel.grid.minor = element_blank(),  # Remover as grades menores
    panel.grid.major = element_blank(),  # Remover as grades maiores
    legend.position = "right",  # Posicionar a legenda à direita
    axis.text = element_text(color = "black")  # Definir a cor do texto dos eixos
  )

# Exibir o gráfico
p1

# Salvar o gráfico em um arquivo PDF
dev.copy2pdf(file="bar_plot_ORA.pdf", width = 12, height = 7)
