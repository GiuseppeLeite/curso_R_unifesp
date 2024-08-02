# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações 
# Aula: Gráfico do hallmark baseado em resultados filtrados - grupos
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Carregar bibliotecas necessárias
library(ggplot2)  # Para visualização de dados
library(scales)   # Para ajustes de escala

# Ler os dados do arquivo "Grupos_combinados.txt" e criar um data frame
data <- read.delim("Grupos_combinados.txt")

# Definir os grupos como um fator com a ordem especificada
data$groups <- factor(data$groups, levels = c("alive_Control", "Dead_Control"))

# Criar o gráfico de dispersão usando ggplot2
p1 <- ggplot(data, aes(x = groups, y = reorder(pathway, NES), size = -log10(padj), groups = groups)) + 
  geom_point(aes(fill = NES), alpha = 0.95, stroke = 0.5, shape = 21, color = "black") +  # Adicionar pontos de dispersão
  scale_fill_gradient2(low = "#20207d", mid = "white", high = "#7d2020", midpoint = 0) +  # Definir paleta de cores para os pontos
  scale_size(range = c(5, 11)) +  # Ajustar o tamanho dos pontos
  scale_y_discrete(labels = label_wrap(45)) +  # Rotular os eixos y com wrap de texto
  theme_bw() +  # Usar tema preto e branco
  theme(
    axis.title = element_blank(),  # Remover títulos dos eixos
    axis.line = element_line(color = 'black'),  # Definir cor das linhas dos eixos
    plot.background = element_blank(),  # Definir fundo do gráfico como vazio
    panel.grid = element_blank(),  # Remover grid
    legend.position = "right",  # Posição da legenda
    axis.text = element_text(size = 12, color = 'black')  # Ajustar tamanho e cor do texto dos eixos
  )

# Exibir o gráfico
p1

# Salvar o gráfico como um arquivo PDF
dev.copy2pdf(file = "hallmark_grupos.pdf", width = 7, height = 6)
dev.off()  # Finalizar o dispositivo gráfico
