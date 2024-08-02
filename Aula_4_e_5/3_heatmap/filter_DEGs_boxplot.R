# Curso: Análise de Dados Ômicos: Fundamentos e Aplicações
# Aula: Boxplot dos DEGs
# Prof. Giuseppe Leite
# giuseppe.gianini@unifesp.br

# Carregar os pacotes necessários para leitura e manipulação de dados
library(readr)         # Para ler arquivos de texto
library(ggpubr)        # Para criar gráficos de publicação
library(RColorBrewer)  # Para paletas de cores
library(rstatix)       # Para testes estatísticos e manipulação de dados

# Definir o diretório de trabalho onde os arquivos de dados estão localizados
# Altere para o diretório apropriado antes de rodar o script
#setwd("caminho/para/seu/diretorio")

# Leitura e manipulação dos dados ---------------------------------------------
# Ler o arquivo CSV contendo os dados de expressão gênica filtrados
filtered_gene_exp_data <- read.csv("filtered_gene_exp_data_rot_anno.csv", sep=";")

# Modificar a variável 'group' para definir os níveis na ordem desejada
filtered_gene_exp_data <- filtered_gene_exp_data %>%
  mutate(group = factor(group, levels = c("alive", "dead")))

# Verificar os níveis da variável 'group'
levels(filtered_gene_exp_data$group)

# Realizar o teste t de Student para comparar a expressão de NUDT4 entre os grupos
test <- filtered_gene_exp_data |> 
  t_test(NUDT4 ~ group, paired = FALSE)

########################
#test <- filtered_gene_exp_data |> 
  #wilcox_test(NUDT4 ~ group, paired = FALSE) 

# Verificar o resultado do teste
print(test)

# Adicionar a posição x e y ao resultado do teste para facilitar a plotagem
test <- test |> add_xy_position(x = "group")

# Criar o gráfico de boxplot usando ggpubr

p1 <- ggboxplot(filtered_gene_exp_data, x = "group", y = "NUDT4", fill ="group", outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.2, color = "black") +  # Adicionar barras de erro
  geom_boxplot(aes(color = factor(group)),  # Adicionar caixas do boxplot com cores personalizadas
               fill = c("#0000E7", "#e70000"),  # Cores de preenchimento para cada grupo
               outlier.color = NA,  # Remover a cor dos outliers
               color = "black",  # Cor das bordas das caixas
               size = 0.4) +  # Espessura das bordas das caixas
  stat_summary(fun = mean, geom = "point", shape = 21, size = 1, fill = "white", color = "black") +  # Adicionar ponto da média
  stat_pvalue_manual(test, hide.ns = FALSE, y.position = 10.6) +  # Adicionar valores p ao gráfico
  theme_bw() +  # Usar o tema de fundo branco
  theme(
    axis.line = element_line(color='black', linewidth = 0.1),  # Cor e espessura das linhas dos eixos
    plot.background = element_blank(),  # Remover o fundo do gráfico
    axis.text.x = element_text(colour = "black", size = 11),  # Cor e tamanho do texto do eixo x
    axis.text.y = element_text(colour = "black", size = 11),  # Cor e tamanho do texto do eixo y
    panel.grid.minor = element_blank(),  # Remover grades menores
    panel.grid.major = element_blank(),  # Remover grades maiores
    axis.title.x = element_blank(),  # Remover o título do eixo x
    legend.position = "none"  # Remover a legenda
  )

# Exibir o gráfico
print(p1)

# Salvar o gráfico em um arquivo PDF
dev.copy2pdf(file="NUDT4.pdf", width = 1.5, height = 2.5)
