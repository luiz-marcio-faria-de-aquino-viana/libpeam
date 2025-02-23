# libpeam
BIBLIOTECA: LIBPEAM.py

LIBPEAM - Bibliotaca de Funções de Probabilidade e Estatística para Aprendizado de Máquina.

Desenvolvida durante o Curso de Pós-Doutoramento em Engenharia Elétrica com Ênfase em Engenharia de Sistemas e Computação (COPPE/UFRJ em 2020).

Disciplina: PEAM - Probabilidade e Estatística para Aprendizado de Máquina
Professores: 
- Rosa Leão (https://www.cos.ufrj.br/index.php/pt-BR/telefones-do-pesc/details/3/1050-rosa);
- Edmundo de Souza e Silva (https://www.cos.ufrj.br/index.php/pt-BR/pessoas-search/details/18/1013-edmundo); e
- Daniel Menasch

Autor(100%): Luiz Marcio Faria de Aquino Viana, Pós-D.Sc.
Prazo de dezenvolvimento da Biblioteca: 1 dia (= 1 sabado!)

Aplicação: Análise dos Dados Diários sobre a COVID-19.

Objetivo: Para fixar os processo estatísticos apresentados em aula, foi desenvolvido a biblioteca LIBPEAM, implementada em Python, com as funções estatísticas necessárias para a resolução dos problemas da lista de exercícios de aula.

RECURSOS:

1. READDATA(): LEITURA DO DATASET PARA UMA LISTA (DESPREZANDO O CABEÇALHO)

2. QTD(), SOMA(), MIN(), MAX(), MEDIA(), VAR(), STD()

3. OUTLIERS(), WORKSET(), LISTA_FREQ()

4. FUNÇÕES:
- Geométrica e Exponencial (definição interna),
- Normal (importada da biblioteca “statistics”),
- Beta,
- F, e
- Chi-Square (importadas da biblioteca: “scipy”)

5. CÁLCULO DO VALOR DE “X” EM FUNÇÃO DA ÁREA (p-value)
- F e Chi-Square

6. CÁLCULO DOS PARÂMETROS - MLE
- Geométrica, Beta, Exponencial e Normal

7. IMPLEMENTAÇÃO DOS TESTES
- One-Way ANOVA
- Chi-Square
- Kolmogorov-Smirnov

8. FUNÇÕES PARA APRESENTAÇÃO GRÁFICA
* implementadas com a biblioteca “matplotlib”

BIBLIOGRAFIA:

PROBABILITY & STATISTICS WITH RELIABILITY, QUEUING, AND COMPUTER SCIENCE APPLICATIONS
AUTOR: Kishor Shridharbhai Trivedi
Prentice-Hall, Inc., Englewood Cliffs, N.J. 07632
ISBN 0-13-711564-4
