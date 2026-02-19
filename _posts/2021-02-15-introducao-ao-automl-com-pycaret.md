---
layout: post
title: "Uma introdução ao AutoML com PyCaret"
date: 2021-02-15
categories: data-science
tags: [pycaret, machine-learning, automl, data-science, python]
---

AutoML é a parte de aprendizado de máquina automatizado, nela processos e tarefas que não requerem tomada de decisão são automatizados facilitando a criação de modelos e seu deploy em larga escala com eficiência e aumentando a produtividade do desenvolvedor.

Nesse tutorial vamos ver o [PyCaret](https://github.com/pycaret/pycaret) criado por [Moez Ali](https://www.linkedin.com/in/profile-moez/). PyCaret é uma biblioteca aberta com pouco uso de código em Python. Ela permite que você prepare seus dados e compare modelos selecionando o melhor de acordo com uma métrica de interesse em poucos minutos. Ela pode ser usada para problemas de [classificação](https://pycaret.readthedocs.io/en/latest/api/classification.html), [regressão](https://pycaret.readthedocs.io/en/latest/api/regression.html), [clusterização](https://pycaret.readthedocs.io/en/latest/api/clustering.html), [detecção de anomalias](https://pycaret.readthedocs.io/en/latest/api/anomaly.html), [processamento de linguagem natural](https://pycaret.readthedocs.io/en/latest/api/nlp.html), [mineração de regras de associação](https://pycaret.readthedocs.io/en/latest/api/arules.html), etc.

No site deles existem alguns [tutoriais](https://pycaret.readthedocs.io/en/latest/tutorials.html) bem interessantes, vale a pena dar uma olhada pra ver se é possível usar em algum projeto do seu interesse. Nas referências desse post vou deixar alguns materiais que levei em conta ao escrever esse. Aqui nós vamos usar dados bem simples do UCI "[Heart failure clinical records Data Set](https://archive.ics.uci.edu/ml/datasets/Heart+failure+clinical+records)" que já estão prontos (fora da realidade) para serem processados para um problema de classificação.

Eu recomendo você usar o Google Colab para testar essa ferramenta, o meu colab que usei para criar esse post pode ser visto [aqui](https://github.com/geocarvalho/tutoriais/blob/main/ml/pycaret_first_impression.ipynb). Recomendo você criar do zero um novo e ir rodando aos poucos e à medida que tiver dúvidas ir jogando no google e pesquisando na documentação do PyCaret.

Na versão 2.x está disponível o uso de GPU para treinar os modelos, mas pra que usar o PyCaret?

- O pessoal usa os modelos criados pelo PyCaret nos dashboards do PowerBI.
- Testar modelos para projetos que você ache interessante usar aprendizagem de máquina (protótipos).

## Instalando o PyCaret

Como estamos usando o Colab, basta usar o comando pip para instalar a ferramenta:

```python
!pip install pycaret
```

## Abrindo os dados de exemplo

Como eu disse antes vamos trabalhar com dados já bonitinhos para a criação de modelos e isso tá longe da realidade ok?! Esse é só um caso para a gente ver algumas utilidades do PyCaret. Os dados eu encontrei dentro de outro github e estou reaproveitando o [link](https://github.com/lorenzodenisi/Heart-Failure-Clinical-Records/blob/master/heart_failure_clinical_records_dataset.csv) para abrirmos aqui. Primeiro importaremos as bibliotecas necessárias e logo depois usaremos a biblioteca [Pandas](https://github.com/pandas-dev/pandas) (sim, subentendo que você já conheça) para abrir os dados no formato de Dataframe.

```python
from sklearn.model_selection import train_test_split
import pycaret.classification as pc
import pandas as pd

data_url = "https://raw.githubusercontent.com/lorenzodenisi/Heart-Failure-Clinical-Records/master/heart_failure_clinical_records_dataset.csv"
data_df = pd.read_csv(data_url)
data_df.head()
```

Com os dados disponíveis, vemos se está tudo ok vendo o tamanho do Dataframe (299 amostras e 13 colunas) e se existe algum dado faltante (não deveria).

```python
data_df.shape
data_df.isnull().sum()
```

Agora é bom dar uma olhada nas colunas para saber qual delas vamos usar como classe de interesse.

```python
data_df.columns
```

No nosso caso vamos usar a coluna "DEATH_EVENT" como classe alvo, sempre é bom ver se estamos trabalhando com dados desbalanceados (apesar que não vamos levar isso em consideração agora).

```python
data_df["DEATH_EVENT"].value_counts().plot(kind="bar")
```

Bora dividir os dados em treino e teste. Esse teste a gente só vai usar no final pra rodar a predição do modelo.

```python
train, test = train_test_split(data_df, test_size=0.2, random_state=42)
```

## Usando PyCaret, finalmente

Como visto antes importamos o `pycaret.classification` por estarmos trabalhando com um problema de classificação, para outros problemas como regressão por exemplo seria uma outra importação, para mais informações sobre isso dá uma olhada na [documentação](https://pycaret.org/setup/). Para ver os métodos e atributos disponíveis:

```python
dir(pc)
```

Por estarmos usando o Google Colab é preciso fazer uma ativação do ambiente:

```python
from pycaret.utils import enable_colab
enable_colab()
```

Agora podemos iniciar a brincadeira, para isso passamos os dados e a coluna alvo para o método `setup`:

```python
clf = pc.setup(data=train, target="DEATH_EVENT")
```

Aqui você precisa dar enter se os tipos das colunas foram inferidos corretamente, senão é possível modificar após escrever `quit` a partir de parâmetro no `setup`, mais informações [aqui](https://pycaret.org/data-types/). Além disso é possível aqui ignorar alguma coluna que você não queira usar como feature para os modelos, basta:

```python
pc.setup(data=train, target="DEATH_EVENT", ignore_features=["ages", "diabetes"])
```

Não é o nosso caso ignorar colunas, então vamos seguir com o comando anterior.

## Comparando modelos

Esse é o primeiro passo recomendado em qualquer experimento supervisionado. É bem simples e muito útil, ele mostra um relatório com cada modelo e as métricas obtidas. Com o parâmetro `sort` você pode especificar a métrica utilizada para ordenar os modelos, por definição é `accuracy` para classificação (`R2` para regressão) mas vou utilizar `F1` aqui.

Essa função treina todos os modelos disponíveis usando os parâmetros predefinidos e avalia seu desempenho usando validação cruzada, `fold=10`. Ela retorna o melhor modelo de acordo com a métrica estabelecida no `sort`. Alguns algoritmos são excluídos da análise pelo seu tempo de processamento, para evitar sua exclusão basta adicionar `turbo=False`.

```python
model = pc.compare_models(sort="F1", turbo=False)
```

É possível adicionar uma lista com os modelos que você gostaria de excluir no processamento com o parâmetro `exclude`, mas não é o nosso caso. Mais informações na [documentação](https://pycaret.org/compare-models/).

```python
pc.compare_models(exclude=["svm"])
```

Os nomes dos modelos estão na primeira coluna, na versão anterior era necessária usar a biblioteca `neatutils` para saber como ficava as abreviações.

## Criando um modelo, a partir do vencedor

Aqui iremos de fato trabalhar o modelo vencedor da comparação anterior. Para isso usamos o `create_model` passando a abreviação do modelo vencedor ou a variável contendo a comparação que retornou o melhor modelo de acordo com a métrica escolhida, no nosso caso o F1-score.

Como visto anteriormente o modelo vencedor para os nossos dados foi o Random Forest. No caso de classificação essa função recebe apenas o modelo e retorna uma tabela com os k-fold da validação cruzada junto com o modelo treinado. O número de folds pode ser definido pelo parâmetro `fold`, mas por pré-definição teremos `fold=10`. A quantidade de valores decimais também pode ser ajustada pelo `round`, sendo quatro decimais por pré-definição. Aqui ainda é possível usar um parâmetro de `ensemble`, mais detalhes na [documentação](https://pycaret.org/create-model/).

```python
best_model = pc.create_model(model)
```

## Tunando os hiperparâmetros

Com a função `tune_model`, os hiperparâmetros são tunados usando o Random grid search com grids predefinidos que podem ser customizados. Para aprendizado supervisionado, essa retorna uma tabela com k-fold validação cruzada com as métricas de avaliação junto com o modelo treinado.

```python
tuned_model = pc.tune_model(best_model)
```

O número de `fold` pode ser alterado, mas por pré-definição é de 10. A quantidade de decimais também pode ser alterada pelo `round`, sendo de 4 por pré-definição. O número de iterações aleatórias feitas pelo Random grid search pode ser mudado pelo parâmetro `n_iter`, sendo de 10 por pré-definição. Aumentando esse valor o tempo de treinamento também aumenta. Além disso, é possível adicionar a métrica que deve ser otimizada pelo parâmetro `optimize`, para classificação é `accuracy` e regressão `R2`.

```python
tuned_model_opt = pc.tune_model(best_model, optimize="F1", n_iter=25)
```

Se comparado, vemos que sem passar os parâmetros já tínhamos o "melhor" modelo. Seguiremos com ele então.

## Avaliando o modelo

Para avaliar o modelo você pode usar o `evaluate_model` que mostrará todas as métricas e gráficos associados ao modelo em uma interface com abas para seleção ou você pode usar o `plot_model` e ir selecionando cada métrica ou gráfico individualmente.

```python
pc.evaluate_model(tuned_model)
```

Por pré-definição o `plot_model` mostra o gráfico de AUC, mas você pode especificar o plot do seu interesse usando o parâmetro `plot`. Mais detalhes de possibilidades na [documentação](https://pycaret.org/plot-model/).

```python
pc.plot_model(tuned_model)
```

```python
pc.plot_model(tuned_model, plot="error")
```

## Otimizando os limites (FN e FP)

Para problemas de classificação o custo de resultados falso-positivos (FP ou erro tipo 1) é normalmente diferente dos falso-negativos (FN ou erro tipo 2). Assim, se você está trabalhando em um projeto em que os erros tipo 1 e tipo 2 possuem impacto diferentes, é possível otimizar sua classificação para um valor limite (threshold) de probabilidade melhorando a função de perda (loss function) apenas definindo o custo de FP e FN separadamente.

Para isso usamos a função `optimize_threshold` que pega o modelo treinado e a função de perda representada pelos FP, FN, TP (verdadeiro-positivos) e TN (verdadeiro-negativos). Essa função retorna um gráfico interativo onde a função de perda no eixo-y é representada por uma função de probabilidade de diferentes valores limites no eixo-x. Uma linha vertical representa o melhor valor da probabilidade limite para seu classificador. Mais detalhes na [documentação](https://pycaret.org/optimize-threshold/).

```python
opt_prob = pc.optimize_threshold(tuned_model, true_negative=1500, false_negative=5000)
```

Não sei o motivo, mas para mim não apareceu nenhuma imagem. Testa aí e me diz se foi.

## Interpretando o modelo

Aqui ele chama de interpretar o modelo, pra mim nada mais é que quais features o modelo deu mais importância. Talvez precise instalar a biblioteca `shap` que é uma abordagem de teoria dos jogos, mais info [aqui](https://github.com/slundberg/shap).

```python
!pip install shap
pc.interpret_model(tuned_model)
```

## Finalizando o modelo

Esse é o último passo de um experimento supervisionado já que quando iniciamos o projeto usando o `setup` parte dos dados não é utilizado no treinamento. Com ele o modelo é treinado uma última vez com os dados completos, mais detalhes na [documentação](https://pycaret.org/finalize-model/). Essa função recebe o modelo treinado e retorna o modelo treinado nos dados completos.

```python
final_model = pc.finalize_model(tuned_model)
```

## Salvando e abrindo o modelo

É bem simples de salvar seu modelo treinado, basta usar a função `save_model` passando o modelo treinado que será salvo todo o pipeline no formato pickle para usar posteriormente.

```python
pc.save_model(final_model, "rf_saved_model_02152021")
```

Já para abrir o modelo é mais simples ainda, basta usar a função `load_model` com o nome do modelo que você usou.

```python
loaded_rf = pc.load_model("rf_saved_model_02152021")
```

Como saída ele mostra o pipeline completo e os valores dos parâmetros do modelo.

## Usando o modelo para predição

Vamos fazer uma predição usando os dados de teste que separamos no começo da análise. Também usarei o modelo que foi salvo e reaberto.

```python
prediction = pc.predict_model(loaded_rf, data=test)
print(prediction.head())
```

Em Label vemos a predição feita pelo modelo final e na coluna vemos o DEATH_EVENT com a realidade. Um modelo simples que erra em alguns casos, mas acerta a maioria.

É possível melhorar esse modelo? Sim, mas por não ser o foco desse post recomendo olhar as referências e ver o que pessoal fez que pode melhorar os resultados desse projetinho.

Confesso que achei que ele seria útil nos dados que tou usando no meu projeto de mestrado. Porém, depois de 3 dias rodando ele estourou memória tentando rodar os meus dados /risos

```
MemoryError: Unable to allocate 1.04 TiB for an array with shape
(377799, 377799) and data type float64
```

Então voltarei a forma clássica de criar modelos para esses pequenos dados que tenho trabalhado. Deixo aqui o material que utilizei pra brincar um pouco com a biblioteca com o intuito de ser útil para alguém.

## Referências

- [Creating the Whole Machine Learning Pipeline with PyCaret](https://towardsdatascience.com/creating-the-whole-machine-learning-pipeline-with-pycaret-db39a3006840)
- [Machine Learning Made Easy by PyCaret](https://towardsdatascience.com/machine-learning-made-easy-by-pycaret-5be22394b1ac)
- [Machine Learning with PyCaret in Python](https://blog.jcharistech.com/2020/07/03/machine-learning-with-pycaret-in-python/)
- [A Gentle Introduction to PyCaret for Machine Learning](https://machinelearningmastery.com/pycaret-for-machine-learning/)
- [PyCaret 2.x — A biblioteca de aprendizagem de máquinas para quem tem prazo](https://medium.com/ensina-ai/pycaret-a-biblioteca-de-aprendizagem-de-m%C3%A1quinas-para-quem-tem-prazo-1c5b09667763)
- [Talks # 7: Moez Ali: Machine learning with PyCaret](https://youtu.be/jlW5kRBwcb0)
- [Machine Learning with PyCaret and Python](https://youtu.be/cnxOGWtwdv8)

---

*Publicado originalmente no [Medium](https://medium.com/@geocarvalho/uma-introdu%C3%A7%C3%A3o-ao-automl-com-pycaret-82a8d1eec83).*
