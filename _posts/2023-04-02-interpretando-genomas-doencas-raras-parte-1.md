---
layout: post
title: "Interpretando genomas para doenças raras — anotações (parte 1)"
date: 2023-04-02
categories: bioinformatics
tags: [bioinformatics, genomics, rare-disease, sequencing, ngs]
image_dir: /assets/images/posts/sequenciamento-doencas-raras
---

Anotações da primeira aula do curso "Interpreting Genomes for Rare Disease" ministrada por Daniel MacArthur, PhD e adição de outros materiais, mais detalhes nas referências.

No campo do diagnóstico de doenças raras, as tecnologias de sequenciamento desempenharam um papel vital na identificação de mutações genéticas e genes causadores de doenças. No entanto, antes de mergulharmos nos componentes práticos do artigo, é crucial ter uma compreensão básica do processo de seleção de casos para sequenciamento. Este artigo abordará o processo de seleção, tanto para famílias que são mais propensas a fornecer um diagnóstico genético quanto para estudos de pesquisa interessados em identificar novos genes de doenças. Também discutiremos brevemente o processo de sequenciamento, as três diferentes tecnologias de sequenciamento com seus pontos positivos e negativos. O objetivo final é fornecer uma compreensão aprofundada das tecnologias de sequenciamento e seu papel no diagnóstico de doenças raras.

![Doenças que caem no espectro de arquitetura genética]({{ page.image_dir }}/01-espectro-arquitetura-genetica.png)
*Doenças que caem no espectro de arquitetura genética (Daniel MacArthur, 2019).*

Doenças raras são apenas uma ponta do espectro de doenças humanas e existem doenças que se encaixam em todo esse espectro de arquitetura genética. Por um lado, temos as doenças em que tradicionalmente trabalhamos no meu grupo: doenças monogênicas muito graves, que podem ser relativamente leves, mas são monogênicas no sentido de que segregam dentro de famílias específicas de maneira dominante, recessiva ou pelo cromossomo X. Elas são normalmente causadas por um número muito pequeno de variantes por família, ou seja, uma ou duas variantes, tipicamente em um único gene, que contribuem para aquela doença e normalmente seguem um padrão mendeliano. Por exemplo, no caso da família dominante mencionada, você tem um único pai afetado que transmitiu a variante da doença para 50% de seus filhos, mas existem algumas ressalvas a serem abordadas futuramente.

Por outro lado, no outro extremo do espectro de doenças, temos as doenças comuns e complexas, que são as que causam a maior mortalidade e morbidade em toda a população humana, como diabetes tipo 2, doenças cardíacas ou câncer. Essas doenças geralmente, mas nem sempre, são massivamente poligênicas, o que significa que em qualquer caso individual existe uma ampla gama de fatores de risco genéticos normalmente pequenos, bem como fatores de risco ambientais que contribuem para o risco daquela pessoa desenvolver a doença. Claro, em alguns casos extremos dessas condições, também existem formas monogênicas, e, em geral, as variantes causais para essas doenças tendem a ser substancialmente mais comuns, mas têm efeitos muito menores do que o que vemos para doenças raras.

É importante ressaltar que este artigo se concentrará principalmente nas doenças raras e monogênicas. Vamos focar em entender as tecnologias de sequenciamento que são usadas para identificar variantes causais para essas doenças e as melhores práticas para interpretar essas variantes.

![Frequência alélica vs tamanho do efeito]({{ page.image_dir }}/02-frequencia-alelica-efeito.png)
*Manolio et al. (2009), Nature. 461:747–753.*

As variantes genéticas que contribuem para o desenvolvimento de doenças seguem um espectro que pode ser representado em dois eixos. No eixo X, temos a frequência alélica da variante, onde as variantes extremamente raras estão localizadas no final esquerdo do espectro e as variantes mais comuns no final direito. Já no eixo Y, temos o tamanho do efeito, que representa a probabilidade de uma pessoa que carrega essa variante desenvolver a doença. As variantes que possuem grande probabilidade de causar a doença estão localizadas no topo do espectro e convergem para o conceito de penetrância. Já as variantes com menor efeito estão localizadas na parte inferior do espectro, que é o que normalmente se pensa em estudos de associação genômica em larga escala (GWAS).

No entanto, é importante ressaltar que no canto superior direito do espectro, onde estaríamos esperando encontrar variantes com grande efeito, quase não existem variantes comuns. Isso se deve ao fato de que, se houvesse uma variante que aumentasse significativamente o risco de doença, a seleção natural provavelmente a eliminaria da população. Por outro lado, no canto inferior esquerdo do espectro, existem variantes extremamente raras e com efeitos muito pequenos que são impossíveis de detectar com as amostras e tecnologias atuais. Para detectar essas variantes, precisaríamos ter milhões de amostras.

As variantes que podemos detectar e estudar estão localizadas no meio do espectro, onde estão as variantes raras com efeitos enormes associadas às doenças mendelianas no canto superior esquerdo e as variantes comuns com efeitos pequenos que são encontradas em estudos GWAS no canto inferior direito. À medida que aplicamos abordagens de sequenciamento de próxima geração em distúrbios complexos mais comuns, estamos começando a descobrir uma variedade de variantes no meio do espectro.

## Estratégia geral

Ao se deparar com casos de doenças raras, nossa estratégia geral é começar com uma família afetada pela doença. Geralmente, essa família já passou por extensos exames prévios, como painéis de genes ou sequenciamento de genes individuais, em busca de causas conhecidas da doença. No entanto, mesmo após anos de tentativas de diagnóstico sem sucesso, nossa primeira abordagem para a análise genética dessas famílias é sequenciar os exomas de todos os membros da família.

Esse método de sequenciamento de exoma (WES) é quase sempre a nossa primeira linha de abordagem para a descoberta genética. Isso ocorre porque, em muitos casos, cerca de 25% a 40%, dependendo da doença em questão, somos capazes de encontrar um diagnóstico usando essa técnica. Quando não conseguimos encontrar um diagnóstico, continuamos a analisar os dados do exoma periodicamente, geralmente a cada 6 a 12 meses, para verificar se alguma nova descoberta genética altera o diagnóstico para aquela família em particular.

Em casos em que não conseguimos encontrar um diagnóstico após vários anos de análise do exoma ou para casos em que há uma forte suspeita de que a condição seja de origem genética, mas o exoma é completamente silencioso, usamos o sequenciamento completo do genoma (WGS) e sempre que possível também realizamos o sequenciamento de RNA dos tecidos afetados.

## Estratégia para seleção das famílias para provável diagnóstico ou pesquisa de gene

No contexto do diagnóstico, vemos que o WES é realmente nosso teste de frente. Obviamente, painéis e várias outras formas de pré-seleção foram críticos em muitos contextos de diagnóstico, mas ir direto para o exoma ou, possivelmente, para um array, achamos que é agora uma abordagem realmente eficiente para ir o mais rápido possível em direção a um diagnóstico genético. Somos um pouco agnósticos em relação à tecnologia exata que usamos, mas em geral, se estivermos em uma situação em que o dinheiro é limitante — e vamos encarar, isso é típico, embora nem sempre — geralmente somos a favor do WES porque é a abordagem mais eficaz em termos de custo para encontrar o máximo de diagnósticos para um determinado gasto. Mas se estamos limitados pelas amostras — especialmente se temos uma família em que o caso específico é extremamente urgente ou clinicamente urgente em alguns casos, ou se temos uma doença específica em que temos um número muito pequeno de famílias que podemos obter — sabemos que queremos sequenciá-las e analisá-las o mais minuciosamente possível, então às vezes vamos direto para o WGS. Isso ocorre normalmente em uma situação em que as amostras são mais limitantes do que o dinheiro.

No contexto da descoberta de genes, nosso objetivo aqui é realmente obter o maior número possível de famílias afetadas por uma doença específica. Tentamos nos concentrar o máximo possível e aumentar nosso tamanho de amostras sem diluir muito a homogeneidade clínica dos casos. No contexto da descoberta de genes, investimos muito tempo em perguntar a colaboradores e famílias sobre pré-seleção, tentando obter famílias onde muito trabalho foi feito para descartar outras causas de doenças. Vamos nos concentrar em estruturas familiares solucionáveis. Isso é principalmente sobre evitar pequenas famílias dominantes e se concentrar em casos em que achamos que há uma chance razoável de realmente poder obter um diagnóstico com outra estrutura familiar que vemos. Quase sempre favoreceremos WES, pelo menos como uma abordagem de primeira linha na descoberta de genes, porque, novamente, a probabilidade de descobrir um novo gene por dinheiro gasto é substancialmente maior no contexto do WES.

### Tipos de famílias que focamos os estudos

![Tipos de famílias estudadas]({{ page.image_dir }}/03-tipos-familias.png)
*Quais tipos de famílias são estudadas (Daniel MacArthur, 2019).*

Estudo genético de famílias é uma das ferramentas mais poderosas para entender as causas de doenças hereditárias. Existem diferentes tipos de famílias que podem ser estudadas, cada uma com suas próprias vantagens e desafios. Um dos tipos de famílias mais interessantes é aquela com pais não afetados e vários filhos afetados. Isso sugere uma doença recessiva, o que muitas vezes pode ser resolvido com apenas uma única família. Identificar o gene candidato responsável pela doença pode ser relativamente fácil nesses casos. Outra explicação possível é mosaicismo e transmissão. No entanto, isso também pode ser resolvido com uma mutação de novo recorrente presente em ambos os filhos.

No entanto, o tipo de família mais comum é um trio com pais não afetados e um filho afetado. Nesses casos, é menos provável que possamos identificar o gene causal com apenas uma única família. Mas se a doença for causada por uma mutação de novo, temos uma boa chance de identificar o gene causal. Além disso, temos uma chance razoável de identificar causas recessivas.

Para famílias dominantes com doenças de início tardio ou compatíveis com a reprodução, geralmente precisamos de famílias substancialmente maiores para ter alguma chance real de descobrir novos genes da doença. Em geral, tendemos a nos concentrar em famílias dominantes muito maiores sempre que pudermos ter acesso a elas. Nessas situações, tentaríamos sequenciar os dois indivíduos afetados mais distantes geneticamente com WES e também tentaríamos realizar microarray para análise de linkage na família inteira.

O último tipo de design de estudo é a coleta de amostras de pacientes com condições de início tardio, pois muitas vezes os pais não estão disponíveis ou são muito difíceis de se obter acesso. Para doenças como distrofia muscular de Lynn Girdle e algumas doenças de retina, estamos construindo uma coleção de amostras de muitos indivíduos afetados pela mesma condição clínica. Em seguida, comparamos centenas de casos com milhares de controles que foram sequenciados na mesma plataforma e procuramos genes que tenham um excesso de mutações patogênicas raras nos casos em relação aos controles.

Em resumo, a escolha do tipo de família ou design de estudo depende da doença que está sendo estudada e da disponibilidade de famílias maiores e amostras de pacientes. Cada tipo de família tem suas próprias vantagens e desafios, mas todos eles podem ser usados para entender melhor as causas de doenças genéticas.

### Outras considerações ao escolher famílias para estudo

Seleção de famílias para estudos de sequenciamento é uma questão complexa. Normalmente a preferência é por famílias com pais não afetados e múltiplos filhos afetados, sugerindo uma desordem recessiva. Nesses casos, geralmente é possível identificar o gene candidato que causa a doença. O segundo tipo de família mais comum é o trio, com pais não afetados e um filho afetado. Embora seja menos provável identificar a causa exata com uma única família, há uma boa chance de identificar mutações de novo. Para famílias dominantes com doenças de início tardio ou compatíveis com a reprodução, são necessárias famílias maiores para descobrir novos genes. Famílias pequenas com apenas uma transmissão tendem a ter pouca sorte em descobrir novos genes. Para doenças de início tardio, está se tornando mais comum coletar grandes grupos de pacientes para análises.

Além disso, a regulamentação é crítica. As famílias precisam ter dado consentimento apropriado para pesquisa e compartilhamento de dados. A disponibilidade de DNA de outros membros da família também é importante em alguns casos. A apresentação clínica incomum pode ser um critério para selecionar uma família. Para WGS é preferível ter acesso a amostras de tecido relevantes para doença para análises de RNA, que podem fornecer informações importantes sobre o impacto funcional das variantes.

Uma grande consideração em nosso contexto de pesquisa, porque a rede é extremamente rigorosa em garantir que as amostras que passam por aqui tenham o consentimento apropriado tanto para sequenciamento quanto para compartilhamento de dados, é a questão regulatória. Portanto, temos as famílias sido adequadamente consentidas para pesquisa e, importantemente, outras consentidas para depósito de seus dados em bancos de dados de acesso controlado como [dbGaP](https://www.ncbi.nlm.nih.gov/gap/), para que outros pesquisadores possam realmente acessar seus dados.

Olhamos muito para a disponibilidade de DNA de outros membros da família e, para alguns projetos de família, isso é absolutamente crítico. Há um critério um pouco vago e incerto aqui, que é a singularidade da apresentação clínica e, tipicamente, descobrimos que, nos casos em que há uma apresentação fortemente sindrômica ou outras características incomuns desse caso em particular, temos muito mais chances de encontrar um novo gene de doença do que se parecer um caso padrão, por exemplo de distrofias musculares.

E o último ponto, para o qual voltarei mais tarde, é que preferimos muito, especialmente para WGS, amostras onde é possível obter acesso a amostras de tecido relevantes para a doença que possamos usar para RNA-seq e a razão para isso é que descobrimos que o RNA-seq nos dá ideias sobre o impacto funcional direto de variantes que afetam a expressão em splicing, o que muitas vezes pode ser muito difícil de prever apenas a partir da sequência de DNA.

> **Lição 1:** Pequenas análises baseadas em famílias geralmente são insuficientes para a descoberta de doenças geneticamente complexas. Isso ocorre porque, assim que você sai da causa monogênica muito estrita da doença, sua descoberta cai precipitadamente. Precisamos olhar para centenas ou possivelmente milhares de amostras para ter o poder de rejeitar essas causas mais complexas geneticamente da doença.

> **Lição 2:** Temos que ser muito cautelosos em relação a relatórios de triagem prévia. A precisão das triagens prévias para testes de painel, por exemplo, depende muito de quem fez a triagem, ou seja, qual laboratório estava envolvido na realização e execução desse teste de painel específico, qual tecnologia foi usada e quando foi feito.

## Como funciona o sequenciamento

![Sequenciamento por corrida]({{ page.image_dir }}/04-sequenciamento-corrida.png)
*Sequenciar bilhões de sequências por corrida de sequenciamento (Daniel MacArthur, 2019).*

Para entender o processo de sequenciamento de terceira geração (trabalhando aqui particularmente com short-read da Illumina), é importante saber que ele basicamente envolve a coleta de uma amostra de DNA de uma pessoa, a quebra dessa amostra em pequenos fragmentos ou leituras individuais de DNA, a leitura de cada um desses segmentos e, em seguida, a união desses segmentos para formar uma imagem da sequência subjacente do genoma ou exoma dessa pessoa. Na prática, se falando de short-read sequencing a maioria das leituras que obtemos pode ter de 100 a 150 pares de bases de comprimento, tanto para WES ou WGS. A partir dessas leituras, podemos formar um conjunto de dados que é representado por um mar de leituras.

![Comparação com genoma de referência]({{ page.image_dir }}/05-genoma-referencia.png)
*Comparar as sequências com o genoma de referência (Daniel MacArthur, 2019).*

Outro passo importante desse processo é o uso de uma sequência de referência. Essa sequência é definida como um padrão de referência usado universalmente. Atualmente, existem duas sequências de referência amplamente utilizadas na comunidade: a build 37 ([hg19 ou GRCh37](https://www.biostars.org/p/123767/)) e a build 38 ([GRCh38 ou hg38](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19)). Muitos estão fazendo a transição da build 37 para a build 38, que já está sendo utilizada em bancos de dados de referência como o [gnomAD](https://gnomad.broadinstitute.org/).

![NGS múltiplas leituras]({{ page.image_dir }}/06-ngs-multiplas-vezes.png)
*NGS nos permite sequenciar uma mesma posição múltiplas vezes (Daniel MacArthur, 2019).*

Para iniciar esse processo, é necessário utilizar uma sequência de referência padrão que representa uma média da sequência da população humana. A seguir, são utilizados algoritmos de alinhamento, como o [BWA](https://github.com/lh3/bwa), para alinhar todas as leituras individuais de DNA contra essa sequência de referência. Embora esse processo de alinhamento não seja simples e possa apresentar erros, quando realizado adequadamente, é possível identificar variantes de sequência comuns ou raras diretamente a partir dos dados de leitura. O número de leituras de DNA em diferentes partes da sequência de referência pode variar, e em algumas áreas é possível identificar uma base que difere pelo menos em algumas leituras daquela do indivíduo de referência. Isso é chamado de variantes ou polimorfismos e a medida de confiança de que essa variante é real é dada pela relação entre o número de bases de referência e as bases não referenciadas, chamada frequência alélica.

Mapear dados de leituras curtas em relação à sequência de referência é um grande desafio. Normalmente, as leituras possuem cerca de 150 pares de bases, o que pode resultar em ambiguidade na hora de alinhar as leituras com a sequência de referência, especialmente em regiões altamente repetitivas. Além disso, a cobertura do genoma pode variar bastante em diferentes partes do genoma, o que afeta a confiança na chamada de variante de um determinado local. A qualidade da posição tende a ser substancialmente pior para inserções e deleções (INDELs) em comparação com variantes de nucleotídeo único (SNPs). A sensibilidade para INDELs ainda está longe da perfeição, mesmo com dados de sequenciamento de boa qualidade, e a taxa de erro para INDELs é certamente muito maior do que para SNPs. Para piorar, a detecção de variantes estruturais (SVs), como deleções, duplicações (CNVs) e outros rearranjos de DNA em larga escala, é muito difícil.

![Tipos de arquivos]({{ page.image_dir }}/07-tipos-arquivos.png)
*Tipos de arquivos — FASTQ, BAM/CRAM, VCF (Daniel MacArthur, 2019).*

Em relação à interpretação de dados de sequenciamento, existem três tipos de arquivos importantes a serem considerados: [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map)/[CRAM](https://en.wikipedia.org/wiki/CRAM_(file_format)) e [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format). O arquivo FASTQ contém os dados brutos de leitura que são gerados pela máquina de sequenciamento. Uma vez que essas leituras são alinhadas à sequência de referência, obtemos um arquivo BAM ou CRAM, que consiste em uma lista de leituras ordenadas pela sua posição em relação ao genoma de referência. O último tipo de arquivo é o VCF, que contém todas as variantes identificadas no arquivo BAM ou CRAM.

![Comparação com bancos de dados]({{ page.image_dir }}/08-comparacao-bancos-dados.png)
*Fazer um genoma ter sentido requer a comparação com muitos outros (Daniel MacArthur, 2019).*

Para compreender as variações que surgem de cada um desses arquivos VCF, é necessário compará-los com bancos de dados populacionais como o gnomAD. Esses conjuntos de dados consistem em indivíduos de referência que foram agregados onde é possível verificar para qualquer variante encontrada no paciente, a frequência com que ela foi observada em uma coleção de dezenas de milhares de indivíduos — 140.000 indivíduos no caso do gnomAD — e o quão comum ela é em diferentes populações.

## Tecnologias disponíveis

### 1. Microarrays

![Genotipagem por Microarrays]({{ page.image_dir }}/09-microarrays.png)
*Genotipagem por Microarrays, região de interesse em vermelho (Daniel MacArthur, 2019).*

Para explorar as variações presentes no genoma de um indivíduo, existem várias tecnologias disponíveis. Uma dessas tecnologias é a genotipagem por meio de microarray. Esta tecnologia permite a investigação de posições específicas no genoma, geralmente posições que contêm variações comuns na população em geral. Os microarrays são extremamente baratos, com centenas de milhares de marcadores podendo ser testados por apenas ~$20. Podendo ser extremamente útil para análises de linkage em famílias. No entanto, os microarrays têm uma limitação importante: eles não detectam a maioria das variações raras, pois suas posições não incluem a maioria das variantes que são geralmente causadoras de doenças raras.

### 2. WES

![Sequenciamento de Exoma]({{ page.image_dir }}/10-wes.png)
*Sequenciamento de Exoma, região de interesse em vermelho (Daniel MacArthur, 2019).*

Uma tecnologia muito mais útil para o nosso contexto é WES que é uma tecnologia focada na parte codificante do genoma humano que conhecidamente dá origem a proteínas, representando cerca de 2% do total do genoma, mas contém praticamente todas as partes que realmente entendemos a função biológica. A grande vantagem do WES é seu custo, o preço para sequenciar o exoma pode chegar a menos de $200. Isso permite que você obtenha sequenciamento para todos os tipos de variantes, incluindo variantes muito raras presentes nas regiões codificantes de proteínas do genoma. O lado negativo do WES é que ele só fornece o que você está procurando, ou seja, apenas as regiões codificantes do genoma, sendo todo o resto descartado.

![Target sequencing]({{ page.image_dir }}/11-target-sequencing.png)
*Sequenciamento com alvo (target sequencing), região de interesse em vermelho (Daniel MacArthur, 2019).*

A ideia básica do WES é usar beads ou alguma outra sonda que contenha sequências individuais que se hibridizam apenas com as partes de sequência que estamos interessados em sequenciar, neste caso, as regiões codificantes de proteínas.

![Cobertura do exoma]({{ page.image_dir }}/12-cobertura-exoma.png)
*A cobertura do exoma é alta mas possui alta variabilidade (Daniel MacArthur, 2019).*

Um dos desafios na análise de dados de sequenciamento de exoma é a grande diferença na cobertura de um exon para outro. Podemos ver que a cobertura média dessas regiões é de 50X a 60X, mas há uma enorme variação. Alguns exons têm cobertura de até 150X, enquanto outros têm cobertura muito baixa ou nenhuma cobertura em uma amostra específica.

![Cobertura do genoma]({{ page.image_dir }}/13-cobertura-genoma.png)
*A cobertura do genoma é bem uniforme quando comparado ao exoma (Daniel MacArthur, 2019).*

Para comparar, é assim que a cobertura de sequenciamento se parece nos dados de sequenciamento do genoma. A cobertura, em média, é menor do que nas sequências do exoma, mas a uniformidade é maravilhosa. Há uma faixa muito estreita de cobertura, onde quase todos os éxons no genoma tendem a ser cobertos dentro de uma precisão bastante apertada.

#### Por que WES como primeira alternativa?

Então, vale a pena perguntar por que fazemos sequenciamento de exoma. A resposta óbvia é que é muito barato e podemos fazer um exoma familiar de trio por menos do que o custo de um único WGS hoje em dia. A cobertura mais alta dos segmentos principais significa que temos, talvez, 60X ou às vezes uma cobertura ainda maior em exons individuais. Mas a grande vantagem do sequenciamento de exoma é que ele nos fornece cobertura custo-efetiva das partes do genoma que podemos realmente interpretar. E quase todas as mutações conhecidas de doenças raras são encontradas nessas regiões codificadoras de proteínas. Achamos que é provável que mais de 85% de todas as variantes causais verdadeiras para condições mendelianas caiam dentro do exoma.

### 3. WGS

![Sequenciamento de genoma completo]({{ page.image_dir }}/14-wgs.png)
*Sequenciamento de genoma completo, região de interesse em vermelho (Daniel MacArthur, 2019).*

O sequenciamento de todo o genoma de uma pessoa é substancialmente mais caro do que o sequenciamento do exoma. Ele cobre quase todas as 3 bilhões de bases do genoma humano, com exceção das partes altamente repetitivas que são desafiadoras de serem mapeadas. No entanto, o sequenciamento do genoma inteiro tem muitas vantagens, como a possibilidade de identificar facilmente variações no número de cópias (CNVs) e variações estruturais (SVs), além da capacidade de descobrir variantes não codificantes.

#### Exemplo de casos resolvidos por WGS

![Pacientes com DMD diagnosticados por WGS]({{ page.image_dir }}/15-dmd-pacientes.png)
*Pacientes com distrofia muscular de Duchenne que foram diagnosticados graças à análise de WGS (Daniel MacArthur, 2019).*

Vou dar um exemplo de um caso que só poderia ser resolvido usando o sequenciamento completo do genoma. Esse caso foi descoberto em um estudo sobre distrofia muscular de Duchenne, que é a forma mais comum de distrofia muscular. Ela afeta principalmente meninos, exceto em casos raros. A doença é causada por mutações no gene DMD, encontrado no cromossomo X.

![Inversão de éxon no gene DMD]({{ page.image_dir }}/16-dmd-inversao-exon.png)
*Caso sem diagnóstico com inversão de éxon do gene DMD (Daniel MacArthur, 2019).*

Na imagem acima, temos parte do gene DMD, que é um gene enorme no cenário genômico. Podemos observar duas coisas. Primeiro, há uma seção de leituras que está completamente ausente aqui. Então, neste menino, parece haver uma deleção, mas ela está no íntron do gene. No entanto, há também esse interessante grupo de leituras que estão se alinhando na orientação oposta ao que esperaríamos. Isso é característico de uma inversão.

![Resolução da inversão no DMD]({{ page.image_dir }}/17-dmd-inversao-resolvida.png)
*Caso de inversão do éxon no gene DMD (Daniel MacArthur, 2019).*

Uma vez que resolvida a posição das leituras individuais, fica claro o que aconteceu: pequenas deleções de um lado e uma maior em verde. Além disso, o éxon encontrado nessa região foi invertido, resultando na sua exclusão da transcrição. Isso resultou em uma mudança no quadro de leitura, levando à perda completa do RNA e da proteína no paciente em questão. Apenas com o WGS é possível detectar rearranjos como este.

#### WGS não resolve tudo

Com o tempo descobrimos que mesmo sequenciando todo o genoma, nem sempre conseguimos solucionar todos os casos que sequenciamos. O WGS completo fornece apenas um aumento de 5% a 10% no poder de diagnóstico para casos que já tiveram WES. Para a grande maioria dos casos que foram negativos para o exoma, o sequenciamento do genoma não fornecerá a resposta.

### 4. RNA-seq

![Sequenciamento de RNA]({{ page.image_dir }}/18-rna-seq.png)
*Sequenciamento de RNA, região de interesse em vermelho (Daniel MacArthur, 2019).*

A tecnologia de sequenciamento de RNA, embora ainda não amplamente utilizada em diagnósticos, é muito poderosa e esperamos que seja cada vez mais adotada. O sequenciamento de RNA é diferente das outras tecnologias mencionadas, pois se concentra no RNA mensageiro em vez dos componentes de DNA do genoma. Ele sequencia todas as regiões expressas do genoma e fornece informações diretas sobre o impacto de variantes genéticas na expressão gênica e splicing. Sua capacidade de detectar interrupções no splicing tem sido poderosa para identificar variantes de splicing não anotadas.

#### Dois casos não diagnosticados de miopatia relacionada ao colágeno

![Paciente com miopatia do colágeno VI]({{ page.image_dir }}/19-col6-paciente.png)
*Paciente com miopatia relacionada ao colágeno VI (Daniel MacArthur, 2019).*

A miopatia do colágeno é uma doença rara em que os distúrbios nos genes do colágeno seis (COL6), geralmente dominantes, mas às vezes recessivos, resultam em uma variedade de fenótipos musculares relativamente graves. O médico da paciente na imagem tinha certeza do diagnóstico de miopatia do colágeno 6, mas não conseguiram identificar uma variante causal, apesar de terem feito o sequenciamento de todas as três cadeias do colágeno 6, WES e WGS.

![Gráfico de sashimi COL6A1]({{ page.image_dir }}/20-sashimi-col6a1.png)
*Gráfico de sashimi onde RNA-seq identifica um ganho de splicing intrônico no gene COL6A1 (Daniel MacArthur, 2019).*

Foram obtidas amostras de RNA muscular desses dois pacientes e sequenciado o RNA. Na última fileira do gráfico de sashimi, vemos o RNA-seq de um controle em uma região do genoma que contém três éxons do gene COL6A1. O que chamou a atenção em ambos os casos diagnosticados com colagenopatia 6 foi que ambos tinham um novo pseudo éxon que aparecia no meio de um intron e que não foi visto em nenhum outro indivíduo normal para quem o sequenciamento de RNA muscular foi feito — cerca de 180 pessoas.

![Variante GC>GT]({{ page.image_dir }}/21-variante-gc-gt.png)
*A variante GC>GT gera um sítio doador de splicing (Daniel MacArthur, 2019).*

Ao analisarmos a sequência do genoma nessa região específica, encontramos em ambos os pacientes uma mutação de GC para GT, resultando em uma variante doador de splice GT no meio de um íntron. Essa mutação GT resultou na criação de um novo local de splice e, portanto, na criação de um novo exon no meio da sequência intrônica. Curiosamente, essa inclusão de exon parece ser específica para o músculo e essa variante se tornou a mutação mais comum nos casos não diagnosticados de miopatia do colágeno 6.

### Quando considerar RNA-seq?

O sequenciamento de RNA é mais útil em casos onde há um gene candidato forte, mas a causa genética ainda não está clara. Especialmente em casos de famílias recessivas, onde há apenas uma mutação em um gene recessivo que se encaixa perfeitamente no fenótipo, mas não há uma segunda variante codificadora que explique o caso daquela família em particular.

Um segundo fator chave é a disponibilidade de tecido. É realmente importante sequenciar o tecido certo. Outra coisa a considerar é que muitas das mutações causais descobertas com RNA-seq afetarão o splicing, muitas vezes fora das regiões de splicing canônicas. Portanto, é útil ter dados de WGS em paralelo com o RNA-seq.

## Alguns pontos a serem discutidos

### Será que o WGS vai substituir WES como teste de entrada?

WES é mais eficiente em termos de custo do que o WGS, vale a pena pensar sobre o momento em que WGS se tornará a tecnologia mais eficiente em termos de custo. Acreditamos que isso provavelmente acontecerá em algum momento nos próximos anos, embora vale mencionar que à medida que o preço do WGS diminui, o preço do WES continua a acompanhar.

### Quais classes de genes relacionados a doenças podem ser descobertas apenas por WGS?

Exemplos disso poderiam ser casos em que todas as mutações causais para aquela doença específica são em regiões do genoma que não codificam proteínas por algum motivo, ou podem estar em alguma região que é inacessível por WES.

### Como lidar com a descoberta de genes em casos de doenças em que a causa genética subjacente não é monogênica?

Especialmente em condições genéticas complexas, onde é provável que haja vários genes que contribuam para esse fenótipo. Nesse caso, a única resposta é o tamanho da amostra. Só com grupos de mais de 1000 indivíduos é que poderíamos começar a obter o poder necessário.

### Quais tipos de tecidos são uma boa opção para RNA-seq?

Um caso a se considerar é a utilização de RNA-seq de alta profundidade do sangue. Além disso, sempre que possível, obter fibroblastos seria muito útil. Eles parecem ser um tipo de tecido mais útil em geral do que o sangue. Acredita-se que isso seja especialmente verdadeiro para condições neurológicas, em que fibroblastos parecem ser um melhor substituto para o cérebro do que o sangue.

## Referências

1. Manolio, Teri A., et al. "Finding the missing heritability of complex diseases." Nature 461.7265 (2009): 747–753.
2. [Interpreting Genomes for Rare Disease: Intro to Next Generation Sequencing — Daniel MacArthur, PhD, 2019](https://youtu.be/57McshGrllQ)
3. [Interpreting genomes for rare disease: variant and gene interpretation, 2019](https://cmg.broadinstitute.org/course-offering)

---

*Publicado originalmente no [Medium](https://geocarvalho.medium.com/introdu%C3%A7%C3%A3o-ao-sequenciamento-de-terceira-gera%C3%A7%C3%A3o-para-doen%C3%A7as-raras-44fd8df18989).*
