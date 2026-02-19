---
layout: post
title: "Interpretando genomas para doenças raras — anotações (parte 2)"
date: 2023-04-08
categories: bioinformatics
tags: [bioinformatics, genomics, rare-disease, sequencing, ngs]
image_dir: /assets/images/posts/genomas-doencas-raras-parte-2
---

Anotações da segunda aula do curso "Interpreting Genomes for Rare Disease" ministrada por Anne O'Donnell-Luria, MD, PhD e adição de outros materiais, mais detalhes nas referências.

Nesse artigo vamos falar sobre o que é a anotação de dados advindos de sequenciamento, o problema dos transcritos que são muito importantes biologicamente, mas em termos de interpretação tornam as coisas mais desafiadoras. Sobre os bancos de dados de referência genômica, como o [gnomAD](https://gnomad.broadinstitute.org/about) e como os utilizamos. Também falaremos um pouco sobre previsões in silico. No final, vamos discutir outros bancos de dados que podem ser úteis na análise de pacientes com doenças raras.

Quando fazemos o sequenciamento de painéis de genes, exoma ou genoma o resultado final mais "bruto" dos pipelines de bioinformática normalmente é um arquivo de chamada de variantes, que é apenas uma lista de todos os locais onde a sequência difere do sequência de referência, que é a sequência genômica humana de referência que foi determinada como padrão. Existem muitos tipos diferentes de variantes que podemos encontrar. Essas variantes vão cair em diferentes transcritos e terão frequências diferentes na população em geral e podemos prever que terão diferentes consequências com base em previsões de algoritmos chamados [in silico](https://www.merriam-webster.com/dictionary/in%20silico).

![Tipos de variantes]({{ page.image_dir }}/01-tipos-variantes.png)
*Tipos de variantes, onde em rosa temos o mRNA e em roxo a proteína resultante (Elise Valkanas, 2019).*

Só para lhe dar um exemplo do que eu quero dizer com tipos de variantes, a tabela acima se restringe a variantes do tipo substituição, onde apenas um nucleotídeo é substituído por outro nucleotídeo. Essa é o tipo de variante mais comum, sendo aproximadamente 80% das variantes reportadas nos bancos de variantes públicos. Voltando a imagem, a sequência de referência é listada no topo. Dentro das substituições temos a variante missense é algo que muda o DNA, o RNA e o aminoácido. A variante sinônima também deve mudar o DNA e o RNA, mas não muda a proteína. Em seguida, a variante frameshift insere uma base no DNA e RNA, todos os aminoácidos seguintes mudem, e uma variante nonsense muda um aminoácido ou uma base que causa um códon de parada (stop codon) e trunca a proteína.

Falando um pouco mais sobre tipos de variantes ainda temos deleção, inserção, duplicação, deleção-inserção, inversão e variantes estruturais.

- As inserções (AAC-GTT > AACAGTT) são variantes onde um ou mais nucleotídeos são inseridos na sequência. Quando a sequência inserida é uma repetição (tandem copy) da sequência original de DNA, chamamos de duplicação. Mesma ideia, quando uma duplicação afeta mais de 50 pb nós nos referimos como CNVs.
- As deleções (AACGTT > AAC-TT) são variantes onde um ou mais nucleotídeos são deletados da sequência, sendo o segundo tipo de variante mais comum. Quando uma deleção afeta mais de 50 pares de bases (pb), nós referimos a ela como variante no número de cópias (Copy Number Variants, CNVs).

Ambas duplicação e deleções ocorrem frequentemente em pequenas regiões repetitivas do DNA.

- As variantes estruturais acabam sendo um termo comumente utilizado para grandes mudanças no cromossomo como translocações e transposições. É importante ressaltar aqui que esse tipo de variante é extremamente difícil de ser detectada utilizando sequenciamento de sequências curtas (short reads, tipo Illumina ou MGI), mesmo que existam algoritmos para isso o tamanho das sequências por si já são um fator limitante para detecção de grandes variantes estruturais e até mesmo as pequenas. Se a variante estrutural for muito grande é possível visualizar por meio de tecnologias de mapeamento óptico ou microscopia (cariótipo).
- As inversões (AACGTT > ACGTTT) ocorrem quando um pedaço da sequência se inverte, a nova sequência é exatamente o reverso-complemento da sequência deletada. Elas possuem um tamanho mínimo de dois nucleotídeos, sendo os casos de um uma simples substituição não é mesmo?!
- As deleção-inserção (AACGTT > AATATT) é a combinação de deleções seguidas de inserções na mesma localização do DNA, parecido com a substituição, um ou mais nucleotídeos são substituídos por um ou mais nucleotídeos.

Se é do seu interesse ver mais estudos sobre variantes estruturais eu tenho um repositório no GitHub com o que acabei achando de material sobre essas variantes e como analisá-las usando sequenciamento genético, só clicar [aqui](https://github.com/geocarvalho/sv-cnv-studies).

## Nomenclatura HGVS

Existem dois tipos de nomenclaturas sendo a primeira a que vemos nos arquivos VCF, mas a mais comum clinicamente é a nomenclatura da Human Genome Variation Society ([HGVS](https://www.hgvs.org/)). Ela utiliza a letra "c." para explicar a mudança no DNA e a letra "p." para a mudança no aminoácido. A nomenclatura HGVS e a VCF não se encaixam perfeitamente. A HGVS faz anotações em relação ao transcrito, então é muito importante saber qual transcrito você está se referindo. Já a VCF faz anotações baseadas nas coordenadas genômicas no genoma de referência utilizado e sempre se refere à fita positiva ou superior (fita 5') do DNA.

- As variantes em sítios de splicing são escritas como a posição "c." mais próxima e se é negativo ou positivo (c.345–2A>T).
- As variantes frameshift são escritas indicando quantas bases estão após a mudança até encontrar um códon de parada (p.Glu54fs ou p.Glu54Glyfs ou Glu54GlyfsX5).
- As variantes nonsense ou stop loss são representadas de três maneiras diferentes, mas todas indicam a introdução de um códon de parada (p.Arg24\* ou p.Arg24X ou p.Arg24Ter).
- Para uma variante INDEL descrevemos as bases que são deletadas com a notação "c." e a mudança real para a proteína. É importante mencionar que INDELs em frame são quando você deleta três bases, ou seja, um aminoácido, mas mantém os aminoácidos seguintes mantendo o quadro de leitura da proteína (c.1115\_1117delCCT, p.Ser372del).
- Para variantes missense ou não-sinônimas, se apresenta a posição do códon, como por exemplo a troca de uma isoleucina para valina (p.Ile123Val).
- Para variantes sinônimas ou silenciosas usamos principalmente a notação "p.", por exemplo, p.= ou p.Ile123= ou p.Ile123Ile para uma isoleucina que se mantém uma isoleucina mesmo com a troca do aminoácido.

![Comparação da anotação HGVS vs VCF]({{ page.image_dir }}/02-hgvs-vs-vcf.png)
*Comparação da anotação para a variante Phe508del no gene CFTR (Anne O'Donnell-Luria, 2019).*

- HGVS: NM\_000492.3(CFTR) c.1521\_1523delCTT ou chr7 g.117199646\_117199648del
- VCF: chr7 117199644 ATCT A

Novamente, a HGVS se refere ao transcrito e alinha tudo à direita (linha azul na imagem acima), enquanto a notação no VCF se alinha tudo à esquerda (linha vermelha na imagem acima), o que pode causar algumas discrepâncias entre os dois sistemas. Como exemplo, podemos usar a variante Phe508del no gene CFTR que é a variante mais comum para fibrose cística. Na anotação presente no VCF, ancoramos à esquerda na letra A, que não é deletada, e depois deletamos as letras CTT (linha vermelha na imagem acima). Já na nomenclatura HGVS, ancoramos à direita, e como sobram as duas letras T, não podemos dizer quais letras foram deletadas (linha azul na imagem acima). Por isso, na notação HGVS, é relatado que as letras CTT foram deletadas. Isso é importante porque se você estiver procurando em um banco de dados que usa HGVS e o banco de dados estiver usando a notação VCF, pode parecer que você está procurando variantes diferentes. Portanto, é crucial que seus dados estejam anotados de uma maneira que se alinhe com outros dados, para obter a anotação correta em cada variante. Assim, é comum termos ambas as anotações quando criamos bancos ou planilhas para análise.

## Genomas de referência

Os genomas de referência são outra complicação, onde temos uma sequência de referência do genoma humano consensual e pesquisadores trabalham na melhoria e preenchimento dessa sequência. Em sua maior parte, os genes não mudam, mas o espaço entre eles está sendo preenchido com detalhes adicionais. Além disso, haplótipos — regiões que variam bastante entre diferentes populações — estão sendo preenchidos, o que torna os genomas cada vez mais complexos. Recentemente foi publicado mais uma versão do genoma humano chamado telomere-to-telomere (T2T), resolvendo várias regiões complexas do genoma que ficaram indisponíveis por anos, aumentando o número de promessas de descobertas e de uso clinicamente [4–6].

Apenas para mostrar um exemplo, esta é a mesma variante no gene CFTR que falamos anteriormente:

- hg37: chr7 117199644 ATCT A
- hg38: chr7 117559590 ATCT A

HG37 é um genoma de referência que foi usado por muito tempo e ainda existem bancos de dados que são baseados apenas nele, enquanto que o HG38 é a referência atualmente recomendada clinicamente. No exemplo acima você verá que as alterações de base são as mesmas, mas as posições são diferentes. Elas estão, na verdade, separadas por 350 kilobases (KB). Se você olhar para uma variante no genoma de referência errado ela pode parecer estar em um gene ou localização totalmente diferente do que você esperava, por isso é muito importante saber em qual construção do genoma você está trabalhando.

Se você tiver o arquivo BAM ou CRAM é possível identificar o genoma que foi utilizado para o alinhamento olhando o cabeçalho do arquivo por meio de ferramentas de bioinformática como o [SAMtools](http://www.htslib.org/doc/samtools-view.html). Mas se você quer visualizar variantes, basta abrir o arquivo na ferramenta [IGV](https://software.broadinstitute.org/software/igv/) que se for o genoma errado você verá um alinhamento colorido que indica variantes na posição em comparação com o genoma:

![BAM alinhado no IGV com genoma errado]({{ page.image_dir }}/03-igv-genoma-errado.png)
*Exemplo de BAM alinhado no hg38 aberto na ferramenta IGV usando o genoma de referência hg19 (George Carvalho, 2023).*

## Métricas de qualidade para variantes

![Métricas de qualidade de variantes]({{ page.image_dir }}/04-metricas-qualidade.png)
*Exemplo de variante com algumas métricas de qualidade comumente usadas (Anne O'Donnell-Luria, 2019).*

Existem vários tipos de métricas de qualidade que são usadas para identificar variantes presentes em um arquivo VCF. No exemplo acima, temos um ID para a amostra (GCGS\_FAM8\_23), o nucleotídeo referência A (que não aparece em negrito) e o G que representa a variante alternativa. O conjunto A/G indica que o paciente é heterozigoto para a variante G na posição onde se espera um A. Se fosse homozigoto para a variante, seria G/G (ambos em negrito), e se fosse homozigoto para a referência, seria A/A.

- **VQSR** (Variant Quality Score Recalibration): Não está presente nesse exemplo, mas é uma métrica bem sofisticada baseada em modelos Gaussianos por mistura para avaliar se uma variante é mais parecida com uma variante real ou com um artefato. Quando olhamos para os dados, é possível encontrar uma lista de pontuações de VQSR para cada variante [7].
- **AB** (Allele Balance ou equilíbrio alélico): é outra característica importante e representa o número de leituras com a variante alternativa em relação ao número total de leituras cobrindo a posição onde se encontra. Neste exemplo, 31% das leituras apresentam a variante alternativa, o que está dentro do esperado para algo heterozigoto de linha germinativa — cerca de 50%. No entanto, vale lembrar que sempre há uma distribuição de leituras e, quanto menor a porcentagem, maior a chance de que essa variante seja somática neste paciente.
- **DP** (Depth ou profundidade): É o número de leituras que cobrem uma determinada posição onde a variante se encontra, aqui com 54 leituras nessa posição. Para WGS consideramos DP=10 e WES consideramos DP=100 como valores mínimos para se reportar uma variante.
- **GQ** (Genotype Quality ou qualidade do genótipo): Essa métrica varia de 0 a 99, sendo que 99 representa a qualidade máxima. Geralmente, olhamos para as variantes com GQ>20. Porém, é importante ressaltar que é possível ter variantes reais com pontuações mais baixas e artefatos com pontuações mais altas.

## Versões dos transcritos de um gene

![Versões de transcritos no gnomAD]({{ page.image_dir }}/05-transcritos-gnomad.png)
*Versões de transcritos para o gene CFTR no site gnomAD (George, 2023).*

Os transcritos são importantes na interpretação das variantes graças a sua localização no gene podendo ter diferentes impactos. É importante mencionar que existem dois bancos de dados principais para transcritos:

- [Ensembl](https://useast.ensembl.org/index.html), que é o EMBL-EBI ou GENCODE e mais voltado para anotações europeias.
- [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/), que é baseado no NCBI americano.

O RefSeq tende a ter menos opções de transcritos, enquanto Ensembl tende a ter mais. Geralmente, o transcrito canônico é o que todos nós falamos e nos referimos como o principal transcrito, mas isso geralmente é apenas porque é o maior e inclui todos os éxons. Isso não significa necessariamente que ele seja a transcrição clinicamente mais relevante. É importante lembrar que existem grupos trabalhando na definição dos transcritos clinicamente mais relevantes, como o projeto [MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/). Ele é um projeto colaborativo que tem como objetivo convergir na anotação de genes e transcritos humanos e definir um conjunto representativo de transcritos e proteínas correspondentes para genes codificantes de proteínas em todo o genoma humano [8]. Cada transcrito MANE representa uma correspondência exata nas regiões exônicas entre um transcrito no RefSeq e sua contraparte na anotação Ensembl/GENCODE, de modo que os dois identificadores possam ser usados simultaneamente. Além disso, um transcrito MANE corresponde perfeitamente ao genoma de referência GRCh38 e é escolhido com base em critérios biologicamente relevantes, como níveis de expressão do transcrito e a conservação das regiões codificadoras.

## Ferramentas de predição de impacto

![Ferramentas de predição de impacto]({{ page.image_dir }}/06-ferramentas-predicao.png)
*Exemplo de ferramentas de predição de impacto de variantes (Anne O'Donnell-Luria, 2019).*

No geral, muitas dessas ferramentas de predição de impacto (análise [in silico](https://www.news-medical.net/life-sciences/What-is-in-Silico.aspx)) usam os mesmos dados de treinamento e métricas semelhantes para os modelos, o que significa que há muita circularidade em usar vários preditores in silico. Por outro lado, é melhor ter mais dados do que menos. Por isso, tendemos a olhar para muitas ferramentas, mas não os consideramos de forma aditiva. Apenas pensamos um pouco sobre o consenso e alguns deles confiamos um pouco mais do que outros. No entanto, para um efeito específico de uma variante, é realmente difícil saber qual é o mais importante. Não vou entrar em detalhes sobre nenhum deles, mas direi que alguns são apenas para variantes missense e outros para todos os tipos de variantes e isso é útil apenas para ter uma pontuação para cada variante. Além disso, existe a ferramenta SpliceAI que foi desenvolvida pela Illumina para predição de junções de splicing, permitindo a predição do impacto de variantes não codificantes como as sinônimas e intrônicas que podem alterar o mecanismo de splicing [9].

## Conservação

*(Seção a ser expandida em futuras atualizações.)*

## Referências

1. [Jaganathan, Kishore, et al. "Predicting splicing from primary sequence with deep learning." Cell 176.3 (2019): 535–548](https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf)
2. [Morales, Joannella, et al. "A joint NCBI and EMBL-EBI transcript set for clinical genomics and research." Nature 604.7905 (2022): 310–315](https://www.nature.com/articles/s41586-022-04558-8)
3. [GATK technical documentation — Variant Quality Score Recalibration (VQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-)
4. [Aganezov, Sergey, et al. "A complete reference genome improves analysis of human genetic variation." Science 376.6588 (2022): eabl3533](https://www.science.org/doi/full/10.1126/science.abl3533)
5. [Alkan, Can, et al. "Implications of the first complete human genome assembly." Genome Research (2022)](https://qmro.qmul.ac.uk/xmlui/bitstream/handle/123456789/78958/Rowe%20Implications%20of%20the%20first%20complete%20human%20genome%20assembly%202022%20Published.pdf?sequence=2)
6. [Miga, Karen H., and Beth A. Sullivan. "Expanding studies of chromosome structure and function in the era of T2T genomics." Human Molecular Genetics 30.R2 (2021): R198-R205](https://www.science.org/doi/full/10.1126/science.abj6987)
7. [Patrinos, George P. Clinical DNA Variant Interpretation: Theory and Practice. Academic Press, 2021](https://www.elsevier.com/books/clinical-dna-variant-interpretation/patrinos/978-0-12-820519-8)
8. [Interpreting Genomes for Rare Disease: Variant and Gene Interpretation — Anne O'Donnell-Luria, MD, PhD, 2019](https://youtu.be/MmAYOgdouJc)
9. [Interpreting genomes for rare disease: variant and gene interpretation, 2019](https://cmg.broadinstitute.org/course-offering)

---

*Publicado originalmente no [Medium](https://geocarvalho.medium.com/interpretando-genomas-para-doen%C3%A7as-raras-anota%C3%A7%C5%8Des-parte-2-b7e85cc724b4).*
