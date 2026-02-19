---
layout: post
title: "Introdução ao WDL — Parte 1"
date: 2021-02-28
categories: bioinformatics
tags: [wdl, bioinformatics, workflow, cloud]
---

O WDL (Workflow Description Language) é uma linguagem de script open-source que permite você especificar workflows de processos de uma forma fácil e com uma sintaxe simples (~palavras deles~).

Foi originalmente desenvolvida para a área de Genômica pelo instituto BROAD, mas pode ser usado em diversas outras áreas que trabalhem com dados.

WDL se lê "Widdle".

Não vou entrar em comparações entre linguagens de workflow, que são várias. Para quem tem interesse, existe essa [lista](https://github.com/common-workflow-language/common-workflow-language/wiki/Existing-Workflow-systems) feita pelo pessoal do CWL (Common Workflow Language). Particularmente, eu tenho mais interesse por [Nextflow](https://github.com/nextflow-io/nextflow) (pretendo escrever tutoriais usando isso posteriormente). Mas comecei usando [CWL](https://github.com/common-workflow-language/common-workflow-language), depois Nextflow, dei uma olhada em [Snakemake](https://github.com/snakemake/snakemake) e hoje tenho estudado WDL.

Isso depende muito do grupo que você está trabalhando e parecido com linguagens de programação você deve entender onde melhor aplicadas elas seriam. WDL por exemplo tem uma relação muito forte com o pessoal da ferramenta [GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360035889771-Pipelining-GATK-with-WDL-and-Cromwell) e o [Terra](https://support.terra.bio/hc/en-us/articles/360037117492-Getting-Started-with-WDL) (baseado no Google Cloud), assim já existem diversos scripts que você pode reutilizar, por exemplo no [BioWDL](https://github.com/biowdl) você já encontra alguns.

## Aprendendo WDL

Iremos trabalhar os tutoriais do "[learn-wdl](https://github.com/openwdl/learn-wdl)", para isso faça download do repositório e abra a pasta com seu editor de código.

Para o Visual Studio Code recomendo baixar a extensão "WDL DevTools".

```bash
git clone https://github.com/openwdl/learn-wdl.git
```

### Ambiente do WDL

Usei como referência esse material do learn-wdl. Eu estou trabalhando em cima do Ubuntu Xenial (16.04). Primeiramente você precisa ter os seguintes pré-requisitos:

**Java:**

```bash
java -version
# openjdk version "1.8.0_275"
# OpenJDK Runtime Environment (build 1.8.0_275-8u275-b01-0ubuntu1~16.04-b01)
# OpenJDK 64-Bit Server VM (build 25.275-b01, mixed mode)
```

Para instalar você pode usar os seguintes comandos:

```bash
sudo apt install openjre-8-headless
# ou
sudo apt update
sudo apt install default-jdk
```

**Docker:**

```bash
docker --version
# Docker version 17.06.1-ce, build 874a737
```

Para instalar o Docker:

```bash
curl -sS https://get.docker.com/ | sh
sudo usermod -aG docker $USER
```

Caso não funcione pra você procure na [documentação](https://docs.docker.com/engine/install/ubuntu/) do Docker.

**Cromwell e Womtool:**

```bash
sudo java -jar cromwell-56.jar --version
# cromwell 56

sudo java -jar womtool-56.jar --version
# womtool 56
```

Para instalar:

```bash
mkdir cromwell; cd cromwell
curl -L -o cromwell-56.jar https://github.com/broadinstitute/cromwell/releases/download/56/cromwell-56.jar
curl -L -o womtool-56.jar https://github.com/broadinstitute/cromwell/releases/download/56/womtool-56.jar
```

É bom copiar o caminho desses `.jar` para usarmos depois pra rodar os scripts.

Também é bom dar uma olhada no [github](https://github.com/broadinstitute/cromwell/releases) para obter a release mais recente das duas ferramentas.

## 1. Como fazer um "hello world"

Volte a pasta do learn-wdl. No caso eu fui para:

```bash
cd ../learn-wdl/1_script_examples/1_hello_worlds/1_hello
```

Vamos analisar o código `hello.wdl` que temos as três principais partes do WDL:

```wdl
version 1.0

workflow HelloWorld {
  call WriteGreeting
}

task WriteGreeting {
  command {
     echo "Hello World"
  }
  output {
     File output_greeting = stdout()
  }
}
```

- Linhas 7–15 nós temos a task que é um comando simples em bash para printar um "Hello" que vai ser um output com nome `output_greeting`.
- Linhas 3–5 temos a definição do workflow que simplesmente chama uma task definida.

## 2. Executando o script

- Cromwell é um serviço em Java que organiza a execução (job scheduler) e pode ser configurado para rodar em diferentes serviços (AWS, GCP, etc).
- WDL é o script que descreve o workflow.

Para começar vamos usar o Cromwell no modo Run que é normalmente utilizado para rodar protótipos, rodando workflows individuais. Podemos usar em uma máquina local (como fizemos anteriormente) ou em uma máquina virtual na nuvem.

Para trabalhar em produção normalmente se trabalha com o modo Server. Usando uma máquina na nuvem que vai ser utilizada como o organizador de execuções (Java scheduler), assim poderemos rodar múltiplos workflows.

### Então como rodamos o script?

```bash
sudo java -jar /caminho/para/cromwell-56.jar run hello.wdl
```

Verás que o output é bem verboso, que no meio em algum momento ele printou nosso "Hello World" com sucesso!

## 3. Sobre o modo Server

Quando vamos rodar nossas análises, nós normalmente utilizamos o modo Server. Você pode utilizar do jeito mais comum que é na "nuvem pública" ou no seu próprio HPC (High Performance Computing):

- **Terra.bio** — Cromwell já disponível no ambiente da GCP.
- **Azure** — conjuntos de VMs com Azure batch.
- **AWS** — conjuntos de VMs com AWS batch.
- **GCP** — conjunto de VMs (máquinas virtuais) com life sciences APIs.

A ideia é você usar grupos de VMs que são ativadas dinamicamente. Os criadores indicam usar o Terra.bio por ter sido criado pela BROAD. Para mais informações veja o tutorial do learn-wdl no YouTube [1].

## 4. Lidando com erros

Vamos analisar agora com o script `x_hello-error.wdl` que possui um erro.

- Para o modo Run que estamos usando aqui a ferramenta `womtool.jar` vai ajudar a detectar possíveis erros no nosso script.
- Primeiramente tenha certeza de estar usando a versão mais atual do Cromwell. Assim você terá certeza que as especificações de linguagem e features são suportadas.

Agora vamos tentar rodar esse script e ver o que ele retorna.

```bash
cd ../2_hello_docker
sudo java -jar /caminho/para/cromwell-56.jar run x_hello-error.wdl
# ...
# Workflow f7b32c3a-6c99-4fdc-b161-4781fde97a2d transitioned to state Failed
```

Verás que ao final do processo uma mensagem indicando um Failed. Mas se você for para o meio do processo verá:

```
Failed to process workflow definition 'HelloWorld' (reason 1 of 1):
Cannot generate outputs for 'call WriteGreeting'.
No such callable exists in [WriteGreetings]
```

Com essa informação já podemos resolver o problema do script, mas vamos ver o que o `womtool` nos mostrará.

```bash
sudo java -jar /caminho/para/womtool-56.jar validate x_hello-error.wdl
# Failed to process workflow definition 'HelloWorld' (reason 1 of 1):
# Cannot generate outputs for 'call WriteGreeting'.
# No such callable exists in [WriteGreetings]
```

Ele nos mostra o mesmo resultado, porém não precisamos ficar procurando no meio dos processos do comando do Cromwell, facilitando nossa vida. Existem outros comandos que podem ser vistos na [documentação](https://cromwell.readthedocs.io/en/stable/WOMtool/) ou usando o `--help`.

Abrindo o arquivo para resolver o bug:

```wdl
version 1.0

workflow HelloWorld {
  call WriteGreeting
}

task WriteGreetings {
  command {
    echo "Hello World"
  }
  output {
    File output_greeting = stdout()
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
```

Vemos a adição de um "s" na `task WriteGreetings` que deveria ser `WriteGreeting`. Resolvendo isso e rodando novamente obtemos sucesso no processo.

```bash
sudo java -jar /caminho/para/womtool-56.jar validate x_hello-error.wdl
# Success!
```

Para o modo Server no Terra.bio ao tentar fazer upload de um script no FireCloud ele vai mostrar o mesmo output do `womtool`.

Dependendo do ambiente na nuvem é preciso ter no script a descrição de `runtime` com os atributos `docker` que na primeira versão não tínhamos.

## 5. Variáveis no script

O script exemplo agora é o `~/learn-wdl/1_script_examples/1_hello_worlds/3_hello_input_task/input-task.wdl`

```wdl
version 1.0

workflow HelloInput {
  call WriteGreeting
}

task WriteGreeting {
  input {
    String name
  }
  command {
    echo 'hello ${name}!'
  }
  output {
    File response = stdout()
  }
  runtime {
    docker: 'ubuntu:latest'
  }
}
```

Primeiramente precisamos validar o script usando o `womtool-56.jar` como visto anteriormente:

```bash
sudo java -jar /caminho/para/womtool-56.jar validate input-task.wdl
# Success!
```

Feito isso, vamos tentar rodar ele usando o Cromwell:

```bash
sudo java -jar /caminho/para/cromwell-56.jar run input-task.wdl
# ...
# Workflow b13559c0-223d-4b08-9548-0a3a2ed1f111 transitioned to state Failed
```

Como vemos que falhou e subindo mais um pouco acharemos o motivo:

```
Required workflow input 'HelloInput.WriteGreeting.name' not specified
```

Então claramente ele quebrou porque não recebeu a variável esperada, nesse caso o `womtool` pode nos ajudar. Nele existe um modo `inputs` que indica se o script possui variáveis a serem dadas de entrada.

```bash
sudo java -jar /caminho/para/womtool-56.jar inputs input-task.wdl
# {
#   "HelloInput.WriteGreeting.name": "String"
# }
```

Esse comando nos retorna o que precisamos dar de entrada, no caso deveria ser um `inputs.json` no formato indicado pelo resultado do comando. Para facilitar nossa vida, podemos jogar isso para um arquivo e manipular `"String"` com o que queremos e então rodar o script novamente passando esse arquivo dessa vez.

```bash
sudo java -jar /caminho/para/womtool-56.jar inputs input-task.wdl > inputs.json
vim inputs.json

cat inputs.json
# {
#   "HelloInput.WriteGreeting.name": "Input"
# }

sudo java -jar /caminho/para/cromwell-56.jar run input-task.wdl --inputs inputs.json
# ...
# echo 'hello Input!'
```

Pensando no modo Server existe esse [artigo](https://support.terra.bio/hc/en-us/articles/360037485511-Add-Variables) no Terra super bem exemplificado de como usar variáveis lá. Além disso, esse [artigo](https://support.terra.bio/hc/en-us/articles/360037484851-Variable-Types-in-WDL) fala sobre os tipos de variáveis disponíveis que é bem útil na hora de você implementar seu pipeline mais na frente.

Tente fazer o mesmo com o script em `~/learn-wdl/1_script_examples/1_hello_worlds/4_hello_input_workflow/input-workflow.wdl`.

## 6. Arquivos e variáveis de ambiente

Agora trabalharemos em cima do script `~/learn-wdl/1_script_examples/1_hello_worlds/5_hello_file/hello-file.wdl`.

```wdl
version 1.0

workflow HelloFile {
    input {
        File file_input
        Int mem_gb
    }
    call InputFile {
        input:
            file_input=file_input,
            mem_gb=mem_gb
    }
}

task InputFile {
    input {
        File file_input
        Int mem_gb
    }
    command {
        bash echo 'The file is ${file_input}!' ${mem_gb}
    }
    output {
        File result = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
        memory: mem_gb + "GB"
        continueOnReturnCode: 126
    }
}
```

Aqui vamos usar o conceito anterior de passar um `.json` como entrada, no caso o arquivo que está na mesma pasta.

```bash
sudo java -jar /caminho/para/cromwell-56.jar run hello-file.wdl --inputs hello-file.json
# ...
# [warn] Local: Key/s [memory] is/are not supported by backend.
# Unsupported attributes will not be part of job executions.
```

Apesar de rodar com sucesso, temos essa observação de `warn` a ser analisada em meio a todos os processos printados na tela. Isso se chama memory assignment, só para mostrar que os parâmetros dependem do ambiente de execução.

Só pra lembrar que não precisa estar na pasta dos scripts para rodá-los, você pode passar o caminho+script que ele será rodado normalmente.

```bash
sudo java -jar /caminho/para/cromwell-56.jar run /caminho/para/hello-file.wdl \
    --inputs /caminho/para/hello-file.json
```

É comum usarmos arquivos que estejam salvos em buckets na nuvem, então é possível por esses arquivos no seu `.json` de entrada (com a devida permissão de acesso ao arquivo no bucket) para o script em WDL. Por exemplo:

```json
{
    "HelloFile.file_input": "https://storage.googleapis.com/wdl-quickstart/wdl-data/input.txt",
    "HelloFile.mem_gb": 5
}
```

Isso é interessante uma vez que vamos rodar pipelines em VMs que estarão disponíveis apenas para gerar os resultados do pipeline, então arquivos que devem persistir para uso no pipeline podem ser incorporados assim.

## 7. Scripts mais legíveis

Abrindo o script em `~/learn-wdl/1_script_examples/1_hello_worlds/6_alias_task/hello-again.wdl` vemos:

```wdl
version 1.0

workflow HelloWorldWithDocker {
    call WriteGreeting as wg1
    call WriteGreeting as wg2
}

task WriteGreeting {
    command {
        echo "Hello Docker"
    }
    output {
        File output_greeting = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}
```

Aqui é apenas para observar a possibilidade de criar alias com as `task` que usamos, no exemplo está repetida a `call` mas com alias diferentes que podemos ver no output quando rodamos o script.

```bash
sudo java -jar /caminho/para/cromwell-56.jar run /caminho/para/hello-again.wdl
# ...
# {
#   "outputs": {
#     "HelloWorldWithDocker.wg1.output_greeting": ".../call-wg1/execution/stdout",
#     "HelloWorldWithDocker.wg2.output_greeting": ".../call-wg2/execution/stdout"
#   }
# }
```

Existe outro exemplo em `~/learn-wdl/1_script_examples/1_hello_worlds/6_alias_task/hello-task.wdl` caso queira refazer o processo. Mas para esse script nós vamos rodar um script que o importará para fazer a chamada das `task` contidas, isso se chama subworkflow. O script maior é o `~/learn-wdl/1_script_examples/1_hello_worlds/6_alias_task/hello-workflow.wdl`

```wdl
version 1.0

import "hello-task.wdl" as HelloTask

workflow HelloWorldWithDocker {
    call HelloTask.WriteGreeting
    call HelloTask.WriteGreeting as wg
}
```

Rodamos esse script que tem duas `call` para a mesma `task` só para exemplificar que sem o `alias` essa chamada daria erro por estar repetida. Além de introduzir o conceito de `import` e poder reutilizar scripts menores em processos maiores.

```bash
sudo java -jar /caminho/para/cromwell-56.jar run /caminho/para/hello-workflow.wdl
# ...
# {
#   "outputs": {
#     "HelloWorldWithDocker.WriteGreeting.output_greeting": ".../call-WriteGreeting/execution/stdout",
#     "HelloWorldWithDocker.wg.output_greeting": ".../call-wg/execution/stdout"
#   }
# }
```

Use os scripts em `~/learn-wdl/1_script_examples/1_hello_worlds/7_subworkflow/hello-workflow.wdl` para treinar tudo que vimos sobre rodar workflow e inputs.

Aqui eu termino a primeira parte da introdução, de acordo com o curso do learn-wdl olhamos todos os arquivos da pasta `~/learn-wdl/1_script_examples/1_hello_worlds/`. Na parte dois que estou escrevendo olharemos a pasta `~/learn-wdl/1_script_examples/2_language_patterns/`, pelo menos essa é a ideia. Esse tutorial pode mudar a medida que eu for estudando, então se ligue nas referências.

## Referências

1. [What do I need to set up to write and execute WDL workflows?](https://support.terra.bio/hc/en-us/articles/360037493871-What-do-I-need-to-set-up-to-write-and-execute-WDL-workflows-)
2. [How to set up GCE VM for WDL](https://github.com/openwdl/learn-wdl/blob/master/5_reference_material/2_WDL-dev-env.md#wdl-dev-env)
3. [Learn WDL Course playlist](https://www.youtube.com/watch?v=RtcW2Zdn_28&list=PL4Q4HssKcxYv5syJKUKRrD8Fbd-_CnxTM)

---

*Publicado originalmente no [Medium](https://geocarvalho.medium.com/introdu%C3%A7%C3%A3o-ao-wdl-parte-1-ba01cf179db2).*
