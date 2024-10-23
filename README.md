# Trabalho de Ajustes de Curvas (Algoritmos Numéricos)

## Introdução

Este projeto faz parte da disciplina de Algoritmos Numéricos, com foco em Ajustes de Curvas. O objetivo do trabalho é implementar três modelos de ajuste de curvas para dados fornecidos e determinar qual deles é o melhor usando o coeficiente de determinação r².

## Objetivo

Dado um conjunto de dados sobre a quantidade de carbono-14 em amostras de obras históricas, o grupo deverá implementar os seguintes modelos de ajuste:

1. Ajuste Linear: N = β₀ + β₁ * t
2. Ajuste Quadrático: N = β₀ + β₁ * t + β₂ * t²
3. Ajuste Exponencial: N = β₀ * e^(β₁ * t)

Onde:
- N é a quantidade de carbono-14 na amostra
- t é a idade da amostra em anos

Com base no coeficiente de determinação r², o modelo que apresentar o melhor ajuste será determinado.

## Descrição do Problema

Os historiadores encontraram uma série de obras em uma biblioteca na Polônia e desejam usar dados de carbono-14 para estimar a idade das obras. Eles realizaram medições de carbono-14 em amostras de obras datadas e geraram uma tabela de dados que será usada como referência.

## Dados Utilizados

Os dados fornecidos consistem em valores de idade da amostra (em anos) e a quantidade de carbono-14 encontrada em cada uma delas. Esses dados serão usados para gerar os ajustes dos três modelos.

## Regressão Utilizada

Este projeto implementa três modelos de regressão:
1. **Ajuste Linear**: N = β₀ + β₁ * t
2. **Ajuste Quadrático**: N = β₀ + β₁ * t + β₂ * t²
3. **Ajuste Exponencial**: N = β₀ * e^(β₁ * t)

A implementação foi feita sem o uso de bibliotecas que realizem cálculos de regressão ou sistemas lineares automaticamente. Todo o cálculo foi implementado pelo grupo.

## Equipe
Este trabalho foi realizado em grupo. Os integrantes são:

1. Bernardo Krause Rodrigues
2. Caio Freire Silva
3. Gabriel Tetzner Menegueti