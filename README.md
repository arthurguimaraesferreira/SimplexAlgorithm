# Trabalho Prático 2 - Pesquisa Operacional (DCC035) - TZ 2025/1  
 **Departamento de Ciência da Computação (DCC)  
 Universidade Federal de Minas Gerais (UFMG)  
 Aluno: Arthur Guimarães Ferreira  
 Professor: Marcio Costa Santos**

Este projeto implementa o **Algoritmo Simplex** utilizando o método **BIG-M** em linguagem C, com foco em desempenho para resolver problemas de Programação Linear no formato:

max c^t x  
st. Ax = b  
    x >= 0  

Alguns exemplos de execução se encontram na pasta 'Examples' e o enunciado do trabalho prático está no PDF.

Seguindo as orientações do enunciado, o programa deve ser compilado e executado com:  
*gcc -m64 -O3 FastSimplexBIGM.c -o FastSimplexBIGM*  
*./FastSimplexBIGM < ./Examples/example01.txt*  


Obs: a primeira versão do Algoritmo Simplex BIG-M foi implementada como no arquivo SlowSimplexBIGM.c.
Essa versão possui algumas redundâncias e cálculos desnecessários, que foram identificados e otimizados na versão seguinte.
O arquivo FastSimplexBIGM.c foi implementado a partir do SlowSimplexBIGM.c, eliminando operações repetitivas, melhorando o uso de memória 
e acelerando etapas críticas do algoritmo, resultando em uma execução significativamente mais rápida.
