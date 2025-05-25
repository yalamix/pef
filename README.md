# Calculadora de Resistência dos Materiais

Um app simples feito com FastAPI e HTMX para resolver problemas da disciplina PEF3307 - Resistência dos Materiais.
Quando finalizado, em teoria será capaz de resolver todos os problemas vistos na disciplina, desde a P1 até a P3.
O conteúdo da P3 ainda tenho que revisar, adicionarei na lista de funcionalidades em breve.
Quando finalizá-lo, irei disponibilizá-lo na internet para qualquer um poder usar sem instalar nada.

## Funcionalidades
- Adicionar vínculos.
- Adicionar forças e momentos.
- Adicionar variáveis.
- Visualização gráfica do problema.
- Pontos de interesse.
- Rotacionar o sistema de coordenadas (em construção).
- Detectar automaticamente as condições de contorno e restrições (em construção).
- Determinar as forças, momentos, deformação e deflexão e gerar seus respectivos gráficos (em construção).
- Dimensionamento de viga (a fazer).
- Dimensionamento de seção circular (a fazer).
- Exportar o problema em arquivo JSON (a fazer).
- Importar o problema a partir de um arquivo JSON (a fazer).

## Como instalar e usar

Instale o Python versão 3.10.6 ou mais recente.
Clone usando git em um diretório de sua escolha ou baixe direto no Github em zip e extraia. 
Em seguida entre no diretório e execute os comandos em sequência:

```
python -m venv venv
```
```
venv\Scripts\activate.bat
```
```
pip install -r requirements.txt
```
```
cd app
```
```
fastapi dev
```

O app estará disponível no seu browser em http://localhost:8000/. Para fechar, basta fechar a janela do terminal. Para usar novamente, basta entrar na pasta app com seu terminal e executar fastapi dev.