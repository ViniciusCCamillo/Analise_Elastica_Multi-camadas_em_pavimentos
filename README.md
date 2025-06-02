# Multilayer Elastic Analysis with Multiple Loads in Pavements | Análise Elástica Multicamadas com Múltiplas Cargas em pavimentos

#  Análise Elástica Multicamadas com Múltiplas Cargas

Este projeto realiza o cálculo de **tensões** e **deformações** em sistemas multicamadas de pavimentos utilizando o modelo elástico linear, conforme abordado no livro **Pavement Analysis and Design** de **Yang H. Huang**.

Suporta múlticas cargas circulares atuando sobre a superfície e retorna os resultados convertidos para coordenadas cartesianas.

---

## ⚙️ Funcionalidades

- Cálculo de tensões normais (`σx`, `σy`, `σz`) e cisalhantes (`τxy`).
- Cálculo de deformações (`εx`, `εy`, `εz`) e distorções (`γxy`).
- Conversão entre coordenadas cartesianas e cilíndricas.
- Cálculo de tensões principais (`s1`, `s2`, `s3`).
- Tensão octaédrica normal e de cisalhamento.
- Suporte a múltiplas cargas em posições diferentes.
- Exportação dos resultados para `.xlsx`.

---

## 📁 Estrutura do Projeto

```text
.
├── OneLoadCalc.py               # Cálculo de tensões para uma única carga
├── MultilayerElasticAnalysis.py # Combinação de múltiplas cargas
├── PosCalculo.py                # Conversão de coordenadas e tensões principais
├── TensoesRadiais.py            # Núcleo de cálculo radial
├── main.py                      # Script de execução com exemplo do livro
├── teste.xlsx                   # Arquivo gerado com os resultados
└── README.md                    # Este arquivo
```

---

## 🧪 Exemplo de Uso

## 📌 Definição da estrutura
```
n = 3  # Camadas (incluindo subleito)
E = np.array([3447.38, 137.90, 34.47])  # MPa
Thickness = np.array([152.4, 304.8])  # mm
v = np.array([0.35, 0.3, 0.45])
isbonded = True   # Não aderido ainda não está funcional
```

## 📍 Pontos de análise
```
main_points = [
    {'x': 0, 'y': 0, 'z': 152.39},
    {'x': 0, 'y': 0, 'z': 457.45}
]
```

## ⚡ Cargas aplicadas
```
load_points = [
    {'x': -190.5, 'y': 0, 'raio': 101.6, 'carga': -0.69},
    {'x': 190.5, 'y': 0, 'raio': 101.6, 'carga': -0.69}
]
```

## ▶️ Execução
```
df = main(main_points, n, Thickness, E, v, isbonded, load_points)
df.to_excel("resultado_analise_elastica.xlsx", index=False)
```

---
## 📦 Requisitos

- Python 3.8+
- Bibliotecas: numpy, pandas e openpyxl (para exportação Excel)

## 📚 Referência
Huang, Y. H. (2004). Pavement Analysis and Design. 2ª edição. Pearson Prentice Hall.

## 🧑‍💻 Autor
Eng. Vinicius Coldebella Camillo
![FOTO PERFIL ](https://github.com/user-attachments/assets/2c8ca216-759a-4547-a602-29fb4535af99)

