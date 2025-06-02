import numpy as np
import pandas as pd
from TensoesRadiais import TensoesRadiais
from PosCalculo import PosCalculo


class OneLoadCalc:
    """
    Classe responsável por calcular tensões e deformações em coordenadas radiais
    para um sistema multicamadas a partir de pontos definidos em coordenadas cartesianas,
    mas retornando os resultados em coordenadas cartesianas.

    OBS.: CALCULA APENAS UMA CARGA CIRCULAR

    """

    def calcular_theta(self, x, y, graus=False):
        """
        Calcula o ângulo theta entre o vetor (x, y) e o eixo x, no sentido anti-horário.

        Args:
            x (float): Coordenada x do ponto.
            y (float): Coordenada y do ponto.
            graus (bool): Se True, retorna o ângulo em graus. Caso contrário, em radianos.

        Returns:
            float: Ângulo theta em graus ou radianos.
        """

        theta = np.arctan2(y, x)
        if theta < 0:
            theta += 2 * np.pi  # Converte para intervalo [0, 2π)
        return np.degrees(theta) if graus else theta

    def preparacao_pontos(self, Thickness, main_points):
        """
        Converte os pontos de coordenadas cartesianas para coordenadas radiais
        e determina a camada a que cada ponto pertence.

        Args:
            Thickness (list or np.ndarray): Espessura das camadas (exceto subleito).
            main_points (list): Lista de dicionários contendo 'x', 'y' e 'z' dos pontos.

        Returns:
            list: Lista de dicionários contendo 'r', 'z', 'layer' e 'angulo' de cada ponto.
        """

        # Calculo da profundidade
        z_coord = [0] + [sum(Thickness[:index+1]) for index, val in enumerate(Thickness)] + [sum(Thickness)]

        # Preencher dicionario com informações necessárias para o cálculo
        Points = []
        for point in main_points:
            # Cálculo dos graus em relação ao eixo X no sentido anti-horário
            graus = self.calcular_theta(x = point['x'], y = point['y'], graus=True)

            # Determinação de qual camada pertence o ponto
            for i in range(1, len(z_coord)):
                if point['z'] >= z_coord[i-1] and point['z'] < z_coord[i]:
                    layer = i
                elif point['z'] >= z_coord[i] and i >= len(z_coord)-1:
                    layer = i

            # Preenchimento do dicionário de coordenadas radiais
            Points.append(
                {'r': (point['x']**2 + point['y']**2)**0.5,
                'z': point['z'],
                'layer': layer, 
                'angulo': graus})
            
        return Points

    def calc(self, main_points, n, Thickness, E, v, isbonded, a, q):
        """
        Executa o cálculo completo das tensões e deformações para pontos dados.

        Args:
            main_points (list): Lista de dicionários com coordenadas cartesianas dos pontos.
            n (int): Número de camadas (incluindo o subleito).
            Thickness (np.ndarray): Espessura das camadas (sem o subleito).
            E (np.ndarray): Módulo de elasticidade de cada camada (MPa).
            v (np.ndarray): Coeficiente de Poisson de cada camada.
            isbonded (bool): Define se as camadas são aderidas.
            a (float): Raio da carga aplicada (mm).
            q (float): Intensidade da carga (MPa, negativo para compressão).

        Returns:
            list: Lista de resultados convertidos para coordenadas cartesianas.
        """

        # Validações dos dados de entrada
        assert len(Thickness) == n - 1, "Thickness deve ter n-1 camadas (sem subleito)"
        assert len(E) == n and len(v) == n, "E e v devem ter comprimento n"

        # Preparação dos pontos - Cartesianas para radiais
        Points = self.preparacao_pontos(Thickness, main_points)

        # Cálculo das tensões e deformações em coordenadas radiais
        Results = TensoesRadiais().R(n, Thickness, E, v, isbonded, Points, a, q)

        # Convertendo coordenadas - Radiais para cartesianas
        Results_converted = PosCalculo().radiais_para_cartesianas(Results, Points)

        return Results_converted


class MultilayerElasticAnalysis:
    """
    Unifica o cálculo de tensõe e deformações para diversas cargas atuando em um mesmo ponto.
    """

    def call(self, main_points, n, Thickness, E, v, isbonded, load_points):
        temp_result = []
        for load in load_points:
            a = load['raio']
            q = load['carga']

            temp_points = []
            for point in main_points:
                temp_points.append(
                    {'x': point['x'] - load['x'],
                    'y': point['y'] - load['y'],
                    'z': point['z']})

            Results_converted = OneLoadCalc().calc(temp_points, n, Thickness, E, v, isbonded, a, q)

            if len(temp_result) == 0:
                temp_result = Results_converted.copy()
            else:
                for point in range(len(temp_result)):
                    for key in temp_result[point].keys():
                        temp_result[point][key] = temp_result[point][key] + Results_converted[point][key]

        for point in range(len(temp_result)):
            temp_dict, dir = PosCalculo().tensoes_principais_cilindricas(temp_result[point])
            temp_result[point]['s1'] = temp_dict['s1']
            temp_result[point]['s2'] = temp_dict['s2']
            temp_result[point]['s3'] = temp_dict['s3']

            # Preparação para cálculo das tensões octaédricas
            s1 = temp_result[point]['s1']
            s2 = temp_result[point]['s2']
            s3 = temp_result[point]['s3']
            txy = temp_result[point]['tau_xy']
            tyz = temp_result[point]['tau_yz']
            txz = temp_result[point]['tau_xz']

            # Cálculo da tensão octaédrica de cisalhamento
            #tensao_octa_cisalhamento = ((((s1 - s2)**2) + ((s2 - s3)**2) + ((s3 - s1)**2))**0.5)/3
            tensao_octa_cisalhamento = ((((s1 - s2)**2) + ((s2 - s3)**2) + ((s3 - s1)**2) + (((txy**2) + (tyz**2) + (txz**2))/6))**0.5)/3
            temp_result[point]['ocata_cisalhamento'] = tensao_octa_cisalhamento

            # Cáluclo da tensão octaédrica normal
            tensao_octa_normal = (s1 + s2 + s3)/3
            temp_result[point]['ocata_normal'] = tensao_octa_normal

        return temp_result


def main(main_points, n, Thickness, E, v, isbonded, load_points):
    final_results = MultilayerElasticAnalysis().call(main_points, n, Thickness, E, v, isbonded, load_points)
    df = pd.DataFrame(final_results)
    df['dels_vertical (mm 10^-2)'] = df['dels_vertical'] / (10**(-2))
    df = df[['sigma_x', 'sigma_y', 'sigma_z', 'tau_xy', 'epsilon_x', 'epsilon_y', 'epsilon_z', 'gamma_xy', 'dels_vertical (mm 10^-2)', 's1', 's2', 's3', 'ocata_cisalhamento', 'ocata_normal']]

    return df


if __name__ == "__main__":
    # Configurações da estrutura - Exemplo do livro do Yang H. Huang (Pavement Analysis and Design)
    n = 3                                             # Número de camadas (incluindo o subleito)
    E = np.array([3447.38, 137.90, 34.47])            # Módulo das camadas, incluindo subleito (MPa)
    Thickness = np.array([152.4, 304.8])              # Espessura das camadas, sem subleito (mm)
    v = np.array([0.35, 0.3, 0.45])                   # Poisson das camadas, incluindo subleito
    isbonded = True
    
    # Configuração dos pontos de análise
    main_points = [                                   # Coordenadas dos pontos a serem analizados (mm)
        {'x': 0, 'y': 0, 'z': 0},
        {'x': 0, 'y': 200, 'z': 0},
        {'x': 0, 'y': 300, 'z': 0},
        {'x': 0, 'y': 450, 'z': 0},
        {'x': 0, 'y': 600, 'z': 0},
        {'x': 0, 'y': 900, 'z': 0},
        {'x': 0, 'y': 1200, 'z': 0}
    ]

    main_points = [                                   # Coordenadas dos pontos a serem analizados (mm)
        {'x': -190.5, 'y': 0, 'z': 152.39},
        {'x': -190.5, 'y': 0, 'z': 457.45},
        {'x': 0, 'y': 0, 'z': 152.39},
        {'x': 0, 'y': 0, 'z': 457.45},
        {'x': -190.5, 'y': 254, 'z': 152.39},
        {'x': -190.5, 'y': 254, 'z': 457.45},
        {'x': 0, 'y': 254, 'z': 152.39},
        {'x': 0, 'y': 254, 'z': 457.45}
    ]

    # Configurações dos pontos de aplicação de carga
    load_points = [
        {'x': -381/2, 'y': -508,                      # Coordenadas dos pontos de aplicação da carga (mm)
         'raio': 101.6,                               # Raio de aplicação da carga (mm)
         'carga': -0.69},                             # Carga aplicada (MPa)
        {'x': -381/2, 'y': 0, 
         'raio': 101.6,
         'carga': -0.69},
        {'x': -381/2, 'y': 508, 
         'raio': 101.6,
         'carga': -0.69},

        {'x': 381/2, 'y': -508, 
         'raio': 101.6,
         'carga': -0.69},
        {'x': 381/2, 'y': 0, 
         'raio': 101.6,
         'carga': -0.69},
        {'x': 381/2, 'y': 508, 
         'raio': 101.6,
         'carga': -0.69}
    ]

    df = main(main_points, n, Thickness, E, v, isbonded, load_points)

    df.to_excel('teste.xlsx')

