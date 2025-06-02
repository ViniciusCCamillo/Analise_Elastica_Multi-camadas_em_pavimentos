import numpy as np


# -------------------------- PONTOS EM ABERTO --------------------------
# 1 - Cisalhamento perto da superfície esta impreciso

# ----------------------------------------------------------------------

class PosCalculo:

    def convert_to_cartesian(self, stress_strain_cyl, theta):
        """
        Converte tensões e deformações de coordenadas cilíndricas (r, θ, z) para cartesianas (x, y, z).
        
        Parâmetros:
            stress_strain_cyl (dict): Tensões e deformações no sistema cilíndrico.
            theta (float): Ângulo em radianos.
        
        Retorna:
            dict: Tensões e deformações no sistema cartesiano.
        """
        # Cálculo de funções trigonométricas necessárias para a conversão
        cos_t = np.cos(theta)           # cosseno de θ
        sin_t = np.sin(theta)           # seno de θ
        cos2_t = cos_t**2               # cosseno ao quadrado
        sin2_t = sin_t**2               # seno ao quadrado
        sin2t = sin_t * cos_t           # seno vezes cosseno (usado na conversão de cisalhamento)

        # Extração dos componentes de tensão no sistema cilíndrico
        sigma_r = stress_strain_cyl['sigma_r']                       # tensão radial
        sigma_t = stress_strain_cyl['sigma_t']                       # tensão circunferencial (theta)
        sigma_z = stress_strain_cyl['sigma_z']                       # tensão axial (z)
        tau_rz = stress_strain_cyl['tau_rz']                         # tensão de cisalhamento entre r e z
        tau_rtheta = stress_strain_cyl.get('tau_rtheta', 0.0)        # cisalhamento entre r e θ (não fornecido)
        tau_thetaz = stress_strain_cyl['tau_thetaz']                 # cisalhamento entre θ e z (não fornecido)

        # Extração dos componentes de deformação no sistema cilíndrico
        epsilon_r = stress_strain_cyl['epsilon_r']                   # deformação radial
        epsilon_t = stress_strain_cyl['epsilon_t']                   # deformação circunferencial
        epsilon_z = stress_strain_cyl['epsilon_z']                   # deformação axial
        gamma_rz = stress_strain_cyl['gamma_rz']                     # Deformação radial de cisalhamento
        gamma_thetaz = stress_strain_cyl['gamma_thetaz']             # (não fornecido)

        # Conversão das tensões para coordenadas cartesianas
        sigma_x = sigma_r * cos2_t + sigma_t * sin2_t                # tensão normal na direção x
        sigma_y = sigma_r * sin2_t + sigma_t * cos2_t                # tensão normal na direção y
        tau_xy = (sigma_r - sigma_t) * sin2t                         # tensão de cisalhamento no plano xy
        tau_xz = tau_rz * cos_t - tau_thetaz * sin_t                 # tensão de cisalhamento no plano xz
        tau_yz = tau_rz * sin_t + tau_thetaz * cos_t                 # tensão de cisalhamento no plano yz

        # Conversão das deformações para coordenadas cartesianas
        epsilon_x = epsilon_r * cos2_t + epsilon_t * sin2_t          # deformação normal na direção x
        epsilon_y = epsilon_r * sin2_t + epsilon_t * cos2_t          # deformação normal na direção y
        gamma_xy = 2 * (epsilon_r - epsilon_t) * sin2t               # deformação de cisalhamento no plano xy
        gamma_xz = 2 * (gamma_rz * cos_t + gamma_thetaz * -sin_t)    # deformação de cisalhamento no plano xz (aproximação)
        gamma_yz = 2 * (gamma_rz * sin_t + gamma_thetaz * cos_t)     # deformação de cisalhamento no plano yz (aproximação)


        tau_xy = (sigma_r - sigma_t) * sin_t * cos_t                         # tensão de cisalhamento no plano xy
        gamma_xy = 2 * (epsilon_r - epsilon_t) * sin_t * cos_t               # deformação de cisalhamento no plano xy
        gamma_xz = gamma_rz * cos_t - gamma_thetaz * sin_t                   # deformação de cisalhamento no plano xz (aproximação)
        gamma_yz = gamma_rz * sin_t + gamma_thetaz * cos_t                   # deformação de cisalhamento no plano yz (aproximação)

        # Retorna um dicionário com os componentes convertidos para coordenadas cartesianas
        return {
            'sigma_x': sigma_x,
            'sigma_y': sigma_y,
            'sigma_z': sigma_z,
            'tau_xy': tau_xy,
            'tau_xz': tau_xz,
            'tau_yz': tau_yz,
            'epsilon_x': epsilon_x,
            'epsilon_y': epsilon_y,
            'epsilon_z': epsilon_z,
            'gamma_xy': gamma_xy,
            'gamma_xz_estimado': gamma_xz,
            'gamma_yz_estimado': gamma_yz,
            'dels_vertical': stress_strain_cyl['w'],
            'dels_radial': stress_strain_cyl['u']
        }

    def tensoes_principais_cilindricas(self, tensoes_cart):
        """
        Converte tensões cilíndricas para cartesianas e calcula as tensões principais.

        Parâmetros:
            stress_strain_cyl (dict): Tensões e deformações no sistema cilíndrico.
            theta (float): Ângulo em radianos.

        Retorna:
            tuple: (autovalores ordenados, autovetores correspondentes)
        """

        # Monta o tensor de tensões cartesiano simétrico
        sigma_tensor = np.array([
            [tensoes_cart['sigma_x'], tensoes_cart['tau_xy'], tensoes_cart['tau_xz']],
            [tensoes_cart['tau_xy'], tensoes_cart['sigma_y'], tensoes_cart['tau_yz']],
            [tensoes_cart['tau_xz'], tensoes_cart['tau_yz'], tensoes_cart['sigma_z']]
        ])

        # Calcula autovalores (tensões principais) e autovetores (direções principais)
        tensoes_principais, direcoes = np.linalg.eigh(sigma_tensor)

        # Ordena decrescentemente
        idx = np.argsort(tensoes_principais)[::-1]
        tensoes_principais = tensoes_principais[idx]
        direcoes = direcoes[:, idx]

        # Transforma a lista de tensões pricipais em dicionário
        tensoes_principais_dict = {
            's1': tensoes_principais[0],
            's2': tensoes_principais[1],
            's3': tensoes_principais[2]
        }

        return tensoes_principais_dict, direcoes

    def radiais_para_cartesianas(self, Results, Points):
        Results_converted = []
        # Transformação para coordendas cartesianas
        for line in range(len(Results)):
            # Cálculo do ângulo em radianos
            theta = Points[line]['angulo'] * (np.pi / 180)

            # Convertendo coordenadas radiais em cartesiandas
            stress_strain_cart = self.convert_to_cartesian(Results[line], theta)

            # Cálculo das tensões principais
            tensoes_principais_dict, direcoes = self.tensoes_principais_cilindricas(stress_strain_cart)

            # Unificando dicionários de tensões e deformações com o de tensões principais
            main_stress_strain = stress_strain_cart | tensoes_principais_dict

            # # Preparação para cálculo das tensões octaédricas
            # s1 = tensoes_principais_dict['s1']
            # s2 = tensoes_principais_dict['s2']
            # s3 = tensoes_principais_dict['s3']
            # txy = stress_strain_cart['tau_xy']
            # tyz = stress_strain_cart['tau_yz']
            # txz = stress_strain_cart['tau_xz']

            # # Cálculo da tensão octaédrica de cisalhamento
            # #tensao_octa_cisalhamento = ((((s1 - s2)**2) + ((s2 - s3)**2) + ((s3 - s1)**2))**0.5)/3
            # tensao_octa_cisalhamento = ((((s1 - s2)**2) + ((s2 - s3)**2) + ((s3 - s1)**2) + (((txy**2) + (tyz**2) + (txz**2))/6))**0.5)/3
            # main_stress_strain['ocata_cisalhamento'] = tensao_octa_cisalhamento

            # # Cáluclo da tensão octaédrica normal
            # tensao_octa_normal = (s1 + s2 + s3)/3
            # main_stress_strain['ocata_normal'] = tensao_octa_normal

            Results_converted.append(main_stress_strain)

        return Results_converted

