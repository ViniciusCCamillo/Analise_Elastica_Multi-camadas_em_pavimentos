import numpy as np
from scipy.special import jv as besselj


class TensoesRadiais:

    def gerar_pontos_cisalhamento(self, ponto_central, Thickness, delta_r=0.01, delta_z=0.01):
        """
        Gera 9 pontos no entorno do ponto central (3x3 em r e z) para cálculo de cisalhamento gamma_rz.
        """
        r0 = ponto_central['r']
        z0 = ponto_central['z']

        # Calculo da profundidade
        z_coord = [0] + [sum(Thickness[:index+1]) for index, val in enumerate(Thickness)] + [sum(Thickness)]

        pontos = []
        for dr in [0, -delta_r, delta_r]:
            for dz in [0, -delta_z, delta_z]:
                new_z = z0 + dz if z0 + dz >= 0 else 0
                # Determinação de qual camada pertence o ponto
                for i in range(1, len(z_coord)):
                    if new_z >= z_coord[i-1] and new_z < z_coord[i]:
                        layer = i
                    elif new_z >= z_coord[i] and i >= len(z_coord)-1:
                        layer = i
                pontos.append({'r': r0 + dr, 'z': new_z, 'layer': layer})
        
        return pontos

    def R(self, n, Thickness, E, v, isbonded, Points, a, q):
        # Calculate the cumulative thicknesses z
        z = np.cumsum(Thickness)
        H = z[-1]  # Last value in z is the total depth (H)
        
        alpha = a / H
        dl = 0.1 / alpha  # Small increment for integration

        # Initialize the results dictionary
        Results = [{'sigma_z': 0, 'sigma_r': 0, 'sigma_t': 0, 'tau_rz': 0, 'w': 0, 'u': 0, 'tau_thetaz': 0} for _ in range(len(Points))]

        # Loop for integration
        for l in np.arange(0, 20, dl):  # Integration points
            for gauss_point in range(1, 5):  # 4-point Gauss quadrature
                # Gauss point logic
                if gauss_point == 1:
                    m = (l + dl / 2) - 0.86114 * (dl / 2)
                    fc = 0.34786 * (dl / 2)
                elif gauss_point == 2:
                    m = (l + dl / 2) - 0.33998 * (dl / 2)
                    fc = 0.65215 * (dl / 2)
                elif gauss_point == 3:
                    m = (l + dl / 2) + 0.33998 * (dl / 2)
                    fc = 0.65215 * (dl / 2)
                elif gauss_point == 4:
                    m = (l + dl / 2) + 0.86114 * (dl / 2)
                    fc = 0.34786 * (dl / 2)

                # Call R_hat function (replace this with your actual implementation)
                Result = self.R_hat(n, Thickness, E, v, isbonded, m, Points)

                if m != 0:
                    for j in range(len(Points)):
                        # Update Results with the calculated values
                        Results[j]['sigma_z'] += fc * (q * alpha * Result[j]['sigma_z'] / m * besselj(1, m * alpha))
                        Results[j]['sigma_r'] += fc * (q * alpha * Result[j]['sigma_r'] / m * besselj(1, m * alpha))
                        Results[j]['sigma_t'] += fc * (q * alpha * Result[j]['sigma_t'] / m * besselj(1, m * alpha))
                        Results[j]['tau_rz'] += fc * (q * alpha * Result[j]['tau_rz'] / m * besselj(1, m * alpha))
                        Results[j]['tau_thetaz'] += fc * (q * alpha * Result[j]['tau_thetaz'] / m * besselj(1, m * alpha))
                        Results[j]['w'] += fc * (q * alpha * Result[j]['w'] / m * besselj(1, m * alpha))
                        Results[j]['u'] += fc * (q * alpha * Result[j]['u'] / m * besselj(1, m * alpha))

        # Calculate strain values for each point
        for j in range(len(Points)):
            # ii = Points[j]['layer']
            ii = Points[j]['layer'] - 1
            Results[j]['epsilon_z'] = (Results[j]['sigma_z'] - v[ii] * (Results[j]['sigma_r'] + Results[j]['sigma_t'])) / E[ii]
            Results[j]['epsilon_r'] = (Results[j]['sigma_r'] - v[ii] * (Results[j]['sigma_z'] + Results[j]['sigma_t'])) / E[ii]
            Results[j]['epsilon_t'] = (Results[j]['sigma_t'] - v[ii] * (Results[j]['sigma_r'] + Results[j]['sigma_z'])) / E[ii]
            
            # Módulo de cisalhamento
            G = E[ii] / (2 * (1 + v[ii]))

            # Deformações de cisalhamento
            Results[j]['gamma_rz'] = Results[j]['tau_rz'] / G
            Results[j]['gamma_thetaz'] = Results[j]['tau_thetaz'] / G

        return Results

    def R_hat(self, n, Thickness, E, v, isbonded, m, Points):
        # Inicialização das variáveis
        z = np.cumsum(Thickness)  # soma cumulativa para a espessura das camadas
        H = z[-1]  # profundidade total (última camada)

        lambda_ = np.append(z / H, np.inf)

        F = np.zeros(n)
        R = np.zeros(n)

        F[0] = np.exp(-m * (lambda_[0] - 0))
        R[0] = (E[0] / E[1]) * (1 + v[1]) / (1 + v[0])

        for i in range(1, n - 1):
            F[i] = np.exp(-m * (lambda_[i] - lambda_[i - 1]))
            R[i] = (E[i] / E[i + 1]) * (1 + v[i + 1]) / (1 + v[i])

        F[n - 1] = np.exp(-m * (lambda_[n - 1] - lambda_[n - 2]))

        M = []
        for i in range(n - 1):
            M1 = np.zeros((4, 4))
            M2 = np.zeros((4, 4))
            
            if isbonded:
                # Construção de M1 e M2 para o caso de camadas coladas
                M1[0, 0] = 1
                M1[1, 0] = 1
                M1[2, 0] = 1
                M1[3, 0] = 1
                M1[0, 1] = F[i]
                M1[1, 1] = -F[i]
                M1[2, 1] = F[i]
                M1[3, 1] = -F[i]
                M1[0, 2] = -(1 - 2 * v[i] - m * lambda_[i])
                M1[1, 2] = 2 * v[i] + m * lambda_[i]
                M1[2, 2] = 1 + m * lambda_[i]
                M1[3, 2] = -(2 - 4 * v[i] - m * lambda_[i])
                M1[0, 3] = (1 - 2 * v[i] + m * lambda_[i]) * F[i]
                M1[1, 3] = (2 * v[i] - m * lambda_[i]) * F[i]
                M1[2, 3] = -(1 - m * lambda_[i]) * F[i]
                M1[3, 3] = -(2 - 4 * v[i] + m * lambda_[i]) * F[i]

                M2[0, 0] = F[i + 1]
                M2[1, 0] = F[i + 1]
                M2[2, 0] = R[i] * F[i + 1]
                M2[3, 0] = R[i] * F[i + 1]
                M2[0, 1] = 1
                M2[1, 1] = -1
                M2[2, 1] = R[i]
                M2[3, 1] = -R[i]
                M2[0, 2] = -(1 - 2 * v[i + 1] - m * lambda_[i]) * F[i + 1]
                M2[1, 2] = (2 * v[i + 1] + m * lambda_[i]) * F[i + 1]
                M2[2, 2] = (1 + m * lambda_[i]) * R[i] * F[i + 1]
                M2[3, 2] = -(2 - 4 * v[i + 1] - m * lambda_[i]) * R[i] * F[i + 1]
                M2[0, 3] = 1 - 2 * v[i + 1] + m * lambda_[i]
                M2[1, 3] = (2 * v[i + 1] - m * lambda_[i])
                M2[2, 3] = -(1 - m * lambda_[i]) * R[i]
                M2[3, 3] = -(2 - 4 * v[i + 1] + m * lambda_[i]) * R[i]
            else:
                # Construção de M1 e M2 para o caso de camadas deslizantes (sem aderência)
                M1[0, 0] = 1
                M1[1, 0] = 1
                M1[2, 0] = 1
                M1[3, 0] = 0
                M1[0, 1] = F[i]
                M1[1, 1] = F[i]
                M1[2, 1] = -F[i]
                M1[3, 1] = 0
                M1[0, 2] = -(1 - 2 * v[i] - m * lambda_[i])
                M1[1, 2] = 1 + m * lambda_[i]
                M1[2, 2] = 2 * v[i] + m * lambda_[i]
                M1[3, 2] = 0
                M1[0, 3] = (1 - 2 * v[i] + m * lambda_[i]) * F[i]
                M1[1, 3] = -(1 - m * lambda_[i]) * F[i]
                M1[2, 3] = (2 * v[i] - m * lambda_[i]) * F[i]
                M1[3, 3] = 0

                M2[0, 0] = F[i + 1]
                M2[1, 0] = R[i] * F[i + 1]
                M2[2, 0] = 0
                M2[3, 0] = F[i + 1]
                M2[0, 1] = 1
                M2[1, 1] = R[i]
                M2[2, 1] = 0
                M2[3, 1] = -1
                M2[0, 2] = -(1 - 2 * v[i + 1] - m * lambda_[i]) * F[i + 1]
                M2[1, 2] = (1 + m * lambda_[i]) * R[i] * F[i + 1]
                M2[2, 2] = 0
                M2[3, 2] = (2 * v[i + 1] + m * lambda_[i]) * F[i + 1]
                M2[0, 3] = 1 - 2 * v[i + 1] + m * lambda_[i]
                M2[1, 3] = -(1 - m * lambda_[i]) * R[i]
                M2[2, 3] = 0
                M2[3, 3] = 2 * v[i + 1] - m * lambda_[i]
            
            # ESTA PARTE NÃO EXISTIA NO CÓDIGO ORIGINAL MAS É NECESSARIAS PARA EVITAR PROBLEMAS COM MATRIZES SINGULARES ------------------------------
            if np.linalg.det(M1) == 0:
                # Parte nova

                # Tentativa 1
                Mi = np.linalg.lstsq(M1, M2, rcond=None)[0]

                # Tentativa 2
                # Mi = np.linalg.pinv(M1) @ M2

                # Tentativa 3
                # ruido = np.random.normal(0, 1e-2, M1.shape)
                # M1_com_ruido = M1 + ruido
                # Mi = np.linalg.solve(M1_com_ruido, M2)

            else:
                # Código original
                Mi = np.linalg.solve(M1, M2)

            M.append(Mi)
            # M.append(np.linalg.lstsq(M1, M2, rcond=None)[0])
            # ----------------------------------------------------------------------------------------------------------------------------------------

        # Multiplicação das matrizes M
        MM = np.eye(4)
        for i in range(n - 1):
            MM = np.dot(MM, M[i])
        MM = MM[:, [1, 3]]

        # Cálculo de A, B, C, D
        b11 = np.exp(-lambda_[0] * m)
        b21 = np.exp(-lambda_[0] * m)
        b12 = 1
        b22 = -1
        c11 = -(1 - 2 * v[0]) * np.exp(-m * lambda_[0])
        c21 = 2 * v[0] * np.exp(-m * lambda_[0])
        c12 = 1 - 2 * v[0]
        c22 = 2 * v[0]

        k11 = b11 * MM[0, 0] + b12 * MM[1, 0] + c11 * MM[2, 0] + c12 * MM[3, 0]
        k12 = b11 * MM[0, 1] + b12 * MM[1, 1] + c11 * MM[2, 1] + c12 * MM[3, 1]
        k21 = b21 * MM[0, 0] + b22 * MM[1, 0] + c21 * MM[2, 0] + c22 * MM[3, 0]
        k22 = b21 * MM[0, 1] + b22 * MM[1, 1] + c21 * MM[2, 1] + c22 * MM[3, 1]

        A = np.zeros(n)
        B = np.zeros(n)
        C = np.zeros(n)
        D = np.zeros(n)

        A[n - 1] = 0
        B[n - 1] = k22 / (k11 * k22 - k12 * k21)
        C[n - 1] = 0
        D[n - 1] = 1 / (k12 - k22 * k11 / k21)

        for i in range(n - 2, -1, -1):
            XX = np.dot(M[i], np.array([A[i + 1], B[i + 1], C[i + 1], D[i + 1]]))
            A[i], B[i], C[i], D[i] = XX

        # Cálculos de sigma_t, tau_rz, w, e u
        Results = []
        for point in Points:
            rho = point['r'] / H
            lmm = point['z'] / H
            ii = point['layer']-1

            if ii != 0:
                sigma_z = -m * besselj(0, m * rho) * (
                    (A[ii] - C[ii] * (1 - 2 * v[ii] - m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                    (B[ii] + D[ii] * (1 - 2 * v[ii] + m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1])))

                if rho == 0:
                    sigma_r = (m * besselj(0, m * rho) - m / 2) * (
                        (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                        (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1]))) + 2 * v[ii] * m * besselj(0, m * rho) * (C[ii] * np.exp(-m * (lambda_[ii] - lmm)) - D[ii] * np.exp(-m * (lmm - lambda_[ii-1])))
                    
                    sigma_t = (m / 2) * (
                        (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                        (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1]))) + 2 * v[ii] * m * besselj(0, m * rho) * (C[ii] * np.exp(-m * (lambda_[ii] - lmm)) - D[ii] * np.exp(-m * (lmm - lambda_[ii-1])))
                
                else:
                    sigma_r = (m * besselj(0, m * rho) - besselj(1, m * rho) / rho) * (
                        (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                        (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1]))) + 2 * v[ii] * m * besselj(0, m * rho) * (C[ii] * np.exp(-m * (lambda_[ii] - lmm)) - D[ii] * np.exp(-m * (lmm - lambda_[ii-1])))

                    sigma_t = besselj(1, m * rho) / rho * (
                        (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                        (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1]))) + 2 * v[ii] * m * besselj(0, m * rho) * (C[ii] * np.exp(-m * (lambda_[ii] - lmm)) - D[ii] * np.exp(-m * (lmm - lambda_[ii-1])))

                tau_rz = m * besselj(1, m * rho) * (
                    (A[ii] + C[ii] * (2 * v[ii] + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) -
                    (B[ii] - D[ii] * (2 * v[ii] - m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1])))
                
                w = -H * (1 + v[ii]) / E[ii] * besselj(0, m * rho) * (
                    (A[ii] - C[ii] * (2 - 4 * v[ii] - m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) -
                    (B[ii] + D[ii] * (2 - 4 * v[ii] + m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1])))

                u = H * (1 + v[ii]) / E[ii] * besselj(1, m * rho) * (
                    (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                    (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1])))

                # ------------------- TAU_THETAZ GERADO POR IA -------------------
                tau_thetaz = m * besselj(1, m * rho) * (
                    (A[ii] + C[ii] * (m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) -
                    (B[ii] - D[ii] * (-m * lmm)) * np.exp(-m * (lmm - lambda_[ii-1])))
                
            else:
                sigma_z = -m * besselj(0, m * rho) * (
                    (A[ii] - C[ii] * (1 - 2 * v[ii] - m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                    (B[ii] + D[ii] * (1 - 2 * v[ii] + m * lmm)) * np.exp(-m * (lmm - 0)))

                if rho == 0:
                    sigma_r = (m * besselj(0, m * rho) - m / 2) * (
                        (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                        (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - 0))) + 2 * v[ii] * m * besselj(0, m * rho) * (C[ii] * np.exp(-m * (lambda_[ii] - lmm)) - D[ii] * np.exp(-m * (lmm - 0)))
                    
                    sigma_t = (m / 2) * (
                        (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                        (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - 0))) + 2 * v[ii] * m * besselj(0, m * rho) * (C[ii] * np.exp(-m * (lambda_[ii] - lmm)) - D[ii] * np.exp(-m * (lmm - 0)))

                else:
                    sigma_r = (m * besselj(0, m * rho) - besselj(1, m * rho) / rho) * (
                        (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                        (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - 0))) + 2 * v[ii] * m * besselj(0, m * rho) * (C[ii] * np.exp(-m * (lambda_[ii] - lmm)) - D[ii] * np.exp(-m * (lmm - 0)))
                    
                    sigma_t = besselj(1, m * rho) / rho * (
                        (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                        (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - 0))) + 2 * v[ii] * m * besselj(0, m * rho) * (C[ii] * np.exp(-m * (lambda_[ii] - lmm)) - D[ii] * np.exp(-m * (lmm - 0)))

                tau_rz = m * besselj(1, m * rho) * (
                    (A[ii] + C[ii] * (2 * v[ii] + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) -
                    (B[ii] - D[ii] * (2 * v[ii] - m * lmm)) * np.exp(-m * (lmm - 0)))
                
                w = -H * (1 + v[ii]) / E[ii] * besselj(0, m * rho) * (
                    (A[ii] - C[ii] * (2 - 4 * v[ii] - m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) -
                    (B[ii] + D[ii] * (2 - 4 * v[ii] + m * lmm)) * np.exp(-m * (lmm - 0)))
                
                u = H * (1 + v[ii]) / E[ii] * besselj(1, m * rho) * (
                    (A[ii] + C[ii] * (1 + m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) +
                    (B[ii] - D[ii] * (1 - m * lmm)) * np.exp(-m * (lmm - 0)))

                # ------------------- TAU_THETAZ GERADO POR IA -------------------
                tau_thetaz = m * besselj(1, m * rho) * (
                    (A[ii] + C[ii] * (m * lmm)) * np.exp(-m * (lambda_[ii] - lmm)) -
                    (B[ii] - D[ii] * (-m * lmm)) * np.exp(-m * (lmm - 0)))

            Results.append({'sigma_z': sigma_z, 'sigma_r': sigma_r, 'sigma_t': sigma_t, 'tau_rz': tau_rz, 'w': w, 'u': u, 'tau_thetaz': 0})

        return Results

