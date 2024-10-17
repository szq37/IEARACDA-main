import numpy as np
def IE(a):
    RM = np.loadtxt('../data/circRNA_disease_asso.txt')
    rs = a**2
    rn = -a
    rw = 1
    def compute_circRNA_similarity_matrix(RM, rs, rn, rw):
        n = RM.shape[0]
        circRNA_similarity_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                col1 = RM[i,:]
                col2 = RM[j,:]
                similarity = sum(compute_similarity_reward_matrix(col1, col2, rs, rn, rw) for col1, col2 in zip(col1, col2))
                circRNA_similarity_matrix[i, j] = similarity
        return circRNA_similarity_matrix
    def compute_disease_similarity_matrix(RM, rs, rn, rw):
        n = RM.shape[1]
        # print(n)
        disease_similarity_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                col1 = RM[:,i]
                col2 = RM[:,j]
                similarity = sum(compute_similarity_reward_matrix(col1, col2, rs, rn, rw) for col1, col2 in zip(col1, col2))
                disease_similarity_matrix[i, j] = similarity
        return disease_similarity_matrix
    def calculate_similarity_matrix(RM):
        n = len(RM)
        sim_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                min_val = np.min(RM)
                max_val = np.max(RM)
                if i==j:
                    sim_matrix[i,j]=1
                else:
                    sim_matrix[i, j] = (RM[i, j] - min_val) / (max_val - min_val)
        return sim_matrix