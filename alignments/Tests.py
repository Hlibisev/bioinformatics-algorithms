from local_alignment import local_alignment
from global_alignment import global_alignment
import pandas as pd


if __name__ == "__main__":
    protein1 = "GGGACTGAG"
    protein2 = "GACT"

    matrix = pd.read_csv("/Users/hlibisev/Documents/GitHub/bioinformatics-algorithms/alignments/BLOSUM")
    print("Локальное выравнивание с BLOSUM: ", local_alignment(protein1, protein2, -4, matrix))
    # Вывод: Локальное выравнивание с BLOSUM:  ('GGGACTGAG_', '__GACT____')

    print("Глобальное выравнивание с BLOSUM: ", global_alignment(protein1, protein2, -4, matrix))
    # Вывод: Глобальное выравнивание с BLOSUM:  ('GGGACTGAG', 'G__ACT___')
