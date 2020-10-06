from local_alignment import local_alignment
from global_alignment import global_alignment
import pandas as pd


if __name__ == "__main__":
    protein1 = "AAGACTGAG"
    protein2 = "GAG"

    matrix = pd.read_csv("/Users/hlibisev/Documents/GitHub/bioinformatics-algorithms/alignments/BLOSUM")
    print("Локальное выравнивание с BLOSUM: ", local_alignment(protein1, protein2, -4, matrix))
    # Вывод: Локальное выравнивание с BLOSUM:  ('AAGACTGAG', '______GAG')

    print("Глобальное выравнивание с BLOSUM: ", global_alignment(protein1, protein2, -4, matrix))
    # Вывод: Глобальное выравнивание с BLOSUM:  ('AAGACTGAG', '__GA__G__')
