#!/usr/bin/env python3

"""
서열 정렬 알고리즘 구현
- Needleman-Wunsch (전역 정렬)
- Smith-Waterman (국소 정렬)
"""

import numpy as np

class SequenceAlignment:
    def __init__(self, match=1, mismatch=-1, gap=-2):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap

    def _score(self, a, b):
        return self.match if a == b else self.mismatch

    def needleman_wunsch(self, seq1, seq2):
        """Needleman-Wunsch 알고리즘을 이용한 전역 정렬"""
        # 행렬 초기화
        m, n = len(seq1), len(seq2)
        score_matrix = np.zeros((m + 1, n + 1))
        traceback = np.zeros((m + 1, n + 1), dtype=str)

        # 첫 행과 열 초기화
        for i in range(m + 1):
            score_matrix[i, 0] = i * self.gap
            traceback[i, 0] = 'up'
        for j in range(n + 1):
            score_matrix[0, j] = j * self.gap
            traceback[0, j] = 'left'

        # 점수 행렬 채우기
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = score_matrix[i-1, j-1] + self._score(seq1[i-1], seq2[j-1])
                delete = score_matrix[i-1, j] + self.gap
                insert = score_matrix[i, j-1] + self.gap
                
                score_matrix[i, j] = max(match, delete, insert)
                
                if score_matrix[i, j] == match:
                    traceback[i, j] = 'diag'
                elif score_matrix[i, j] == delete:
                    traceback[i, j] = 'up'
                else:
                    traceback[i, j] = 'left'

        # Traceback
        align1, align2 = [], []
        i, j = m, n
        while i > 0 or j > 0:
            if traceback[i, j] == 'diag':
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif traceback[i, j] == 'up':
                align1.append(seq1[i-1])
                align2.append('-')
                i -= 1
            else:
                align1.append('-')
                align2.append(seq2[j-1])
                j -= 1

        return ''.join(reversed(align1)), ''.join(reversed(align2))

    def smith_waterman(self, seq1, seq2):
        """Smith-Waterman 알고리즘을 이용한 국소 정렬"""
        # 행렬 초기화
        m, n = len(seq1), len(seq2)
        score_matrix = np.zeros((m + 1, n + 1))
        traceback = np.zeros((m + 1, n + 1), dtype=str)

        # 최대 점수 위치 추적
        max_score = 0
        max_pos = (0, 0)

        # 점수 행렬 채우기
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = score_matrix[i-1, j-1] + self._score(seq1[i-1], seq2[j-1])
                delete = score_matrix[i-1, j] + self.gap
                insert = score_matrix[i, j-1] + self.gap
                
                score_matrix[i, j] = max(0, match, delete, insert)
                
                if score_matrix[i, j] == match:
                    traceback[i, j] = 'diag'
                elif score_matrix[i, j] == delete:
                    traceback[i, j] = 'up'
                elif score_matrix[i, j] == insert:
                    traceback[i, j] = 'left'
                
                if score_matrix[i, j] > max_score:
                    max_score = score_matrix[i, j]
                    max_pos = (i, j)

        # Traceback
        align1, align2 = [], []
        i, j = max_pos
        
        while i > 0 and j > 0 and score_matrix[i, j] > 0:
            if traceback[i, j] == 'diag':
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
            elif traceback[i, j] == 'up':
                align1.append(seq1[i-1])
                align2.append('-')
                i -= 1
            else:
                align1.append('-')
                align2.append(seq2[j-1])
                j -= 1

        return ''.join(reversed(align1)), ''.join(reversed(align2))

def main():
    # 예시 서열
    seq1 = "GCATGCU"
    seq2 = "GATTACA"

    # 정렬 객체 생성
    aligner = SequenceAlignment()

    # Needleman-Wunsch 전역 정렬
    print("Needleman-Wunsch Global Alignment:")
    align1, align2 = aligner.needleman_wunsch(seq1, seq2)
    print(align1)
    print(align2)
    print()

    # Smith-Waterman 국소 정렬
    print("Smith-Waterman Local Alignment:")
    align1, align2 = aligner.smith_waterman(seq1, seq2)
    print(align1)
    print(align2)

if __name__ == "__main__":
    main()