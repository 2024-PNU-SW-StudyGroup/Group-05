from .utils import calc_distances  # setup.py에 npmodule.utils

def print_distances(a, b):
    ed, md = calc_distances(a, b)
    print(f"Euclidean distance: {ed}")
    print(f"Manhattan distance: {md}")