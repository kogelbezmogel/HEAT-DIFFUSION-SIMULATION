import imageio
import os
from functools import cmp_to_key


def c_comp(file1, file2):
    return int(file1.split('_')[1].split(".")[0]) - int(file2.split("_")[1].split(".")[0])


if __name__ == "__main__":

    result_path = "./PLOTS"
    filenames = [os.path.join(result_path, path) for path in os.listdir(result_path)]
    filenames.sort(key=cmp_to_key(c_comp))
    images = []

    with imageio.get_writer('./GIFS/heat2D_0.gif', mode='I', fps=30) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)