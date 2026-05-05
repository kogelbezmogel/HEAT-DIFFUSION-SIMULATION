from PIL import Image
import os
from functools import cmp_to_key


def c_comp(file1, file2):
    return int(file1.split('_')[1].split(".")[0]) - int(file2.split("_")[1].split(".")[0])


if __name__ == "__main__":

    result_path = "./PLOTS"
    filenames = [os.path.join(result_path, path) for path in os.listdir(result_path)]
    filenames.sort(key=cmp_to_key(c_comp))
    images = [ Image.open(path) for path in filenames]

    images[0].save("./GIFS/heat_2D.gif", save_all=True, append_images=images[1:], duration=100, loop=0)


    