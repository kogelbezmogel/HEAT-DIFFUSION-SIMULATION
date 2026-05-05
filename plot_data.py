import imageio
import pandas as pd
import matplotlib.pyplot as plt
import os



if __name__ == "__main__":

    base_path = "./RESULT"
    result_path = "./PLOTS"
    files = [path for path in os.listdir(base_path)]

    for path in files: 
        
        data = pd.read_csv(os.path.join(base_path, path), delimiter=";", header=None)
        plt.imshow(data, cmap="plasma", interpolation="nearest", vmin=0, vmax=100)
        plt.savefig(os.path.join(result_path, path.split(".")[0] + ".png"))
        plt.close()

