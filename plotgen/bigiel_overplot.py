""" This is the initial test version for the automatic bigiel plot (previously
    I have organised this by hand which is **NOT** ideal """

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

original = Image.open("bigiel_overplot.png")

plt.imshow(original)
plt.show()

