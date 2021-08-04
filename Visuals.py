import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


class Plotter3D:
    #use the scroll wheel to change the third dimension of an image
    def __init__(self, ax, image):
        self.ax = ax
        ax.set_title('Scroll through MRI slices')

        self.image = image
        rows, cols, self.slices = image.shape
        self.ind = self.slices//2

        self.im = ax.imshow(self.image[:, :, self.ind])
        self.update()

    def on_scroll(self, event):
        print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = min((self.ind + 1), self.slices-1)
        else:
            self.ind = max((self.ind - 1) , 0)
        self.update()

    def update(self):
        self.im.set_data(self.image[:, :, self.ind])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()


# fig, ax = plt.subplots(1, 1)

# X = np.random.rand(20, 20, 40)

# tracker = Plotter3D(ax, X)


# fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
# plt.show()
def SliderPlot(image, specialChannelNames):
    global fig, ax, sliderT, sliderO, array, img, specialChannels, num_bVals, previous_j
    num_bVals = image.shape[3]
    specialChannels = specialChannelNames
    array = image 
    #this function takes a 4d array with the 3rd index being the slice and the 4th index being the b value 
    fig, ax = plt.subplots()
    img = ax.imshow(array[:,:,0,0], origin = 'upper')

    axT = fig.add_axes([0.2, 0.95, 0.65, 0.03])
    axO = fig.add_axes([0.2, 0.90, 0.65, 0.03])

    sliderT = Slider(axT, 'Slice', 0, array.shape[2]-1, valinit=0, valfmt='%i')
    sliderO = Slider(axO, 'b Val', 0, array.shape[3]-1, valinit=0, valfmt='%i')
    previous_j = 0

    sliderT.on_changed(update)
    sliderO.on_changed(update)
    plt.show()

def update(val):
    i = int(sliderT.val)
    j = int(sliderO.val)
    im = array[:,:,i,j]

    img = img = ax.imshow(array[:,:,i,j], origin = 'upper')


    #Only want to reload the colour bar for special images like D, D*, f. So we must clear the axis and make a new image for these types of images. Then when we go back to regular signal images need to recreate
    #specialImageRange = range(num_bVals - 3, num_bVals - 1)

    fig.canvas.draw_idle()
    


  