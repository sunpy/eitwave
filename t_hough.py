from scikits.image.transform import hough, probabilistic_hough
from scikits.image.filter import canny
from scikits.image import data

import numpy as np
import matplotlib.pyplot as plt

def htLine(distance,angle,img):
    shape = img.shape
    nx = shape[0]
    ny = shape[1]
    eps = 1.0/float(ny)

    if abs(np.sin(angle)) > eps:
        gradient = - np.cos(angle) / np.sin(angle)
        constant = distance / np.sin(angle)
        for x in range(0,nx):
            y = gradient*x + constant
            if y <= ny-1 and y >= 0:
                img[x,y] = 255
    else:
        img[distance,:] = 255
    

    return np.transpose(img)


# Construct test image

image = np.zeros((100, 100))


# Classic straight-line Hough transform

idx = np.arange(25, 75)
image[idx[::-1], idx] = 255
image[idx, idx] = 255
image[:,90] = 255

h, theta, d = hough(image)

plt.figure(figsize=(12, 5))

plt.subplot(131)
plt.imshow(image, cmap=plt.cm.gray)
plt.title('Input image')

plt.subplot(132)
plt.imshow(np.log(1 + h),
           extent=[np.rad2deg(theta[-1]), np.rad2deg(theta[0]),
                   d[-1], d[0]],
           cmap=plt.cm.gray, aspect=1/1.5)
plt.title('Hough transform')
plt.xlabel('Angles (degrees)')
plt.ylabel('Distance (pixels)')

votethresh = 20
indices =  (h >votethresh).nonzero()

distances = d[indices[0]]
theta = theta[indices[1]]
n =len(indices[1])
invTransform = np.zeros((100, 100))
for i in range(0,n):
    nextLine = htLine( distances[i],theta[i], np.zeros((100, 100))   )
    invTransform = invTransform + nextLine

plt.subplot(133)
plt.imshow(invTransform)
plt.show()

# Line finding, using the Probabilistic Hough Transform

image = data.camera()
edges = canny(image, 2, 1, 25)
lines = probabilistic_hough(edges, threshold=10, line_length=5, line_gap=3)

plt.figure(figsize=(12, 4))

plt.subplot(131)
plt.imshow(image, cmap=plt.cm.gray)
plt.title('Input image')

plt.subplot(132)
plt.imshow(edges, cmap=plt.cm.gray)
plt.title('Sobel edges')

plt.subplot(133)
plt.imshow(edges * 0)

for line in lines:
    p0, p1 = line
    plt.plot((p0[0], p1[0]), (p0[1], p1[1]))

plt.title('Lines found with PHT')
plt.axis('image')
plt.show()
