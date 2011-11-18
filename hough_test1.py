import numpy as np
import scikits.image.transform.hough as scikits_hough
import datetime

img = N.zeros((1024,1024),dtype=float)

start = datetime.datetime.utcnow()

end = datetime.datetime.utcnow()

print 'Time taken for scikits_hough '
print end - start


start = datetime.datetime.utcnow()
out,angles,d = houghtf(img)
end = datetime.datetime.utcnow()

print 'Time taken for houghtf '
print end - start





