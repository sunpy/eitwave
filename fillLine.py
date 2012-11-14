def fillLine(pos1,pos2,img):
    shape=img.shape
    ny = shape[0]
    nx = shape[1]
    if pos2[0] == pos1[0]:
        m = 9999
    else:
        m = (pos2[1] - pos1[1]) / (pos2[0] - pos1[0])
        
    constant = (pos2[1] - m*pos2[0])
    
    for x in range(pos1[0],pos2[0]):
        y = m*x + constant
        if y <= ny-1 and y>= 0:
            img[y,x] = 255

    return img
        
