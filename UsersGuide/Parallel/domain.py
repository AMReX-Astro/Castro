import matplotlib.pyplot as plt

def drawBox(ll, uu, nx=1, ny=1, gridColor="0.5", ng=0):

    # draw the frame
    plt.plot([ll[0], ll[0], uu[0], uu[0], ll[0]],
               [ll[1], uu[1], uu[1], ll[1], ll[1]], color="k", lw=2)

    if nx > 1:
        # draw the x grid lines
        dx = (uu[0] - ll[0])/nx
        for n in range(nx):
            plt.plot([ll[0]+n*dx, ll[0]+n*dx],
                     [ll[1], uu[1]], color=gridColor, ls=":")

    if ny > 1:
        # draw the y grid lines
        dy = (uu[1] - ll[1])/ny
        for n in range(ny):
            plt.plot([ll[0], uu[0]],
                     [ll[1]+n*dy, ll[1]+n*dy], color=gridColor, ls=":")
    
    # ghostcells?  
    if (ng > 0):
        print "here"
        xmin = ll[0]-ng*dx
        xmax = uu[0]+ng*dx
        ymin = ll[1]-ng*dy
        ymax = uu[1]+ng*dy
        plt.plot([xmin, xmin, xmax, xmax, xmin], 
                   [ymin, ymax, ymax, ymin, ymin], 
                   ls="--", color="r")


if __name__ == "__main__":

    # define our boxes
    boxes = []
    boxes.append([(0, 20), (20,40)])
    boxes.append([(10,13), (32,20)])
    boxes.append([(20,20), (40,28)])
    boxes.append([(32,0),  (52,20)])

    dl = 3
    
    plt.clf()
    
    for e, b in enumerate(boxes):
        lo, hi = b        
        drawBox([lo[0], lo[1]], [hi[0], hi[1]])
        plt.text(lo[0]+dl, hi[1]-dl, "{}".format(e))
        
    a = plt.gca()
    a.set_aspect("equal", "datalim")
    plt.axis("off")

    plt.ylim(-5, 45)
    
    f = plt.gcf()
    f.set_size_inches(6.0,4.25)

    plt.savefig("domain.png")
    plt.savefig("domain.eps", bbox_inches="tight")


    # now with tiling
    plt.clf()
    
    ntile = 8
    
    for e, b in enumerate(boxes):
        lo, hi = b

        nx = max(1,(hi[0]-lo[0])/ntile)
        ny = max(1,(hi[1]-lo[1])/ntile)

        print nx, ny
        drawBox([lo[0], lo[1]], [hi[0], hi[1]], nx=nx, ny=ny)
        plt.text(lo[0]+dl, hi[1]-dl, "{}".format(e))
        

    a = plt.gca()
    a.set_aspect("equal", "datalim")
    plt.axis("off")

    plt.xlim(-1,54)
    plt.ylim(-1, 41)
    
    f = plt.gcf()
    f.set_size_inches(6.0,4.25)

    plt.savefig("domain-tile.png")
    plt.savefig("domain-tile.eps", bbox_inches="tight")
        
